#!/usr/bin/env python3
"""
Extract domain models from AF2 PDBs based on DPAM-AI predictions.

Converts full-length AlphaFold Database (AFDB) PDB models into a domain-centric
library by splitting structures according to DPAM-AI domain predictions. Each
predicted domain is extracted, capped with ACE/NME termini, renumbered to match
UniProt canonical numbering, and written as a standalone PDB file with metadata.
"""

import argparse
import io
import os
import random
import statistics
import subprocess
import sys
import tempfile
import time

import pandas as pd
import requests
from Bio import AlignIO
from Bio.PDB import Chain, Model, PDBIO, PDBParser, Structure
from MDAnalysis import Universe

# Add script directory to path for add_caps import
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)
from add_caps import add_caps

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

AA_NAME_MAPPING = {
    "VAL": "V",
    "ILE": "I",
    "LEU": "L",
    "GLU": "E",
    "GLN": "Q",
    "ASP": "D",
    "ASN": "N",
    "HIS": "H",
    "TRP": "W",
    "PHE": "F",
    "TYR": "Y",
    "ARG": "R",
    "LYS": "K",
    "SER": "S",
    "THR": "T",
    "MET": "M",
    "ALA": "A",
    "GLY": "G",
    "PRO": "P",
    "CYS": "C",
    "ACE": "X",
    "NME": "X",
}

UNIPROT_DELAY_SEC = 0.15

# -----------------------------------------------------------------------------
# Biopython / MDAnalysis helpers
# -----------------------------------------------------------------------------

_PARSER = PDBParser()


def _struct_to_universe(struct, pdb_io):
    """Convert a biopython structure to a MDAnalysis universe."""
    with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w+", delete=True) as tmp:
        pdb_io.set_structure(struct)
        pdb_io.save(tmp.name)
        tmp.flush()
        return Universe(tmp.name)


def _universe_to_struct(u):
    """Convert an MDAnalysis universe to a biopython structure."""
    with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w+", delete=True) as tmp:
        u.atoms.write(tmp.name)
        tmp.flush()
        return _PARSER.get_structure("converted", tmp.name)


def parse_range_to_residue_ids(range_str):
    """
    Parse Range string into set of residue IDs.
    '1-100' -> {1,2,...,100}
    '1-100,200-300' -> {1,2,...,100,200,...,300}
    """
    result = set()
    for part in str(range_str).split(","):
        part = part.strip()
        start, end = map(int, part.split("-"))
        result.update(range(start, end + 1))
    return result


def truncate_structure(struct, residue_ids, pdb_io):
    """
    Truncate a biopython structure to specified residue IDs (single chain),
    then add ACE and NME capping groups. Supports discontinuous ranges.
    """
    residue_ids = set(residue_ids)
    trunc_struct = Structure.Structure(struct.id)
    for model in struct:
        new_model = Model.Model(model.id)
        for chain in model:
            new_chain = Chain.Chain(chain.id)
            for res in chain:
                resnum = res.id[1]
                if resnum in residue_ids:
                    new_chain.add(res.copy())
            if len(new_chain):
                new_model.add(new_chain)
        if len(new_model):
            trunc_struct.add(new_model)
        break

    u = _struct_to_universe(trunc_struct, pdb_io)
    u = add_caps(u)
    trunc_struct = _universe_to_struct(u)
    return trunc_struct


def pdb_to_sequence(structure, chain="A"):
    """Extract sequence and residue numbers from a structure."""
    residues = list(structure[0][chain].get_residues())
    resnums = [r.id[1] for r in residues]
    if len(resnums) != len(set(resnums)):
        raise ValueError("Residue numbers should be unique")
    seq_list = []
    if resnums:
        prev_resnum = resnums[0] - 1
        for residue, resnum in zip(residues, resnums):
            gap_len = resnum - prev_resnum - 1
            if gap_len > 0:
                seq_list.extend(["-"] * gap_len)
            seq_list.append(AA_NAME_MAPPING.get(residue.resname, f"[{residue.resname}]"))
            prev_resnum = resnum
        seq = "".join(seq_list)
    else:
        seq = ""
    return seq, resnums


# -----------------------------------------------------------------------------
# UniProt helpers
# -----------------------------------------------------------------------------


def get_uniprot_entry(uniprot_id, max_retries=3, uniprot_delay=UNIPROT_DELAY_SEC):
    """Fetch full UniProt entry from REST API with retries and throttling."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    last_error = None
    for attempt in range(max_retries):
        time.sleep(uniprot_delay)
        try:
            response = requests.get(url, timeout=30)
            if response.ok:
                return response.json()
            if response.status_code in (429, 500, 502, 503):
                last_error = f"HTTP {response.status_code}"
                if attempt < max_retries - 1:
                    wait = (2 ** attempt) + random.uniform(0, 1)
                    time.sleep(wait)
                    continue
            else:
                print(
                    f"[ERROR] UniProt fetch for {uniprot_id}: {response.status_code}",
                    file=sys.stderr,
                )
                return None
        except requests.exceptions.RequestException as e:
            last_error = e
            if attempt < max_retries - 1:
                wait = (2 ** attempt) + random.uniform(0, 1)
                time.sleep(wait)
            else:
                break
    print(
        f"[ERROR] UniProt fetch for {uniprot_id} failed after {max_retries} attempts: {last_error}",
        file=sys.stderr,
    )
    return None


def _get_uniprot_sequence(entry):
    """Extract sequence from UniProt entry."""
    return entry["sequence"]["value"] if entry and "sequence" in entry else ""

# Function to get sequence from UniParc (As Fallback for UniProt)
def _get_uniparc_sequence(entry):
    uniParcId = entry['extraAttributes']['uniParcId']
    url = f"https://rest.uniprot.org/uniparc/%7Bupi%7D?upi={uniParcId}&fields=Sequence"
    response = requests.get(url)
    if not response.ok:
        print("[ERROR]", response.content)
        return None
    sequence = response.json()['sequence']['value']
    return sequence


# -----------------------------------------------------------------------------
# Alignment and renumbering
# -----------------------------------------------------------------------------


def needle_alignment(seq1, seq2):
    """Run EMBOSS needle alignment. Returns alignment string or None on failure."""
    f1_path = f2_path = outf_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fa"
        ) as f1, tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fa"
        ) as f2, tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".aln"
        ) as outf:
            f1.write(f">seq1\n{seq1}\n")
            f1.flush()
            f2.write(f">seq2\n{seq2}\n")
            f2.flush()
            f1_path, f2_path, outf_path = f1.name, f2.name, outf.name

        cmd = [
            "needle",
            "-asequence", f1_path,
            "-bsequence", f2_path,
            "-sprotein1",
            "-sprotein2",
            "-gapopen", "10",
            "-gapextend", "0.5",
            "-outfile", outf_path,
            "-auto",
        ]
        subprocess.run(cmd, check=True, capture_output=True)

        with open(outf_path) as f:
            return f.read()
    except subprocess.CalledProcessError as e:
        print(
            f"Error running EMBOSS needle: {e.stderr.decode() if e.stderr else e}",
            file=sys.stderr,
        )
        return None
    finally:
        for p in (f1_path, f2_path, outf_path):
            if p and os.path.exists(p):
                try:
                    os.unlink(p)
                except OSError:
                    pass


def map_pdb_resnums_to_uniprot(pdb_seq, pdb_resnums, uniprot_seq):
    """
    Align PDB sequence to UniProt sequence and return mapping of
    (pdb_resnum, uniprot_resnum) pairs (1-based).
    """
    alignment = needle_alignment(pdb_seq, uniprot_seq)
    if alignment is None:
        return []

    aln = AlignIO.read(io.StringIO(alignment), "emboss")
    seq1_aligned = str(aln[0].seq)
    seq2_aligned = str(aln[1].seq)

    mapping = []
    pdb_idx = 0
    uniprot_resnum = 0

    for i in range(len(seq1_aligned)):
        c1, c2 = seq1_aligned[i], seq2_aligned[i]
        if c1 == "-":
            if c2 != "-":
                uniprot_resnum += 1
            continue
        if c2 == "-":
            pdb_idx += 1
            continue
        uniprot_resnum += 1
        if c1 == "X":
            pdb_idx += 1
            continue
        mapping.append((pdb_resnums[pdb_idx], uniprot_resnum))
        pdb_idx += 1

    return mapping


def mode_offset(mapping):
    """Compute the mode of (uniprot_resnum - pdb_resnum) from mapping."""
    if not mapping:
        return None
    offsets = [u - p for p, u in mapping]
    return statistics.mode(offsets)


def renumber_structure(struct, offset):
    """Create a copy with residues renumbered by offset. None treated as 0."""
    offset = 0 if offset is None else offset
    new_struct = struct.copy()
    for model in new_struct:
        for chain in model:
            for res in chain:
                res.id = (res.id[0], res.id[1] + offset, res.id[2])
    return new_struct


def _max_residue_number(struct):
    """Return the maximum residue number in the structure (any chain)."""
    max_resid = 0
    for model in struct:
        for chain in model:
            for res in chain:
                max_resid = max(max_resid, res.id[1])
    return max_resid


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------


def run(domains_path, in_pdb_dir, out_pdb_dir, out_csv_path, max_retries=3, uniprot_delay=UNIPROT_DELAY_SEC):
    """
    Process one subset: read domains, extract domain PDBs, write CSV.
    """
    os.makedirs(out_pdb_dir, exist_ok=True)
    os.makedirs(os.path.dirname(out_csv_path) or ".", exist_ok=True)

    df_pred = pd.read_csv(domains_path, sep="\t")
    df_pred.rename(columns={"Protein": "uniprot_accn"}, inplace=True)

    df_pred["id"] = (
        df_pred["uniprot_accn"].astype(str)
        + "-F1-"
        + df_pred["Domain"].astype(str)
        + "_A"
    )
    df_pred.insert(0, "id", df_pred.pop("id"))

    header = df_pred.head(0).copy()
    header["offset"] = None
    header["pdb_seq"] = None
    header.to_csv(out_csv_path, index=False)

    parser = PDBParser()
    pdb_io = PDBIO()
    unique_uniprot = df_pred["uniprot_accn"].unique()

    for uniprot_accn in unique_uniprot:
        print(f"Splitting {uniprot_accn} into domains")

        df_pred_uniprot = df_pred[df_pred["uniprot_accn"] == uniprot_accn]
        af2_pdb_path = os.path.join(in_pdb_dir, f"{uniprot_accn}.pdb")

        if not os.path.isfile(af2_pdb_path):
            print(f"[WARN] Missing PDB: {af2_pdb_path}", file=sys.stderr)
            continue

        try:
            af2_struct = parser.get_structure(str(uniprot_accn), af2_pdb_path)
        except Exception as e:
            print(f"[ERROR] Failed to parse {af2_pdb_path}: {e}", file=sys.stderr)
            continue

        uniprot_entry = get_uniprot_entry(uniprot_accn, max_retries=max_retries, uniprot_delay=uniprot_delay)
        uniprot_seq = _get_uniprot_sequence(uniprot_entry)
        if not uniprot_seq:
            try:
                uniprot_seq = _get_uniparc_sequence(uniprot_entry)
            except Exception as e:
                print(
                    f"[ERROR] Could not fetch UniProt sequence for {uniprot_accn}: {e}",
                    file=sys.stderr,
                )
                continue

        for _, row in df_pred_uniprot.iterrows():
            domain_id = row["id"]
            domain = row["Domain"]
            domain_range = row["Range"]

            # Parse Range (supports "1-100" or "1-100,200-300")
            residue_ids = parse_range_to_residue_ids(domain_range)
            print(f"    Truncating {domain} from {domain_range}")

            try:
                domain_struct = truncate_structure(
                    af2_struct, residue_ids, pdb_io
                )
            except Exception as e:
                print(
                    f"[ERROR] Truncate failed for {domain_id}: {e}",
                    file=sys.stderr,
                )
                continue

            pdb_seq, pdb_resnums = pdb_to_sequence(domain_struct)
            row["pdb_seq"] = pdb_seq

            mapping = map_pdb_resnums_to_uniprot(pdb_seq, pdb_resnums, uniprot_seq)
            offset = mode_offset(mapping) if mapping else None

            domain_struct_renumbered = renumber_structure(domain_struct, offset)
            max_resid = _max_residue_number(domain_struct_renumbered)

            if max_resid > 9999:
                struct_to_save = domain_struct
                row["offset"] = "NA"
                row["Range"] = domain_range
            else:
                struct_to_save = domain_struct_renumbered
                effective_offset = 0 if offset is None else offset
                ranges = []
                for part in str(domain_range).split(","):
                    start, end = map(int, part.strip().split("-"))
                    ranges.append((start + effective_offset, end + effective_offset))
                row["Range"] = ",".join(f"{s}-{e}" for s, e in ranges)
                row["offset"] = offset

            domain_pdb_path = os.path.join(out_pdb_dir, f"{domain_id}.pdb")
            pdb_io.set_structure(struct_to_save)
            pdb_io.save(domain_pdb_path)

            with open(out_csv_path, "a") as f:
                row_df = row.to_frame().T
                row_df = row_df[header.columns]
                row_df.to_csv(f, header=False, index=False)


def main():
    ap = argparse.ArgumentParser(
        description="Extract domain PDBs from AF2 models using DPAM-AI predictions."
    )
    ap.add_argument(
        "--domains",
        required=True,
        help="Path to {subset_id}_domains tab-delimited file",
    )
    ap.add_argument(
        "--in_pdb_dir",
        required=True,
        help="Directory with full-length AF2 PDBs ({uniprot_id}.pdb)",
    )
    ap.add_argument(
        "--out_pdb_dir",
        required=True,
        help="Output directory for domain PDB files",
    )
    ap.add_argument(
        "--out_csv",
        required=True,
        help="Output path for domains_info.csv",
    )
    ap.add_argument(
        "--max-retries",
        type=int,
        default=3,
        help="Max retries for UniProt fetch (default: 3)",
    )
    ap.add_argument(
        "--uniprot-delay",
        type=float,
        default=UNIPROT_DELAY_SEC,
        help="Delay in seconds between UniProt requests (default: 0.15)",
    )
    args = ap.parse_args()

    if not os.path.isfile(args.domains):
        print(f"[ERROR] Domains file not found: {args.domains}", file=sys.stderr)
        sys.exit(1)

    run(
        args.domains,
        args.in_pdb_dir,
        args.out_pdb_dir,
        args.out_csv,
        max_retries=args.max_retries,
        uniprot_delay=args.uniprot_delay,
    )


if __name__ == "__main__":
    main()
