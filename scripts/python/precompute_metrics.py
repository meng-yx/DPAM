#!/usr/bin/env python3
"""
Precompute domain metrics for DPAM-AI AFDB domainome.

Reads an input CSV, computes domain metrics in a single pass with shared caches
for UniProt entries and PDB structures, and writes the enriched CSV to output.
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import mdtraj as md
from tqdm import tqdm
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# Add script directory to path for imports
_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.python.handle_gaps import (
    _parse_range_fragment,
    handle_gaps_pdb_seq,
    UnlinkableGapError,
)
from scripts.python.DPAM_domains_to_pdb import (
    get_uniprot_entry,
    _get_uniprot_sequence,
    _get_uniparc_sequence,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# Default paths
# -----------------------------------------------------------------------------

DEFAULT_DEEPTMHMM_ROOT = "/work/upthomae/Meng/hp_list_DeepTMHMM/out/"
DEFAULT_MEDIAN_LINKER_CSV = _REPO_ROOT / "notebooks" / "median_linker_length.csv"

# -----------------------------------------------------------------------------
# Helper: parse_range_to_resi
# -----------------------------------------------------------------------------


def parse_range_to_resi(range_str: str) -> List[int]:
    """Parse Range string (e.g. '86-400' or '1-90,141-180,206-220' or '-16-213') to list of residue IDs."""
    resi = []
    for part in str(range_str).split(","):
        part = part.strip()
        if not part:
            continue
        try:
            start, end = _parse_range_fragment(part)
            resi.extend(range(start, end + 1))
        except ValueError:
            resi.append(int(part))
    return resi


# -----------------------------------------------------------------------------
# Step 2: DeepTMHMM
# -----------------------------------------------------------------------------


def parse_deeptmhmm(
    row: pd.Series,
    deeptmhmm_root: Path,
    placeholder: str = "?",
) -> Optional[str]:
    """
    Per-residue labels: S (signal), I (inside), M (alpha mem), B (beta mem), P (periplasm), O (outside).
    Returns annotation for domain residues (from Range) only.
    """
    file = deeptmhmm_root / row["uniprot_accn"] / "predicted_topologies.3line"
    if not file.exists():
        return None

    with open(file, "r") as f:
        header, sequence, annotation = f.read().strip().splitlines()
        sequence = np.array(list(sequence))
        annotation = np.array(list(annotation))

    if len(sequence) != len(annotation):
        return None

    domain_indices = np.array(parse_range_to_resi(row["Range"]))
    global_index = domain_indices - 1  # 0-based
    if np.any(global_index < 0) or np.any(global_index >= len(annotation)):
        return None
    domain_labels = annotation[global_index]
    if not np.isin(domain_labels, ["B", "P", "M", "I", "O", "S", placeholder]).all():
        return None
    return "".join(domain_labels)


def intracellular_fraction_from_annotation(annotation: Optional[str]) -> Optional[float]:
    """Fraction of 'I' (inside cell) in deeptmhmm annotation."""
    if annotation is None or not isinstance(annotation, str):
        return None
    if len(annotation) == 0:
        return None
    return annotation.count("I") / len(annotation)


# -----------------------------------------------------------------------------
# Step 3: mdtraj DSSP (identical to notebook)
# -----------------------------------------------------------------------------


def get_mdtraj_annotation(pdb_path: str) -> list:
    """
    Annotations:
    H: Helix (H, G or I in the original DSSP code)
    E: Strand (E or B in the original DSSP code)
    C: Coil (T or S in the original DSSP code)
    NA: Assigned if it's not a protein residue
    """
    traj = md.load(pdb_path)
    dssp_codes = md.compute_dssp(traj, simplified=True)
    dssp_codes = dssp_codes[0]
    return dssp_codes


def find_sse(resi: List[int], annotation: List[str]) -> List[Dict[str, Any]]:
    """Find contiguous segments with the same secondary structure."""
    segments = []
    for i, label in zip(resi, annotation):
        if (len(segments) > 0) and (
            i == segments[-1]["end"] + 1 and label == segments[-1]["label"]
        ):
            segments[-1]["end"] = i
        else:
            segments.append({"start": i, "end": i, "label": label})
    return segments


def mdtraj_dssp_string(row: pd.Series, pdb_dir: Path) -> Optional[str]:
    """Compute per-residue DSSP annotation (H, E, C) for domain PDB. Excludes ACE/NME caps."""
    pdb_path = pdb_dir / f"{row['id']}.pdb"
    try:
        labels = get_mdtraj_annotation(str(pdb_path))
        labels_str = "".join([x if x in {"H", "E", "C"} else "?" for x in labels])
        return labels_str
    except Exception:
        return None


def secondary_structure_metrics(row) -> Tuple[Optional[int], Optional[int], Optional[int], Optional[float], Optional[float]]:
    """
    Compute total_n_sse, n_helix, n_strand, helix_frac, strand_frac from mdtraj_dssp_annotation.
    DSSP: H=helix, E=strand, C=coil.
    """
    if row.get("mdtraj_dssp_annotation") is None or pd.isna(row.get("mdtraj_dssp_annotation")):
        return None, None, None, None, None

    dssp_labels = list(row["mdtraj_dssp_annotation"])
    residue_ids = parse_range_to_resi(row["Range"])

    segments = find_sse(residue_ids, dssp_labels)

    # count number of characters excluding ?
    length = len(dssp_labels) - dssp_labels.count("?")
    if length == 0:
        return None, None, None, None, None

    total_n_sse = int(len(segments))
    n_helix = int(sum(s["label"] == "H" for s in segments))
    n_strand = int(sum(s["label"] == "E" for s in segments))
    helix_frac = np.sum([x == "H" for x in dssp_labels]) / length
    strand_frac = np.sum([x == "E" for x in dssp_labels]) / length

    return total_n_sse, n_helix, n_strand, helix_frac, strand_frac


# -----------------------------------------------------------------------------
# Step 4: SASA
# -----------------------------------------------------------------------------


def get_sasa(structure, sub_surface_key=None) -> float:
    """Compute solvent-accessible surface area using ShrakeRupley."""
    computed = False
    for res in structure.get_atoms():
        computed = hasattr(res, "sasa")
        break
    if not computed:
        ShrakeRupley().compute(structure, level="A")
    if sub_surface_key is not None:
        structure = structure[sub_surface_key]
    return sum(a.sasa for a in structure.get_atoms())


# -----------------------------------------------------------------------------
# Step 5: protein_name, gene_name
# -----------------------------------------------------------------------------


def _get_protein_name(entry: Optional[Dict]) -> Optional[str]:
    if entry is None:
        return None
    try:
        if "recommendedName" not in entry["proteinDescription"]:
            names = entry["proteinDescription"].get("submissionNames", [])
            if len(names) == 1:
                return names[0]["fullName"]["value"]
            return None
        return entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except (KeyError, TypeError):
        return None


def _get_gene_name(entry: Optional[Dict]) -> Optional[str]:
    if entry is None:
        return None
    try:
        return entry["genes"][0]["geneName"]["value"]
    except (KeyError, TypeError, IndexError):
        return None


# -----------------------------------------------------------------------------
# Caches
# -----------------------------------------------------------------------------


def _get_uniprot_sequence_from_entry(entry: Optional[Dict]) -> Optional[str]:
    """Get sequence from UniProt entry, with UniParc fallback."""
    if entry is None:
        return None
    seq = _get_uniprot_sequence(entry)
    if not seq:
        try:
            seq = _get_uniparc_sequence(entry)
        except Exception:
            return None
    return seq if seq else None


# -----------------------------------------------------------------------------
# Main processing
# -----------------------------------------------------------------------------


def process_row(
    row: pd.Series,
    pdb_dir: Path,
    df_median_linker_length: pd.DataFrame,
    uniprot_cache: Dict[str, Optional[Dict]],
    pdb_structure_cache: Dict[str, Any],
    deeptmhmm_root: Optional[Path],
    min_gap_size: int,
    uniprot_timeout: int = 30,
    uniprot_max_retries: int = 3,
) -> Dict[str, Any]:
    """
    Process a single row through all 5 steps. Returns dict of new column values.
    """
    result = {}
    domain_id = row["id"]
    uniprot_accn = row["uniprot_accn"]
    pdb_path = pdb_dir / f"{domain_id}.pdb"
    pdb_path_str = str(pdb_path)

    try:
        # Resolve caches
        if uniprot_accn not in uniprot_cache:
            uniprot_cache[uniprot_accn] = get_uniprot_entry(
                uniprot_accn,
                max_retries=uniprot_max_retries,
                timeout=uniprot_timeout,
            )
        uniprot_entry = uniprot_cache[uniprot_accn]

        if pdb_path_str not in pdb_structure_cache:
            if not pdb_path.exists():
                pdb_structure_cache[pdb_path_str] = None
            else:
                pdb_structure_cache[pdb_path_str] = PDBParser(QUIET=True).get_structure(
                    "protein", pdb_path_str
                )
        pdb_structure = pdb_structure_cache[pdb_path_str]

        # Step 1: domain_seq
        full_sequence = _get_uniprot_sequence_from_entry(uniprot_entry)
        if not full_sequence:
            result.update({
                "full_sequence": None,
                "num_gaps": None,
                "linker_types": None,
                "linker_seqs": None,
                "domain_seq": None,
            })
        else:
            try:
                step1 = handle_gaps_pdb_seq(
                    pdb_seq=row["pdb_seq"],
                    full_sequence=full_sequence,
                    range_str=row["Range"],
                    domain_id=domain_id,
                    domain_pdb_path=pdb_path_str,
                    df_median_linker_length=df_median_linker_length,
                    min_gap_size=min_gap_size,
                    seed=hash(domain_id) % (2**32),
                    structure=pdb_structure,
                )
                result.update(step1)
            except (UnlinkableGapError, FileNotFoundError, ValueError) as e:
                logger.warning("Step 1 failed for %s: %s", domain_id, e)
                result.update({
                    "full_sequence": full_sequence,
                    "num_gaps": None,
                    "linker_types": None,
                    "linker_seqs": None,
                    "domain_seq": None,
                })

        # Step 2: deeptmhmm
        if deeptmhmm_root is not None:
            ann = parse_deeptmhmm(row, deeptmhmm_root)
            result["deeptmhmm_annotation"] = ann
            result["intracellular_fraction"] = intracellular_fraction_from_annotation(ann)
        else:
            result["deeptmhmm_annotation"] = None
            result["intracellular_fraction"] = None

        # Step 3: DSSP (uses get_mdtraj_annotation/mdtraj_dssp_string/secondary_structure_metrics from notebook)
        result["mdtraj_dssp_annotation"] = mdtraj_dssp_string(row, pdb_dir)
        row_with_dssp = {**row.to_dict(), **result}
        sse = secondary_structure_metrics(row_with_dssp)
        result["total_n_sse"], result["n_helix"], result["n_strand"], result["helix_frac"], result["strand_frac"] = sse

        # Step 4: SASA
        if pdb_structure is not None:
            try:
                sasa = get_sasa(pdb_structure)
                num_residues = sum(1 for _ in pdb_structure.get_residues())
                norm_sasa = sasa / num_residues if num_residues > 0 else None
                result["sasa"] = sasa
                result["norm_sasa"] = norm_sasa
            except Exception as e:
                logger.warning("Step 4 (SASA) failed for %s: %s", domain_id, e)
                result["sasa"] = result["norm_sasa"] = None
        else:
            result["sasa"] = result["norm_sasa"] = None

        # Step 5: protein_name, gene_name
        result["protein_name"] = _get_protein_name(uniprot_entry)
        result["gene_name"] = _get_gene_name(uniprot_entry)

    except Exception as e:
        logger.warning("Row %s failed: %s", domain_id, e)
        # Fill with None for any missing keys
        for key in [
            "full_sequence", "num_gaps", "linker_types", "linker_seqs", "domain_seq",
            "deeptmhmm_annotation", "intracellular_fraction",
            "mdtraj_dssp_annotation", "total_n_sse", "n_helix", "n_strand", "helix_frac", "strand_frac",
            "sasa", "norm_sasa",
            "protein_name", "gene_name",
        ]:
            if key not in result:
                result[key] = None

    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Precompute domain metrics for DPAM-AI AFDB domainome."
    )
    parser.add_argument("--input", required=True, type=Path, help="Path to input CSV")
    parser.add_argument("--output", required=True, type=Path, help="Path to output CSV")
    parser.add_argument(
        "--pdb-dir",
        required=True,
        type=Path,
        help="Directory containing {id}.pdb files",
    )
    parser.add_argument(
        "--deeptmhmm-root",
        type=Path,
        default=None,
        help=f"Path to DeepTMHMM output (default: {DEFAULT_DEEPTMHMM_ROOT} or skip Step 2 if not found)",
    )
    parser.add_argument(
        "--median-linker-csv",
        type=Path,
        default=DEFAULT_MEDIAN_LINKER_CSV,
        help=f"Path to median_linker_length.csv (default: {DEFAULT_MEDIAN_LINKER_CSV})",
    )
    parser.add_argument(
        "--min-gap-size",
        type=int,
        default=20,
        help="Min gap size for handle_gaps (default: 20)",
    )
    parser.add_argument(
        "--uniprot-timeout",
        type=int,
        default=60,
        help="UniProt API request timeout in seconds (default: 60)",
    )
    parser.add_argument(
        "--uniprot-max-retries",
        type=int,
        default=5,
        help="Max retries for UniProt fetch (default: 5)",
    )
    args = parser.parse_args()

    # Resolve deeptmhmm root
    deeptmhmm_root = args.deeptmhmm_root
    if deeptmhmm_root is None and Path(DEFAULT_DEEPTMHMM_ROOT).exists():
        deeptmhmm_root = Path(DEFAULT_DEEPTMHMM_ROOT)
    elif deeptmhmm_root is not None and not deeptmhmm_root.exists():
        logger.warning("DeepTMHMM root %s does not exist, skipping Step 2", deeptmhmm_root)
        deeptmhmm_root = None

    if not args.input.exists():
        logger.error("Input file not found: %s", args.input)
        sys.exit(1)

    if not args.median_linker_csv.exists():
        logger.error("Median linker CSV not found: %s", args.median_linker_csv)
        sys.exit(1)

    df = pd.read_csv(args.input)
    required_cols = {"id", "uniprot_accn", "Range", "pdb_seq"}
    missing = required_cols - set(df.columns)
    if missing:
        logger.error("Input CSV missing required columns: %s", missing)
        sys.exit(1)

    df_median_linker_length = pd.read_csv(args.median_linker_csv)

    uniprot_cache: Dict[str, Optional[Dict]] = {}
    pdb_structure_cache: Dict[str, Any] = {}

    rows_data = []
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing"):
        row_result = process_row(
            row,
            pdb_dir=args.pdb_dir,
            df_median_linker_length=df_median_linker_length,
            uniprot_cache=uniprot_cache,
            pdb_structure_cache=pdb_structure_cache,
            deeptmhmm_root=deeptmhmm_root,
            min_gap_size=args.min_gap_size,
            uniprot_timeout=args.uniprot_timeout,
            uniprot_max_retries=args.uniprot_max_retries,
        )
        rows_data.append(row_result)

    result_df = pd.DataFrame(rows_data)
    for col in result_df.columns:
        df[col] = result_df[col].values

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    logger.info("Wrote %d rows to %s", len(df), args.output)


if __name__ == "__main__":
    main()
