"""
Reusable gap-handling logic for DPAM-AI generated domains.

Handles gaps in protein sequences by either filling with UniProt sequence
(for short gaps) or grafting G/S/D linkers (for long gaps). Linker length
is determined by CA-CA distance lookup in median_linker_length.csv.
"""

import os
import re
import random
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1


# -----------------------------------------------------------------------------
# From random_loop_seq.ipynb
# -----------------------------------------------------------------------------


def generate_loop_sequence(
    length: int,
    seed: Optional[int] = None,
    s_percent: Tuple[float, float] = (0.2, 0.4),
    ds_ratio: Tuple[float, float] = (0.1, 0.3),
    max_consecutive_g: int = 4,
    max_consecutive_s: int = 2,
) -> str:
    """
    Generate a flexible loop sequence satisfying specific criteria.

    Args:
        length: Length of the sequence (must be >= 1)
        seed: Random seed for reproducibility. If None, generation is stochastic.
        s_percent: Tuple (min, max) for percent of S residues (e.g., (0.2, 0.4)).
        ds_ratio: Tuple (min, max) for D/S mutation ratio (e.g., (0.1, 0.3)).
        max_consecutive_g: Maximum allowed consecutive G residues.
        max_consecutive_s: Maximum allowed consecutive S residues.

    Returns:
        A protein sequence string consisting of G, S, and D residues.
    """
    if length < 1:
        raise ValueError("Length must be at least 1")

    rng = random.Random(seed)

    if not (isinstance(s_percent, (list, tuple)) and len(s_percent) == 2):
        raise ValueError("s_percent must be a list or tuple of length 2")
    s_percent_min, s_percent_max = s_percent
    s_percentage = rng.uniform(s_percent_min, s_percent_max)
    target_s_count = round(length * s_percentage)

    sequence = []
    for i in range(length):
        if i == 0 or i == length - 1:
            sequence.append("G")
        else:
            consecutive_g_before = 0
            for j in range(i - 1, -1, -1):
                if sequence[j] == "G":
                    consecutive_g_before += 1
                else:
                    break

            consecutive_s_before = 0
            for j in range(i - 1, -1, -1):
                if sequence[j] == "S":
                    consecutive_s_before += 1
                else:
                    break

            must_place_s = consecutive_g_before >= max_consecutive_g
            can_place_s = consecutive_s_before < max_consecutive_s
            current_s_count = sequence.count("S")
            remaining_positions = length - i - 1
            s_still_needed = target_s_count - current_s_count

            if must_place_s and can_place_s:
                sequence.append("S")
            elif must_place_s and not can_place_s:
                sequence.append("G")
            elif s_still_needed > 0 and can_place_s:
                s_probability = (
                    s_still_needed / (remaining_positions + 1)
                    if remaining_positions >= 0
                    else 0
                )
                if rng.random() < s_probability:
                    sequence.append("S")
                else:
                    sequence.append("G")
            else:
                sequence.append("G")

    s_positions = [i for i, aa in enumerate(sequence) if aa == "S"]
    total_s = len(s_positions)

    if total_s > 0:
        if not (isinstance(ds_ratio, (list, tuple)) and len(ds_ratio) == 2):
            raise ValueError("ds_ratio must be a list or tuple of length 2")
        d_prob_min, d_prob_max = ds_ratio
        target_d_s_ratio = rng.uniform(d_prob_min, d_prob_max)
        num_d = round(total_s * target_d_s_ratio / (1 + target_d_s_ratio))
        num_d = max(0, min(num_d, total_s))
        if num_d > 0:
            positions_to_mutate = rng.sample(s_positions, num_d)
            for pos in positions_to_mutate:
                sequence[pos] = "D"

    return "".join(sequence)


def check_loop_sequence(sequence: str) -> Tuple[float, float, int, int]:
    """
    Compute the S content, D/S ratio, and max consecutive G and S residues.
    """
    g_count = sequence.count("G")
    s_count = sequence.count("S")
    d_count = sequence.count("D")

    consecutive_g = 0
    max_consecutive_g = 0
    for char in sequence:
        if char == "G":
            consecutive_g += 1
            max_consecutive_g = max(max_consecutive_g, consecutive_g)
        else:
            consecutive_g = 0

    consecutive_s = 0
    max_consecutive_s = 0
    for char in sequence:
        if char == "S":
            consecutive_s += 1
            max_consecutive_s = max(max_consecutive_s, consecutive_s)
        else:
            consecutive_s = 0

    s_content = (s_count + d_count) / len(sequence)
    d_s_ratio = d_count / s_count if s_count > 0 else 0
    return s_content, d_s_ratio, max_consecutive_g, max_consecutive_s


# -----------------------------------------------------------------------------
# From domains_to_library.ipynb - Range-based utilities
# -----------------------------------------------------------------------------


def num_gaps(range_str: str) -> int:
    """Return the number of gaps (count of commas) in a range string."""
    return range_str.count(",")


def domain_length_gap(range_str: str) -> int:
    """Calculate total length by summing fragment lengths."""
    ranges = [tuple(map(int, frag.split("-"))) for frag in range_str.split(",")]
    return sum(end - start + 1 for start, end in ranges)


def domain_length_continuous(range_str: str) -> int:
    """Calculate continuous length from min start to max end."""
    fragments = range_str.split(",")
    starts = []
    ends = []
    for frag in fragments:
        frag = frag.strip()
        if "-" in frag:
            start_str, end_str = frag.split("-")
            starts.append(int(start_str))
            ends.append(int(end_str))
        else:
            starts.append(int(frag))
            ends.append(int(frag))
    return max(ends) - min(starts) + 1


class UnlinkableGapError(Exception):
    """Raised when a gap cannot be bridged (CA distance exceeds table maximum)."""

    def __init__(self, gap_indices, distances, max_linkable_dist):
        self.gap_indices = gap_indices
        self.distances = distances
        self.max_linkable_dist = max_linkable_dist
        message = (
            f"Cannot link gaps at indices {gap_indices} with distances "
            f"{[f'{d:.2f}' for d in distances]} Å. "
            f"Maximum linkable distance is {max_linkable_dist:.2f} Å."
        )
        super().__init__(message)


def get_CA_coord(struct, resi: int, chain=None):
    """Get CA coordinate of a residue in a structure."""
    if chain is None:
        chain = list(struct.get_chains())[0]
    for residue in chain:
        if residue.id[1] == resi:
            return residue["CA"].get_coord()
    return None


def _lookup_linker_length(
    ca_distance: float, df_median_linker_length: pd.DataFrame
) -> Optional[int]:
    """
    Look up linker length (seq_dist) from median_linker_length table based on CA-CA distance.
    Returns None if distance exceeds table maximum.
    """
    max_calpha_dist = df_median_linker_length["calpha_dist"].max()
    if ca_distance > max_calpha_dist:
        return None
    filtered_df = df_median_linker_length[
        df_median_linker_length["calpha_dist"] > ca_distance
    ]
    if filtered_df.empty:
        return None
    return int(filtered_df.iloc[0]["seq_dist"])


def get_protein_sequence(struct, resi_start: int, resi_end: int, chain_id=None) -> str:
    """Extract protein sequence from structure between residue indices."""
    model = struct[0]
    chains = list(model.get_chains())
    chain = chains[0] if chain_id is None else model[chain_id]
    residues_in_range = [
        res for res in chain if res.id[1] >= resi_start and res.id[1] <= resi_end
    ]
    seq = ""
    for res in residues_in_range:
        if "CA" in res and res.resname in protein_letters_3to1:
            seq += protein_letters_3to1[res.resname]
    return seq


# -----------------------------------------------------------------------------
# pdb_seq workflow for database_summary
# -----------------------------------------------------------------------------


def _parse_range_fragment(frag: str) -> Tuple[int, int]:
    """
    Parse a range fragment like '65-99' or '-205--176' into (start, end).
    Handles both positive and negative residue indices.
    """
    match = re.match(r"^(-?\d+)-(-?\d+)$", frag.strip())
    if not match:
        raise ValueError(f"Invalid range fragment: {frag!r}")
    return int(match.group(1)), int(match.group(2))


def parse_pdb_seq_gaps(pdb_seq: str) -> Tuple[list, list]:
    """
    Parse pdb_seq (domain sequence with X caps and - gaps) into segments and gap lengths.

    Format: X(seg1)X(gap_dashes)X(seg2)X(gap_dashes)X(seg3)X - X marks segment boundaries.
    Strips outer X and boundary X between segments and gaps.

    Args:
        pdb_seq: Sequence with X for ACE/NME caps and - for gaps (e.g. XABC---DEFX)

    Returns:
        (segments, gap_lengths)
        - segments: list of sequence strings between gaps (X boundary markers removed)
        - gap_lengths: list of int, length of each gap in residues (count of dashes)
    """
    stripped = pdb_seq.strip()
    if stripped.startswith("X"):
        stripped = stripped[1:]
    if stripped.endswith("X"):
        stripped = stripped[:-1]

    segments = []
    gap_lengths = []
    i = 0
    while i < len(stripped):
        if stripped[i] == "-":
            gap_count = 0
            while i < len(stripped) and stripped[i] == "-":
                gap_count += 1
                i += 1
            gap_lengths.append(gap_count)
        else:
            seg_start = i
            while i < len(stripped) and stripped[i] != "-":
                i += 1
            seg = stripped[seg_start:i]
            if seg.endswith("X"):
                seg = seg[:-1]
            if seg.startswith("X"):
                seg = seg[1:]
            if seg:
                segments.append(seg)

    return segments, gap_lengths


def handle_gaps_pdb_seq(
    pdb_seq: str,
    full_sequence: str,
    range_str: str,
    domain_id: str,
    domain_pdb_path: str,
    df_median_linker_length: pd.DataFrame,
    min_gap_size: int = 20,
    seed: Optional[int] = None,
    structure=None,
) -> dict:
    """
    Handle gaps in pdb_seq and return domain_seq with grafted linkers.

    Linker length is always determined by CA-based distance from domain PDB
    and lookup in median_linker_length.csv.

    Args:
        pdb_seq: Domain sequence with X caps and - gaps
        full_sequence: Full UniProt sequence for this entry
        range_str: Comma-separated ranges (e.g. "1-90,141-180,206-220")
        domain_id: Domain identifier (for error messages)
        domain_pdb_path: Path to domain PDB file (required for CA distance)
        df_median_linker_length: Linker reference DataFrame
        min_gap_size: Gaps smaller than this use uniprot sequence
        seed: Random seed for linker generation (generate_loop_sequence)
        structure: Optional pre-loaded Bio.PDB structure (avoids reload when cached)

    Returns:
        dict with keys: full_sequence, num_gaps, linker_types, linker_seqs, domain_seq
    """
    segments, gap_lengths = parse_pdb_seq_gaps(pdb_seq)
    fragments = [f.strip() for f in range_str.split(",")]

    if len(segments) != len(fragments):
        raise ValueError(
            f"Segment count ({len(segments)}) does not match fragment count ({len(fragments)}) "
            f"for domain {domain_id}"
        )

    if len(gap_lengths) == 0:
        domain_seq = "".join(segments)
        return {
            "full_sequence": full_sequence,
            "num_gaps": 0,
            "linker_types": None,
            "linker_seqs": None,
            "domain_seq": domain_seq,
        }

    if structure is not None:
        struct = structure
    else:
        parser = PDBParser(QUIET=True)
        if not os.path.exists(domain_pdb_path):
            raise FileNotFoundError(f"Domain PDB not found: {domain_pdb_path}")
        struct = parser.get_structure("protein", domain_pdb_path)

    linker_types_list = []
    linker_seqs_list = []
    domain_parts = []

    for i in range(len(segments)):
        domain_parts.append(segments[i])

        if i < len(gap_lengths):
            frag_i = fragments[i]
            frag_next = fragments[i + 1]
            frag_i_start, frag_i_end = _parse_range_fragment(frag_i)
            frag_next_start, frag_next_end = _parse_range_fragment(frag_next)
            gap_length = frag_next_start - frag_i_end - 1

            if gap_length < min_gap_size:
                gap_resi_start = frag_i_end + 1
                gap_resi_end = frag_next_start - 1
                if gap_resi_start <= gap_resi_end:
                    linker_seq = full_sequence[gap_resi_start - 1 : gap_resi_end]
                else:
                    linker_seq = ""
                linker_type = "uniprot"
            else:
                ca1 = get_CA_coord(struct, frag_i_end)
                ca2 = get_CA_coord(struct, frag_next_start)
                if ca1 is None or ca2 is None:
                    raise ValueError(
                        f"Cannot get CA coordinates for residues {frag_i_end} and "
                        f"{frag_next_start} in domain {domain_id}"
                    )
                ca_distance = float(np.linalg.norm(np.array(ca1) - np.array(ca2)))
                linker_length = _lookup_linker_length(ca_distance, df_median_linker_length)
                if linker_length is None:
                    linker_length = int(df_median_linker_length["seq_dist"].max())
                linker_seq = generate_loop_sequence(length=linker_length, seed=seed)

                if gap_length > 0 and len(linker_seq) / gap_length > 0.5:
                    gap_resi_start = frag_i_end + 1
                    gap_resi_end = frag_next_start - 1
                    if gap_resi_start <= gap_resi_end:
                        linker_seq = full_sequence[gap_resi_start - 1 : gap_resi_end]
                    else:
                        linker_seq = ""
                    linker_type = "uniprot"
                else:
                    linker_type = "GS"

            linker_types_list.append(linker_type)
            linker_seqs_list.append(linker_seq)
            domain_parts.append(linker_seq)

    domain_seq = "".join(domain_parts)
    return {
        "full_sequence": full_sequence,
        "num_gaps": len(gap_lengths),
        "linker_types": ",".join(linker_types_list),
        "linker_seqs": ",".join(linker_seqs_list),
        "domain_seq": domain_seq,
    }


