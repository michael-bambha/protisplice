"""
File: feature_eng.py
Description: Module for feature engineering and data augmentation.
Includes functions for injecting or removing canonical motifs from
splice sites to create stronger decoys.
"""

from typing import Optional
import random
from protisplice import JunctionType, JunctionData


def inject_consensus(
    junc: JunctionData,
    junc_type=JunctionType,
    idx: Optional[int] = None,
) -> str:
    """Replace the center indices of a sequence window with the donor or acceptor
    consensus sequence. Can be used for creating higher quality / harder decoy
    sequences.

    Args:
        seq (str): Sequence of nucleotides
        junc_type (str): Donor or acceptor
        center_idx (int): Index of the junction where you would expect
        the consensus sequence to be located in a true splice site.
        For example...

    Returns:
        str: Modified sequence with donor or acceptor consensus
    """
    seq = junc.sequence
    if idx is None:
        idx = _get_consensus_index(junc, junc_type)
    if idx < 0 or idx + 1 >= len(seq):
        raise IndexError(
            f"Consensus start idx {idx} out of bounds for sequence of length {len(seq)}"
        )
    consensus = _motif_for(junc_type)

    return seq[:idx] + consensus + seq[idx + 2:]


def remove_consensus(
    junc: JunctionData,
    junc_type: JunctionType,
    idx: Optional[int] = None,
    *,
    only_if_present: bool = True,
    replacement: Optional[str] = None,
) -> str:
    """
    Destroy/remove the canonical dinucleotide at the junction to form a negative.
      - If only_if_present=True, will no-op unless the canonical motif is present.
      - If replacement is provided, it must be a 2-mer != canonical motif.
      - Otherwise picks a random non-canonical 2-mer.
    """
    seq = junc.sequence
    mutant = "NN"
    if idx is None:
        idx = _get_consensus_index(junc, junc_type)
    if idx < 0 or idx + 1 >= len(seq):
        raise IndexError(
            f"Consensus start idx {idx} out of bounds for sequence of length {len(seq)}"
        )

    canonical = _motif_for(junc_type)
    current = seq[idx: idx + 2]

    if only_if_present and current != canonical:
        return seq  # leave as-is if not canonical here

    if replacement is not None:
        if len(replacement) != 2 or replacement == canonical:
            raise ValueError(
                "replacement must be a 2-mer different from the canonical motif"
            )
        mutant = replacement
    else:
        _pick_noncanon_dinuc(canonical)

    return seq[:idx] + mutant + seq[idx + 2:]


def _pick_noncanon_dinuc(canonical: str):
    bases = "ACGT"
    while True:
        mutant = random.choice(bases) + random.choice(bases)
        if mutant != canonical:
            return mutant


def _get_consensus_index(
    junc: JunctionData, junc_type: Optional[JunctionType] = None
) -> int:
    """Helper function for getting correct coordinates of the
    motif based on junction.coord. If no junction type is passed in,
    the function will look for the junction type of the passed in JunctionData object.


    Args:
        junc (JunctionData): JunctionData object

    Returns:
        int: Start coordinate of the AG or GT dinucleotide
    """
    exon_idx = junc.junction.coord - junc.window_start
    if not junc_type:
        junc_type = junc.junction.junction_type
    if junc.junction.junction_type == JunctionType.DONOR:
        # Donor site = GT starting at exon_idx (first 2 bases of intron)
        return exon_idx
    if junc.junction.junction_type == JunctionType.ACCEPTOR:
        # Acceptor site = AG ending right before exon_idx
        return exon_idx - 2

    raise ValueError(
        f"Junction type must be either donor or acceptor, not {junc_type}."
    )


def _motif_for(junc_type: JunctionType) -> str:
    return "GT" if junc_type == JunctionType.DONOR else "AG"
