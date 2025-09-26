"""
File: feature_eng.py
Description: Module for feature engineering and data augmentation.
Includes functions for injecting or removing canonical motifs from
splice sites to create stronger decoys.
"""

from typing import Optional
from copy import deepcopy
from .data_models import JunctionType, JunctionData, ExtractionParams
from .utils import pick_noncanon_dinuc, get_consensus_index, motif_for


def inject_consensus(
    junc: JunctionData,
    extraction_params: ExtractionParams,
    junc_type: JunctionType,
    idx: Optional[int] = None,
) -> JunctionData:
    """Replace the center indices of a sequence window with the donor or acceptor
    consensus sequence. Can be used for creating higher quality / harder decoy
    sequences.

    Args:
        junc_type (str): Donor or acceptor
        extraction_params (ExtractionParams): parameters that were passed in to the extractor.
        n_exon of extraction_params will be used to calculate determine the proper
        location to inject the sequence.
        idx (Optional[int]): If you prefer to supply an index to inject the sequence, this
        value will be used instead of the position calculated from ExtractionParams.n_exon.
    Returns:
        JunctionData: New JunctionData with modified sequence with donor or acceptor consensus
    """
    new_junc = deepcopy(junc)  # copy to not modify object in-place
    seq = new_junc.sequence
    if idx is None:
        idx = get_consensus_index(extraction_params, junc, junc_type)
    if idx < 0 or idx + 1 >= len(seq):
        raise IndexError(
            f"Consensus start idx {idx} out of bounds for sequence of length {len(seq)}"
        )
    consensus = motif_for(junc_type)
    new_seq = seq[:idx] + consensus + seq[idx + 2:]
    new_junc.sequence = new_seq
    return new_junc


def remove_consensus(
    junc: JunctionData,
    extraction_params: Optional[ExtractionParams] = None,
    idx: Optional[int] = None,
    *,
    only_if_present: bool = True,
    replacement: Optional[str] = None,
) -> JunctionData:
    """Destroy/remove the canonical dinucleotide at the junction to form a negative.
      - If only_if_present=True, will no-op unless the canonical motif is present.
      - If replacement is provided, it must be a 2-mer != canonical motif.
      - Otherwise picks a random non-canonical 2-mer.

    Args:
        junc (JunctionData): JunctionData object containing the sequence to inject
        idx (Optional[int], optional): Starting index of the consensus sequence
        within the string. Defaults to None.
        only_if_present (bool, optional): Only will remove the dinucleotide
        if a canonical motif is found. Defaults to True.
        replacement (Optional[str], optional): Dinucleotide to replace
        the consensus. Defaults to None.
        seed (Optional[int]): Random state. Defaults to 100.

    Raises:
        IndexError: If the index is less than zero or greater than the length
        of the sequence
        ValueError: If both extraction_params and idx are None
        ValueError: If the replacement is not a string of length 2, or
        if the replacement is the same as the canonical motif.

    Returns:
        JunctionData: New JunctionData with consensus motif mutated
    """
    new_junc = deepcopy(junc)
    seq = new_junc.sequence
    junc_type = junc.junction.junction_type
    mutant = "NN"

    if idx is None:
        if extraction_params is None:
            raise ValueError("extraction params required when idx is None")
        idx = get_consensus_index(extraction_params, junc, junc_type)
    if idx < 0 or idx + 1 >= len(seq):
        raise IndexError(
            f"Consensus start idx {idx} out of bounds for sequence of length {len(seq)}"
        )

    canonical = motif_for(junc_type)
    current = seq[idx: idx+2]

    if only_if_present and current != canonical:
        return junc  # leave as-is if not canonical here

    if replacement is not None:
        if len(replacement) != 2 or replacement == canonical:
            raise ValueError(
                "replacement must be a 2-mer different from the canonical motif"
            )
        mutant = replacement
    else:
        mutant = pick_noncanon_dinuc(canonical)

    new_seq = (
        seq[:idx] + mutant + seq[idx + 2:]
    )  # replace dinuc with mutant at correct pos
    new_junc.sequence = new_seq
    return new_junc
