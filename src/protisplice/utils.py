"""
File: utils.py
Description: Utility function module
"""

import random
from typing import Optional
from Bio.Seq import Seq
from protisplice.data_models import (
    StrandType,
    JunctionData,
    JunctionType,
    ExtractionParams,
)
from protisplice.motif_scoring import calculate_shannon_entropy


def apply_strand(seq: str, strand: StrandType) -> str:
    """Handle reverse complementing when strand type is (-)

    Args:
        seq (str): Sequence to RC if strand is (-)
        strand (StrandType): StrandType of the sequence (+) or (-)

    Returns:
        str: Original sequence if strand is (+) or reverse complement
        if strand is (-)
    """
    return str(Seq(seq).reverse_complement()) if strand == StrandType.NEGATIVE else seq


def pick_noncanon_dinuc(canonical: str) -> str:
    """Pick a random dinucleotide != the
    passed in canonical motif

    Args:
        canonical (str): Canonical dinucleotide, either "AG" or "GT"
        seed (int): Random state.
    Returns:
        str: Any random dinculeotide combination not equivalent
        to what was passed in
    """
    bases = "ACGT"
    while True:
        mutant = random.choice(bases) + random.choice(bases)
        if mutant != canonical:
            return mutant


def get_consensus_index(
    extraction_params: ExtractionParams,
    junc: Optional[JunctionData] = None,
    junc_type: Optional[JunctionType] = None,
) -> int:
    """Helper function for getting correct coordinates of the
    motif based on junction.coord. If no junction type is passed in,
    the function will look for the junction type of the passed in JunctionData object.


    Args:
        junc (JunctionData): JunctionData object
        extraction_params (ExtractionParams): Parameters passed into the extractor.
        ExtractionParams.n_exon will be used to calculate the proper location of
        the consensus motif depending whether it is acceptor or donor.
        junc_type (Optional[JunctionType]): Junction type of the sequence to
        find the consensus motif index. If no argument is supplied, the
        junction type will be obtained from junc.junction.junction_type.

    Raises:
        ValueError: If the junction type is not donor or acceptor.

    Returns:
        int: Start coordinate of the AG or GT dinucleotide
    """
    if junc is None and junc_type is None:
        raise ValueError(
            "One of junc or junc_type must be passed in,"
            "as the motif consensus is ambiguous."
        )
    exon_idx = extraction_params.n_exon
    if not junc_type:
        junc_type = junc.junction.junction_type
    if junc_type == JunctionType.DONOR:
        # Donor site = GT starting at exon_idx (first 2 bases of intron)
        return exon_idx
    if junc_type == JunctionType.ACCEPTOR:
        # Acceptor site = AG ending right before exon_idx
        return exon_idx - 2

    raise ValueError(
        f"Junction type must be either donor or acceptor, not {junc_type}."
    )


def motif_for(junc_type: JunctionType) -> str:
    """Return GT for donors, AG for acceptors

    Args:
        junc_type (JunctionType): Junction type (Donor or Acceptor)

    Returns:
        str: "GT" for donor, "AG" for acceptor
    """
    if junc_type not in (JunctionType.ACCEPTOR, JunctionType.DONOR):
        raise ValueError(
            f"Junction type must be either donor or acceptor, not {junc_type}."
        )
    return "GT" if junc_type == JunctionType.DONOR else "AG"

