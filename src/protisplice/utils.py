"""
File: utils.py
Description: Utility function module
"""
from Bio.Seq import Seq
from protisplice.data_models import StrandType


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
