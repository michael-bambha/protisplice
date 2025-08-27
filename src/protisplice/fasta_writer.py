"""
File: fasta_writer.py
Description: FASTA writing functionality
"""

from typing import List
from .data_models import JunctionData


def write_sequences(sequences: List[JunctionData], output_path: str) -> int:
    """Write the identified splice junctions out in FASTA format.

    Args:
        sequences (List[JunctionData]): List of JunctionData objects containing the
        SpliceJunction,
        window start, window end, and sequence.
        output_path (str): Path for the output FASTA

    Returns:
        int: Number of sequences written to the file
    """
    count = 0
    with open(output_path, "w", encoding="utf-8") as f:
        for data in sequences:
            if data.sequence:
                f.write(f"{data.to_fasta_header()}\n{data.sequence}\n")
                count += 1
    return count
