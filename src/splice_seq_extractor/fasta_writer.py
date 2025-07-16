"""
File: fasta_writer.py
Description: FASTA writing functionality
"""

from typing import List, Optional
from .data_models import JunctionData, ExtractionResults


class FastaWriter:
    """Handles writing sequence data to FASTA files"""

    @staticmethod
    def write_sequences(sequences: List[JunctionData], output_path: str) -> int:
        """Write sequence data to FASTA file"""
        count = 0
        with open(output_path, "w", encoding="utf-8") as f:
            for data in sequences:
                if data.sequence:
                    f.write(f"{data.to_fasta_header()}\n{data.sequence}\n")
                    count += 1
        return count

    @staticmethod
    def write_results(
        results: ExtractionResults,
        positive_output: Optional[str] = None,
        negative_output: Optional[str] = None,
    ) -> tuple:
        """Write extraction results to files"""
        pos_count = 0
        neg_count = 0

        if results.positive_sequences and positive_output:
            pos_count = FastaWriter.write_sequences(
                results.positive_sequences, positive_output
            )

        if results.negative_sequences and negative_output:
            neg_count = FastaWriter.write_sequences(
                results.negative_sequences, negative_output
            )

        return pos_count, neg_count
