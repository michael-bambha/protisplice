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

    @staticmethod
    def write_results(
        results: ExtractionResults,
        positive_output: Optional[str] = None,
        negative_output: Optional[str] = None,
    ) -> tuple:
        """Combination method for writing both true and false splice sites out to two
        separate files

        Args:
            results (ExtractionResults): ExtractionResults object, which contains two lists
            of JunctionData objects for both true and false splice sites.
            positive_output (Optional[str], optional): Output path for true splice sites.
            Defaults to None.
            negative_output (Optional[str], optional): Output path for false splice sites.
            Defaults to None.

        Returns:
            tuple: Tuple of (pos_count, neg_count) representing the number of sequences for true
            and false splice sites, respectively, written out to the files.
        """
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
