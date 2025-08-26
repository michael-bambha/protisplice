"""
File: expression.py
Author: Michael Bambha
Contact: bambha.m@northeastern.edu
Description: Functionality of parsing and filtering transcripts from expression data.
"""

import statistics
from typing import Dict


class ExpressionFilter:
    """
    Class for filtering gene expression based on a threshold
    """

    def __init__(self, data: Dict[str, float], threshold: float):
        self.data = data
        self.threshold = threshold

    def passes_threshold(self, transcript_id: str) -> bool:
        """Determine whether an expression value is >= threshold

        Args:
            transcript_id (str): ID of transcript to check

        Returns:
            bool: T/F if passes threshold or not
        """
        expression_val = self._find_expression_value(transcript_id)
        return expression_val >= self.threshold

    @staticmethod
    def normalize_transcript_id(transcript_id: str) -> str:
        """Keep transcript IDs the same across different formats.
        Note that GENCODE transcripts often have version numbers trailing the ID.
        e.g. TRANS00001.3 -- Ensembl annotations do not have these.

        This will also remove leading "transcript:".

        Args:
            transcript_id (str): ID of the transcript

        Returns:
            str: Transcript ID without leading strings or trailing version numbers.
        """
        if "." in transcript_id:
            transcript_id = transcript_id.split(".")[0]
        if transcript_id.startswith("transcript:"):
            transcript_id = transcript_id.replace("transcript:", "")
        return transcript_id

    def _find_expression_value(self, transcript_id: str) -> float:
        """Returns the expression for a transcript ID

        Args:
            transcript_id (str): ID of the transcript

        Returns:
            float: expression value for that transcript
        """
        normalized_id = self.normalize_transcript_id(transcript_id)
        if normalized_id in self.data:
            return self.data[normalized_id]
        return 0.0


class ExpressionParser:
    """
    Class for parsing various gene expression formats
    """

    @staticmethod
    def parse_file(file_path: str, format_type: str) -> Dict[str, float]:
        """Organizes parsing of various expression file formats

        Args:
            file_path (str): Path of expression file
            format_type (str): Format of expression file (Currently supported: kallisto, gtex)

        Raises:
            ValueError: If unsupported file format is supplied

        Returns:
            Dict[str, float]:  Dict of {ID: expr_val} for all IDs
        """
        if format_type == "kallisto":
            return ExpressionParser._parse_kallisto(file_path)
        if format_type == "gtex":
            return ExpressionParser._parse_gtex(file_path)
        raise ValueError(f"Unsupported format type for {format_type}.")

    @staticmethod
    def _parse_kallisto(file_path: str) -> Dict[str, float]:
        """Parses Kallisto's abundance.tsv format

        Args:
            file_path (str): Path to kallisto abundance file

        Returns:
            Dict[str, float]: Dict of {ID: expr_val} for all IDs
        """
        expression_data = {}
        with open(file_path, "r", encoding="utf-8") as f:
            header = next(f).strip().split("\t")
            tpm_idx = header.index("tpm")
            target_idx = header.index("target_id")

            for line in f:
                fields = line.strip().split("\t")
                transcript_id = fields[target_idx]
                tpm = float(fields[tpm_idx])
                expression_data[transcript_id] = tpm
        return expression_data

    @staticmethod
    def _parse_gtex(file_path: str, min_expression: float = 0.0) -> Dict[str, float]:
        """Parse GTEx format with streaming

        Args:
            file_path (str): Path to GTEx expression file
            min_expression (float, optional): Minimum expression to include. Defaults to 0.0.

        Returns:
            Dict[str, float]: Dict of {ID: expr_val} for all IDs
        """
        gene_expressions = {}

        with open(
            file_path, "r", encoding="utf-8"
        ) as f:  # 3 header lines in gtex format
            next(f)
            next(f)
            next(f)

            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 4:
                    continue

                gene_id = fields[1]
                expr_values = [
                    float(val_str)
                    for val_str in fields[3:]
                    if val_str and float(val_str) >= min_expression
                ]

                if expr_values:
                    avg_expr = statistics.mean(expr_values)
                    gene_expressions[gene_id] = avg_expr

        return gene_expressions
