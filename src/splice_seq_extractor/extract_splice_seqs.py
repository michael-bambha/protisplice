"""
Main SpliceSeqExtractor class for splice sequence extraction.

Author: Michael Bambha
Contact: bambha.m@northeastern.edu
Description: A Python class for obtaining true positive and false positive
sequences around splice sites to be used for downstream model training.
"""

import logging
from typing import Dict, List, Optional, Tuple

from .data_models import (
    ExtractionParams,
    ExtractionResults,
    TranscriptFilter,
    SpliceJunction,
    JunctionData,
    Transcript,
)
from .expression import ExpressionParser, ExpressionFilter
from .gff_parser import GFFParser
from .sequence_extraction import SequenceExtractor
from .splice_junction_extractor import SpliceJunctionExtractor
from .fasta_writer import FastaWriter


class SpliceSeqExtractor:
    """Main class for splice sequence extraction"""

    def __init__(
        self,
        gff_path: str,
        fasta_path: str,
        params: Optional[ExtractionParams] = None,
        transcript_filter: TranscriptFilter = TranscriptFilter.PROTEIN_CODING,
        expression_file: Optional[str] = None,
        min_expression: float = 1.0,
        expression_format: str = "kallisto",
    ):
        """
        Initialize the SpliceSeqExtractor

        Args:
            gff_path: Path to GFF3 file
            fasta_path: Path to FASTA file (must be indexed)
            params: Extraction parameters
            transcript_filter: Filter for transcript types
            expression_file: Optional path to expression file
            min_expression: Minimum expression threshold
            expression_format: Format of expression file
        """
        self.gff_path = gff_path
        self.fasta_path = fasta_path
        self.params = params or ExtractionParams()
        self.transcript_filter = transcript_filter

        logging.basicConfig(
            level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
        )
        self.logger = logging.getLogger(__name__)

        self.expression_filter = None
        if expression_file:
            expression_data = ExpressionParser.parse_file(
                expression_file, expression_format
            )
            self.expression_filter = ExpressionFilter(expression_data, min_expression)
            self.logger.info(
                "Loaded expression data for %d transcripts", len(expression_data)
            )
        self.gff_parser = GFFParser(gff_path, transcript_filter, self.expression_filter)
        self.sequence_extractor = SequenceExtractor(fasta_path, self.params)
        self.junction_extractor = SpliceJunctionExtractor()

        # cache
        self._transcripts = None
        self._junctions = None

    @property
    def transcripts(self) -> Dict[str, Transcript]:
        """Get transcripts (lazy loaded)"""
        if self._transcripts is None:
            self.logger.info("Parsing transcripts from GFF file...")
            self._transcripts = self.gff_parser.parse_transcripts()
            self.logger.info("Found %d transcripts", len(self._transcripts))
        return self._transcripts

    @property
    def junctions(self) -> List[SpliceJunction]:
        """Get splice junctions (lazy loaded)"""
        if self._junctions is None:
            self.logger.info("Extracting splice junctions...")
            self._junctions = self.junction_extractor.identify_splice_junctions(
                self.transcripts
            )
            self.logger.info("Found %d splice junctions", len(self._junctions))
        return self._junctions

    def extract_positive_sequences(self) -> List[JunctionData]:
        """Extract positive splice site sequences"""
        self.logger.info("Extracting positive sequences...")
        return self.sequence_extractor.extract_splice_sites(self.junctions)

    def extract_negative_sequences(
        self, target_count: Optional[int] = None
    ) -> List[JunctionData]:
        """Extract negative sequences from intronic regions"""
        if target_count is None:
            target_count = len(self.junctions)

        self.logger.info("Extracting negative sequences...")
        return self.sequence_extractor.sample_introns(self.transcripts, target_count)

    def extract_all_sequences(
        self, negative_count: Optional[int] = None
    ) -> ExtractionResults:
        """Extract both positive and negative sequences"""
        positive_sequences = self.extract_positive_sequences()

        if negative_count is None:
            negative_count = len(positive_sequences)

        negative_sequences = self.extract_negative_sequences(negative_count)

        results = ExtractionResults(
            positive_sequences=positive_sequences, negative_sequences=negative_sequences
        )

        self.logger.info(
            "Extraction complete -- found %d positive and %d negative sequences",
            results.positive_count,
            results.negative_count,
        )

        return results

    def extract_splice_sites(self) -> ExtractionResults:
        """Extract only positive sequences"""
        positive_sequences = self.extract_positive_sequences()

        results = ExtractionResults(positive_sequences=positive_sequences)

        self.logger.info(
            "Positive extraction complete -- found %d sequences", results.positive_count
        )

        return results

    def sample_introns(self, target_count: int) -> ExtractionResults:
        """Extract only negative sequences"""
        negative_sequences = self.extract_negative_sequences(target_count)

        results = ExtractionResults(negative_sequences=negative_sequences)

        self.logger.info(
            "Negative extraction complete -- found %d sequences", results.negative_count
        )

        return results

    def write_sequences_to_fasta(
        self,
        results: ExtractionResults,
        positive_output: Optional[str] = None,
        negative_output: Optional[str] = None,
    ) -> Optional[Tuple]:
        """Write sequences to FASTA files"""
        return FastaWriter.write_results(results, positive_output, negative_output)

    def extract_and_write(
        self,
        positive_output: Optional[str] = None,
        negative_output: Optional[str] = None,
        negative_count: Optional[int] = None,
    ) -> Optional[Tuple]:
        """Extract sequences and write to FASTA files in one step"""
        results = None
        if positive_output and negative_output:
            results = self.extract_all_sequences(negative_count)
        elif positive_output:
            results = self.extract_splice_sites()
        elif negative_output:
            if negative_count is None:
                negative_count = len(self.junctions)
            results = self.sample_introns(negative_count)
        else:
            raise ValueError("Must provide at least one output file path")

        return self.write_sequences_to_fasta(results, positive_output, negative_output)

    def get_stats(self) -> Dict[str, int]:
        """Get statistics about the extraction parameters and data"""
        return {
            "window_size": self.params.window_size,
            "exon_bases": self.params.n_exon,
            "intron_bases": self.params.n_intron,
            "buffer_size": self.params.buffer_size,
            "transcript_count": len(self.transcripts),
            "junction_count": len(self.junctions),
        }
