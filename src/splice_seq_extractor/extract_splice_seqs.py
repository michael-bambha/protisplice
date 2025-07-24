"""
File: extract_splice_seqs.py
Author: Michael Bambha
Contact: bambha.m@northeastern.edu
Description: Main SpliceSeqExtractor class for splice sequence extraction.
Users can define a window of bases into the intron and exon regions of an
identified splice junction and extract out the sequence. Additionally,
methods for randomly sampling introns of transcripts are provided.
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
        """Lazy loaded transcripts

        Returns:
            Dict[str, Transcript]: Dict of ID: transcript
        """
        if self._transcripts is None:
            self.logger.info("Parsing transcripts from GFF file...")
            self._transcripts = self.gff_parser.parse_transcripts()
            self.logger.info("Found %d transcripts", len(self._transcripts))
        return self._transcripts

    @property
    def junctions(self) -> List[SpliceJunction]:
        """Lazy loaded junctions

        Returns:
            List[SpliceJunction]: List of {id: transcriptid_junctype_num, seqid: chr#,
            coord: 1-based coord of first exon base, strand: + or -, junc_type: donor or acceptor}
        """
        if self._junctions is None:
            self.logger.info("Extracting splice junctions...")
            self._junctions = self.junction_extractor.identify_splice_junctions(
                self.transcripts
            )
            self.logger.info("Found %d splice junctions", len(self._junctions))
        return self._junctions

    def extract_positive_sequences(self) -> List[JunctionData]:
        """Finds the sequences for identified true splice sites.

        Returns:
            List[JunctionData]: List of {junction: SpliceJunction, win_start: win_start,
            win_end: win_end, sequence: sequence}
        """
        self.logger.info("Extracting positive sequences...")
        return self.sequence_extractor.extract_splice_sites(self.junctions)

    def extract_negative_sequences(
        self, target_count: Optional[int] = None, seed: int = 100
    ) -> List[JunctionData]:
        """Finds the sequence for randomly determined start position of an intron.

        Args:
            target_count (Optional[int], optional): Max number of sequences to sample.
            Defaults to None.
            seed (int): Random seed. Defaults to 100.

        Returns:
            List[JunctionData]: List of {junction: SpliceJunction, win_start: win_start,
            win_end: win_end, sequence: sequence}
        """
        if target_count is None:
            target_count = len(self.junctions)

        self.logger.info("Extracting negative sequences...")
        return self.sequence_extractor.sample_introns(
            self.transcripts, target_count, seed
        )

    def extract_splice_sites(self) -> List[JunctionData]:
        """Find sequences around true splice sites.

        Returns:
            List[JunctionData]:  List of JunctionData objects with splice site sequences
        """
        return self.sequence_extractor.extract_splice_sites(self.junctions)

    def sample_introns(
        self, target_count: int, random_seed: Optional[int] = None
    ) -> List[JunctionData]:
        """Sample introns of extracted transcripts. Finds a

        Args:
            target_count (int): Maximum number of introns to sample
            random_seed (Optional[int], optional): _description_. Defaults to None.

        Returns:
            List[JunctionData]: _description_
        """
        return self.sequence_extractor.sample_introns(
            self.transcripts, target_count, random_seed
        )

    def write_sequences_to_fasta(
        self,
        results: ExtractionResults,
        positive_output: Optional[str] = None,
        negative_output: Optional[str] = None,
    ) -> Optional[Tuple]:
        """Write sequences to FASTA files"""
        return FastaWriter.write_results(results, positive_output, negative_output)

    def get_info(self) -> Dict[str, int]:
        """Returns parameters and info about the extractor. Includes exon bases,
        intron bases, window size, buffer size, transcript count, and splice junction count.

        Returns:
            Dict[str, int]: Dict of exon bases,
        intron bases, window size, buffer size, transcript count, and splice junction count.
        """
        return {
            "window_size": self.params.window_size,
            "exon_bases": self.params.n_exon,
            "intron_bases": self.params.n_intron,
            "buffer_size": self.params.buffer_size,
            "transcript_count": len(self.transcripts),
            "junction_count": len(self.junctions),
        }
