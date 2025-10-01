"""
File: extract_splice_seqs.py
Author: Michael Bambha
Contact: bambha.m@northeastern.edu
Description: Main SpliceSeqExtractor class for splice sequence extraction.
Users can define a window of bases into the intron and exon regions of an
identified splice junction and extract out the sequence. Additionally,
methods for randomly sampling introns, exons, and intergenic regions
are included.
"""

from typing import Dict, List, Optional

from .data_models import (
    ExtractionParams,
    SamplingParams,
    TranscriptFilter,
    SpliceJunction,
    JunctionData,
    Transcript,
    Gene,
)
from .expression import ExpressionParser, ExpressionFilter
from .gff_parser import GFFParser
from .sequence_extraction import SequenceExtractor
from .splice_junction_extractor import SpliceJunctionExtractor
from .regional_sampler import RegionalSampler


class SpliceSeqExtractor:
    """Main class for splice sequence extraction"""

    def __init__(
        self,
        gff_path: str,
        fasta_path: str,
        extraction_params: Optional[ExtractionParams] = None,
        sampling_params: Optional[SamplingParams] = None,
        transcript_filter: TranscriptFilter = TranscriptFilter.ALL,
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
        self.extraction_params = extraction_params or ExtractionParams()
        self.sampling_params = sampling_params or SamplingParams()
        self.transcript_filter = transcript_filter

        self.expression_filter = None
        if expression_file:
            expression_data = ExpressionParser.parse_file(
                expression_file, expression_format
            )
            self.expression_filter = ExpressionFilter(expression_data, min_expression)
        self.gff_parser = GFFParser(gff_path, transcript_filter, self.expression_filter)
        self.sequence_extractor = SequenceExtractor(fasta_path, self.extraction_params)
        self.junction_extractor = SpliceJunctionExtractor()
        self.regional_sampler = RegionalSampler(fasta_path, self.sampling_params)

        # cache
        self._transcripts = None
        self._junctions = None
        self._genes = None
        self._chromosome_lengths = None

    @property
    def transcripts(self) -> Dict[str, Transcript]:
        """Lazy loaded transcripts

        Returns:
            Dict[str, Transcript]: Dict of ID: transcript
        """
        if self._transcripts is None:
            self._transcripts = self.gff_parser.parse_transcripts()
        return self._transcripts

    @property
    def junctions(self) -> List[SpliceJunction]:
        """Lazy loaded junctions

        Returns:
            List[SpliceJunction]: List of {id: transcriptid_junctype_num, seqid: chr#,
            coord: 1-based coord of first exon base, strand: + or -, junc_type: donor or acceptor}
        """
        if self._junctions is None:
            self._junctions = self.junction_extractor.identify_splice_junctions(
                self.transcripts
            )
        return self._junctions

    @property
    def unique_junctions(self) -> List[SpliceJunction]:
        """Obtain only unique junctions by genomic coordinates. Junctions are parsed
        per-transcript, so many duplicates will be present in the original junctions
        property due to alternative splicing.

        Returns:
            List[SpliceJunction]: List of {id: transcriptid_junctype_num, seqid: chr#,
            coord: 1-based coord of first exon base, strand: + or -, junc_type: donor or acceptor}
            for only unique seqid, coord, strand, and junc_type
        """
        return list(set(self.junctions))

    @property
    def genes(self) -> Dict[str, Gene]:
        """Lazy-loaded genes

        Returns:
            Dict[str, Gene]: Dict of ID: gene
        """
        if self._genes is None:
            self._genes = self.gff_parser.parse_genes()
        return self._genes

    @property
    def chromosome_lengths(self) -> Dict[str, int]:
        """Lazy-loaded chromosome lengths

        Returns:
            Dict[str, int]: Dict of chromosome: length
        """
        if self._chromosome_lengths is None:
            self._chromosome_lengths = self.regional_sampler.get_chromosome_lengths()
        return self._chromosome_lengths

    def extract_splice_sites(self) -> List[JunctionData]:
        """Find sequences around true splice sites.

        Returns:
            List[JunctionData]:  List of {junction: SpliceJunction, win_start: win_start,
            win_end: win_end, sequence: sequence}
        """
        return self.sequence_extractor.extract_splice_sites(self.junctions)

    def extract_intronic_regions(
        self, target_count: int, random_seed: Optional[int] = None
    ) -> List[JunctionData]:
        """Sample introns of extracted transcripts. Finds a

        Args:
            target_count (int): Maximum number of introns to sample
            random_seed (Optional[int], optional): Random state. Defaults to None.

        Returns:
            List[JunctionData]: List of {junction: SpliceJunction, win_start: win_start,
            win_end: win_end, sequence: sequence}
        """
        return self.regional_sampler.sample_introns(
            self.transcripts, target_count, random_seed
        )

    def extract_exonic_sequences(
        self, target_count: int, seed: int = 100
    ) -> List[JunctionData]:
        """Sample exonic regions within a defined buffer region.

        Args:
            target_count (int): Maximum number of exons to sample
            seed (int, optional): Random state. Defaults to 100.

        Returns:
            List[JunctionData]: List of {junction: SpliceJunction, win_start: win_start,
            win_end: win_end, sequence: sequence}
        """
        return self.regional_sampler.sample_exonic_regions(
            self.transcripts, target_count, seed
        )

    def extract_intergenic_sequences(
        self, target_count: int, seed: int = 100
    ) -> List[JunctionData]:
        """Sample intergenic regions of a defined length.

        Args:
            target_count (int): Number of sequences to sample
            seed (int, optional): Random state. Defaults to 100.

        Returns:
            List[JunctionData]: List of {junction: SpliceJunction, win_start: win_start,
            win_end: win_end, sequence: sequence}
        """
        return self.regional_sampler.sample_intergenic_regions(
            self.genes, self.chromosome_lengths, target_count, seed
        )

    def get_info(self) -> Dict[str, int]:
        """Returns parameters and info about the extractor. Includes exon bases,
        intron bases, window size, buffer size, transcript count, and splice junction count.

        Returns:
            Dict[str, int]: Dict of exon bases,
        intron bases, window size, buffer size, transcript count, and splice junction count.
        """
        return {
            "window_size": self.sampling_params.window_size,
            "exon_bases": self.extraction_params.n_exon,
            "intron_bases": self.extraction_params.n_intron,
            "buffer_size": self.sampling_params.buffer_size,
            "transcript_count": len(self.transcripts),
            "junction_count": len(self.junctions),
        }
