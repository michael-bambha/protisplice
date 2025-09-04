"""
File: test_extract_splice_seqs.py
Description: Updated unit tests for main SpliceSeqExtractor class
"""

# pylint:disable=protected-access

from unittest.mock import Mock, patch
from protisplice import (
    SpliceSeqExtractor,
    ExtractionParams,
    SamplingParams,
    TranscriptFilter,
    JunctionData,
)


class TestSpliceSeqExtractor:
    """Test main SpliceSeqExtractor class"""

    def test_init_basic(self, test_gff_file, test_fasta_file):
        """Test basic initialization"""
        extractor = SpliceSeqExtractor(str(test_gff_file), str(test_fasta_file))

        assert extractor.gff_path == str(test_gff_file)
        assert extractor.fasta_path == str(test_fasta_file)
        assert isinstance(extractor.extraction_params, ExtractionParams)
        assert isinstance(extractor.sampling_params, SamplingParams)
        assert extractor.transcript_filter == TranscriptFilter.ALL
        assert extractor.expression_filter is None

    def test_init_with_custom_params(self, test_gff_file, test_fasta_file):
        """Test initialization with custom parameters"""
        extraction_params = ExtractionParams(n_exon=50, n_intron=100, buffer_size=75)
        sampling_params = SamplingParams(window_size=150, buffer_size=60)

        extractor = SpliceSeqExtractor(
            str(test_gff_file),
            str(test_fasta_file),
            extraction_params=extraction_params,
            sampling_params=sampling_params,
            transcript_filter=TranscriptFilter.PROTEIN_CODING,
        )

        assert extractor.extraction_params == extraction_params
        assert extractor.sampling_params == sampling_params
        assert extractor.transcript_filter == TranscriptFilter.PROTEIN_CODING

    def test_init_with_expression_file(
        self, test_gff_file, test_fasta_file, test_kallisto_file
    ):
        """Test initialization with expression filtering"""
        extractor = SpliceSeqExtractor(
            str(test_gff_file),
            str(test_fasta_file),
            expression_file=str(test_kallisto_file),
            min_expression=2.0,
            expression_format="kallisto",
        )

        assert extractor.expression_filter is not None
        assert extractor.expression_filter.threshold == 2.0

    @patch("protisplice.extract_splice_seqs.GFFParser")
    def test_transcripts_lazy_loading(
        self, mock_gff_parser, test_gff_file, test_fasta_file
    ):
        """Test lazy loading of transcripts"""
        mock_parser_instance = Mock()
        mock_gff_parser.return_value = mock_parser_instance
        mock_parser_instance.parse_transcripts.return_value = {"transcript1": Mock()}

        extractor = SpliceSeqExtractor(str(test_gff_file), str(test_fasta_file))

        # First access should trigger parsing
        transcripts = extractor.transcripts
        mock_parser_instance.parse_transcripts.assert_called_once()

        # Second access should use cached result
        transcripts2 = extractor.transcripts
        assert mock_parser_instance.parse_transcripts.call_count == 1
        assert transcripts is transcripts2

    @patch("protisplice.extract_splice_seqs.SpliceJunctionExtractor")
    def test_junctions_lazy_loading(
        self, mock_junction_extractor, test_gff_file, test_fasta_file
    ):
        """Test lazy loading of junctions"""
        mock_extractor_instance = Mock()
        mock_junction_extractor.return_value = mock_extractor_instance
        mock_extractor_instance.identify_splice_junctions.return_value = [Mock()]

        extractor = SpliceSeqExtractor(str(test_gff_file), str(test_fasta_file))
        extractor._transcripts = {"transcript1": Mock()}  # Set cached transcripts

        # First access should trigger junction identification
        junctions = extractor.junctions
        mock_extractor_instance.identify_splice_junctions.assert_called_once()

        # Second access should use cached result
        junctions2 = extractor.junctions
        assert mock_extractor_instance.identify_splice_junctions.call_count == 1
        assert junctions is junctions2

    @patch("protisplice.extract_splice_seqs.RegionalSampler")
    def test_extract_intronic_regions(
        self, mock_sampler, test_gff_file, test_fasta_file
    ):
        """Test extracting intronic regions"""
        mock_sampler_instance = Mock()
        mock_sampler.return_value = mock_sampler_instance
        mock_junction_data = [Mock(spec=JunctionData)]
        mock_sampler_instance.sample_introns.return_value = mock_junction_data

        extractor = SpliceSeqExtractor(str(test_gff_file), str(test_fasta_file))
        extractor._transcripts = {"transcript1": Mock()}

        result = extractor.extract_intronic_regions(target_count=10, random_seed=123)

        mock_sampler_instance.sample_introns.assert_called_once_with(
            extractor._transcripts, 10, 123
        )
        assert result == mock_junction_data

    @patch("protisplice.extract_splice_seqs.RegionalSampler")
    def test_extract_exonic_sequences(
        self, mock_sampler, test_gff_file, test_fasta_file
    ):
        """Test extracting exonic sequences"""
        mock_sampler_instance = Mock()
        mock_sampler.return_value = mock_sampler_instance
        mock_junction_data = [Mock(spec=JunctionData)]
        mock_sampler_instance.sample_exonic_regions.return_value = mock_junction_data

        extractor = SpliceSeqExtractor(str(test_gff_file), str(test_fasta_file))
        extractor._transcripts = {"transcript1": Mock()}

        result = extractor.extract_exonic_sequences(target_count=15, seed=456)

        mock_sampler_instance.sample_exonic_regions.assert_called_once_with(
            extractor._transcripts, 15, 456
        )
        assert result == mock_junction_data

    @patch("protisplice.extract_splice_seqs.RegionalSampler")
    def test_extract_intergenic_sequences(
        self, mock_sampler, test_gff_file, test_fasta_file
    ):
        """Test extracting intergenic sequences"""
        mock_sampler_instance = Mock()
        mock_sampler.return_value = mock_sampler_instance
        mock_junction_data = [Mock(spec=JunctionData)]
        mock_sampler_instance.sample_intergenic_regions.return_value = (
            mock_junction_data
        )

        extractor = SpliceSeqExtractor(str(test_gff_file), str(test_fasta_file))
        extractor._genes = {"gene1": Mock()}
        extractor._chromosome_lengths = {"chr1": 1000}

        result = extractor.extract_intergenic_sequences(target_count=20, seed=789)

        mock_sampler_instance.sample_intergenic_regions.assert_called_once_with(
            extractor._genes, extractor._chromosome_lengths, 20, 789
        )
        assert result == mock_junction_data

    def test_get_info(self, test_gff_file, test_fasta_file):
        """Test getting extractor information"""
        extraction_params = ExtractionParams(n_exon=30, n_intron=70, buffer_size=40)
        sampling_params = SamplingParams(window_size=100, buffer_size=25)

        extractor = SpliceSeqExtractor(
            str(test_gff_file),
            str(test_fasta_file),
            extraction_params=extraction_params,
            sampling_params=sampling_params,
        )

        # Mock cached data
        extractor._transcripts = {"t1": Mock(), "t2": Mock()}
        extractor._junctions = [Mock(), Mock(), Mock()]

        info = extractor.get_info()

        expected = {
            "window_size": 100,
            "exon_bases": 30,
            "intron_bases": 70,
            "buffer_size": 25,
            "transcript_count": 2,
            "junction_count": 3,
        }

        assert info == expected

    def test_extract_splice_sites_alias(self, test_gff_file, test_fasta_file):
        """Test that extract_splice_sites is an alias for extract_positive_sequences"""
        with patch(
            "protisplice.extract_splice_seqs.SequenceExtractor"
        ) as mock_seq_extractor:
            mock_extractor_instance = Mock()
            mock_seq_extractor.return_value = mock_extractor_instance
            mock_junction_data = [Mock(spec=JunctionData)]
            mock_extractor_instance.extract_splice_sites.return_value = (
                mock_junction_data
            )

            extractor = SpliceSeqExtractor(str(test_gff_file), str(test_fasta_file))
            extractor._junctions = [Mock()]

            result = extractor.extract_splice_sites()

            # Should call the same method as extract_positive_sequences
            mock_extractor_instance.extract_splice_sites.assert_called_once()
            assert result == mock_junction_data


class TestSpliceSeqExtractorIntegration:
    """Integration tests for SpliceSeqExtractor"""

    def test_full_workflow_basic(self, test_gff_file, test_fasta_file):
        """Test basic end-to-end workflow"""
        extractor = SpliceSeqExtractor(str(test_gff_file), str(test_fasta_file))

        # These should work without errors
        transcripts = extractor.transcripts
        assert isinstance(transcripts, dict)

        junctions = extractor.junctions
        assert isinstance(junctions, list)

        info = extractor.get_info()
        assert isinstance(info, dict)
        assert all(
            key in info
            for key in [
                "window_size",
                "exon_bases",
                "intron_bases",
                "buffer_size",
                "transcript_count",
                "junction_count",
            ]
        )

    def test_full_workflow_with_expression(
        self, test_gff_file, test_fasta_file, test_kallisto_file
    ):
        """Test workflow with expression filtering"""
        extractor = SpliceSeqExtractor(
            str(test_gff_file),
            str(test_fasta_file),
            expression_file=str(test_kallisto_file),
            min_expression=1.0,
        )

        transcripts = extractor.transcripts
        # Should have filtered transcripts based on expression
        assert isinstance(transcripts, dict)

        info = extractor.get_info()
        assert info["transcript_count"] >= 0  # May be reduced due to filtering
