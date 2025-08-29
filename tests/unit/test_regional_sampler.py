"""
File: test_regional_sampler.py
Description: Unit tests for regional sampling module
"""
# pylint:disable=protected-access

from unittest.mock import patch
import pytest
from protisplice.regional_sampler import RegionalSampler
from protisplice import (
    SamplingParams,
    Transcript,
    TranscriptInfo,
    Gene,
    StrandType,
    JunctionType,
)


class TestRegionalSampler:
    """Test RegionalSampler class"""

    def test_init_valid_fasta(self, test_fasta_file):
        """Test initialization with valid FASTA file"""
        params = SamplingParams(window_size=120, buffer_size=50)
        sampler = RegionalSampler(str(test_fasta_file), params)
        assert sampler.fasta_path == test_fasta_file
        assert sampler.params == params

    def test_init_missing_fasta(self):
        """Test initialization with missing FASTA raises FileNotFoundError"""
        params = SamplingParams()
        with pytest.raises(FileNotFoundError, match="FASTA file not found"):
            RegionalSampler("nonexistent.fasta", params)

    def test_get_chromosome_lengths(self, test_fasta_file):
        """Test getting chromosome lengths"""
        params = SamplingParams()
        sampler = RegionalSampler(str(test_fasta_file), params)
        lengths = sampler.get_chromosome_lengths()

        assert isinstance(lengths, dict)
        assert "chr1" in lengths
        assert "chr2" in lengths
        assert lengths["chr1"] > 0
        assert lengths["chr2"] > 0

    @patch("protisplice.regional_sampler.add_intron_coords")
    def test_sample_introns(self, mock_add_introns, test_fasta_file, sample_transcript):
        """Test intron sampling"""
        # Setup mock
        transcript_with_introns = Transcript(
            info=TranscriptInfo(seqid="chr1", strand=StrandType.POSITIVE),
            exons=[(100, 200), (300, 400)],
            introns=[(201, 299)],
        )
        mock_add_introns.return_value = {"transcript1": transcript_with_introns}

        params = SamplingParams(window_size=50, buffer_size=20)
        sampler = RegionalSampler(str(test_fasta_file), params)

        transcripts = {"transcript1": sample_transcript}
        result = sampler.sample_introns(transcripts, target_count=5, seed=42)

        assert isinstance(result, list)
        mock_add_introns.assert_called_once_with(transcripts)

    def test_sample_exonic_regions(self, test_fasta_file, sample_transcript):
        """Test exonic region sampling"""
        params = SamplingParams(window_size=50, buffer_size=10)
        sampler = RegionalSampler(str(test_fasta_file), params)

        transcripts = {"transcript1": sample_transcript}
        result = sampler.sample_exonic_regions(transcripts, target_count=2, seed=42)

        assert isinstance(result, list)
        # Check that junction data has correct junction type
        if result:
            assert all(jd.junction.junction_type == JunctionType.EXON for jd in result)

    def test_sample_intergenic_regions(self, test_fasta_file):
        """Test intergenic region sampling"""
        params = SamplingParams(window_size=50, buffer_size=10)
        sampler = RegionalSampler(str(test_fasta_file), params)

        # Create sample genes
        genes = {
            "gene1": Gene(
                gene_id="gene1",
                seq_id="chr1",
                start=100,
                end=200,
                strand=StrandType.POSITIVE,
            ),
            "gene2": Gene(
                gene_id="gene2",
                seq_id="chr1",
                start=300,
                end=400,
                strand=StrandType.POSITIVE,
            ),
        }

        chromosome_lengths = {"chr1": 1000, "chr2": 500}

        result = sampler.sample_intergenic_regions(
            genes, chromosome_lengths, target_count=2, seed=42
        )

        assert isinstance(result, list)
        if result:
            assert all(
                jd.junction.junction_type == JunctionType.INTERGENIC for jd in result
            )

    def test_identify_intergenic_regions(self, test_fasta_file):
        """Test identification of intergenic regions"""
        params = SamplingParams()
        sampler = RegionalSampler(str(test_fasta_file), params)

        genes = {
            "gene1": Gene(
                gene_id="gene1",
                seq_id="chr1",
                start=100,
                end=200,
                strand=StrandType.POSITIVE,
            ),
            "gene2": Gene(
                gene_id="gene2",
                seq_id="chr1",
                start=300,
                end=400,
                strand=StrandType.POSITIVE,
            ),
        }

        chromosome_lengths = {"chr1": 500}

        regions = sampler._identify_intergenic_regions(genes, chromosome_lengths)

        assert "chr1" in regions
        expected_regions = [
            (1, 99),  # Before gene1
            (201, 299),  # Between genes
            (401, 500),  # After gene2
        ]
        assert regions["chr1"] == expected_regions

    def test_sample_from_intron_too_small(self, test_fasta_file):
        """Test that small introns are skipped"""
        params = SamplingParams(window_size=100, buffer_size=50)  # Large requirements
        sampler = RegionalSampler(str(test_fasta_file), params)

        transcript = Transcript(
            info=TranscriptInfo(seqid="chr1", strand=StrandType.POSITIVE),
            exons=[(100, 200), (250, 300)],  # Small intron: 201-249 (49bp)
        )

        with patch("pysam.FastaFile") as mock_fasta:
            result = sampler._sample_from_intron(
                mock_fasta, "transcript1", transcript, 201, 249, 0
            )
            assert result is None

    def test_sample_from_exon_sufficient_size(self, test_fasta_file):
        """Test sampling from exon with sufficient size"""
        params = SamplingParams(window_size=50, buffer_size=10)
        sampler = RegionalSampler(str(test_fasta_file), params)

        exon_data = {
            "transcript_id": "transcript1",
            "exon_index": 0,
            "start": 100,
            "end": 200,  # 101bp exon
            "seqid": "chr1",
            "strand": StrandType.POSITIVE,
        }

        with patch("pysam.FastaFile") as mock_fasta, patch(
            "protisplice.regional_sampler.extract_sequence", return_value="A" * 50
        ):

            result = sampler._sample_from_exon(mock_fasta, exon_data, 0)

            assert result is not None
            assert result.junction.junction_type == JunctionType.EXON
            assert len(result.sequence) == 50

    def test_sample_from_intergenic_region(self, test_fasta_file):
        """Test sampling from intergenic region"""
        params = SamplingParams(window_size=50, buffer_size=10)
        sampler = RegionalSampler(str(test_fasta_file), params)

        region = (100, 200)  # 101bp region

        with patch("pysam.FastaFile") as mock_fasta, patch(
            "protisplice.regional_sampler.extract_sequence", return_value="A" * 50
        ):

            result = sampler._sample_from_intergenic_region(
                mock_fasta, "chr1", region, 0
            )

            assert result is not None
            assert result.junction.junction_type == JunctionType.INTERGENIC
            assert result.junction.seqid == "chr1"
            assert len(result.sequence) == 50


class TestSamplingParams:
    """Test SamplingParams data class"""

    def test_valid_params(self):
        """Test valid sampling parameters"""
        params = SamplingParams(window_size=120, buffer_size=50)
        assert params.window_size == 120
        assert params.buffer_size == 50

    def test_default_params(self):
        """Test default parameters"""
        params = SamplingParams()
        assert params.window_size == 120
        assert params.buffer_size == 50

    def test_invalid_window_size(self):
        """Test that invalid window size raises ValueError"""
        with pytest.raises(ValueError, match="Window size must be positive"):
            SamplingParams(window_size=0)

        with pytest.raises(ValueError, match="Window size must be positive"):
            SamplingParams(window_size=-10)

    def test_invalid_buffer_size(self):
        """Test that invalid buffer size raises ValueError"""
        with pytest.raises(ValueError, match="Buffer size must be non-negative"):
            SamplingParams(buffer_size=-5)
