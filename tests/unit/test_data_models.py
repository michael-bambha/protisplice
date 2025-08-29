"""
File: test_data_models.py
Description: Updated unit tests for data structures
"""

import pytest

from protisplice import (
    ExtractionParams,
    SamplingParams,
    JunctionData,
    SpliceJunction,
    JunctionType,
    StrandType,
    Transcript,
    TranscriptInfo,
    Gene,
    TranscriptFilter,
)


class TestExtractionParams:
    """Test extraction parameters"""

    def test_valid_params(self):
        """Test with default params"""
        params = ExtractionParams(n_exon=40, n_intron=80, buffer_size=50)
        assert params.n_exon == 40
        assert params.n_intron == 80
        assert params.buffer_size == 50
        assert params.window_size == 120

    def test_default_params(self):
        """Test default parameter values"""
        params = ExtractionParams()
        assert params.n_exon == 40
        assert params.n_intron == 80
        assert params.buffer_size == 50
        assert params.window_size == 120

    def test_neg_params_raise_error(self):
        """Negative number in params should raise ValueError"""
        with pytest.raises(ValueError, match="All parameters must be non-negative!"):
            ExtractionParams(n_exon=-1, n_intron=80, buffer_size=50)

        with pytest.raises(ValueError, match="All parameters must be non-negative!"):
            ExtractionParams(n_exon=40, n_intron=-1, buffer_size=50)

        with pytest.raises(ValueError, match="All parameters must be non-negative!"):
            ExtractionParams(n_exon=40, n_intron=80, buffer_size=-1)

    def test_non_int_params_raise_error(self):
        """Non-integer types in params should raise ValueError"""
        with pytest.raises(ValueError, match="All parameters must be integers!"):
            ExtractionParams(n_exon=40.5, n_intron=80, buffer_size=50)

        with pytest.raises(ValueError, match="All parameters must be integers!"):
            ExtractionParams(n_exon=40, n_intron=80.5, buffer_size=50)

    def test_both_zero_raises_error(self):
        """ValueError should be raised if n_intron + n_exon = 0"""
        with pytest.raises(
            ValueError, match="At least one of n_exon or n_intron must be > 0"
        ):
            ExtractionParams(n_exon=0, n_intron=0, buffer_size=50)

    def test_one_zero_allowed(self):
        """One of n_exon or n_intron can be zero"""
        params1 = ExtractionParams(n_exon=0, n_intron=80, buffer_size=50)
        assert params1.window_size == 80

        params2 = ExtractionParams(n_exon=40, n_intron=0, buffer_size=50)
        assert params2.window_size == 40

    def test_win_size(self):
        """Test the window size calculation"""
        params = ExtractionParams(n_exon=30, n_intron=70, buffer_size=50)
        assert params.window_size == 100


class TestSamplingParams:
    """Test sampling parameters"""

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

    def test_zero_buffer_size_allowed(self):
        """Test that zero buffer size is allowed"""
        params = SamplingParams(window_size=100, buffer_size=0)
        assert params.buffer_size == 0


class TestJunctionData:
    """Test JunctionData dataclass"""

    def test_fasta_header_generation(self):
        """Test the FASTA headers generation from junctions"""
        junction = SpliceJunction(
            id="transcript1_donor_0",
            seqid="chr1",
            coord=200,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction, window_start=160, window_end=280, sequence="ATGC"
        )
        expected_header = ">chr1_donor_+_160_280"
        assert expected_header == junction_data.to_fasta_header()

    def test_fasta_header_negative_strand(self):
        """Test FASTA header for negative strand"""
        junction = SpliceJunction(
            id="transcript1_acceptor_1",
            seqid="chr2",
            coord=500,
            strand=StrandType.NEGATIVE,
            junction_type=JunctionType.ACCEPTOR,
        )
        junction_data = JunctionData(
            junction=junction, window_start=450, window_end=550, sequence="GCTA"
        )
        expected_header = ">chr2_acceptor_-_450_550"
        assert expected_header == junction_data.to_fasta_header()

    def test_fasta_header_different_junction_types(self):
        """Test FASTA headers for different junction types"""
        junction_types = [
            JunctionType.DONOR,
            JunctionType.ACCEPTOR,
            JunctionType.INTRON,
            JunctionType.EXON,
            JunctionType.INTERGENIC,
        ]

        for jtype in junction_types:
            junction = SpliceJunction(
                id=f"test_{jtype.value}",
                seqid="chr1",
                coord=100,
                strand=StrandType.POSITIVE,
                junction_type=jtype,
            )
            junction_data = JunctionData(
                junction=junction, window_start=50, window_end=150
            )
            header = junction_data.to_fasta_header()
            assert jtype.value in header

    def test_junction_data_with_none_sequence(self):
        """Test JunctionData with None sequence"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=100,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction, window_start=50, window_end=150, sequence=None
        )
        assert junction_data.sequence is None
        # Header should still work
        header = junction_data.to_fasta_header()
        assert ">chr1_donor_+_50_150" == header


class TestSpliceJunction:
    """Test SpliceJunction dataclass"""

    def test_splice_junction_creation(self):
        """Test creating splice junction"""
        junction = SpliceJunction(
            id="transcript1_donor_0",
            seqid="chr1",
            coord=200,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )

        assert junction.id == "transcript1_donor_0"
        assert junction.seqid == "chr1"
        assert junction.coord == 200
        assert junction.strand == StrandType.POSITIVE
        assert junction.junction_type == JunctionType.DONOR

    def test_splice_junction_different_types(self):
        """Test splice junctions with different types"""
        for jtype in JunctionType:
            junction = SpliceJunction(
                id=f"test_{jtype.value}",
                seqid="chr1",
                coord=100,
                strand=StrandType.POSITIVE,
                junction_type=jtype,
            )
            assert junction.junction_type == jtype


class TestTranscript:
    """Test Transcript dataclass"""

    def test_transcript_creation(self):
        """Test creating transcript"""
        info = TranscriptInfo(seqid="chr1", strand=StrandType.POSITIVE)
        transcript = Transcript(info=info, exons=[(100, 200), (300, 400)])

        assert transcript.info.seqid == "chr1"
        assert transcript.info.strand == StrandType.POSITIVE
        assert len(transcript.exons) == 2
        assert transcript.exons[0] == (100, 200)
        assert transcript.introns is None

    def test_transcript_with_introns(self):
        """Test transcript with introns"""
        info = TranscriptInfo(seqid="chr1", strand=StrandType.POSITIVE)
        transcript = Transcript(
            info=info, exons=[(100, 200), (300, 400)], introns=[(201, 299)]
        )

        assert len(transcript.introns) == 1
        assert transcript.introns[0] == (201, 299)


class TestGene:
    """Test Gene dataclass"""

    def test_gene_creation(self):
        """Test creating gene"""
        gene = Gene(
            gene_id="ENSG001",
            seq_id="chr1",
            start=1000,
            end=5000,
            strand=StrandType.POSITIVE,
        )

        assert gene.gene_id == "ENSG001"
        assert gene.seq_id == "chr1"
        assert gene.start == 1000
        assert gene.end == 5000
        assert gene.strand == StrandType.POSITIVE

    def test_gene_negative_strand(self):
        """Test gene on negative strand"""
        gene = Gene(
            gene_id="ENSG002",
            seq_id="chr2",
            start=2000,
            end=8000,
            strand=StrandType.NEGATIVE,
        )

        assert gene.strand == StrandType.NEGATIVE


class TestEnums:
    """Test enum classes"""

    def test_strand_type_enum(self):
        """Test StrandType enum"""
        assert StrandType.POSITIVE.value == "+"
        assert StrandType.NEGATIVE.value == "-"

    def test_junction_type_enum(self):
        """Test JunctionType enum"""
        assert JunctionType.DONOR.value == "donor"
        assert JunctionType.ACCEPTOR.value == "acceptor"
        assert JunctionType.INTRON.value == "intron"
        assert JunctionType.EXON.value == "exon"
        assert JunctionType.INTERGENIC.value == "intergenic"

    def test_transcript_filter_enum(self):
        """Test TranscriptFilter enum"""
        assert TranscriptFilter.ALL.value == "all"
        assert TranscriptFilter.PROTEIN_CODING.value == "protein_coding"
        assert TranscriptFilter.EXPRESSED.value == "expressed"


class TestTranscriptInfo:
    """Test TranscriptInfo dataclass"""

    def test_transcript_info_creation(self):
        """Test creating transcript info"""
        info = TranscriptInfo(seqid="chr1", strand=StrandType.POSITIVE)
        assert info.seqid == "chr1"
        assert info.strand == StrandType.POSITIVE

    def test_transcript_info_negative_strand(self):
        """Test transcript info on negative strand"""
        info = TranscriptInfo(seqid="chr2", strand=StrandType.NEGATIVE)
        assert info.seqid == "chr2"
        assert info.strand == StrandType.NEGATIVE
