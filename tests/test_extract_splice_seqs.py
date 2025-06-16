"""
Unit tests for extract_splice_seqs.py
Run with: pytest tests/
"""
# pylint:disable=protected-access, redefined-outer-name, unused-argument, no-value-for-parameter

from pathlib import Path
import tempfile
import os
from unittest.mock import Mock, patch, mock_open
import pytest
from extract_splice_seqs import (
    SpliceSeqExtractor,
    StrandType,
    JunctionType,
    SpliceJunction,
    TranscriptInfo,
    Transcript,
    ExtractionParams,
    ExpressionFilter,
    ExpressionParser,
)

# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def sample_extraction_params():
    """Standard extraction parameters for testing"""
    return ExtractionParams(n_exon=40, n_intron=80, buffer_size=50)


@pytest.fixture
def sample_transcript():
    """Sample transcript with multiple exons"""
    info = TranscriptInfo(seqid="chr1", strand=StrandType.POSITIVE)
    exons = [(100, 200), (300, 400), (500, 600)]
    return Transcript(info=info, exons=exons)


@pytest.fixture
def sample_negative_strand_transcript():
    """Sample transcript on negative strand"""
    info = TranscriptInfo(seqid="chr1", strand=StrandType.NEGATIVE)
    exons = [(100, 200), (300, 400), (500, 600)]
    return Transcript(info=info, exons=exons)


@pytest.fixture
def sample_splice_junction():
    """Sample splice junction for testing"""
    return SpliceJunction(
        id="test_donor_1",
        seqid="chr1",
        coord=200,
        strand=StrandType.POSITIVE,
        junction_type=JunctionType.DONOR,
    )


@pytest.fixture
def sample_gff_content():
    """Sample GFF3 content for testing"""
    return """##gff-version 3
chr1	test	mRNA	100	600	.	+	.	ID=transcript:TRANS001;biotype=protein_coding
chr1	test	exon	100	200	.	+	.	Parent=transcript:TRANS001
chr1	test	exon	300	400	.	+	.	Parent=transcript:TRANS001
chr1	test	exon	500	600	.	+	.	Parent=transcript:TRANS001
chr2	test	mRNA	1000	1500	.	-	.	ID=transcript:TRANS002;biotype=lncRNA
chr2	test	exon	1000	1200	.	-	.	Parent=transcript:TRANS002
chr2	test	exon	1300	1500	.	-	.	Parent=transcript:TRANS002"""


@pytest.fixture
def sample_kallisto_content():
    """Sample Kallisto expression data"""
    return """target_id	length	eff_length	est_counts	tpm
TRANS001	500	350	100	5.0
TRANS002	400	250	50	2.5
TRANS003	600	450	200	10.0"""


@pytest.fixture
def mock_fasta_file():
    """Mock pysam.FastaFile for testing"""
    mock_fasta = Mock()
    mock_fasta.references = ["chr1", "chr2"]
    mock_fasta.get_reference_length.return_value = 10000
    mock_fasta.fetch.return_value = "ATCGATCGATCGATCG" * 10  # 160bp sequence
    return mock_fasta


@pytest.fixture
def temp_files():
    """Create temporary files for testing"""
    with tempfile.TemporaryDirectory() as temp_dir:
        gff_path = Path(temp_dir) / "test.gff3"
        fasta_path = Path(temp_dir) / "test.fa"
        fai_path = Path(temp_dir) / "test.fa.fai"

        # Create empty files
        gff_path.touch()
        fasta_path.touch()
        fai_path.touch()

        yield {
            "gff": str(gff_path),
            "fasta": str(fasta_path),
            "fai": str(fai_path),
            "dir": temp_dir,
        }


# ============================================================================
# DATA CLASS TESTS
# ============================================================================


class TestDataClasses:
    """Test data classes and enums"""

    def test_strand_type_enum(self):
        """Test StrandType enum values"""
        assert StrandType.POSITIVE.value == "+"
        assert StrandType.NEGATIVE.value == "-"

    def test_junction_type_enum(self):
        """Test JunctionType enum values"""
        assert JunctionType.DONOR.value == "donor"
        assert JunctionType.ACCEPTOR.value == "acceptor"
        assert JunctionType.INTRON.value == "intron"

    def test_splice_junction_creation(self):
        """Test SpliceJunction data class"""
        junction = SpliceJunction(
            id="test_junction",
            seqid="chr1",
            coord=100,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        assert junction.id == "test_junction"
        assert junction.seqid == "chr1"
        assert junction.coord == 100
        assert junction.strand == StrandType.POSITIVE
        assert junction.junction_type == JunctionType.DONOR

    def test_extraction_params_defaults(self):
        """Test ExtractionParams default values"""
        params = ExtractionParams()
        assert params.n_exon == 40
        assert params.n_intron == 80
        assert params.buffer_size == 50

    def test_transcript_with_introns(self, sample_transcript):
        """Test transcript with intron coordinates"""
        sample_transcript.introns = [(201, 299), (401, 499)]
        assert len(sample_transcript.introns) == 2
        assert sample_transcript.introns[0] == (201, 299)


# ============================================================================
# EXPRESSION FILTER TESTS
# ============================================================================


class TestExpressionFilter:
    """Test expression filtering functionality"""

    def test_expression_filter_creation(self):
        """Test ExpressionFilter creation"""
        data = {"TRANS001": 5.0, "TRANS002": 2.5}
        filter_obj = ExpressionFilter(data=data, threshold=3.0)
        assert filter_obj.threshold == 3.0
        assert len(filter_obj.data) == 2

    def test_passes_threshold_true(self):
        """Test expression passes threshold"""
        data = {"TRANS001": 5.0}
        filter_obj = ExpressionFilter(data=data, threshold=3.0)
        assert filter_obj.passes_threshold("TRANS001") is True

    def test_passes_threshold_false(self):
        """Test expression fails threshold"""
        data = {"TRANS001": 2.0}
        filter_obj = ExpressionFilter(data=data, threshold=3.0)
        assert filter_obj.passes_threshold("TRANS001") is False

    def test_missing_transcript_returns_zero(self):
        """Test missing transcript returns 0.0 expression"""
        data = {"TRANS001": 5.0}
        filter_obj = ExpressionFilter(data=data, threshold=3.0)
        assert filter_obj._find_expression_value("MISSING") == 0.0
        assert filter_obj.passes_threshold("MISSING") is False

    def test_normalize_transcript_id(self):
        """Test transcript ID normalization"""
        data = {"TRANS001.1": 5.0}
        filter_obj = ExpressionFilter(data=data, threshold=3.0)

        # Test with version number
        assert filter_obj.passes_threshold("TRANS001") is True
        assert filter_obj.passes_threshold("TRANS001.1") is True

        # Test with transcript prefix
        data = {"transcript:TRANS001": 5.0}
        filter_obj = ExpressionFilter(data=data, threshold=3.0)
        assert filter_obj.passes_threshold("TRANS001") is True


# ============================================================================
# EXPRESSION PARSER TESTS
# ============================================================================


class TestExpressionParser:
    """Test expression file parsing"""

    def test_parse_kallisto_file(self, sample_kallisto_content):
        """Test parsing Kallisto expression file"""
        with patch("builtins.open", mock_open(read_data=sample_kallisto_content)):
            result = ExpressionParser.parse_file("test.tsv", "kallisto")

            assert len(result) == 3
            assert result["TRANS001"] == 5.0
            assert result["TRANS002"] == 2.5
            assert result["TRANS003"] == 10.0

    def test_parse_unsupported_format(self):
        """Test parsing unsupported format returns None"""
        result = ExpressionParser.parse_file("test.tsv", "unsupported")
        assert result is None


# ============================================================================
# SPLICE SEQ EXTRACTOR TESTS
# ============================================================================


class TestSpliceSeqExtractor:
    """Test main SpliceSeqExtractor class"""

    def test_initialization(self, temp_files, sample_extraction_params):
        """Test extractor initialization"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        assert extractor.gff_path == Path(temp_files["gff"])
        assert extractor.fasta_path == Path(temp_files["fasta"])
        assert extractor.params == sample_extraction_params
        assert extractor.transcript_filter == "all"

    def test_initialization_with_expression(
        self, temp_files, sample_extraction_params, sample_kallisto_content
    ):
        """Test extractor initialization with expression data"""
        expr_file = Path(temp_files["dir"]) / "expression.tsv"
        with open(expr_file, "w", encoding="utf-8") as f:
            f.write(sample_kallisto_content)

        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
            expression_file=str(expr_file),
            min_expression=3.0,
        )

        assert extractor.expression_filter is not None
        assert extractor.expression_filter.threshold == 3.0

    def test_validation_missing_fasta(self, temp_files, sample_extraction_params):
        """Test validation fails for missing FASTA file"""
        os.remove(temp_files["fasta"])

        with pytest.raises(FileNotFoundError, match="FASTA file not found"):
            SpliceSeqExtractor(
                gff_path=temp_files["gff"],
                fasta_path=temp_files["fasta"],
                params=sample_extraction_params,
            )

    def test_validation_missing_fai(self, temp_files, sample_extraction_params):
        """Test validation fails for missing FASTA index"""
        os.remove(temp_files["fai"])

        with pytest.raises(FileNotFoundError, match="Fasta index.*not found"):
            SpliceSeqExtractor(
                gff_path=temp_files["gff"],
                fasta_path=temp_files["fasta"],
                params=sample_extraction_params,
            )

    def test_validation_negative_params(self, temp_files):
        """Test validation fails for negative parameters"""
        with pytest.raises(ValueError, match="must be positive integers"):
            SpliceSeqExtractor(
                gff_path=temp_files["gff"],
                fasta_path=temp_files["fasta"],
                params=ExtractionParams(n_exon=-1, n_intron=80),
            )

    def test_get_sequence_stats(self, temp_files, sample_extraction_params):
        """Test sequence statistics calculation"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        stats = extractor.get_sequence_stats()
        assert stats["window_size"] == 120  # 40 + 80
        assert stats["exon_bases"] == 40
        assert stats["intron_bases"] == 80
        assert stats["buffer_size"] == 50


# ============================================================================
# GFF PARSING TESTS
# ============================================================================


class TestGFFParsing:
    """Test GFF3 file parsing functionality"""

    def test_parse_gff_line_valid_exon(self, temp_files, sample_extraction_params):
        """Test parsing valid exon line"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        line = "chr1\ttest\texon\t100\t200\t.\t+\t.\tParent=transcript:TRANS001"
        result = extractor._parse_gff_line(line, 1, {})

        assert result is not None
        transcript_id, seqid, start, end, strand = result
        assert transcript_id == "TRANS001"
        assert seqid == "chr1"
        assert start == 100
        assert end == 200
        assert strand == "+"

    def test_parse_gff_line_comment(self, temp_files, sample_extraction_params):
        """Test parsing comment line returns None"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        line = "##gff-version 3"
        result = extractor._parse_gff_line(line, 1, {})
        assert result is None

    def test_parse_gff_line_non_exon(self, temp_files, sample_extraction_params):
        """Test parsing non-exon line returns None"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        line = "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=GENE001"
        result = extractor._parse_gff_line(line, 1, {})
        assert result is None

    def test_extract_transcript_attributes(self, temp_files, sample_extraction_params):
        """Test extracting attributes from GFF3 field"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        attrs_str = "ID=transcript:TRANS001;Parent=GENE001;biotype=protein_coding"
        attrs = extractor._extract_transcript_attributes(attrs_str)

        assert attrs["ID"] == "transcript:TRANS001"
        assert attrs["Parent"] == "GENE001"
        assert attrs["biotype"] == "protein_coding"

    @patch("builtins.open", new_callable=mock_open)
    def test_parse_transcripts(
        self, mock_file, temp_files, sample_extraction_params, sample_gff_content
    ):
        """Test parsing complete transcript information"""
        mock_file.return_value.read_data = sample_gff_content
        mock_file.return_value.__iter__ = lambda self: iter(
            sample_gff_content.splitlines()
        )

        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        transcripts = extractor._parse_transcripts()

        assert len(transcripts) == 2
        assert "TRANS001" in transcripts
        assert "TRANS002" in transcripts

        trans1 = transcripts["TRANS001"]
        assert trans1.info.seqid == "chr1"
        assert trans1.info.strand == StrandType.POSITIVE
        assert len(trans1.exons) == 3


# ============================================================================
# JUNCTION EXTRACTION TESTS
# ============================================================================


class TestJunctionExtraction:
    """Test splice junction identification"""

    def test_get_splice_junctions_positive_strand(self, sample_transcript):
        """Test junction extraction for positive strand transcript"""
        extractor = Mock()
        extractor._create_junctions_for_transcript = (
            SpliceSeqExtractor._create_junctions_for_transcript.__get__(extractor)
        )
        junctions = extractor._create_junctions_for_transcript(
            "TRANS001", sample_transcript, sample_transcript.exons
        )

        # Should have 2 acceptors and 2 donors for 3 exons
        assert len(junctions) == 4

        # Check junction types and coordinates
        junction_types = [j.junction_type for j in junctions]
        assert junction_types.count(JunctionType.ACCEPTOR) == 2
        assert junction_types.count(JunctionType.DONOR) == 2

    def test_get_splice_junctions_negative_strand(
        self, sample_negative_strand_transcript
    ):
        """Test junction extraction for negative strand transcript"""
        extractor = Mock()
        extractor._create_junctions_for_transcript = (
            SpliceSeqExtractor._create_junctions_for_transcript.__get__(extractor)
        )

        junctions = extractor._create_junctions_for_transcript(
            "TRANS001",
            sample_negative_strand_transcript,
            sample_negative_strand_transcript.exons,
        )

        assert len(junctions) == 4

        # For negative strand, coordinates should be different
        acceptor_coords = [
            j.coord for j in junctions if j.junction_type == JunctionType.ACCEPTOR
        ]
        donor_coords = [
            j.coord for j in junctions if j.junction_type == JunctionType.DONOR
        ]

        assert len(acceptor_coords) == 2
        assert len(donor_coords) == 2

    def test_single_exon_transcript_no_junctions(self):
        """Test that single exon transcripts produce no junctions"""
        info = TranscriptInfo(seqid="chr1", strand=StrandType.POSITIVE)
        single_exon_transcript = Transcript(info=info, exons=[(100, 200)])

        extractor = Mock()
        extractor._create_junctions_for_transcript = (
            SpliceSeqExtractor._create_junctions_for_transcript.__get__(extractor)
        )

        junctions = extractor._create_junctions_for_transcript(
            "TRANS001", single_exon_transcript, single_exon_transcript.exons
        )

        assert len(junctions) == 0


# ============================================================================
# COORDINATE CALCULATION TESTS
# ============================================================================


class TestCoordinateCalculation:
    """Test window coordinate calculations"""

    def test_get_window_coords_positive_donor(
        self, temp_files, sample_extraction_params
    ):
        """Test window coordinates for positive strand donor"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=200,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )

        start, end = extractor._get_window_coords(junction)
        # Donor positive: start = coord - n_exon + 1, end = coord + n_intron
        assert start == 161  # 200 - 40 + 1
        assert end == 280  # 200 + 80

    def test_get_window_coords_positive_acceptor(
        self, temp_files, sample_extraction_params
    ):
        """Test window coordinates for positive strand acceptor"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=300,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.ACCEPTOR,
        )

        start, end = extractor._get_window_coords(junction)
        # Acceptor positive: start = coord - n_intron, end = coord + n_exon - 1
        assert start == 220  # 300 - 80
        assert end == 339  # 300 + 40 - 1

    def test_get_window_coords_negative_donor(
        self, temp_files, sample_extraction_params
    ):
        """Test window coordinates for negative strand donor"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=200,
            strand=StrandType.NEGATIVE,
            junction_type=JunctionType.DONOR,
        )

        start, end = extractor._get_window_coords(junction)
        # Donor negative: start = coord - n_intron, end = coord + n_exon - 1
        assert start == 120  # 200 - 80
        assert end == 239  # 200 + 40 - 1


# ============================================================================
# SEQUENCE EXTRACTION TESTS
# ============================================================================


class TestSequenceExtraction:
    """Test sequence extraction from FASTA"""

    def test_extract_sequence_valid(
        self, temp_files, sample_extraction_params, mock_fasta_file
    ):
        """Test successful sequence extraction"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        seq = extractor._extract_sequence(mock_fasta_file, "chr1", 100, 200)
        assert seq is not None
        assert len(seq) > 0
        mock_fasta_file.fetch.assert_called_once()

    def test_extract_sequence_invalid_coords(
        self, temp_files, sample_extraction_params, mock_fasta_file
    ):
        """Test sequence extraction with invalid coordinates"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        # Invalid coordinates (start > end)
        seq = extractor._extract_sequence(mock_fasta_file, "chr1", 200, 100)
        assert seq is None

        # Invalid coordinates (start < 1)
        seq = extractor._extract_sequence(mock_fasta_file, "chr1", 0, 100)
        assert seq is None

    def test_extract_sequence_missing_seqid(self, temp_files, sample_extraction_params):
        """Test sequence extraction with missing sequence ID"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        mock_fasta = Mock()
        mock_fasta.references = ["chr1"]  # chr2 not present

        seq = extractor._extract_sequence(mock_fasta, "chr2", 100, 200)
        assert seq is None


# ============================================================================
# INTRON COORDINATE TESTS
# ============================================================================


class TestIntronCoordinates:
    """Test intron coordinate calculation"""

    def test_add_intron_coords(
        self, sample_transcript, temp_files, sample_extraction_params
    ):
        """Test adding intron coordinates to transcripts"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        transcripts = {"TRANS001": sample_transcript}
        transcripts_with_introns = extractor._add_intron_coords(transcripts)

        transcript = transcripts_with_introns["TRANS001"]
        assert transcript.introns is not None
        assert len(transcript.introns) == 2

        # Check intron coordinates
        # Exons: (100,200), (300,400), (500,600)
        # Introns: (201,299), (401,499)
        assert transcript.introns[0] == (201, 299)
        assert transcript.introns[1] == (401, 499)

    def test_add_intron_coords_single_exon(self, temp_files, sample_extraction_params):
        """Test intron coordinates for single exon transcript"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        info = TranscriptInfo(seqid="chr1", strand=StrandType.POSITIVE)
        single_exon = Transcript(info=info, exons=[(100, 200)])
        transcripts = {"TRANS001": single_exon}

        transcripts_with_introns = extractor._add_intron_coords(transcripts)
        transcript = transcripts_with_introns["TRANS001"]

        assert transcript.introns == []


# ============================================================================
# FILTERING TESTS
# ============================================================================


class TestFiltering:
    """Test transcript filtering functionality"""

    def test_passes_filter_all(self, temp_files, sample_extraction_params):
        """Test 'all' filter passes everything"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
            transcript_filter="all",
        )

        # Should pass regardless of biotype
        assert extractor._passes_filter("TRANS001", {"TRANS001": "lncRNA"}) is True
        assert (
            extractor._passes_filter("TRANS002", {"TRANS002": "protein_coding"}) is True
        )

    def test_passes_filter_protein_coding(self, temp_files, sample_extraction_params):
        """Test protein_coding filter"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
            transcript_filter="protein_coding",
        )

        biotypes = {"TRANS001": "protein_coding", "TRANS002": "lncRNA"}

        assert extractor._passes_filter("TRANS001", biotypes) is True
        assert extractor._passes_filter("TRANS002", biotypes) is False

    def test_passes_filter_with_expression(
        self, temp_files, sample_extraction_params, sample_kallisto_content
    ):
        """Test filtering with expression threshold"""
        expr_file = Path(temp_files["dir"]) / "expression.tsv"
        with open(expr_file, "w", encoding="utf-8") as f:
            f.write(sample_kallisto_content)

        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
            expression_file=str(expr_file),
            min_expression=3.0,
        )

        # TRANS001 has TPM=5.0, TRANS002 has TPM=2.5
        assert extractor._passes_filter("TRANS001", {}) is True
        assert extractor._passes_filter("TRANS002", {}) is False


# ============================================================================
# INTEGRATION TESTS
# ============================================================================


class TestIntegration:
    """Integration tests for complete workflows"""

    @patch("pysam.FastaFile")
    @patch("builtins.open", new_callable=mock_open)
    def test_extract_sequences_workflow(
        self,
        mock_file,
        mock_fasta_class,
        temp_files,
        sample_extraction_params,
        sample_gff_content,
    ):
        """Test complete sequence extraction workflow"""
        # Setup mocks for file reading
        mock_file.return_value.__iter__ = lambda self: iter(
            sample_gff_content.splitlines()
        )

        # Setup mock FASTA file
        mock_fasta = Mock()
        mock_fasta.references = ["chr1", "chr2"]
        mock_fasta.get_reference_length.return_value = 10000
        mock_fasta.fetch.return_value = "ATCGATCGATCGATCG" * 10  # 160bp sequence
        mock_fasta_class.return_value.__enter__.return_value = mock_fasta

        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        with tempfile.NamedTemporaryFile(
            mode="w", delete=False
        ) as pos_file, tempfile.NamedTemporaryFile(mode="w", delete=False) as neg_file:

            pos_count, neg_count = extractor.extract_sequences(
                pos_file.name, neg_file.name
            )

            assert pos_count > 0
            # Note: neg_count might be 0 if introns are too small for the buffer size
            # This is expected behavior, so we test that the function completes successfully
            assert neg_count >= 0

            # Clean up
            os.unlink(pos_file.name)
            os.unlink(neg_file.name)

    @patch("pysam.FastaFile")
    @patch("builtins.open", new_callable=mock_open)
    def test_extract_sequences_with_large_introns(
        self, mock_file, mock_fasta_class, temp_files
    ):
        """Test extraction with large introns that can accommodate negative samples"""
        # Create GFF content with large introns
        large_intron_gff = """##gff-version 3
chr1\ttest\tmRNA\t100\t2000\t.\t+\t.\tID=transcript:TRANS001;biotype=protein_coding
chr1\ttest\texon\t100\t200\t.\t+\t.\tParent=transcript:TRANS001
chr1\ttest\texon\t1800\t2000\t.\t+\t.\tParent=transcript:TRANS001"""

        # Setup mocks
        mock_file.return_value.__iter__ = lambda self: iter(
            large_intron_gff.splitlines()
        )

        mock_fasta = Mock()
        mock_fasta.references = ["chr1"]
        mock_fasta.get_reference_length.return_value = 10000
        mock_fasta.fetch.return_value = "ATCGATCGATCGATCG" * 20  # Longer sequence
        mock_fasta_class.return_value.__enter__.return_value = mock_fasta

        # Use smaller parameters to ensure negative samples can be extracted
        params = ExtractionParams(n_exon=20, n_intron=40, buffer_size=25)

        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"], fasta_path=temp_files["fasta"], params=params
        )

        with tempfile.NamedTemporaryFile(
            mode="w", delete=False
        ) as pos_file, tempfile.NamedTemporaryFile(mode="w", delete=False) as neg_file:

            pos_count, neg_count = extractor.extract_sequences(
                pos_file.name, neg_file.name
            )

            assert pos_count > 0
            assert neg_count > 0  # Should have negative samples with large introns

            # Clean up
            os.unlink(pos_file.name)
            os.unlink(neg_file.name)


# ============================================================================
# ERROR HANDLING TESTS
# ============================================================================


class TestErrorHandling:
    """Test error handling and edge cases"""

    def test_malformed_gff_line_returns_none(
        self, temp_files, sample_extraction_params
    ):
        """Test that malformed GFF lines return None"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        # Test with insufficient fields
        malformed_line = "chr1\ttest\texon\t100"  # Missing fields
        result = extractor._parse_gff_line(malformed_line, 1, {})
        assert result is None

        # Test with empty line
        empty_line = ""
        result = extractor._parse_gff_line(empty_line, 1, {})
        assert result is None

        # Test with comment line
        comment_line = "# This is a comment"
        result = extractor._parse_gff_line(comment_line, 1, {})
        assert result is None

    def test_invalid_coordinates_in_gff_returns_none(
        self, temp_files, sample_extraction_params
    ):
        """Test that invalid coordinates return None"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        # Test with non-numeric coordinates (this should return None)
        invalid_line = "chr1\ttest\texon\tXXX\t200\t.\t+\t.\tParent=transcript:TRANS001"
        result = extractor._parse_gff_line(invalid_line, 1, {})
        assert result is None

        # Test with negative coordinates (this actually parses successfully!)
        # The coordinate validation happens later in _extract_sequence, not here
        negative_line = (
            "chr1\ttest\texon\t-100\t200\t.\t+\t.\tParent=transcript:TRANS001"
        )
        result = extractor._parse_gff_line(negative_line, 1, {})
        assert result is not None  # _parse_gff_line doesn't validate coordinate values
        _, _, start, _, _ = result
        assert start == -100  # The negative coordinate is parsed as-is

    def test_gff_line_without_parent_returns_none(
        self, temp_files, sample_extraction_params
    ):
        """Test that exon lines without Parent attribute return None"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        # Test with missing Parent attribute
        line_no_parent = "chr1\ttest\texon\t100\t200\t.\t+\t.\tID=exon001"
        result = extractor._parse_gff_line(line_no_parent, 1, {})
        assert result is None

        # Test with empty attributes
        line_empty_attrs = "chr1\ttest\texon\t100\t200\t.\t+\t."
        result = extractor._parse_gff_line(line_empty_attrs, 1, {})
        assert result is None

    def test_error_handling_graceful_degradation(
        self, temp_files, sample_extraction_params
    ):
        """Test that all types of malformed lines are handled gracefully"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        # Test comprehensive list of malformed lines
        malformed_lines = [
            "chr1\ttest\texon\t100",  # Too few fields
            "chr1\ttest\texon\tINVALID\t200\t.\t+\t.\tParent=transcript:TRANS001",  # Invalid start
            "chr1\ttest\texon\t100\tINVALID\t.\t+\t.\tParent=transcript:TRANS001",  # Invalid end
            "chr1\ttest\tgene\t100\t200\t.\t+\t.\tParent=transcript:TRANS001",  # Not an exon
            "",  # Empty line
            "# Comment line",  # Comment
            "chr1\ttest\texon\t100\t200\t.\t+\t.\tID=exon001",  # No Parent
            "chr1\ttest\texon\t100\t200\t.\t+\t.",  # Empty attributes
        ]

        # All malformed lines should return None (graceful error handling)
        for i, line in enumerate(malformed_lines):
            result = extractor._parse_gff_line(line, i + 1, {})
            assert result is None, f"Expected None for malformed line {i+1}: {line}"

        # But a valid line should work
        valid_line = "chr1\ttest\texon\t100\t200\t.\t+\t.\tParent=transcript:TRANS001"
        result = extractor._parse_gff_line(valid_line, 1, {})
        assert result is not None, "Valid line should parse successfully"

        # Verify the parsed data
        transcript_id, seqid, start, end, strand = result
        assert transcript_id == "TRANS001"
        assert seqid == "chr1"
        assert start == 100
        assert end == 200
        assert strand == "+"

    def test_coordinate_validation_in_sequence_extraction(
        self, temp_files, sample_extraction_params
    ):
        """Test that coordinate validation happens during sequence extraction"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        mock_fasta = Mock()
        mock_fasta.references = ["chr1"]
        mock_fasta.get_reference_length.return_value = 1000

        # Test that negative coordinates are rejected during sequence extraction
        seq = extractor._extract_sequence(mock_fasta, "chr1", -100, 200)
        assert seq is None

        # Test that start > end is rejected
        seq = extractor._extract_sequence(mock_fasta, "chr1", 200, 100)
        assert seq is None

        # Test that valid coordinates work
        mock_fasta.fetch.return_value = "ATCGATCGATCGATCG"
        seq = extractor._extract_sequence(mock_fasta, "chr1", 100, 200)
        assert seq is not None

    def test_empty_gff_file_handling(self, temp_files, sample_extraction_params):
        """Test handling of empty or comment-only GFF files"""
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        with patch(
            "builtins.open",
            mock_open(read_data="##gff-version 3\n# This is a comment\n"),
        ):
            transcripts = extractor._parse_transcripts()
            assert len(transcripts) == 0


# ============================================================================
# PARAMETERIZED TESTS
# ============================================================================


class TestParameterized:
    """Parameterized tests for different scenarios"""

    @pytest.mark.parametrize(
        "strand,junction_type,expected_start_offset,expected_end_offset",
        [
            (
                StrandType.POSITIVE,
                JunctionType.DONOR,
                -39,
                80,
            ),  # coord - n_exon + 1, coord + n_intron
            (
                StrandType.POSITIVE,
                JunctionType.ACCEPTOR,
                -80,
                39,
            ),  # coord - n_intron, coord + n_exon - 1
            (
                StrandType.NEGATIVE,
                JunctionType.DONOR,
                -80,
                39,
            ),  # coord - n_intron, coord + n_exon - 1
            (
                StrandType.NEGATIVE,
                JunctionType.ACCEPTOR,
                -39,
                80,
            ),  # coord - n_exon + 1, coord + n_intron
        ],
    )
    def test_window_coordinates_all_combinations(
        self,
        temp_files,
        strand,
        junction_type,
        expected_start_offset,
        expected_end_offset,
    ):
        """Test window coordinate calculation for all strand/junction combinations"""
        params = ExtractionParams(n_exon=40, n_intron=80)
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"], fasta_path=temp_files["fasta"], params=params
        )

        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=1000,
            strand=strand,
            junction_type=junction_type,
        )

        start, end = extractor._get_window_coords(junction)

        assert start == 1000 + expected_start_offset
        assert end == 1000 + expected_end_offset

    @pytest.mark.parametrize(
        "transcript_filter,biotype,expected",
        [
            ("all", "protein_coding", True),
            ("all", "lncRNA", True),
            ("all", "pseudogene", True),
            ("protein_coding", "protein_coding", True),
            ("protein_coding", "lncRNA", False),
            ("protein_coding", "pseudogene", False),
        ],
    )
    def test_transcript_filtering(
        self, temp_files, transcript_filter, biotype, expected
    ):
        """Test transcript filtering with different biotypes"""
        params = ExtractionParams()
        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=params,
            transcript_filter=transcript_filter,
        )

        biotypes = {"TEST_TRANS": biotype}
        result = extractor._passes_filter("TEST_TRANS", biotypes)

        assert result == expected


# ============================================================================
# PERFORMANCE TESTS
# ============================================================================


class TestPerformance:
    """Tests for performance-critical functionality"""

    def test_large_transcript_processing(self, temp_files, sample_extraction_params):
        """Test processing transcript with many exons"""
        # Create transcript with 100 exons
        info = TranscriptInfo(seqid="chr1", strand=StrandType.POSITIVE)
        exons = [(i * 200, i * 200 + 100) for i in range(100)]
        large_transcript = Transcript(info=info, exons=exons)

        extractor = SpliceSeqExtractor(
            gff_path=temp_files["gff"],
            fasta_path=temp_files["fasta"],
            params=sample_extraction_params,
        )

        junctions = extractor._create_junctions_for_transcript(
            "LARGE_TRANS", large_transcript, large_transcript.exons
        )

        # Should have 99 acceptors and 99 donors for 100 exons
        assert len(junctions) == 198

        # Test intron coordinate calculation
        transcripts = {"LARGE_TRANS": large_transcript}
        transcripts_with_introns = extractor._add_intron_coords(transcripts)

        # Should have 99 introns
        assert len(transcripts_with_introns["LARGE_TRANS"].introns) == 99


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
