"""
File: test_sequence_extraction.py
Description: Tests for sequence extraction logic
"""
# pylint: disable=protected-access
import pytest
from protisplice import SequenceExtractor
from protisplice import JunctionType, StrandType, SpliceJunction


class TestSequenceExtractor:
    """Test for the SequenceExtractor class"""

    def test_init_with_files(self, test_fasta_file, extraction_params):
        """Test that the init works correctly with valid files"""
        extractor = SequenceExtractor(str(test_fasta_file), extraction_params)
        assert extractor.fasta_path == test_fasta_file
        assert extractor.params == extraction_params

    def test_init_missing_fasta(self, extraction_params):
        """Test that FileNotFoundError is raised if FASTA is not found"""
        with pytest.raises(FileNotFoundError, match="FASTA file not found"):
            SequenceExtractor("doesnotexist.fasta", extraction_params)

    def test_init_missing_index(self, temp_dir, extraction_params):
        """Test that missing .fai creates new one"""
        fasta_path = temp_dir / "test.fasta"
        fai_path = fasta_path.with_suffix(fasta_path.suffix + ".fai")
        fasta_path.write_text(">chr1\nATCG\n")
        assert not fai_path.exists()
        SequenceExtractor(str(fasta_path), extraction_params)
        assert fai_path.exists()
        fai_path.unlink()

    def test_extract_splice_sites(
        self, test_fasta_file, extraction_params, sample_junctions
    ):
        """Test splice site sequence extraction"""
        extractor = SequenceExtractor(str(test_fasta_file), extraction_params)
        sequences = extractor.extract_splice_sites(sample_junctions)
        assert len(sequences) > 0
        for seq_data in sequences:
            assert seq_data.sequence is not None
            assert len(seq_data.sequence) == extraction_params.window_size
            assert seq_data.window_start is not None
            assert seq_data.window_end is not None

    def test_get_window_coords_positive_donor(self, extraction_params):
        """Test window coordinate calculation for positive strand donor"""
        extractor = SequenceExtractor.__new__(SequenceExtractor)
        extractor.params = extraction_params
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=200,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )

        start, end = extractor._get_window_coords(junction)

        # For positive donor: start = coord - n_exon + 1, end = coord + n_intron
        expected_start = 200 - 40 + 1  # 161
        expected_end = 200 + 80  # 280

        assert start == expected_start
        assert end == expected_end
