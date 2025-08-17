"""
File: test_data_models.py
Description: Unit tests for data structures
"""

import pytest

from protisplice import (
    ExtractionParams,
    JunctionData,
    ExtractionResults,
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

    def test_neg_params_raise_error(self):
        """Negative number in params should raise ValueError"""
        with pytest.raises(ValueError, match="All parameters must be non-negative!"):
            ExtractionParams(n_exon=-1, n_intron=80, buffer_size=50)

    def test_non_int_params_raise_error(self):
        """Non-integer types in params should raise ValueError"""
        with pytest.raises(ValueError, match="All parameters must be integers!"):
            ExtractionParams(n_exon=40.5, n_intron=80, buffer_size=50)

    def test_both_zero_raises_error(self):
        """ValueError should be raised if n_intron + n_exon = 0"""
        with pytest.raises(
            ValueError, match="At least one of n_exon or n_intron must be > 0"
        ):
            ExtractionParams(n_exon=0, n_intron=0, buffer_size=50)

    def test_win_size(self):
        """Test the window size calculation"""
        params = ExtractionParams(n_exon=30, n_intron=70, buffer_size=50)
        assert params.window_size == 100


class TestJunctionData:
    """Test JunctionData dataclass"""

    def test_fasta_header_generation(self, sample_junctions):
        """Test the FASTA headers generation from junctions"""
        junction = sample_junctions[0]
        junction_data = JunctionData(
            junction=junction, window_start=160, window_end=280, sequence="ATGC"
        )
        expected_header = ">chr1_donor_+_160_280"
        assert expected_header == junction_data.to_fasta_header()


class TestExtractionResults:
    """Test ExtractionResults dataclass"""

    def test_empty_results(self):
        """Test ExtractionResults with no input"""
        results = ExtractionResults()
        assert results.positive_count == 0
        assert results.negative_count == 0
        assert results.total_count == 0

    def test_results_with_data(self, sample_junctions):
        """Test ExtractionResults with mock data"""
        junction_data = [
            JunctionData(junction=sample_junctions[0], window_start=1, window_end=100),
            JunctionData(junction=sample_junctions[1], window_start=1, window_end=100),
        ]
        results = ExtractionResults(
            positive_sequences=junction_data, negative_sequences=junction_data[:1]
        )
        assert results.positive_count == 2
        assert results.negative_count == 1
        assert results.total_count == 3

    def test_get_stats(self, sample_junctions):
        """Test the info retriever"""
        junction_data = [
            JunctionData(junction=sample_junctions[0], window_start=1, window_end=100)
        ]
        results = ExtractionResults(positive_sequences=junction_data)
        expected_stats = {
            "positive_count": 1,
            "negative_count": 0,
            "total_count": 1,
            "has_positive": True,
            "has_negative": False,
        }
        assert expected_stats == results.get_stats()
