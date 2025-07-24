"""
File: test_expression.py
Description: Test for gene expression filtering module
"""

import pytest
from splice_seq_extractor.expression import ExpressionFilter, ExpressionParser


class TestExpressionFilter:
    """Test the ExpressionFilter class"""

    def test_passes_threshold(self, test_expression_data):
        """Test threshold checker"""
        filt = ExpressionFilter(test_expression_data, threshold=1.0)
        assert filt.passes_threshold("transcript1") is True  # 5.0 >= 1.0
        assert filt.passes_threshold("transcript2") is False  # 0.5 < 1.0
        assert filt.passes_threshold("transcript3") is True  # 2.0 >= 1.0
        assert filt.passes_threshold("nonexistent") is False  # 0.0 < 1.0

    def test_normalize_transcript_id(self):
        """Test transcript ID reformatting (normalization)"""
        # Test version num removal
        assert (
            ExpressionFilter.normalize_transcript_id("ENST00000123456.3")
            == "ENST00000123456"
        )
        # Prefix removal
        assert (
            ExpressionFilter.normalize_transcript_id("transcript:ENST00000123456")
            == "ENST00000123456"
        )
        # Version num and prefix
        assert (
            ExpressionFilter.normalize_transcript_id("transcript:ENST00000123456.3")
            == "ENST00000123456"
        )
        # No change
        assert (
            ExpressionFilter.normalize_transcript_id("ENST00000123456")
            == "ENST00000123456"
        )


class TestExpressionParser:
    """Test ExpressionParser class"""

    def test__parse_kallisto(self, test_kallisto_file):
        """Test Kallisto format parsing"""
        data = ExpressionParser.parse_file(str(test_kallisto_file), "kallisto")

        expected = {"transcript1": 5.0, "transcript2": 0.5, "transcript3": 2.0}

        assert data == expected

    def test_unsupported_format(self, test_kallisto_file):
        """Test unsupported formats raise ValueError"""
        with pytest.raises(ValueError, match="Unsupported format type for unsupported"):
            ExpressionParser.parse_file(str(test_kallisto_file), "unsupported")
