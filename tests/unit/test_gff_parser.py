"""
File: test_gff_parser.py
Description: Test GFF parsing module
"""
# pylint:disable=W0212

import pytest
from protisplice import GFFParser
from protisplice import TranscriptFilter, StrandType
from protisplice import ExpressionFilter


class TestGFFParser:
    """Test GFFParser class"""

    def test_parse_transcripts_basic(self, test_gff_file):
        """Test basic GFF parsing functionality"""
        parser = GFFParser(str(test_gff_file), TranscriptFilter.ALL)
        transcripts = parser.parse_transcripts()
        assert len(transcripts) == 2
        assert "transcript1" in transcripts
        assert "transcript2" in transcripts

        t1 = transcripts["transcript1"]
        assert t1.info.strand == StrandType.POSITIVE
        assert t1.info.seqid == "chr1"
        assert len(t1.exons) == 3
        assert t1.exons == [(100, 200), (300, 400), (500, 600)]

    def test_parse_transcripts_with_filter(self, test_gff_file, test_expression_data):
        """Test GFF parsing with expression filtering"""
        expression_filter = ExpressionFilter(test_expression_data, threshold=1.0)
        parser = GFFParser(str(test_gff_file), TranscriptFilter.ALL, expression_filter)
        transcripts = parser.parse_transcripts()
        assert "transcript1" in transcripts
        assert "transcript2" not in transcripts  # expression of 0.5 < threshold

    def test_file_not_found(self):
        """Non-existent file should raise FileNotFoundError"""
        with pytest.raises(FileNotFoundError, match="GFF3 file .* not found"):
            GFFParser("doesnotexist.gff3", TranscriptFilter.ALL)

    def test_extract_gff_attrs(self, test_gff_file):
        """Test GFF attribute parsing"""
        parser = GFFParser(str(test_gff_file), TranscriptFilter.ALL)
        attrs_str = "ID=transcript1;Parent=gene1;biotype=protein_coding"
        attrs = parser._extract_gff_attrs(attrs_str)
        expected = {"ID": "transcript1", "Parent": "gene1", "biotype": "protein_coding"}
        assert attrs == expected
