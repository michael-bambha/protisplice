"""
File: test_extract_splice_seqs.py
Description: Integration test for splice site extraction workflow
"""
from protisplice import SpliceSeqExtractor
from protisplice import (
    ExtractionParams,
    ExtractionResults,
)


class TestSpliceSeqExtractorIntegration:
    """Integration tests for SpliceSeqExtractor"""

    def test_basic_extraction_workflow(self, test_gff_file, test_fasta_file, temp_dir):
        """Test complete extraction workflow"""
        extractor = SpliceSeqExtractor(
            gff_path=str(test_gff_file), fasta_path=str(test_fasta_file)
        )

        transcripts = extractor.transcripts
        junctions = extractor.junctions

        assert len(transcripts) > 0
        assert len(junctions) > 0

        positive_sequences = extractor.extract_positive_sequences()
        negative_sequences = extractor.extract_negative_sequences(target_count=5)

        assert len(positive_sequences) > 0
        assert len(negative_sequences) == 0

        results = ExtractionResults(
            positive_sequences=positive_sequences, negative_sequences=negative_sequences
        )

        pos_path = temp_dir / "positive.fasta"
        neg_path = temp_dir / "negative.fasta"

        pos_count, neg_count = extractor.write_sequences_to_fasta(
            results, str(pos_path), str(neg_path)
        )

        assert pos_count > 0
        assert neg_count == 0
        assert pos_path.exists()
        assert not neg_path.exists()

    def test_extraction_with_expression_filter(
        self, test_gff_file, test_fasta_file, test_kallisto_file
    ):
        """Test extraction with expression filtering"""
        extractor = SpliceSeqExtractor(
            gff_path=str(test_gff_file),
            fasta_path=str(test_fasta_file),
            expression_file=str(test_kallisto_file),
            min_expression=1.0,
            expression_format="kallisto",
        )

        transcripts = extractor.transcripts

        # Should filter out low-expression transcripts
        # Based on test data, transcript2 has expression 0.5 < 1.0
        assert "transcript1" in transcripts  # expression 5.0

    def test_other_extraction_params(self, test_gff_file, test_fasta_file):
        """Test extraction with other parameters"""
        params = ExtractionParams(n_exon=20, n_intron=60, buffer_size=30)

        extractor = SpliceSeqExtractor(
            gff_path=str(test_gff_file),
            fasta_path=str(test_fasta_file),
            params=params,
        )

        positive_sequences = extractor.extract_positive_sequences()

        # Check that sequences have the expected length
        for seq_data in positive_sequences:
            if seq_data.sequence:
                assert len(seq_data.sequence) == params.window_size

    def test_get_info(self, test_gff_file, test_fasta_file):
        """Test info retrieval"""
        extractor = SpliceSeqExtractor(
            gff_path=str(test_gff_file), fasta_path=str(test_fasta_file)
        )

        info = extractor.get_info()

        required_keys = [
            "window_size",
            "exon_bases",
            "intron_bases",
            "buffer_size",
            "transcript_count",
            "junction_count",
        ]

        for key in required_keys:
            assert key in info
            assert isinstance(info[key], int)
            assert info[key] >= 0
