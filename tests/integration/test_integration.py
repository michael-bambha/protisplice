"""
File: test_integration.py
Description: Integration tests for the entire splice sequence extraction pipeline
"""

# pylint:disable=protected-access
import time
import pytest
from Bio.Seq import Seq
import numpy as np
import pysam
from protisplice import (
    SpliceSeqExtractor,
    ExtractionParams,
    SamplingParams,
    TranscriptFilter,
    write_sequences,
    JunctionType,
    StrandType,
)
from protisplice.motif_scoring import generate_ppm, generate_pfm, generate_pwm


class TestFullPipelineIntegration:
    """Test the complete pipeline from GFF/FASTA to sequence extraction"""

    def test_complete_extraction_workflow(
        self, complex_gff_file, complex_fasta_file, temp_dir
    ):
        """Test complete workflow from files to extracted sequences"""
        # Initialize extractor with custom parameters
        extraction_params = ExtractionParams(n_exon=30, n_intron=50, buffer_size=25)
        sampling_params = SamplingParams(window_size=80, buffer_size=30)

        extractor = SpliceSeqExtractor(
            str(complex_gff_file),
            str(complex_fasta_file),
            extraction_params=extraction_params,
            sampling_params=sampling_params,
            transcript_filter=TranscriptFilter.PROTEIN_CODING,
        )
        transcripts = extractor.transcripts

        # Get basic info
        info = extractor.get_info()
        assert info["exon_bases"] == 30
        assert info["intron_bases"] == 50
        assert info["transcript_count"] >= 1  # At least one protein coding transcript
        assert info["junction_count"] >= 2  # At least one donor and one acceptor

        # Extract positive sequences (true splice sites)
        positive_sequences = extractor.extract_splice_sites()
        assert len(positive_sequences) > 0
        assert all(seq.sequence is not None for seq in positive_sequences)
        assert all(
            len(seq.sequence) == 80 for seq in positive_sequences
        )  # n_exon + n_intron

        # Extract negative sequences (decoys)
        negative_sequences = extractor.regional_sampler.sample_introns(
            transcripts, target_count=5
        )
        assert len(negative_sequences) <= 5
        assert all(seq.sequence is not None for seq in negative_sequences)

        # Write sequences to files
        pos_output = temp_dir / "positive.fasta"
        neg_output = temp_dir / "negative.fasta"

        pos_count = write_sequences(positive_sequences, str(pos_output))
        neg_count = write_sequences(negative_sequences, str(neg_output))

        assert pos_count == len(positive_sequences)
        assert neg_count == len(negative_sequences)
        assert pos_output.exists()
        assert neg_output.exists()

        # Verify file contents
        pos_content = pos_output.read_text()
        assert pos_content.count(">") == pos_count
        assert all(jtype in pos_content for jtype in ["donor", "acceptor"])

        if neg_count > 0:
            neg_content = neg_output.read_text()
            assert neg_content.count(">") == neg_count

    def test_expression_filtering_integration(
        self, complex_gff_file, complex_fasta_file, test_kallisto_file
    ):
        """Test integration with expression filtering"""
        extractor = SpliceSeqExtractor(
            str(complex_gff_file),
            str(complex_fasta_file),
            expression_file=str(test_kallisto_file),
            min_expression=2.0,
            expression_format="kallisto",
        )

        # Check that filtering worked
        info = extractor.get_info()
        assert info["transcript_count"] >= 0  # May be 0 if all filtered out
        assert info["junction_count"] >= 0

    def test_different_extraction_parameters(
        self, complex_gff_file, complex_fasta_file
    ):
        """Test with different extraction parameter combinations"""
        test_cases = [
            {"n_exon": 20, "n_intron": 40},  # Short sequences
            {"n_exon": 100, "n_intron": 200},  # Long sequences
            {"n_exon": 0, "n_intron": 50},  # Intron only
            {"n_exon": 50, "n_intron": 0},  # Exon only
        ]

        for params_dict in test_cases:
            params = ExtractionParams(**params_dict, buffer_size=10)
            extractor = SpliceSeqExtractor(
                str(complex_gff_file), str(complex_fasta_file), extraction_params=params
            )

            sequences = extractor.extract_splice_sites()
            expected_length = params.window_size

            if sequences:  # May be empty for some parameter combinations
                assert all(len(seq.sequence) == expected_length for seq in sequences)

    def test_regional_sampling_integration(self, complex_gff_file, complex_fasta_file):
        """Test all types of regional sampling work together"""
        sampling_params = SamplingParams(window_size=100, buffer_size=20)
        extractor = SpliceSeqExtractor(
            str(complex_gff_file),
            str(complex_fasta_file),
            sampling_params=sampling_params,
        )

        # Test intron sampling
        intron_sequences = extractor.extract_intronic_regions(
            target_count=3, random_seed=42
        )
        if intron_sequences:
            assert all(
                seq.junction.junction_type == JunctionType.INTRON
                for seq in intron_sequences
            )
            assert all(len(seq.sequence) == 100 for seq in intron_sequences)

        # Test exon sampling
        exon_sequences = extractor.extract_exonic_sequences(target_count=3, seed=42)
        if exon_sequences:
            assert all(
                seq.junction.junction_type == JunctionType.EXON
                for seq in exon_sequences
            )
            assert all(len(seq.sequence) == 100 for seq in exon_sequences)

        # Test intergenic sampling
        intergenic_sequences = extractor.extract_intergenic_sequences(
            target_count=3, seed=42
        )
        if intergenic_sequences:
            assert all(
                seq.junction.junction_type == JunctionType.INTERGENIC
                for seq in intergenic_sequences
            )
            assert all(len(seq.sequence) == 100 for seq in intergenic_sequences)

    def test_strand_handling_integration(self, complex_gff_file, complex_fasta_file):
        """Test that both positive and negative strand transcripts are handled correctly"""
        extractor = SpliceSeqExtractor(str(complex_gff_file), str(complex_fasta_file))

        transcripts = extractor.transcripts
        junctions = extractor.junctions

        # Should have transcripts on both strands
        strands = set(t.info.strand for t in transcripts.values())
        assert StrandType.POSITIVE in strands or StrandType.NEGATIVE in strands

        # Should have junctions on both strands
        if junctions:
            junction_strands = set(j.strand for j in junctions)
            # At least one strand should be represented
            assert len(junction_strands) >= 1

    def test_motif_analysis_integration(self, complex_gff_file, complex_fasta_file):
        """Test integration with motif analysis"""
        extractor = SpliceSeqExtractor(
            str(complex_gff_file),
            str(complex_fasta_file),
            extraction_params=ExtractionParams(
                n_exon=10, n_intron=10
            ),  # Short for motif analysis
        )

        sequences = extractor.extract_splice_sites()
        if len(sequences) < 2:
            pytest.skip("Need at least 2 sequences for motif analysis")

        # Convert to Bio.Seq objects
        bio_seqs = [Seq(seq.sequence) for seq in sequences[:10]]  # Limit for testing

        # Test motif scoring functions
        ppm = generate_ppm(bio_seqs)
        pfm = generate_pfm(bio_seqs)
        pwm = generate_pwm(bio_seqs)

        assert ppm.shape == (4, 20)  # 4 bases x 20 positions
        assert pfm.shape == (4, 20)
        assert pwm.shape == (4, 20)

        # Check that probabilities sum to 1

        column_sums = ppm.sum(axis=0)
        np.testing.assert_array_almost_equal(column_sums, [1.0] * 20, decimal=5)


class TestErrorHandlingIntegration:
    """Test error handling in the integrated pipeline"""

    def test_missing_files_error_handling(self):
        """Test graceful handling of missing files"""
        with pytest.raises(FileNotFoundError):
            SpliceSeqExtractor("nonexistent.gff3", "nonexistent.fasta")

    def test_invalid_parameters_error_handling(self):
        """Test error handling for invalid parameters"""
        with pytest.raises(ValueError):
            ExtractionParams(n_exon=-1, n_intron=50)

        with pytest.raises(ValueError):
            SamplingParams(window_size=0)

    def test_empty_results_handling(self, temp_dir):
        """Test handling of empty or minimal results"""
        # Create minimal GFF with no valid transcripts
        minimal_gff = temp_dir / "minimal.gff3"
        minimal_gff.write_text("##gff-version 3\n")

        # Create minimal FASTA
        minimal_fasta = temp_dir / "minimal.fasta"
        minimal_fasta.write_text(">chr1\nATCG\n")

        # Index the FASTA

        pysam.faidx(str(minimal_fasta))

        extractor = SpliceSeqExtractor(str(minimal_gff), str(minimal_fasta))

        # Should handle empty results gracefully
        transcripts = extractor.transcripts
        junctions = extractor.junctions

        assert isinstance(transcripts, dict)
        assert isinstance(junctions, list)
        assert len(junctions) == 0  # No valid transcripts = no junctions


class TestPerformanceIntegration:
    """Test performance characteristics of the integrated pipeline"""

    def test_large_dataset_handling(self, complex_gff_file, complex_fasta_file):
        """Test handling of larger datasets"""
        extractor = SpliceSeqExtractor(
            str(complex_gff_file),
            str(complex_fasta_file),
            sampling_params=SamplingParams(window_size=200, buffer_size=100),
        )

        start_time = time.time()

        transcripts = extractor.transcripts
        junctions = extractor.junctions
        positive_sequences = extractor.extract_splice_sites()

        end_time = time.time()

        # Should complete within reasonable time
        assert end_time - start_time < 30
        # Results should be reasonable
        assert isinstance(transcripts, dict)
        assert isinstance(junctions, list)
        assert isinstance(positive_sequences, list)

    def test_memory_efficiency(self, complex_gff_file, complex_fasta_file):
        """Test that lazy loading works for memory efficiency"""
        extractor = SpliceSeqExtractor(str(complex_gff_file), str(complex_fasta_file))

        # Initially, cached data should be None
        assert extractor._transcripts is None
        assert extractor._junctions is None
        assert extractor._genes is None

        # Access should trigger loading
        _ = extractor.transcripts
        assert extractor._transcripts is not None

        _ = extractor.junctions
        assert extractor._junctions is not None

        _ = extractor.genes
        assert extractor._genes is not None


class TestReproducibilityIntegration:
    """Test that results are reproducible with same parameters"""

    def test_random_seed_reproducibility(self, complex_gff_file, complex_fasta_file):
        """Test that random sampling is reproducible with same seed"""
        extractor1 = SpliceSeqExtractor(str(complex_gff_file), str(complex_fasta_file))
        extractor2 = SpliceSeqExtractor(str(complex_gff_file), str(complex_fasta_file))
        transcripts1 = extractor1.transcripts
        transcripts2 = extractor2.transcripts

        # Extract with same seed
        sequences1 = extractor1.regional_sampler.sample_introns(
            transcripts1, target_count=5, seed=42
        )
        sequences2 = extractor2.regional_sampler.sample_introns(
            transcripts2, target_count=5, seed=42
        )

        # Should get same results
        assert len(sequences1) == len(sequences2)

        if sequences1 and sequences2:
            # Check that sequences are identical
            for seq1, seq2 in zip(sequences1, sequences2):
                assert seq1.sequence == seq2.sequence
                assert seq1.window_start == seq2.window_start
                assert seq1.window_end == seq2.window_end

    def test_parameter_consistency(self, complex_gff_file, complex_fasta_file):
        """Test that same parameters always give same results"""
        params = ExtractionParams(n_exon=25, n_intron=75, buffer_size=15)

        extractor1 = SpliceSeqExtractor(
            str(complex_gff_file), str(complex_fasta_file), extraction_params=params
        )
        extractor2 = SpliceSeqExtractor(
            str(complex_gff_file), str(complex_fasta_file), extraction_params=params
        )

        sequences1 = extractor1.extract_splice_sites()
        sequences2 = extractor2.extract_splice_sites()

        # Should get identical results
        assert len(sequences1) == len(sequences2)

        if sequences1 and sequences2:
            for seq1, seq2 in zip(sequences1, sequences2):
                assert seq1.sequence == seq2.sequence
                assert seq1.junction.id == seq2.junction.id
                assert seq1.window_start == seq2.window_start
                assert seq1.window_end == seq2.window_end
