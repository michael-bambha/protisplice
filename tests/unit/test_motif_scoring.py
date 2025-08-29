"""
File: test_motif_scoring.py
Description: Unit tests for motif scoring module
"""

import pytest
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from protisplice.motif_scoring import generate_ppm, generate_pfm, generate_pwm


class TestMotifScoring:
    """Test motif scoring functions"""

    @pytest.fixture
    def sample_sequences(self):
        """Sample aligned sequences for testing"""
        return [Seq("ATCG"), Seq("ATGG"), Seq("ACCG"), Seq("ATCG")]

    @pytest.fixture
    def donor_sequences(self):
        """Sample donor site sequences"""
        return [
            Seq("AAGTAAGT"),  # GT donor
            Seq("AGGTAAGT"),  # GT donor
            Seq("ACGTAAGT"),  # GT donor
            Seq("ATGTAAGT"),  # GT donor
        ]

    def test_generate_ppm_basic(self, sample_sequences):
        """Test basic PPM generation"""
        ppm = generate_ppm(sample_sequences)

        assert isinstance(ppm, pd.DataFrame)
        assert ppm.shape == (4, 4)  # 4 nucleotides x 4 positions
        assert list(ppm.index) == ["A", "C", "G", "T"]

        # Check that probabilities sum to 1 for each position
        column_sums = ppm.sum(axis=0)
        np.testing.assert_array_almost_equal(column_sums, [1.0] * 4, decimal=5)

        # Check first position (all A's)
        assert ppm.loc["A", 0] == 1.0
        assert ppm.loc["C", 0] == 0.0
        assert ppm.loc["G", 0] == 0.0
        assert ppm.loc["T", 0] == 0.0

    def test_generate_pfm_basic(self, sample_sequences):
        """Test basic PFM generation"""
        pfm = generate_pfm(sample_sequences)

        assert isinstance(pfm, pd.DataFrame)
        assert pfm.shape == (4, 4)

        # Check that counts are integers
        assert all(pfm.dtypes == "int64")

        # Check that column sums equal number of sequences
        column_sums = pfm.sum(axis=0)
        assert all(column_sums == 4)  # 4 sequences

        # Check first position counts (all A's) - use integer indices
        assert pfm.loc["A", 0] == 4
        assert pfm.loc["C", 0] == 0
        assert pfm.loc["G", 0] == 0
        assert pfm.loc["T", 0] == 0

    def test_generate_pwm_default_background(self, sample_sequences):
        """Test PWM generation with default background frequencies"""
        pwm = generate_pwm(sample_sequences)

        assert isinstance(pwm, pd.DataFrame)
        assert pwm.shape == (4, 4)

        # Check that values are log2 ratios
        # For first position (all A's): log2(1.0/0.25) = 2.0
        expected_first_col = np.log2(1.0 / 0.25)
        assert abs(pwm.loc["A", 0] - expected_first_col) < 1e-10

        # Other nucleotides should be very negative (due to epsilon)
        assert pwm.loc["C", 0] < -10
        assert pwm.loc["G", 0] < -10
        assert pwm.loc["T", 0] < -10

    def test_generate_pwm_custom_background(self, sample_sequences):
        """Test PWM generation with custom background frequencies"""
        background = {"A": 0.4, "C": 0.1, "G": 0.1, "T": 0.4}
        pwm = generate_pwm(sample_sequences, background_freq=background)

        assert isinstance(pwm, pd.DataFrame)
        # For first position (all A's): log2(1.0/0.4)
        expected_first_col = np.log2(1.0 / 0.4)
        # Use a more reasonable tolerance for floating point comparison
        assert abs(pwm.loc["A", 0] - expected_first_col) < 1e-6

    def test_donor_site_motif(self, donor_sequences):
        """Test motif generation with donor site sequences"""
        ppm = generate_ppm(donor_sequences)
        pfm = generate_pfm(donor_sequences)
        pwm = generate_pwm(donor_sequences)

        # All sequences should have G at position 2 and T at position 3
        assert ppm.loc["G", 2] == 1.0  # 100% G at donor position
        assert ppm.loc["T", 3] == 1.0  # 100% T at donor position

        assert pfm.loc["G", 2] == 4  # Count of 4 G's
        assert pfm.loc["T", 3] == 4  # Count of 4 T's

        # PWM should show high scores for G and T at donor positions
        assert pwm.loc["G", 2] > 1.5  # High positive score
        assert pwm.loc["T", 3] > 1.5  # High positive score

    def test_empty_sequences_raises_error(self):
        """Test that empty sequence list raises appropriate error"""
        with pytest.raises((ValueError, IndexError)):
            generate_ppm([])

    def test_unequal_length_sequences(self):
        """Test with sequences of unequal length"""
        sequences = [Seq("ATCG"), Seq("AT")]  # Different lengths

        # Should raise an error or handle gracefully
        with pytest.raises(Exception):
            generate_ppm(sequences)

    def test_single_sequence(self):
        """Test with single sequence"""
        sequences = [Seq("ATCG")]
        ppm = generate_ppm(sequences)

        assert ppm.shape == (4, 4)
        # All positions should have probability 1.0 for the observed base
        assert ppm.loc["A", 0] == 1.0
        assert ppm.loc["T", 1] == 1.0
        assert ppm.loc["C", 2] == 1.0
        assert ppm.loc["G", 3] == 1.0

    def test_all_same_sequences(self):
        """Test with identical sequences"""
        sequences = [Seq("AAAA")] * 5
        ppm = generate_ppm(sequences)
        pfm = generate_pfm(sequences)

        # All positions should be 100% A
        for pos in range(4):
            assert ppm.loc["A", pos] == 1.0
            assert pfm.loc["A", pos] == 5
            for base in ["C", "G", "T"]:
                assert ppm.loc[base, pos] == 0.0
                assert pfm.loc[base, pos] == 0

    def test_background_freq_validation(self, sample_sequences):
        """Test that background frequencies are properly handled"""
        # Background frequencies that don't sum to 1
        background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.2}  # Sums to 0.9

        # Should still work (pandas will normalize or handle appropriately)
        pwm = generate_pwm(sample_sequences, background_freq=background)
        assert isinstance(pwm, pd.DataFrame)

    def test_pwm_matrix_orientation(self, sample_sequences):
        """Test that PWM matrix has correct orientation"""
        pwm = generate_pwm(sample_sequences)

        # Rows should be nucleotides, columns should be positions
        assert len(pwm.index) == 4  # A, C, G, T
        assert len(pwm.columns) == 4  # 4 positions in sequence

        # Check that we can access specific nucleotide-position combinations
        assert isinstance(pwm.loc["A", 0], (float, np.float64))
        assert isinstance(pwm.loc["G", 2], (float, np.float64))
