"""
File: test_feature_eng.py
Description: Unit tests for feature engineering module
"""

import pytest
from protisplice import JunctionData, SpliceJunction, JunctionType, StrandType
from protisplice.feature_eng import (
    inject_consensus,
    remove_consensus,
    _get_consensus_index,
    _motif_for,
    _pick_noncanon_dinuc,
)


class TestInjectConsensus:
    """Test consensus injection functionality"""

    def test_inject_donor_consensus(self):
        """Test injecting donor consensus (GT)"""
        junction = SpliceJunction(
            id="test_donor",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction,
            window_start=10,
            window_end=90,
            sequence="A" * 40 + "NN" + "T" * 38,  # 80bp total
        )

        result = inject_consensus(junction_data, JunctionType.DONOR)
        expected = "A" * 40 + "GT" + "T" * 38
        assert result == expected

    def test_inject_acceptor_consensus(self):
        """Test injecting acceptor consensus (AG)"""
        junction = SpliceJunction(
            id="test_acceptor",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.ACCEPTOR,
        )
        junction_data = JunctionData(
            junction=junction,
            window_start=10,
            window_end=90,
            sequence="A" * 38 + "NN" + "T" * 40,  # 80bp total
        )

        result = inject_consensus(junction_data, JunctionType.ACCEPTOR)
        expected = "A" * 38 + "AG" + "T" * 40
        assert result == expected

    def test_inject_with_custom_index(self):
        """Test consensus injection with custom index"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction, window_start=10, window_end=90, sequence="A" * 80
        )

        result = inject_consensus(junction_data, JunctionType.DONOR, idx=20)
        expected = "A" * 20 + "GT" + "A" * 58
        assert result == expected

    def test_inject_out_of_bounds_raises_error(self):
        """Test that out of bounds index raises IndexError"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction,
            window_start=10,
            window_end=90,
            sequence="ATCG",  # Only 4bp
        )

        with pytest.raises(IndexError, match="Consensus start idx .* out of bounds"):
            inject_consensus(junction_data, JunctionType.DONOR, idx=10)


class TestRemoveConsensus:
    """Test consensus removal functionality"""

    def test_remove_donor_consensus_when_present(self):
        """Test removing donor consensus when GT is present"""
        junction = SpliceJunction(
            id="test_donor",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction,
            window_start=10,
            window_end=90,
            sequence="A" * 40 + "GT" + "T" * 38,
        )

        result = remove_consensus(junction_data, JunctionType.DONOR)
        # Should replace GT with some non-GT dinucleotide
        assert result[:40] == "A" * 40
        assert result[42:] == "T" * 38
        assert result[40:42] != "GT"
        assert len(result[40:42]) == 2

    def test_remove_acceptor_consensus_when_present(self):
        """Test removing acceptor consensus when AG is present"""
        junction = SpliceJunction(
            id="test_acceptor",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.ACCEPTOR,
        )
        junction_data = JunctionData(
            junction=junction,
            window_start=10,
            window_end=90,
            sequence="A" * 38 + "AG" + "T" * 40,
        )

        result = remove_consensus(junction_data, JunctionType.ACCEPTOR)
        assert result[:38] == "A" * 38
        assert result[40:] == "T" * 40
        assert result[38:40] != "AG"

    def test_remove_consensus_only_if_present_false(self):
        """Test removing consensus when only_if_present=False"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction,
            window_start=10,
            window_end=90,
            sequence="A" * 40 + "TT" + "T" * 38,  # No GT present
        )

        result = remove_consensus(
            junction_data, JunctionType.DONOR, only_if_present=False
        )
        # Should still replace TT with non-GT dinucleotide
        assert result[:40] == "A" * 40
        assert result[42:] == "T" * 38
        assert result[40:42] != "GT"
        assert result[40:42] != "TT"  # Should be replaced

    def test_remove_consensus_only_if_present_true_no_change(self):
        """Test no change when consensus not present and only_if_present=True"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction,
            window_start=10,
            window_end=90,
            sequence="A" * 40 + "TT" + "T" * 38,  # No GT present
        )

        result = remove_consensus(
            junction_data, JunctionType.DONOR, only_if_present=True
        )
        # Should remain unchanged
        assert result == junction_data.sequence

    def test_remove_consensus_custom_replacement(self):
        """Test consensus removal with custom replacement"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction,
            window_start=10,
            window_end=90,
            sequence="A" * 40 + "GT" + "T" * 38,
        )

        result = remove_consensus(junction_data, JunctionType.DONOR, replacement="CC")
        expected = "A" * 40 + "CC" + "T" * 38
        assert result == expected

    def test_invalid_replacement_raises_error(self):
        """Test that invalid replacement raises ValueError"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction,
            window_start=10,
            window_end=90,
            sequence="A" * 40 + "GT" + "T" * 38,
        )

        with pytest.raises(ValueError, match="replacement must be a 2-mer"):
            remove_consensus(junction_data, JunctionType.DONOR, replacement="GT")

        with pytest.raises(ValueError, match="replacement must be a 2-mer"):
            remove_consensus(junction_data, JunctionType.DONOR, replacement="A")


class TestHelperFunctions:
    """Test helper functions"""

    def test_get_consensus_index_donor(self):
        """Test getting consensus index for donor site"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction, window_start=10, window_end=90, sequence="A" * 80
        )

        idx = _get_consensus_index(junction_data, JunctionType.DONOR)
        expected = junction.coord - junction_data.window_start  # 50 - 10 = 40
        assert idx == expected

    def test_get_consensus_index_acceptor(self):
        """Test getting consensus index for acceptor site"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.ACCEPTOR,
        )
        junction_data = JunctionData(
            junction=junction, window_start=10, window_end=90, sequence="A" * 80
        )

        idx = _get_consensus_index(junction_data, JunctionType.ACCEPTOR)
        expected = junction.coord - junction_data.window_start - 2  # 50 - 10 - 2 = 38
        assert idx == expected

    def test_get_consensus_index_invalid_type_raises_error(self):
        """Test that invalid junction type raises ValueError"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.INTRON,  # Invalid for consensus
        )
        junction_data = JunctionData(
            junction=junction, window_start=10, window_end=90, sequence="A" * 80
        )

        with pytest.raises(
            ValueError, match="Junction type must be either donor or acceptor"
        ):
            _get_consensus_index(junction_data, JunctionType.INTRON)

    def test_motif_for_donor(self):
        """Test getting motif for donor"""
        assert _motif_for(JunctionType.DONOR) == "GT"

    def test_motif_for_acceptor(self):
        """Test getting motif for acceptor"""
        assert _motif_for(JunctionType.ACCEPTOR) == "AG"

    def test_pick_noncanon_dinuc(self):
        """Test picking non-canonical dinucleotide"""
        result = _pick_noncanon_dinuc("GT")
        assert len(result) == 2
        assert result != "GT"
        assert all(base in "ACGT" for base in result)

        result = _pick_noncanon_dinuc("AG")
        assert len(result) == 2
        assert result != "AG"
        assert all(base in "ACGT" for base in result)
