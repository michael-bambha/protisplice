"""
File: test_feature_eng.py
Description: Unit tests for feature engineering module
"""
# pylint:disable=redefined-outer-name

import pytest
from protisplice import (
    JunctionData,
    SpliceJunction,
    JunctionType,
    StrandType,
    ExtractionParams,
)
from protisplice.feature_eng import (
    inject_consensus,
    remove_consensus,
    _get_consensus_index,
    _motif_for,
    _pick_noncanon_dinuc,
)


@pytest.fixture
def extraction_params():
    """Standard extraction parameters for testing"""
    return ExtractionParams(n_exon=40, n_intron=80, buffer_size=50)


class TestInjectConsensus:
    """Test consensus injection functionality"""

    def test_inject_donor_consensus(self, extraction_params):
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
            window_end=130,  # 120bp window
            sequence="A" * 40 + "NN" + "T" * 78,  # 120bp total
        )

        result = inject_consensus(junction_data, extraction_params, JunctionType.DONOR)
        expected = "A" * 40 + "GT" + "T" * 78
        assert result.sequence == expected
        assert isinstance(result, JunctionData)

    def test_inject_acceptor_consensus(self, extraction_params):
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
            window_end=130,  # 120bp window
            sequence="A" * 38 + "NN" + "T" * 80,  # 120bp total
        )

        result = inject_consensus(
            junction_data, extraction_params, JunctionType.ACCEPTOR
        )
        expected = "A" * 38 + "AG" + "T" * 80
        assert result.sequence == expected
        assert isinstance(result, JunctionData)

    def test_inject_with_custom_index(self, extraction_params):
        """Test consensus injection with custom index"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction, window_start=10, window_end=130, sequence="A" * 120
        )

        result = inject_consensus(
            junction_data, extraction_params, JunctionType.DONOR, idx=20
        )
        expected = "A" * 20 + "GT" + "A" * 98
        assert result.sequence == expected
        assert isinstance(result, JunctionData)

    def test_inject_out_of_bounds_raises_error(self, extraction_params):
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
            window_end=14,
            sequence="ATCG",  # Only 4bp
        )

        with pytest.raises(IndexError, match="Consensus start idx .* out of bounds"):
            inject_consensus(
                junction_data, extraction_params, JunctionType.DONOR, idx=10
            )


class TestRemoveConsensus:
    """Test consensus removal functionality"""

    def test_remove_donor_consensus_when_present(self, extraction_params):
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
            window_end=130,
            sequence="A" * 40 + "GT" + "T" * 78,
        )

        result = remove_consensus(junction_data, extraction_params)
        # Should replace GT with some non-GT dinucleotide
        assert result.sequence[:40] == "A" * 40
        assert result.sequence[42:] == "T" * 78
        assert result.sequence[40:42] != "GT"
        assert len(result.sequence[40:42]) == 2
        assert isinstance(result, JunctionData)

    def test_remove_acceptor_consensus_when_present(self, extraction_params):
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
            window_end=130,
            sequence="A" * 38 + "AG" + "T" * 80,
        )

        result = remove_consensus(junction_data, extraction_params)
        assert result.sequence[:38] == "A" * 38
        assert result.sequence[40:] == "T" * 80
        assert result.sequence[38:40] != "AG"
        assert isinstance(result, JunctionData)

    def test_remove_consensus_only_if_present_false(self, extraction_params):
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
            window_end=130,
            sequence="A" * 40 + "TT" + "T" * 78,  # No GT present
        )

        result = remove_consensus(
            junction_data, extraction_params, only_if_present=False
        )
        # Should still replace TT with non-GT dinucleotide
        assert result.sequence[:40] == "A" * 40
        assert result.sequence[42:] == "T" * 78
        assert result.sequence[40:42] != "GT"
        assert result.sequence[40:42] != "TT"  # Should be replaced
        assert isinstance(result, JunctionData)

    def test_remove_consensus_only_if_present_true_no_change(self, extraction_params):
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
            window_end=130,
            sequence="A" * 40 + "TT" + "T" * 78,  # No GT present
        )
        original_sequence = junction_data.sequence

        result = remove_consensus(
            junction_data, extraction_params, only_if_present=True
        )
        # Should remain unchanged
        assert result.sequence == original_sequence
        assert isinstance(result, JunctionData)

    def test_remove_consensus_custom_replacement(self, extraction_params):
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
            window_end=130,
            sequence="A" * 40 + "GT" + "T" * 78,
        )

        result = remove_consensus(junction_data, extraction_params, replacement="CC")
        expected = "A" * 40 + "CC" + "T" * 78
        assert result.sequence == expected
        assert isinstance(result, JunctionData)

    def test_invalid_replacement_raises_error(self, extraction_params):
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
            window_end=130,
            sequence="A" * 40 + "GT" + "T" * 78,
        )

        with pytest.raises(ValueError, match="replacement must be a 2-mer"):
            remove_consensus(junction_data, extraction_params, replacement="GT")

        with pytest.raises(ValueError, match="replacement must be a 2-mer"):
            remove_consensus(junction_data, extraction_params, replacement="A")


class TestHelperFunctions:
    """Test helper functions"""

    def test_get_consensus_index_donor(self, extraction_params):
        """Test getting consensus index for donor site"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction, window_start=10, window_end=130, sequence="A" * 120
        )

        idx = _get_consensus_index(junction_data, extraction_params, JunctionType.DONOR)
        expected = extraction_params.n_exon  # Should be 40 for donor
        assert idx == expected

    def test_get_consensus_index_acceptor(self, extraction_params):
        """Test getting consensus index for acceptor site"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.ACCEPTOR,
        )
        junction_data = JunctionData(
            junction=junction, window_start=10, window_end=130, sequence="A" * 120
        )

        idx = _get_consensus_index(
            junction_data, extraction_params, JunctionType.ACCEPTOR
        )
        expected = extraction_params.n_exon - 2  # Should be 38 for acceptor
        assert idx == expected

    def test_get_consensus_index_invalid_type_raises_error(self, extraction_params):
        """Test that invalid junction type raises ValueError"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.INTRON,  # Invalid for consensus
        )
        junction_data = JunctionData(
            junction=junction, window_start=10, window_end=130, sequence="A" * 120
        )

        with pytest.raises(
            ValueError, match="Junction type must be either donor or acceptor"
        ):
            _get_consensus_index(junction_data, extraction_params, JunctionType.INTRON)

    def test_get_consensus_index_uses_junction_type_when_none_provided(
        self, extraction_params
    ):
        """Test that function uses junction's type when no junc_type parameter provided"""
        junction = SpliceJunction(
            id="test",
            seqid="chr1",
            coord=50,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )
        junction_data = JunctionData(
            junction=junction, window_start=10, window_end=130, sequence="A" * 120
        )

        idx = _get_consensus_index(
            junction_data, extraction_params
        )  # No junc_type provided
        expected = extraction_params.n_exon  # Should use DONOR from junction
        assert idx == expected

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
