"""
File: test_fasta_writer.py
Description: Updated tests for FASTA writing functionality
"""

from protisplice import write_sequences
from protisplice.data_models import (
    JunctionData,
    SpliceJunction,
    JunctionType,
    StrandType,
)


class TestFastaWriter:
    """Test FASTA writing functionality"""

    def test_write_sequences(self, temp_dir):
        """Test writing sequences to FASTA file"""
        # Create sample junction data
        junction1 = SpliceJunction(
            id="transcript1_donor_0",
            seqid="chr1",
            coord=200,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )

        junction2 = SpliceJunction(
            id="transcript1_acceptor_1",
            seqid="chr1",
            coord=300,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.ACCEPTOR,
        )

        sequences = [
            JunctionData(
                junction=junction1,
                window_start=160,
                window_end=280,
                sequence="ATCGATCGATCG" + "A" * 108,  # 120bp total
            ),
            JunctionData(
                junction=junction2,
                window_start=220,
                window_end=380,
                sequence="GCTAGCTAGCTA" + "T" * 108,  # 120bp total
            ),
        ]

        output_path = temp_dir / "test_output.fasta"
        count = write_sequences(sequences, str(output_path))

        assert count == 2
        assert output_path.exists()

        content = output_path.read_text()
        assert ">chr1_donor_+_160_280" in content
        assert "ATCGATCGATCG" in content
        assert ">chr1_acceptor_+_220_380" in content
        assert "GCTAGCTAGCTA" in content

        # Check that sequences are properly formatted
        lines = content.strip().split("\n")
        assert len(lines) == 4  # 2 headers + 2 sequences
        assert lines[0].startswith(">")
        assert lines[2].startswith(">")
        assert not lines[1].startswith(">")
        assert not lines[3].startswith(">")

    def test_write_sequences_with_none_sequence(self, temp_dir):
        """Test writing sequences where some have None sequence"""
        junction = SpliceJunction(
            id="test_junction",
            seqid="chr1",
            coord=100,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )

        sequences = [
            JunctionData(
                junction=junction,
                window_start=50,
                window_end=150,
                sequence="ATCGATCG",
            ),
            JunctionData(
                junction=junction,
                window_start=200,
                window_end=300,
                sequence=None,  # This should be skipped
            ),
            JunctionData(
                junction=junction,
                window_start=350,
                window_end=450,
                sequence="GCTAGCTA",
            ),
        ]

        output_path = temp_dir / "test_output_with_none.fasta"
        count = write_sequences(sequences, str(output_path))

        # Should only write 2 sequences (skipping the None)
        assert count == 2
        assert output_path.exists()

        content = output_path.read_text()
        assert content.count(">") == 2
        assert "ATCGATCG" in content
        assert "GCTAGCTA" in content

    def test_write_sequences_empty_list(self, temp_dir):
        """Test writing empty sequence list"""
        output_path = temp_dir / "empty_output.fasta"
        count = write_sequences([], str(output_path))

        assert count == 0
        assert output_path.exists()
        assert output_path.read_text() == ""

    def test_write_sequences_different_junction_types(self, temp_dir):
        """Test writing sequences with different junction types"""
        junctions = [
            SpliceJunction(
                id="donor_test",
                seqid="chr1",
                coord=100,
                strand=StrandType.POSITIVE,
                junction_type=JunctionType.DONOR,
            ),
            SpliceJunction(
                id="acceptor_test",
                seqid="chr1",
                coord=200,
                strand=StrandType.NEGATIVE,
                junction_type=JunctionType.ACCEPTOR,
            ),
            SpliceJunction(
                id="intron_test",
                seqid="chr2",
                coord=300,
                strand=StrandType.POSITIVE,
                junction_type=JunctionType.INTRON,
            ),
            SpliceJunction(
                id="exon_test",
                seqid="chr2",
                coord=400,
                strand=StrandType.NEGATIVE,
                junction_type=JunctionType.EXON,
            ),
            SpliceJunction(
                id="intergenic_test",
                seqid="chr3",
                coord=500,
                strand=StrandType.POSITIVE,
                junction_type=JunctionType.INTERGENIC,
            ),
        ]

        sequences = [
            JunctionData(
                junction=junction,
                window_start=junction.coord - 50,
                window_end=junction.coord + 50,
                sequence="A" * 100,
            )
            for junction in junctions
        ]

        output_path = temp_dir / "different_types.fasta"
        count = write_sequences(sequences, str(output_path))

        assert count == 5
        content = output_path.read_text()

        # Check that all junction types are present
        assert "donor" in content
        assert "acceptor" in content
        assert "intron" in content
        assert "exon" in content
        assert "intergenic" in content

        # Check that strands are present
        assert "_+_" in content
        assert "_-_" in content

        # Check that different chromosomes are present
        assert "chr1" in content
        assert "chr2" in content
        assert "chr3" in content

    def test_fasta_header_format(self, temp_dir):
        """Test that FASTA headers are correctly formatted"""
        junction = SpliceJunction(
            id="test_junction",
            seqid="chr10",
            coord=12345,
            strand=StrandType.NEGATIVE,
            junction_type=JunctionType.ACCEPTOR,
        )

        sequence = JunctionData(
            junction=junction,
            window_start=12300,
            window_end=12400,
            sequence="ATCG" * 25,  # 100bp
        )

        output_path = temp_dir / "header_test.fasta"
        count = write_sequences([sequence], str(output_path))

        assert count == 1
        content = output_path.read_text()

        # Header format should be: >seqid_junctiontype_strand_start_end
        expected_header = ">chr10_acceptor_-_12300_12400"
        assert expected_header in content

    def test_write_sequences_long_sequences(self, temp_dir):
        """Test writing very long sequences"""
        junction = SpliceJunction(
            id="long_test",
            seqid="chr1",
            coord=1000,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )

        # Create a long sequence (1000bp)
        long_sequence = "ATCG" * 250

        sequence = JunctionData(
            junction=junction, window_start=500, window_end=1500, sequence=long_sequence
        )

        output_path = temp_dir / "long_sequence.fasta"
        count = write_sequences([sequence], str(output_path))

        assert count == 1
        content = output_path.read_text()

        # Should contain the full long sequence
        assert long_sequence in content
        assert len(content.split("\n")[1]) == 1000  # Sequence line should be 1000bp

    def test_write_sequences_special_characters_in_path(self, temp_dir):
        """Test writing to paths with special characters"""
        junction = SpliceJunction(
            id="special_test",
            seqid="chr1",
            coord=100,
            strand=StrandType.POSITIVE,
            junction_type=JunctionType.DONOR,
        )

        sequence = JunctionData(
            junction=junction, window_start=50, window_end=150, sequence="ATCGATCG"
        )

        # Path with spaces and special characters
        output_path = temp_dir / "test output file.fasta"
        count = write_sequences([sequence], str(output_path))

        assert count == 1
        assert output_path.exists()

        content = output_path.read_text()
        assert "ATCGATCG" in content
