"""
File: test_fasta_writer.py
Description: Test FASTA writing functionality
"""

from protisplice import FastaWriter
from protisplice import ExtractionResults, JunctionData


class TestFastaWriter:
    """Test FastaWriter class"""

    def test_write_sequences(self, temp_dir, sample_junctions):
        """Test writing sequences to FASTA file"""

        sequences = [
            JunctionData(
                junction=sample_junctions[0],
                window_start=160,
                window_end=280,
                sequence="ATCGATCGATCG",
            ),
            JunctionData(
                junction=sample_junctions[1],
                window_start=220,
                window_end=380,
                sequence="GCTAGCTAGCTA",
            ),
        ]

        output_path = temp_dir / "test_output.fasta"
        count = FastaWriter.write_sequences(sequences, str(output_path))

        assert count == 2
        assert output_path.exists()

        content = output_path.read_text()
        assert ">chr1_donor_+_160_280" in content
        assert "ATCGATCGATCG" in content
        assert ">chr1_acceptor_+_220_380" in content
        assert "GCTAGCTAGCTA" in content

    def test_write_results(self, temp_dir, sample_junctions):
        """Test writing extraction results"""

        positive_seq = JunctionData(
            junction=sample_junctions[0],
            window_start=160,
            window_end=280,
            sequence="ATCGATCG",
        )

        negative_seq = JunctionData(
            junction=sample_junctions[1],
            window_start=220,
            window_end=380,
            sequence="GCTAGCTA",
        )

        results = ExtractionResults(
            positive_sequences=[positive_seq], negative_sequences=[negative_seq]
        )

        pos_path = temp_dir / "positive.fasta"
        neg_path = temp_dir / "negative.fasta"

        pos_count, neg_count = FastaWriter.write_results(
            results, str(pos_path), str(neg_path)
        )

        assert pos_count == 1
        assert neg_count == 1
        assert pos_path.exists()
        assert neg_path.exists()
