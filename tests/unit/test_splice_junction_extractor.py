"""
File: test_splice_junction_extractor.py
Description: Test splice junction extraction logic
"""

from splice_seq_extractor.splice_junction_extractor import SpliceJunctionExtractor
from splice_seq_extractor.data_models import (
    JunctionType,
    StrandType,
    Transcript,
    TranscriptInfo,
)


class TestSpliceJunctionExtractor:
    """Tests for SpliceJunctionExtractor class"""

    def test_identify_splice_junctions_pos(self, sample_transcript):
        """Test junction identification for transcripts on positive strand"""
        extractor = SpliceJunctionExtractor()
        transcripts = {"transcript1": sample_transcript}
        junctions = extractor.identify_splice_junctions(transcripts)
        assert len(junctions) == 4

        junction_types = [j.junction_type for j in junctions]
        assert junction_types.count(JunctionType.ACCEPTOR) == 2
        assert junction_types.count(JunctionType.DONOR) == 2

        acceptor_coords = [
            j.coord for j in junctions if j.junction_type == JunctionType.ACCEPTOR
        ]
        donor_coords = [
            j.coord for j in junctions if j.junction_type == JunctionType.DONOR
        ]

        assert sorted(acceptor_coords) == [300, 500]
        assert sorted(donor_coords) == [200, 400]

    def test_single_exon_transcript(self):
        """Test that transcripts with one exon return no junctions"""
        transcript = Transcript(
            info=TranscriptInfo(seqid="chr1", strand=StrandType.POSITIVE),
            exons=[(100, 200)],
        )
        extractor = SpliceJunctionExtractor()
        transcripts = {"transcript1": transcript}
        junctions = extractor.identify_splice_junctions(transcripts)
        assert len(junctions) == 0
