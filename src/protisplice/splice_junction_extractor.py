"""
File: splice_junction_extractor.py
Description: Splice junction extraction logic
"""

from typing import Dict, List

from .data_models import SpliceJunction, JunctionType, StrandType, Transcript


class SpliceJunctionExtractor:
    """Extracts splice junctions from transcript data"""

    def identify_splice_junctions(
        self, transcripts: Dict[str, Transcript]
    ) -> List[SpliceJunction]:
        """Get splice junction coordinates from transcripts"""
        junctions = []

        for transcript_id, transcript in transcripts.items():
            if len(transcript.exons) < 2:
                continue

            sorted_exons = sorted(transcript.exons, key=lambda x: x[0])
            junctions.extend(
                self._create_junctions_for_transcript(
                    transcript_id, transcript, sorted_exons
                )
            )

        return junctions

    def _create_junctions_for_transcript(
        self, transcript_id: str, transcript: Transcript, sorted_exons: List[tuple]
    ) -> List[SpliceJunction]:
        """Create junction objects for a single transcript"""
        junctions = []

        for i, (exon_start, exon_end) in enumerate(sorted_exons):
            # Acceptor sites (except for first exon)
            if i > 0:
                coord = (
                    exon_start
                    if transcript.info.strand == StrandType.POSITIVE
                    else exon_end
                )
                junctions.append(
                    SpliceJunction(
                        id=f"{transcript_id}_acceptor_{i}",
                        seqid=transcript.info.seqid,
                        coord=coord,
                        strand=transcript.info.strand,
                        junction_type=JunctionType.ACCEPTOR,
                    )
                )

            # Donor sites (except for last exon)
            if i < len(sorted_exons) - 1:
                coord = (
                    exon_end
                    if transcript.info.strand == StrandType.POSITIVE
                    else exon_start
                )
                junctions.append(
                    SpliceJunction(
                        id=f"{transcript_id}_donor_{i}",
                        seqid=transcript.info.seqid,
                        coord=coord,
                        strand=transcript.info.strand,
                        junction_type=JunctionType.DONOR,
                    )
                )

        return junctions
