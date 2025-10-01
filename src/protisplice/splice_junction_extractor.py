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
        """Takes in a dictionary of transcripts and returns a list
        of SpliceJunction objects, which contains the id of the parent sequence,
        the transcript ID, 1-based coordinate of the first exon base, and
        the strand.

        Args:
            transcripts (Dict[str, Transcript]): Dict of transcript ID: Transcript object.
            Transcript objects contain a list of tuples of exon start/end coordinates,
            list of tuples of intron start/end coordinates, and a dictionary containing
            the transcript ID and strand.


        Returns:
            List[SpliceJunction]: List of SpliceJunction objects, which contains the
            id of the parent sequence, the transcript ID, 1-based coordinate of the first
            exon base, and the strand.
        """
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
        """Create SpliceJunction objects for a single transcript

        Args:
            transcript_id (str): ID of the transcript
            transcript (Transcript): Transcript object
            sorted_exons (List[tuple]): Sorted list of exon start, end coords.
            (1-based)

        Returns:
            List[SpliceJunction]: List of SpliceJunction objects for one transcript, which
            contains the id of the parent sequence, the transcript ID, 1-based coordinate
            of the first exon base, and the strand.
        """
        junctions = []

        for i, (exon_start, exon_end) in enumerate(sorted_exons):
            if i > 0:
                if transcript.info.strand == StrandType.POSITIVE:
                    junction_type = JunctionType.ACCEPTOR
                else:
                    junction_type = JunctionType.DONOR

                junctions.append(
                    SpliceJunction(
                        id=f"{transcript_id}_{junction_type.value}_{i}",
                        seqid=transcript.info.seqid,
                        coord=exon_start,
                        strand=transcript.info.strand,
                        junction_type=junction_type,
                    )
                )

            if i < len(sorted_exons) - 1:
                if transcript.info.strand == StrandType.POSITIVE:
                    junction_type = JunctionType.DONOR
                else:
                    junction_type = JunctionType.ACCEPTOR

                junctions.append(
                    SpliceJunction(
                        id=f"{transcript_id}_{junction_type.value}_{i}",
                        seqid=transcript.info.seqid,
                        coord=exon_end,
                        strand=transcript.info.strand,
                        junction_type=junction_type,
                    )
                )

        return junctions
