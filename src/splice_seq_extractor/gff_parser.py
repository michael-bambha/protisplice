"""
File: gff_parser.py
Description: GFF3 parsing functionality
"""

from collections import defaultdict
from pathlib import Path
from typing import Dict, Tuple, Optional
from .data_models import Transcript, TranscriptInfo, StrandType, TranscriptFilter
from .expression import ExpressionFilter


class GFFParser:
    """
    Class for parsing GFF3 files.
    """

    def __init__(
        self,
        gff_path: str,
        transcript_filter: TranscriptFilter.ALL,
        expression_filter: Optional[ExpressionFilter] = None,
    ):
        self.gff_path = Path(gff_path)
        self.transcript_filter = transcript_filter
        self.expression_filter = expression_filter

        if not self.gff_path.exists():
            raise FileNotFoundError(
                f"GFF3 file {gff_path} not found -- double check the path!"
            )

    def parse_transcripts(self) -> Dict[str, Transcript]:
        """Parse GFF3 file to extract transcript information"""
        transcript_biotypes = {}
        if self.transcript_filter != TranscriptFilter.ALL:
            transcript_biotypes = self._collect_transcript_biotypes()
        transcripts = defaultdict(
            lambda: Transcript(
                info=TranscriptInfo(seqid="", strand=StrandType.POSITIVE), exons=[]
            )
        )

        with open(self.gff_path, "r", encoding="utf-8") as file:
            for _, line in enumerate(file, 1):
                transcript_data = self._parse_gff_line(line, transcript_biotypes)
                if transcript_data:
                    transcript_id, seqid, start, end, strand = transcript_data
                    transcript = transcripts[transcript_id]
                    transcript.exons.append((start, end))

                    if not transcript.info.seqid:
                        transcript.info.seqid = seqid
                        transcript.info.strand = StrandType(strand)

        return {
            tid: transcript
            for tid, transcript in transcripts.items()
            if transcript.exons
        }

    def _parse_gff_line(
        self, line: str, transcript_biotypes: Dict[str, str]
    ) -> Optional[Tuple]:
        """Parse GFF3 line for exons and return transcript data

        Args:
            line (str): Line to parse
            transcript_biotypes (Dict[str, str]): Dict of {ID: biotype}

        Returns:
            Optional[Tuple]: _description_
        """
        try:
            if line.startswith("#") or not line.strip():
                return None

            fields = line.strip().split("\t")

            if len(fields) < 9 or fields[2] != "exon":
                return None

            attrs = self._extract_gff_attrs(fields[8])
            parent = attrs.get("Parent", "")
            transcript_id = None

            if parent.startswith("ENST") or parent.startswith("ENSMU"):
                transcript_id = parent
            elif parent and not parent.startswith("gene:"):
                transcript_id = parent

            if not transcript_id:
                return None

            transcript_id = ExpressionFilter.normalize_transcript_id(transcript_id)

            if not self._passes_filter(transcript_id, transcript_biotypes):
                return None

            seqid, start, end, strand = (
                fields[0],
                int(fields[3]),
                int(fields[4]),
                fields[6],
            )

            return transcript_id, seqid, start, end, strand

        except (ValueError, IndexError):
            return None

    def _collect_transcript_biotypes(self) -> Dict[str, str]:
        """Get the annotation biotype associated with each transcript

        Returns:
            Dict[str, str]: Dict of {transcript_id: biotype}
        """
        transcript_biotypes = {}

        with open(self.gff_path, "r", encoding="utf-8") as file:
            for line in file:
                if line.startswith("#") or not line.strip():
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue

                feature_type = fields[2]
                if feature_type not in [
                    "mRNA",
                    "transcript",
                    "lncRNA",
                    "miRNA",
                    "ncRNA",
                    "rRNA",
                    "snoRNA",
                    "snRNA",
                    "tRNA",
                ]:
                    continue

                attrs = self._extract_gff_attrs(fields[8])
                transcript_id = attrs.get("ID", "")
                transcript_id = ExpressionFilter.normalize_transcript_id(transcript_id)

                biotype = None
                for attr_name in [
                    "biotype",
                    "gene_type",
                    "transcript_type",
                    "gene_biotype",
                    "transcript_biotype",
                ]:
                    if attr_name in attrs:
                        biotype = attrs[attr_name]
                        break

                if biotype:
                    transcript_biotypes[transcript_id] = biotype

        return transcript_biotypes

    def _extract_gff_attrs(self, gff_attributes: str) -> Dict[str, str]:
        """Parse the 9th field of a GFF3 formatted line

        Args:
            gff_attributes (str): 9th field of GFF3 line

        Returns:
            Dict[str, str]: Parsed field into dict
        """
        attrs = {}
        for part in filter(None, gff_attributes.strip().split(";")):
            if "=" in part:
                key, value = part.split("=", 1)
                attrs[key.strip()] = value.strip()
        return attrs

    def _passes_filter(
        self, transcript_id: str, transcript_biotypes: Dict[str, str]
    ) -> bool:
        """Check if a filter passes expression and biotype filters

        Args:
            transcript_id (str): ID of the transcript
            transcript_biotypes (Dict[str, str]): Dict of {transcript_id: biotype}

        Returns:
            bool: T/F if transcript passes filter
        """
        if self.transcript_filter == TranscriptFilter.PROTEIN_CODING:
            biotype = transcript_biotypes.get(transcript_id, "")
            if biotype != "protein_coding":
                return False

        if self.expression_filter and not self.expression_filter.passes_threshold(
            transcript_id
        ):
            return False

        return True
