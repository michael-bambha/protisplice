"""
File: extract_splice_seqs.py
Author: Michael Bambha
Contact: bambha.m@northeastern.edu
Description: A Python script for obtaining true positive and false positive
sequences around splice sites to be used for downstream model training.
"""

# pylint: disable=no-member
import argparse
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Tuple, Optional, TextIO
import logging
from pathlib import Path
import random
import pysam
from Bio.Seq import Seq

RANDOM_SEED = 100


class StrandType(Enum):
    """Enum for strand types (positive or negative)"""

    POSITIVE = "+"
    NEGATIVE = "-"


class JunctionType(Enum):
    """Enum for junction type (donor/acceptor/intron)"""

    DONOR = "donor"
    ACCEPTOR = "acceptor"
    INTRON = "intron"


@dataclass
class SpliceJunction:
    """Data class for splice junction info"""

    id: str
    seqid: str
    coord: int
    strand: StrandType
    junction_type: JunctionType


class TranscriptFilter(Enum):
    """Enum for transcript filter"""

    ALL = "all"
    PROTEIN_CODING = "protein_coding"
    EXPRESSED = "expressed"  # implement RNA-seq integration later


@dataclass
class TranscriptInfo:
    """Data class for transcript info"""

    seqid: str
    strand: StrandType


@dataclass
class Transcript:
    """Data class for transcript with exons and introns"""

    info: TranscriptInfo
    exons: List[Tuple[int, int]]
    introns: Optional[List[Tuple[int, int]]] = None


@dataclass
class ExtractionParams:
    """Parameters for sequence extraction"""

    n_exon: int = 40
    n_intron: int = 80
    buffer_size: int = 50


@dataclass
class NegativeSamplingContext:
    """Context for negative sample extraction"""

    file: TextIO
    fasta: pysam.FastaFile
    window_size: int
    remaining_samples: int
    samples_written: int = 0


@dataclass
class IntronCoordinates:
    """Intron coordinate information"""

    start: int
    end: int


@dataclass
class SequenceWindow:
    """Sequence and its window coordinates"""

    sequence: str
    start: int
    end: int


@dataclass
class ExpressionData:
    """Expression Data dictionary"""

    transcript_id: str
    expression_val: float
    gene_id: Optional[str] = None


@dataclass
class ExpressionFilter:
    """Class for filtering expression data based on a threshold"""

    data: Dict[str, float]  # tid: express_val
    threshold: float

    def passes_threshold(self, transcript_id: str) -> bool:
        """Determine whether or not expression value is >= threshold"""
        expression_val = self._find_expression_value(transcript_id)
        return expression_val >= self.threshold

    def _find_expression_value(self, transcript_id: str) -> float:
        """Get the expression for a particular transcript"""
        normalized_id = self._normalize_transcript_id(transcript_id)
        if normalized_id in self.data:
            return self.data[normalized_id]
        for expr_id in self.data:
            if self._normalize_transcript_id(expr_id) == normalized_id:
                return self.data[expr_id]

        return 0.0

    def _normalize_transcript_id(self, transcript_id: str) -> str:
        """Normalize transcript IDs across different expression formats"""
        if "." in transcript_id:
            transcript_id = transcript_id.split(".")[0]
        if transcript_id.startswith("transcript:"):
            transcript_id = transcript_id.replace("transcript:", "")
        return transcript_id


class ExpressionParser:
    """Class for building parsers for different gene expression formats"""

    @staticmethod
    def parse_file(file_path: str, format_type: str) -> Dict[str, float]:
        """General method for parsing different expression formats"""
        if format_type == "kallisto":
            return ExpressionParser._parse_kallisto(file_path)
        # elif format_type == "salmon":
        #     return ExpressionParser._parse_salmon(file_path)
        # elif format_type == "stringtie":
        #     return ExpressionParser._parse_stringtie(file_path)
        return None

    @staticmethod
    def _parse_kallisto(file_path: str) -> Dict[str, float]:
        expression_data = {}
        with open(file_path, "r", encoding="utf-8") as f:
            header = next(f).strip().split("\t")
            tpm_idx = header.index("tpm")
            target_idx = header.index("target_id")

            for line in f:
                fields = line.strip().split("\t")
                transcript_id = fields[target_idx]
                tpm = float(fields[tpm_idx])
                expression_data[transcript_id] = tpm
        return expression_data


class SpliceSeqExtractor:
    """Main class for splice sequence extraction"""

    def __init__(
        self,
        gff_path: str,
        fasta_path: str,
        params: ExtractionParams,
        transcript_filter: str = "all",
        expression_file: Optional[str] = None,
        min_expression: float = 1.0,
        expression_format: str = "kallisto",
    ):
        """Initialize with file paths and extraction parameters"""
        self.gff_path = Path(gff_path)
        self.fasta_path = Path(fasta_path)
        self.params = params
        self.transcript_filter = transcript_filter

        self._validate_inputs()
        logging.basicConfig(
            level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
        )
        self.logger = logging.getLogger(__name__)
        self.expression_filter = None

        if expression_file:
            expression_data = ExpressionParser.parse_file(
                expression_file, expression_format
            )
            self.expression_filter = ExpressionFilter(expression_data, min_expression)
            self.logger.info(
                "Loaded expression data for %d transcripts", len(expression_data)
            )

    def _validate_inputs(self) -> None:
        """Validate input files and parameters"""
        if not self.fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {self.fasta_path}")

        fai_path = Path(f"{self.fasta_path}.fai")
        if not fai_path.exists():
            raise FileNotFoundError(
                f"Fasta index {fai_path} not found. Run samtools faidx."
            )

        if not self.gff_path.exists():
            raise FileNotFoundError(f"GFF3 file not found: {self.gff_path}")

        if not all(x >= 0 for x in [self.params.n_exon, self.params.n_intron]):
            raise ValueError("n_intron and n_exon must be positive integers!")

        if not self.params.buffer_size >= 0:
            raise ValueError("Buffer size cannot be negative!")

        if not all(
            isinstance(x, int)
            for x in [self.params.n_exon, self.params.n_intron, self.params.buffer_size]
        ):
            raise ValueError("n_intron and n_exon must be integers!")

    def extract_sequences(
        self, positive_output: str, negative_output: str
    ) -> Tuple[int, int]:
        """Extract positive and negative sequences"""
        self.logger.info("Parsing transcripts from GFF file...")
        transcripts = self._parse_transcripts()

        self.logger.info("Identifying splice junctions...")
        junctions = self._get_splice_junctions(transcripts)

        self.logger.info("Extracting positive samples to %s...", positive_output)
        positive_count = self._extract_positive_samples(junctions, positive_output)

        self.logger.info("Finding intron coordinates...")
        transcripts_with_introns = self._add_intron_coords(transcripts)

        self.logger.info("Extracting negative samples to %s...", negative_output)
        negative_count = self._extract_negative_samples(
            transcripts_with_introns, negative_output, positive_count
        )

        self.logger.info(
            "Extraction complete -- found %d positive samples and %d negative samples.",
            positive_count,
            negative_count,
        )
        return positive_count, negative_count

    def get_sequence_stats(self) -> Dict[str, int]:
        """Get statistics about sequences"""
        window_size = self.params.n_exon + self.params.n_intron
        return {
            "window_size": window_size,
            "exon_bases": self.params.n_exon,
            "intron_bases": self.params.n_intron,
            "buffer_size": self.params.buffer_size,
        }

    def _parse_transcripts(self) -> Dict[str, Transcript]:
        """Parse GFF3 file to extract transcript information"""
        transcript_biotypes = {}
        if self.transcript_filter != "all":
            transcript_biotypes = self._collect_transcript_biotypes()
            self.logger.info(
                "Found %d transcripts with biotype info", len(transcript_biotypes)
            )

        transcripts = defaultdict(
            lambda: Transcript(
                info=TranscriptInfo(seqid="", strand=StrandType.POSITIVE), exons=[]
            )
        )

        with open(self.gff_path, "r", encoding="utf-8") as file:
            for line_num, line in enumerate(file, 1):
                transcript_data = self._parse_gff_line(
                    line, line_num, transcript_biotypes
                )
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
        self, line: str, line_num: int, transcript_biotypes: Dict[str, str]
    ) -> Optional[Tuple]:
        """Parse a single GFF line and return transcript data if valid"""
        try:
            if line.startswith("#") or not line.strip():
                return None

            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                return None

            attrs = self._extract_transcript_attributes(fields[8])
            parent = attrs.get("Parent", "")
            transcript_id = None
            if parent.startswith("transcript:"):
                transcript_id = parent.replace("transcript:", "")

            if not transcript_id:
                return None

            if transcript_biotypes is not None and not self._passes_filter(
                transcript_id, transcript_biotypes
            ):
                return None

            seqid, start, end, strand = (
                fields[0],
                int(fields[3]),
                int(fields[4]),
                fields[6],
            )
            return transcript_id, seqid, start, end, strand

        except (IndexError, ValueError) as error:
            self.logger.warning("Skipping malformed line %d: %s", line_num, error)
            return None

    def _collect_transcript_biotypes(self) -> Dict[str, str]:
        """Collect transcript biotype information from mRNA lines"""
        transcript_biotypes = {}

        with open(self.gff_path, "r", encoding="utf-8") as file:
            for line in file:
                if line.startswith("#") or not line.strip():
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 9 or fields[2] != "mRNA":
                    continue

                attrs = self._extract_transcript_attributes(fields[8])
                transcript_id = attrs.get("ID", "")
                if transcript_id.startswith("transcript:"):
                    transcript_id = transcript_id.replace("transcript:", "")
                    biotype = attrs.get("biotype", "")
                    if biotype:
                        transcript_biotypes[transcript_id] = biotype

        return transcript_biotypes

    def _extract_transcript_attributes(self, gff_attributes: str) -> Dict[str, str]:
        """Extract transcript ID from GFF3 attributes field"""
        attrs = {}
        for part in filter(None, gff_attributes.strip().split(";")):
            if "=" in part:
                key, value = part.split("=", 1)
                attrs[key.strip()] = value.strip()
        return attrs

    def _passes_filter(
        self, transcript_id: str, transcript_biotypes: Dict[str, str]
    ) -> bool:
        """Check if the transcript passes filter"""

        if self.transcript_filter == "protein_coding":
            biotype = transcript_biotypes.get(transcript_id, "")
            if biotype != "protein_coding":
                return False

        if self.expression_filter and not self.expression_filter.passes_threshold(
            transcript_id
        ):
            return False

        return True

    def _get_splice_junctions(
        self, transcripts: Dict[str, Transcript]
    ) -> List[SpliceJunction]:
        """Get splice junction coordinates from transcripts"""
        junctions = []
        # Need 2+ exons to have a defined junction
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
        self,
        transcript_id: str,
        transcript: Transcript,
        sorted_exons: List[Tuple[int, int]],
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

    def _get_window_coords(
        self, junction: SpliceJunction
    ) -> Tuple[Optional[int], Optional[int]]:
        strand = junction.strand
        junc_type = junction.junction_type
        coord = junction.coord

        if strand == StrandType.POSITIVE:
            if junc_type == JunctionType.DONOR:
                start = coord - self.params.n_exon + 1
                end = coord + self.params.n_intron
                return start, end

            if junc_type == JunctionType.ACCEPTOR:
                start = coord - self.params.n_intron
                end = coord + self.params.n_exon - 1
                return start, end

        if strand == StrandType.NEGATIVE:
            if junc_type == JunctionType.DONOR:
                start = coord - self.params.n_intron
                end = coord + self.params.n_exon - 1
                return start, end

            if junc_type == JunctionType.ACCEPTOR:
                start = coord - self.params.n_exon + 1
                end = coord + self.params.n_intron
                return start, end

        return None, None

    def _extract_sequence(
        self, fasta: pysam.FastaFile, seq_id: str, win_start: int, win_end: int
    ) -> Optional[str]:
        """Extract sequence from FASTA file with coordinate validation"""
        if win_start > win_end or win_start < 1:
            return None

        if seq_id not in fasta.references:
            return None

        seq_len = fasta.get_reference_length(seq_id)
        win_start_0based = max(0, win_start - 1)  # Convert to 0-based
        win_end_0based = min(
            win_end, seq_len
        )  # Keep as 1-based since pysam end is exclusive

        if win_start_0based >= win_end_0based:
            return None

        extracted = fasta.fetch(seq_id, win_start_0based, win_end_0based)
        return extracted

    def _extract_positive_samples(
        self, junctions: List[SpliceJunction], output_path: str
    ) -> int:
        """Extract positive splice site sequences"""
        count = 0

        with open(output_path, "w", encoding="utf-8") as file:
            with pysam.FastaFile(str(self.fasta_path)) as fasta:
                for junction in junctions:
                    count += self._process_junction(file, fasta, junction)

        return count

    def _process_junction(
        self, file: TextIO, fasta: pysam.FastaFile, junction: SpliceJunction
    ) -> int:
        """Process a single junction and write to file if valid"""
        win_start, win_end = self._get_window_coords(junction)

        if win_start is None or win_end is None:
            return 0

        seq = self._extract_sequence(fasta, junction.seqid, win_start, win_end)
        if not seq:
            return 0

        # Get reverse complement for (-) seqs
        if junction.strand == StrandType.NEGATIVE:
            seq = str(Seq(seq).reverse_complement())

        seq_window = SequenceWindow(sequence=seq, start=win_start, end=win_end)

        self._write_fasta_entry(file, junction, seq_window)
        return 1

    def _add_intron_coords(
        self, transcripts: Dict[str, Transcript]
    ) -> Dict[str, Transcript]:
        """Add intron coordinates to transcript objects"""
        for transcript in transcripts.values():
            introns = []
            sorted_exons = sorted(transcript.exons, key=lambda x: x[0])

            for i in range(len(sorted_exons) - 1):
                intron_start = sorted_exons[i][1] + 1
                intron_end = sorted_exons[i + 1][0] - 1

                if intron_start <= intron_end:
                    introns.append((intron_start, intron_end))

            transcript.introns = introns
        return transcripts

    def _extract_negative_samples(
        self, transcripts: Dict[str, Transcript], output_path: str, target_count: int
    ) -> int:
        """Extract negative samples from intronic regions"""
        random.seed(RANDOM_SEED)
        samples_written = 0
        window_size = self.params.n_exon + self.params.n_intron

        with open(output_path, "w", encoding="utf-8") as file:
            with pysam.FastaFile(str(self.fasta_path)) as fasta:
                transcript_ids = list(transcripts.keys())
                random.shuffle(transcript_ids)

                for transcript_id in transcript_ids:
                    if samples_written >= target_count:
                        break

                    context = NegativeSamplingContext(
                        file=file,
                        fasta=fasta,
                        window_size=window_size,
                        remaining_samples=target_count - samples_written,
                        samples_written=samples_written,
                    )

                    extracted = self._extract_from_transcript(
                        context, transcript_id, transcripts[transcript_id]
                    )
                    samples_written += extracted

        return samples_written

    def _extract_from_transcript(
        self,
        context: NegativeSamplingContext,
        transcript_id: str,
        transcript: Transcript,
    ) -> int:
        """Extract negative samples from a single transcript"""
        if not transcript.introns:
            return 0

        samples_from_transcript = 0

        for intron_start, intron_end in transcript.introns:
            if samples_from_transcript >= context.remaining_samples:
                break

            intron_coords = IntronCoordinates(start=intron_start, end=intron_end)
            sample_extracted = self._extract_from_intron(
                context, transcript_id, transcript, intron_coords
            )
            samples_from_transcript += sample_extracted
            context.samples_written += sample_extracted

        return samples_from_transcript

    def _extract_from_intron(
        self,
        context: NegativeSamplingContext,
        transcript_id: str,
        transcript: Transcript,
        intron_coords: IntronCoordinates,
    ) -> int:
        """Extract one negative sample from an intron region"""
        # Define safe extraction region (avoid splice sites)
        safe_start = intron_coords.start + self.params.buffer_size
        safe_end = intron_coords.end - self.params.buffer_size

        if safe_end - safe_start + 1 < context.window_size:
            return 0

        # Random sampling within safe region
        max_start = safe_end - context.window_size + 1
        if safe_start > max_start:
            return 0

        rand_start = random.randint(safe_start, max_start)
        rand_end = rand_start + context.window_size - 1

        seq = self._extract_sequence(
            context.fasta, transcript.info.seqid, rand_start, rand_end
        )
        if not seq:
            return 0

        # Apply reverse complement for negative strand
        if transcript.info.strand == StrandType.NEGATIVE:
            seq = str(Seq(seq).reverse_complement())

        # Create pseudo-junction for consistent output format
        pseudo_junction = SpliceJunction(
            id=f"{transcript_id}_intron_{context.samples_written}",
            seqid=transcript.info.seqid,
            coord=rand_start,
            strand=transcript.info.strand,
            junction_type=JunctionType.INTRON,
        )

        seq_window = SequenceWindow(sequence=seq, start=rand_start, end=rand_end)
        self._write_fasta_entry(context.file, pseudo_junction, seq_window)
        return 1

    def _write_fasta_entry(
        self, file_handle: TextIO, junction: SpliceJunction, seq_window: SequenceWindow
    ) -> None:
        """Write a FASTA entry to file"""
        header = (
            f">{junction.seqid}_{junction.junction_type.value}_"
            f"{junction.strand.value}_{seq_window.start}_{seq_window.end}"
        )
        file_handle.write(f"{header}\n{seq_window.sequence}\n")


def get_cli_args() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Extract sequences around splice site junctions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-g",
        "--gff",
        required=True,
        type=str,
        metavar="FILE_PATH",
        help="Path to input GFF3 file",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        type=str,
        metavar="FILE_PATH",
        help="Path to input FASTA file (must be indexed with samtools faidx)",
    )
    parser.add_argument(
        "-o1",
        "--out1",
        required=True,
        type=str,
        metavar="FILE_PATH",
        help="Path to output file for positive sequences",
    )
    parser.add_argument(
        "-o2",
        "--out2",
        required=True,
        type=str,
        metavar="FILE_PATH",
        help="Path to output file for negative sequences",
    )
    parser.add_argument(
        "-ne",
        "--n_exon",
        default=40,
        type=int,
        metavar="INT",
        help="Number of bases to include from exon region",
    )
    parser.add_argument(
        "-ni",
        "--n_intron",
        default=80,
        type=int,
        metavar="INT",
        help="Number of bases to include from intron region",
    )
    parser.add_argument(
        "-b",
        "--buffer",
        default=50,
        type=int,
        metavar="INT",
        help="Buffer size around splice sites for negative sampling",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )
    parser.add_argument(
        "--transcript-filter",
        type=str,
        choices=["all", "protein_coding"],
        help="Include all or only protein coding transcripts",
        default="protein_coding",
    )
    parser.add_argument(
        "--expression-file",
        type=str,
        help="Path to expression quantification file (TSV format)",
    )
    parser.add_argument(
        "--min-expression",
        type=float,
        help="Minimum expression to filter on (TPM/FPKM)",
    )
    parser.add_argument(
        "--expression-format",
        choices=["kallisto", "salmon", "stringtie"],
        help="Format of expression file (kallisto / salmon / stringtie)",
    )

    return parser.parse_args()


def main() -> None:
    """Main entry point"""
    args = get_cli_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    params = ExtractionParams(
        n_exon=args.n_exon, n_intron=args.n_intron, buffer_size=args.buffer
    )

    extractor = SpliceSeqExtractor(
        gff_path=args.gff,
        fasta_path=args.fasta,
        params=params,
        transcript_filter=args.transcript_filter,
        expression_file=args.expression_file,
        min_expression=args.min_expression,
        expression_format=args.expression_format,
    )

    try:
        pos_count, neg_count = extractor.extract_sequences(args.out1, args.out2)
        print(
            f"Successfully extracted {pos_count} positive and "
            f"{neg_count} negative sequences"
        )
    except Exception as error:
        logging.error("Extraction failed: %s", error)
        raise


if __name__ == "__main__":
    main()
