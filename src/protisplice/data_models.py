"""
File: data_models.py
Description: Data classes for sequence extraction
"""

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Tuple, Optional, Any, TextIO


class StrandType(Enum):
    """Enum for strand types (positive or negative)"""

    POSITIVE = "+"
    NEGATIVE = "-"


class JunctionType(Enum):
    """Enum for junction type (donor/acceptor/intron)"""

    DONOR = "donor"
    ACCEPTOR = "acceptor"
    INTRON = "intron"


class TranscriptFilter(Enum):
    """Enum for transcript filter"""

    ALL = "all"
    PROTEIN_CODING = "protein_coding"
    EXPRESSED = "expressed"


@dataclass
class SpliceJunction:
    """Data class for splice junction info"""

    id: str
    seqid: str
    coord: int  # 1-based coord of first exon base
    strand: StrandType
    junction_type: JunctionType


@dataclass
class JunctionData:
    """Data class for junction data with sequence and coordinates"""

    junction: SpliceJunction
    window_start: int
    window_end: int
    sequence: Optional[str] = None

    def to_fasta_header(self) -> str:
        """Create FASTA header for the splice junction"""
        return (
            f">{self.junction.seqid}_{self.junction.junction_type.value}_"
            f"{self.junction.strand.value}_{self.window_start}_{self.window_end}"
        )


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

    def __post_init__(self):
        """
        Raises:
            ValueError: if any param is less than 0
            ValueError: if any param is not an integer
            ValueError: if both n_exon and n_intron are 0
        """
        if not all(x >= 0 for x in [self.n_exon, self.n_intron, self.buffer_size]):
            raise ValueError("All parameters must be non-negative!")

        if not all(
            isinstance(x, int) for x in [self.n_exon, self.n_intron, self.buffer_size]
        ):
            raise ValueError("All parameters must be integers!")

        if self.n_exon == 0 and self.n_intron == 0:
            raise ValueError("At least one of n_exon or n_intron must be > 0!")

    @property
    def window_size(self) -> int:
        """Total window size"""
        return self.n_exon + self.n_intron


@dataclass
class ExtractionResults:
    """Results from sequence extraction"""

    positive_sequences: Optional[List[JunctionData]] = None
    negative_sequences: Optional[List[JunctionData]] = None

    @property
    def positive_count(self) -> int:
        """Number of positive sequences"""
        return len(self.positive_sequences) if self.positive_sequences else 0

    @property
    def negative_count(self) -> int:
        """Number of negative sequences"""
        return len(self.negative_sequences) if self.negative_sequences else 0

    @property
    def total_count(self) -> int:
        """Total number of sequences"""
        return self.positive_count + self.negative_count

    def get_stats(self) -> Dict[str, int]:
        """Get extraction statistics"""
        return {
            "positive_count": self.positive_count,
            "negative_count": self.negative_count,
            "total_count": self.total_count,
            "has_positive": self.positive_sequences is not None,
            "has_negative": self.negative_sequences is not None,
        }


@dataclass
class NegativeSamplingContext:
    """Context for negative sample extraction"""

    file: TextIO
    fasta: Any  # pysam.FastaFile
    window_size: int
    remaining_samples: int
    samples_written: int = 0


@dataclass
class SequenceWindow:
    """Sequence and its window coordinates"""

    sequence: str
    start: int
    end: int
