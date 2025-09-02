"""
Package for splice sequence extraction.

This package provides functionality for extracting splice site sequences from a FASTA and
associated GFF3 annotation file.
"""

from .extract_splice_seqs import SpliceSeqExtractor
from .sequence_extraction import SequenceExtractor
from .splice_junction_extractor import SpliceJunctionExtractor
from .data_models import (
    ExtractionParams,
    TranscriptFilter,
    SpliceJunction,
    JunctionData,
    JunctionType,
    StrandType,
    Transcript,
    TranscriptInfo,
    SamplingParams,
    Gene
)
from .expression import ExpressionFilter, ExpressionParser
from .gff_parser import GFFParser
from .fasta_writer import write_sequences
from .motif_scoring import generate_ppm, generate_pfm, generate_pwm
from .feature_eng import inject_consensus, remove_consensus

__all__ = [
    "SpliceSeqExtractor",
    "ExtractionParams",
    "TranscriptFilter",
    "Transcript",
    "TranscriptInfo",
    "SpliceJunction",
    "JunctionData",
    "JunctionType",
    "StrandType",
    "write_sequences",
    "ExpressionFilter",
    "ExpressionParser",
    "GFFParser",
    "Gene",
    "SamplingParams",
    "SequenceExtractor",
    "SpliceJunctionExtractor",
    "generate_ppm",
    "generate_pwm",
    "generate_pfm",
    "inject_consensus",
    "remove_consensus"
]

__version__ = "0.1.0"
__author__ = "Michael Bambha"
__email__ = "bambha.m@northeastern.edu"
