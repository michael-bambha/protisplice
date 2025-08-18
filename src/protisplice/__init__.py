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
    ExtractionResults,
    TranscriptFilter,
    SpliceJunction,
    JunctionData,
    JunctionType,
    StrandType,
    Transcript,
    TranscriptInfo,
)
from .expression import ExpressionFilter, ExpressionParser
from .gff_parser import GFFParser
from .fasta_writer import FastaWriter
from .motif_scoring import generate_ppm, generate_pfm, generate_pwm

__all__ = [
    "SpliceSeqExtractor",
    "ExtractionParams",
    "ExtractionResults",
    "TranscriptFilter",
    "Transcript",
    "TranscriptInfo",
    "SpliceJunction",
    "JunctionData",
    "JunctionType",
    "StrandType",
    "FastaWriter",
    "ExpressionFilter",
    "ExpressionParser",
    "GFFParser",
    "SequenceExtractor",
    "SpliceJunctionExtractor",
    "generate_ppm",
    "generate_pwm",
    "generate_pfm",
]

__version__ = "0.1.0"
__author__ = "Michael Bambha"
__email__ = "bambha.m@northeastern.edu"
