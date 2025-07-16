"""
Package for splice sequence extraction.

This package provides functionality for extracting splice site sequences from a FASTA and
associated GFF3 annotation file.
"""

from .extract_splice_seqs import SpliceSeqExtractor

from .data_models import (
    ExtractionParams,
    ExtractionResults,
    TranscriptFilter,
    SpliceJunction,
    JunctionData,
    JunctionType,
    StrandType,
)
from .fasta_writer import FastaWriter

__all__ = [
    "SpliceSeqExtractor",
    "ExtractionParams",
    "ExtractionResults",
    "TranscriptFilter",
    "SpliceJunction",
    "JunctionData",
    "JunctionType",
    "StrandType",
    "FastaWriter",
]

__version__ = "0.1.0"
__author__ = "Michael Bambha"
__email__ = "bambha.m@northeastern.edu"
