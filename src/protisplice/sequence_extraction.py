"""
File: sequence_extraction.py
Author: Michael Bambha
Contact: bambha.m@northeastern.edu
Description: A Python script for obtaining true positive and false positive
sequences around splice sites to be used for downstream model training.
"""

# pylint:disable=no-member

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pysam
from .data_models import (
    SpliceJunction,
    JunctionData,
    JunctionType,
    StrandType,
    Transcript,
    ExtractionParams,
)
from .utils import apply_strand


def add_intron_coords(transcripts: Dict[str, Transcript]) -> Dict[str, Transcript]:
    """Adds introns to Transcript object.

    Args:
        transcripts (Dict[str, Transcript]): Dict of ID: Transcript. Transcript
        introns are initialized to None

    Returns:
        Dict[str, Transcript]: Dict of transcript_id: Transcript. Transcript object:
        {{TranscriptInfo: seqid, strand}, {exons: List[start, end]},
        {introns: List[start, end] = None}}
    """
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


def extract_sequence(
    fasta: pysam.FastaFile, seq_id: str, win_start: int, win_end: int
) -> Optional[str]:
    """Uses Pysam's FastaFile.fetch() to find the sequence of an indexed
    FASTA file. Takes in 1-based start and end coordinates, converts to
    0-based for Pysam compatability, then fetches the sequence.

    Args:
        fasta (pysam.FastaFile): pysam FastaFile object, built from indexed FASTA
        seq_id (str): Reference ID for the FASTA ID (must match reference)
        win_start (int): 1-based inclusive start of the sequence
        win_end (int): 1-based exclusive end of the sequence

    Returns:
        Optional[str]: String of the fetched sequence, or None if coordinates
        are not valid, or if sequence ID not found in references.
    """
    if win_start > win_end or win_start < 1:
        return None

    if seq_id not in fasta.references:
        return None

    seq_len = fasta.get_reference_length(seq_id)
    win_start_0based = max(0, win_start - 1)  # pysam needs 0-based for fetch
    win_end_0based = min(win_end, seq_len)  # pysam end is exclusive

    if win_start_0based >= win_end_0based:
        return None

    extracted = fasta.fetch(seq_id, win_start_0based, win_end_0based)
    return extracted


class SequenceExtractor:
    """Handles sequence extraction from FASTA files"""

    def __init__(self, fasta_path: str, params: ExtractionParams):
        self.fasta_path = Path(fasta_path)
        self.params = params

        if not self.fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {self.fasta_path}")

        fai_path = Path(f"{self.fasta_path}.fai")
        if not fai_path.exists():
            print(f"FASTA index not found -- indexing {fasta_path}...")
            pysam.faidx(fasta_path)
            print(f"FASTA index {fasta_path}.fai created!")

    def extract_splice_sites(
        self, junctions: List[SpliceJunction]
    ) -> List[JunctionData]:
        """Get the sequences around known splice coordinates.

        Args:
            junctions (List[SpliceJunction]): SpliceJunction object containing
            transcript ID, seqID (chr#), coord (1-based coord of first exon base
            at the junction), strand, and junction type.

        Returns:
            List[JunctionData]: List of Dicts in the format {junction: SpliceJunction,
            win_start: win_start, win_end: win_end, seq: seq}.
        """
        sequences = []

        with pysam.FastaFile(str(self.fasta_path)) as fasta:
            for junction in junctions:
                junction_data = self._process_junction(fasta, junction)
                if junction_data:
                    sequences.append(junction_data)

        return sequences

    def _process_junction(
        self, fasta: pysam.FastaFile, junction: SpliceJunction
    ) -> Optional[JunctionData]:
        """Finds the window coordinates using _get_window_coords, then extracts
        the respective sequence using _extract_sequence to obtain sequences
        around splice sites, as dictated by user's input params. Will return
        the reverse complement for sequences identified as (-) strand.

        Args:
            fasta (pysam.FastaFile): pysam FastaFile object built from user input
            FASTA path
            junction (SpliceJunction): SpliceJunction object containing transcript ID, seqID (chr#),
            coord (1-based coord of first exon base at the junction), strand, and junction type.

        Returns:
            Optional[JunctionData]: List of Dicts in the format {junction: SpliceJunction,
            win_start: win_start, win_end: win_end, seq: seq}.
        """
        win_start, win_end = self._get_window_coords(junction)

        if win_start is None or win_end is None:
            return None

        seq = extract_sequence(fasta, junction.seqid, win_start, win_end)
        if not seq:
            return None

        seq = apply_strand(seq, junction.strand)

        return JunctionData(
            junction=junction, window_start=win_start, window_end=win_end, sequence=seq
        )

    def _get_window_coords(
        self, junction: SpliceJunction
    ) -> Tuple[Optional[int], Optional[int]]:
        """Based on identified junctions, applies the user-defined
        number of bases to include for the intron and exon regions.
        Sequence coordinates are identified based on the supplied coordinate
        adjustments, and function returns (start, end) of the adjusted
        sequence.

        Args:
            junction (SpliceJunction): SpliceJunction object containing transcript ID, seqID (chr#),
            coord (1-based coord of first exon base at the junction), strand, and junction type.


        Returns:
            Tuple[Optional[int], Optional[int]]: Tuple of (start, end), 1-based coordinates
            calculated from the user's input parameters and the junction coordinates or None
            if strand type and/or junction type are not defined.
        """
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
