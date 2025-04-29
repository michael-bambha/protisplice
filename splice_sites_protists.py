"""
File: splice_sites_protists.py
Author: Michael Bambha
Contact: bambha.m@northeastern.edu
Description: A python script for an ML workflow predicting alternative splice sites in protists.
"""

import argparse
from collections import defaultdict
import re
from typing import Dict, List, Tuple, Any
import pysam


def get_cli_args():
    """
    Parse command line args: gtf, fasta, N_exon, N_intron

    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Extract sequences around splice site junctions."
        )
    parser.add_argument(
        "-g", "--gtf",
        required=True,
        type=str,
        metavar="FILE_PATH",
        help="Path to input GTF file."
    )
    parser.add_argument(
        "-f", "--fasta",
        required=True,
        type=str,
        metavar="FILE_PATH",
        help="Path to input FASTA index file (.fai) indexed by samtools faidx."
    )
    parser.add_argument(
        "-ne", "--n_exon",
        required=True,
        metavar="INT",
        help="Number of bases to include in the exon region of the window."
    )
    parser.add_argument(
        "-ni", "--n_intron",
        required=True,
        metavar="INT",
        help="Number of bases to include in the intron region of the window"
    )


def parse_transcript_id(gtf_str: str) -> Dict[str, str]:
    """Parse column 9 of a GTF file into a dictionary by separating k/v pairs
    by semicolons.

    Args:
        gtf_str (str): Path to GTF file

    Returns:
        Dict[str, str]: Dictionary in the structure:
            {"transcript_id": transcript_id, "gene_id": gene_id ...}
    """
    attrs = {}
    for part in filter(None, gtf_str.strip().split(";")):
        match = re.match(r'\s*(\S+)\s+"([^"]+)"\s*', part)
        if match:
            k, v = match.groups()
            attrs[k] = v
    return attrs


def group_exons_by_transcript(gtf: str) -> Dict[str, Dict[str, Any]]:
    """
    Obtain all exons for a given transcript in a GTF annotation file.
    Example line format (tab-delineated):

    FR824046	ena	exon	1244	1570	.	-	.
    transcript_id "transcript:CCA13858";

    Args:
        gtf (str): GTF file of annotations

    Returns:
        transcripts (Dict[str, Dict[str, Any]]): A dictionary where keys are transcript IDs.
        Each value is another dictionary with two keys:
            'info': {'seqid': str, 'strand': str} - Chromosome/contig and strand.
            'exons': List[Tuple[int, int]] - A list of (start, end) tuples
                     for each exon belonging to the transcript. Coordinates are integers.}
    """
    transcripts = defaultdict(lambda: {"info": {}, "exons": []})
    with open(gtf, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            if fields[2] == "exon":
                attrs = parse_transcript_id(fields[8])
                transcript_id = attrs.get("transcript_id")
                seqid, start, end, strand = [fields[i] for i in [0, 3, 4, 6]]
                transcript_entry = transcripts[transcript_id]
                transcript_entry['exons'].append((int(start), int(end)))
                if not transcript_entry['info']:
                    transcript_entry['info']['seqid'] = seqid
                    transcript_entry['info']['strand'] = strand
    return transcripts


def get_splice_junctions(transcripts: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Find coordinates (start/end) of splice junctions and whether they are donor or acceptor regions.

    Args:
        transcripts (Dict[str, Dict[str, Any]]): A dictionary where keys are transcript IDs.
        Each value is another dictionary with two keys:
            'info': {'seqid': str, 'strand': str} - Chromosome/contig and strand.
            'exons': List[Tuple[int, int]] - A list of (start, end) tuples
                     for each exon belonging to the transcript. Coordinates are integers.}

    Returns:
        List[Dict[str, Any]]: List of Dictionaries, where each element of the list is a splice junction.
        Each junction is a dictionary with 5 keys:
            {'id': }
    """
    junctions = []
    for transcript_id in transcripts:
        transcript = transcripts[transcript_id]
        exons = transcript['exons']
        # sort exons by the start coordinates
        exons_sorted = sorted(exons, key=lambda exon: exon[0])
        strand = transcript["info"]["strand"]
        seqid = transcript["info"]["seqid"]

        if len(exons_sorted) > 1:
            for i in range(len(exons_sorted)):
                exon = exons_sorted[i]
                if i > 0:
                    if strand == "+":
                        junctions.append({'id': f"{transcript_id}_acceptor_{i}", 'seqid': seqid, 'coord': exon[0], 'strand': '+', 'type': 'acceptor'})
                    elif strand == "-":
                        junctions.append({'id': f"{transcript_id}_acceptor_{i}", 'seqid': seqid, 'coord': exon[1], 'strand': '-', 'type': 'acceptor'})
                if i < len(exons) - 1:
                    if strand == "+":
                        junctions.append({'id': f"{transcript_id}_donor_{i}", 'seqid': seqid, 'coord': exon[1], 'strand': '+', 'type': 'donor'})
                    elif strand == "-":
                        junctions.append({'id': f"{transcript_id}_donor_{i}", 'seqid': seqid, 'coord': exon[0], 'strand': '-', 'type': 'donor'})
    return junctions


def get_window_coords(strand: str, junc_type: str, coord: int, N_exon: int, N_intron: int) -> Tuple:
    """Obtains the window coordinates for a sequence, given a known splice site junction coordinate
    and pre-defined window lengths in exon and intron regions.

    Args:
        strand (str): + or - strand of DNA
        junc_type (str): acceptor or donor
        coord (int): splice junction coordinate
        N_exon (int): number of bases to include in exon region.
        N_intron (int): number of bases to include in intron region

    Returns:
        Tuple: Tuple of coordinates for the window(start, end)
    """
    win_start, win_end = None, None
    if strand == "+":
        if junc_type == "donor":
            win_start = coord - N_exon + 1
            win_end = coord + N_intron
        else:
            win_start = coord - N_intron
            win_end = coord + N_exon - 1
    elif strand == "-":
        if junc_type == "donor":
            win_start = coord - N_intron
            win_end = coord + N_exon - 1
    else:
        win_start = coord - N_exon + 1
        win_end = coord + N_intron
    return win_start, win_end


def extract_sequence(fasta: str, seq_id: str, win_start: int, win_end: int) -> str:
    """Obtain a sequence from an indexed FASTA file, given window coordinates.

    Args:
        fasta (str): path to indexed fasta file
        seq_id (str): sequence ID for the desired sequence
        win_start (int): start coordinate (inclusive)
        win_end (int): end coordinate (inclusive)

    Returns:
        str: Sequence of the desired window
    """
    with pysam.FastaFile(fasta) as f:
        if seq_id not in f.references:
            return None
        # pysam takes 0-based coords -- need to subtract 1 from start
        win_start = win_start - 1
        seq = f.fetch(seq_id, win_start, win_end)
    return seq

