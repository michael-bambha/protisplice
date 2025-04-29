"""
File: splice_sites_protists.py
Author: Michael Bambha
Contact: bambha.m@northeastern.edu
Description: A python script for an ML workflow predicting alternative splice sites in protists.
"""

import argparse
from collections import defaultdict
from typing import Dict, List, Tuple, Any, Optional, TextIO
import re
import pysam
from Bio.seq import Seq


def main():
    """
    Business logic
    """
    args = get_cli_args()
    transcripts = group_exons_by_transcript(args.gtf)
    junctions = get_splice_junctions(transcripts)
    output_file = "output.txt"
    with open(output_file, "w", encoding='utf-8') as f:
        for junction in junctions:
            seqid = junction['seqid']
            coord = junction['coord']
            strand = junction['strand']
            junc_type = junction['type']
            win_start, win_end = get_window_coords(strand, junc_type, coord, args.n_exon,
                                                   args.n_intron)
            if win_start and win_end:
                seq = extract_sequence(args.fasta, seqid, win_start, win_end)
                if seq and strand == "-":  # take reverse complement on the (-) strand seqs
                    seq = str(Seq(seq).reverse_complement())
                write_seq_to_file(seq, f)


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
        default=40,
        type=int,
        metavar="INT",
        help="Number of bases to include in the exon region of the window."
    )
    parser.add_argument(
        "-ni", "--n_intron",
        default=80,
        type=int,
        metavar="INT",
        help="Number of bases to include in the intron region of the window"
    )
    return parser.parse_args()


def _parse_transcript_id(gtf_str: str) -> Dict[str, str]:
    """Parse column 9 of a GTF file into a dictionary by separating k/v pairs
    by semicolons.

    Args:
        gtf_str (str): 9th column of a GTF file

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
    with open(gtf, "r", encoding='utf-8') as f:
        for line in f:
            fields = line.strip().split("\t")
            if fields[2] == "exon":
                attrs = _parse_transcript_id(fields[8])
                transcript_id = attrs.get("transcript_id")
                seqid, start, end, strand = [fields[i] for i in [0, 3, 4, 6]]
                transcript_entry = transcripts[transcript_id]
                transcript_entry['exons'].append((int(start), int(end)))
                if not transcript_entry['info']:
                    transcript_entry['info']['seqid'] = seqid
                    transcript_entry['info']['strand'] = strand
    return transcripts


def get_splice_junctions(transcripts: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Find coordinates (start/end) of splice junctions and whether
      they are donor or acceptor regions.
    Args:
        transcripts (Dict[str, Dict[str, Any]]): A dictionary where keys are transcript IDs.
        Each value is another dictionary with two keys:
            'info': {'seqid': str, 'strand': str} - Chromosome/contig and strand.
            'exons': List[Tuple[int, int]] - A list of (start, end) tuples
                     for each exon belonging to the transcript. Coordinates are integers.}

    Returns:
        List[Dict[str, Any]]: List of Dictionaries, where each element of the list
          is a splice junction. Each junction is a dictionary with 5 keys:
            {'id': str, 'seqid': str, 'coord': int, 'strand': str, 'type': str}
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
            for i, exon in enumerate(exons_sorted):
                # exon is a tuple as (start, end)
                start = exon[0]
                end = exon[1]
                if i > 0:
                    acceptor_coord = None
                    if strand == "+":
                        # acceptor is at the start of the current exon
                        acceptor_coord = start
                        junctions.append({
                            'id': f"{transcript_id}_acceptor_{i}",  # index refers to preceding exon
                            'seqid': seqid,
                            'coord': acceptor_coord,
                            'strand': '+',
                            'type': 'acceptor'
                        })
                    elif strand == "-":
                        # acceptor is at the end of the current exon
                        acceptor_coord = end
                        junctions.append({
                            'id': f"{transcript_id}_acceptor_{i}",
                            'seqid': seqid,
                            'coord': acceptor_coord,
                            'strand': '-',
                            'type': 'acceptor'
                        })
                if i < len(exons_sorted) - 1:
                    donor_coord = None
                    if strand == "+":
                        # donor is at the end of the current exon
                        donor_coord = end
                        junctions.append({
                            'id': f"{transcript_id}_donor_{i}",  # index refers to following exon
                            'seqid': seqid,
                            'coord': donor_coord,
                            'strand': '+',
                            'type': 'donor'
                        })
                    elif strand == "-":
                        # donor is at the start of the current exon
                        donor_coord = start
                        junctions.append({
                            'id': f"{transcript_id}_donor_{i}",
                            'seqid': seqid,
                            'coord': donor_coord,
                            'strand': '-',
                            'type': 'donor'
                        })
    return junctions


def get_window_coords(strand: str, junc_type: str, coord: int, n_exon: int, n_intron: int) -> Tuple:
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
            win_start = coord - n_exon + 1
            win_end = coord + n_intron
        elif junc_type == "acceptor":
            win_start = coord - n_intron
            win_end = coord + n_exon - 1
    elif strand == "-":
        if junc_type == "donor":
            win_start = coord - n_intron
            win_end = coord + n_exon - 1
        elif junc_type == "acceptor":
            win_start = coord - n_exon + 1
            win_end = coord + n_intron
    return win_start, win_end


def extract_sequence(fasta: str, seq_id: str, win_start: int, win_end: int) -> Optional[str]:
    """Obtain a sequence from an indexed FASTA file, given window coordinates.

    Args:
        fasta (str): path to indexed fasta file
        seq_id (str): sequence ID for the desired sequence
        win_start (int): start coordinate (inclusive)
        win_end (int): end coordinate (inclusive)

    Returns:
        str: Sequence of the desired window
    """
    if win_start > win_end:
        return ""
    if win_start < 1:
        win_start = 1
        if win_start > win_end:
            return ""
    seq = ""
    with pysam.FastaFile(fasta) as f:  # pylint: disable=no-member
        if seq_id not in f.references:
            return ""
        seq_len = f.get_reference_length(seq_id)
        win_start = win_start - 1  # pysam takes 0-based coords
        if win_end > seq_len:  # truncate end if it is larger than seq length
            win_end = seq_len
        if win_start >= win_end:  # check after end is truncated also
            return ""
        seq = f.fetch(seq_id, win_start, win_end)
    return seq


def write_seq_to_file(seq: str, f: TextIO) -> None:
    """Write sequences to a file.

    Args:
        seq (str): sequence extracted from FASTA
        f (TextIO): Output file to write the seq
    """
    if seq:
        f.write(f"{seq}\n")


main()
