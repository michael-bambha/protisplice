"""
File: splice_sites_protists.py
Author: Michael Bambha
Contact: bambha.m@northeastern.edu
Description: A Python script for obtaining true positive and false positive
sequences around splice sites in protists to be used for downstream model training.
"""
# pylint: disable=no-member
import argparse
from collections import defaultdict
from typing import Dict, List, Tuple, Any, Optional, TextIO
import random
import re
import pysam
from Bio.Seq import Seq


def main() -> None:
    """
    Business logic
    """
    args = get_cli_args()
    transcripts = group_exons_by_transcript(args.gtf)
    junctions = get_splice_junctions(transcripts)
    count = extract_positive_samples(junctions,
                                     args.fasta,
                                     "seqs_positive.fa",
                                     args.n_exon,
                                     args.n_intron)
    transcripts_with_introns = get_intron_coords(transcripts)

    extract_negative_samples(transcripts_with_introns,
                             args.fasta,
                             "seqs_negative.fa",
                             args.n_exon + args.n_intron,  # make equal length to (+) seqs
                             args.buffer,
                             count)  # sample the same number of (-) samples


def get_cli_args() -> argparse.Namespace:
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
    parser.add_argument(
        "-b", "--buffer",
        default=50,
        type=int,
        metavar="INT",
        help="Buffer size for intron region"
    )
    return parser.parse_args()


def extract_positive_samples(junctions: List[Dict[str, Any]],
                             fasta_path: str,
                             output_path: str,
                             n_exon: int,
                             n_intron: int) -> int:
    """Logic for obtaining true positives

    Args:
        junctions (List[Dict[str, Any]]): List of Dictionaries, where each element of the list
            is a splice junction. Each junction is a dictionary with 5 keys:
            {'id': str, 'seqid': str, 'coord': int, 'strand': str, 'type': str}
        fasta_path (str): path to FASTA file
        output_path (str): path to output file
        n_exon (int): Number of bases to include in the exon region of the window.
        n_intron (int): Number of bases to include in the intron region of the window.

    Returns:
        int: The number of samples written to the file.
    """
    count = 0
    with open(output_path, "w", encoding='utf-8') as f:
        with pysam.FastaFile(fasta_path) as fasta:
            for junction in junctions:
                seqid = junction['seqid']
                coord = junction['coord']
                strand = junction['strand']
                junc_type = junction['type']
                win_start, win_end = get_window_coords(strand, junc_type, coord, n_exon, n_intron)
                if win_start is not None and win_end is not None:
                    seq = extract_sequence(fasta, seqid, win_start, win_end)
                    if seq:
                        if strand == "-":
                            seq = str(Seq(seq).reverse_complement())  # get RC for (-) strands
                        write_seq_to_file(seq, f)
                        count += 1
    return count


def extract_negative_samples(transcripts: Dict[str, Dict[str, Any]],
                             fasta_path: str,
                             output_path: str,
                             window_size: int,
                             buffer_size: int,
                             num_samples: int) -> None:
    """Logic for obtaining negative samples.

    Args:
        transcripts (Dict[str, Dict[str, Any]]): A dictionary where keys are transcript IDs.
        Each value is another dictionary with three keys:
            'info': {'seqid': str, 'strand': str} - Chromosome/contig and strand.
            'exons': List[Tuple[int, int]] - A list of (start, end) tuples
                     for each exon belonging to the transcript. Coordinates are integers.
            'introns': List[Tuple[int, int]] A list of (start, end) tuples
                     for each intron between subsequent exons. Coordinates are integers.
        fasta_path (str): path to FASTA file
        output_path (str): path to output file
        window_size (int): size of sequence to extract
        buffer_size (int): size of buffer in intron window
        num_samples (int): number of sequences to sample
    """
    samples_written = 0
    with open(output_path, "w", encoding='utf-8') as f:
        with pysam.FastaFile(fasta_path) as fasta:
            transcript_ids = list(transcripts.keys())
            random.shuffle(transcript_ids)
            for transcript_id in transcript_ids:
                if samples_written >= num_samples:
                    break
                transcript_data = transcripts[transcript_id]
                seqid = transcript_data['info']['seqid']
                introns = transcript_data.get('introns', [])
                for intron_start, intron_end in introns:
                    if samples_written >= num_samples:
                        break
                    safe_start = intron_start + buffer_size
                    safe_end = intron_end - buffer_size
                    if safe_start < safe_end and (safe_end - safe_start + 1) >= window_size:
                        max_possible_start = safe_end - window_size + 1
                        if safe_start <= max_possible_start:
                            rand_coord = random.randint(safe_start, max_possible_start)
                            win_start = rand_coord
                            win_end = rand_coord + window_size - 1
                            seq = extract_sequence(fasta, seqid, win_start, win_end)
                            if seq:
                                write_seq_to_file(seq, f)
                                samples_written += 1


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
            try:
                if fields[2] == "exon":
                    attrs = _parse_transcript_id(fields[8])
                    transcript_id = attrs.get("transcript_id")
                    seqid, start, end, strand = [fields[i] for i in [0, 3, 4, 6]]
                    transcript_entry = transcripts[transcript_id]
                    transcript_entry['exons'].append((int(start), int(end)))
                    if not transcript_entry['info']:
                        transcript_entry['info']['seqid'] = seqid
                        transcript_entry['info']['strand'] = strand
            except IndexError:
                raise
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
    def _add_junction_entry(
        transcript_id_val: str,
        seqid_val: str,
        strand_val: str,
        coord_val: int,
        junction_type_val: str,
        exon_index: int
    ) -> None:
        junctions.append({
            'id': f"{transcript_id_val}_{junction_type_val}_{exon_index}",
            'seqid': seqid_val,
            'coord': coord_val,
            'strand': strand_val,  # The strand of the junction is the transcript's strand
            'type': junction_type_val
        })
    junctions = []
    for transcript_id, transcript_data in transcripts.items():
        transcript_info = transcript_data['info']
        exons = transcript_data['exons']
        # sort exons by the start coordinates
        exons_sorted = sorted(exons, key=lambda exon: exon[0])
        if len(exons_sorted) <= 1:
            continue
        seqid = transcript_info['seqid']
        strand = transcript_info['strand']
        for i, exon_coords in enumerate(exons_sorted):
            exon_start, exon_end = exon_coords
            if i > 0:
                acceptor_coord = exon_start if strand == "+" else exon_end
                _add_junction_entry(transcript_id, seqid, strand, acceptor_coord, 'acceptor', i)
            if i < len(exons_sorted) - 1:
                donor_coord = exon_end if strand == "+" else exon_start
                _add_junction_entry(transcript_id, seqid, strand, donor_coord, 'donor', i)
    return junctions


def get_intron_coords(transcripts: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """Find the intron coordinates between all exons in the transcripts.

    Args:
        transcripts (Dict[str, Dict[str, Any]]): A dictionary where keys are transcript IDs.
        Each value is another dictionary with two keys:
            'info': {'seqid': str, 'strand': str} - Chromosome/contig and strand.
            'exons': List[Tuple[int, int]] - A list of (start, end) tuples
                     for each exon belonging to the transcript. Coordinates are integers.}

    Returns:
        Dict[str, Dict[str, Any]]: A dictionary where keys are transcript IDs.
        Each value is another dictionary with three keys:
            'info': {'seqid': str, 'strand': str} - Chromosome/contig and strand.
            'exons': List[Tuple[int, int]] - A list of (start, end) tuples
                     for each exon belonging to the transcript. Coordinates are integers.
            'introns': List[Tuple[int, int]] - A list of (start, end) tuples
                     for each intron between the defined exons. Coordinates are integers.
    """
    for transcript_id in transcripts:
        introns = []
        exons = transcripts[transcript_id]['exons']
        if len(exons) > 1:  # need 2 exons to define the intron between them
            for i in range(len(exons) - 1):
                exon1 = exons[i]
                exon2 = exons[i + 1]
                intron_start = exon1[1] + 1  # start of intron is 1 + end of exon
                intron_end = exon2[0] - 1  # end of intron is the start of next exon - 1
                if intron_start <= intron_end:
                    introns.append((intron_start, intron_end))
        transcripts[transcript_id]['introns'] = introns
    return transcripts


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


def extract_sequence(fasta: pysam.FastaFile, seq_id: str,
                     win_start: int, win_end: int) -> Optional[str]:
    """Obtain a sequence from an indexed FASTA file, given window coordinates.

    Args:
        fasta (pysam.FastaFile): pysam FASTA object
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
    if seq_id not in fasta.references:
        return ""
    seq_len = fasta.get_reference_length(seq_id)
    win_start = win_start - 1  # pysam takes 0-based coords
    win_end = min(win_end, seq_len)
    if win_start >= win_end:  # check after end is truncated also
        return ""
    seq = fasta.fetch(seq_id, win_start, win_end)
    return seq


def write_seq_to_file(seq: str, f: TextIO) -> None:
    """Write sequences to a file.

    Args:
        seq (str): sequence extracted from FASTA
        f (TextIO): Output file to write the seq
    """
    if seq:
        f.write(f"{seq}\n")


if __name__ == "__main__":
    main()
