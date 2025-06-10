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
import os
import random
import pysam
from Bio.Seq import Seq


def main() -> None:
    """
    Business logic
    """
    args = get_cli_args()
    if not os.path.exists(f"{args.fasta}.fai"):
        raise FileNotFoundError(
            f"Missing FASTA index: {args.fasta}.fai. Run samtools faidx."
        )
    if args.n_intron == 0 or args.n_exon == 0:
        raise ValueError("n_intron and n_exon must be nonzero!")
    transcripts = group_exons_by_transcript(args.gff)
    junctions = get_splice_junctions(transcripts)
    count = extract_positive_samples(
        junctions, args.fasta, args.out1, args.n_exon, args.n_intron
    )
    transcripts_with_introns = get_intron_coords(transcripts)
    extract_negative_samples(
        transcripts_with_introns,
        args.fasta,
        args.out2,
        args.n_exon + args.n_intron,  # make equal length to (+) seqs
        args.buffer,
        count,
    )  # sample the same number of (-) samples, if possible


def get_cli_args() -> argparse.Namespace:
    """
    Parse command line args: gff, fasta, n_exon, n_intron

    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Extract sequences around splice site junctions."
    )
    parser.add_argument(
        "-g",
        "--gff",
        required=True,
        type=str,
        metavar="FILE_PATH",
        help="Path to input GFF3 file.",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        type=str,
        metavar="FILE_PATH",
        help="Path to input FASTA file. Corresponding .fai must also exist.",
    )
    parser.add_argument(
        "-o1",
        "--out1",
        required=True,
        type=str,
        metavar="FILE_PATH",
        help="Path to output file for positive seqs.",
    )
    parser.add_argument(
        "-o2",
        "--out2",
        required=True,
        type=str,
        metavar="FILE_PATH",
        help="Path to output file for negative seqs",
    )
    parser.add_argument(
        "-ne",
        "--n_exon",
        default=40,
        type=int,
        metavar="INT",
        help="Number of bases to include in the exon region of the window.",
    )
    parser.add_argument(
        "-ni",
        "--n_intron",
        default=80,
        type=int,
        metavar="INT",
        help="Number of bases to include in the intron region of the window",
    )
    parser.add_argument(
        "-b",
        "--buffer",
        default=50,
        type=int,
        metavar="INT",
        help="Buffer size for intron region",
    )
    return parser.parse_args()


def extract_positive_samples(
    junctions: List[Dict[str, Any]],
    fasta_path: str,
    output_path: str,
    n_exon: int,
    n_intron: int,
) -> int:
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
    with open(output_path, "w", encoding="utf-8") as f:
        with pysam.FastaFile(fasta_path) as fasta:
            for junction in junctions:
                seqid = junction["seqid"]
                coord = junction["coord"]
                strand = junction["strand"]
                junc_type = junction["type"]
                win_start, win_end = get_window_coords(
                    strand, junc_type, coord, n_exon, n_intron
                )
                if win_start is not None and win_end is not None:
                    seq = extract_sequence(fasta, seqid, win_start, win_end)
                    if seq:
                        if strand == "-":
                            seq = str(
                                Seq(seq).reverse_complement()
                            )  # get RC for (-) strands
                        write_seq_to_file(
                            seq, seqid, junc_type, strand, win_start, win_end, f
                        )
                        count += 1
    return count


def extract_negative_samples(
    transcripts: Dict[str, Dict[str, Any]],
    fasta_path: str,
    output_path: str,
    window_size: int,
    buffer_size: int,
    num_samples: int,
) -> None:
    """
    Extracts negative sample sequences from intronic regions.

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
    random.seed(100)
    samples_written = 0
    with open(output_path, "w", encoding="utf-8") as f:
        with pysam.FastaFile(fasta_path) as fasta:
            transcript_ids = list(transcripts.keys())
            random.shuffle(transcript_ids)
            for transcript_id in transcript_ids:
                if samples_written >= num_samples:
                    break
                transcript_data = transcripts[transcript_id]
                seqid = transcript_data["info"]["seqid"]
                strand = transcript_data["info"]["strand"]
                introns = transcript_data.get("introns", [])
                for intron_start, intron_end in introns:
                    if samples_written >= num_samples:
                        break
                    # define a "safe start" region
                    # this makes sure we don't accidentally capture a splice jxn
                    safe_start = intron_start + buffer_size
                    safe_end = intron_end - buffer_size
                    if (
                        safe_start < safe_end
                        and (safe_end - safe_start + 1) >= window_size
                    ):
                        max_possible_start = safe_end - window_size + 1
                        if safe_start <= max_possible_start:
                            rand_coord = random.randint(safe_start, max_possible_start)
                            win_start = rand_coord
                            win_end = rand_coord + window_size - 1
                            seq = extract_sequence(fasta, seqid, win_start, win_end)
                            if seq:
                                if strand == "-":
                                    seq = str(Seq(seq).reverse_complement())
                                write_seq_to_file(
                                    seq, seqid, "intron", strand, win_start, win_end, f
                                )
                                samples_written += 1


def _parse_transcript_id(gff_str: str) -> Dict[str, str]:
    """Parse column 9 of a GFF3 file into a dictionary by separating k/v pairs
    by semicolons.

    Args:
        gff_str (str): 9th column of a GFF3 file

    Returns:
        Dict[str, str]: Dictionary in the structure:
            {"Parent": parent_id, "ID": feature_id, ...}
    """
    attrs = {}
    for part in filter(None, gff_str.strip().split(";")):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def group_exons_by_transcript(gff: str) -> Dict[str, Dict[str, Any]]:
    """
    Obtain all exons for a given transcript in a GFF3 annotation file.
    Example line format (tab-delineated):

    FR824046	ena	exon	1244	1570	.	-	.
    Parent=transcript:CCA13858;Name=CCA13858-1;constitutive=1;ensembl_end_phase=0;
    ensembl_phase=0;exon_id=CCA13858-1;rank=1

    Args:
        gff (str): GFF3 file of annotations

    Returns:
        transcripts (Dict[str, Dict[str, Any]]): A dictionary where keys are transcript IDs.
        Each value is another dictionary with two keys:
            'info': {'seqid': str, 'strand': str} - Chromosome/contig and strand.
            'exons': List[Tuple[int, int]] - A list of (start, end) tuples
                     for each exon belonging to the transcript. Coordinates are integers.
    """
    transcripts = defaultdict(lambda: {"info": {}, "exons": []})
    with open(gff, "r", encoding="utf-8") as f:
        for line in f:
            # Skip comment lines and empty lines
            if line.startswith("#") or not line.strip():
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:  # GFF3 requires 9 columns
                continue

            try:
                if fields[2] == "exon":
                    attrs = _parse_transcript_id(fields[8])
                    parent = attrs.get("Parent", "")
                    if parent.startswith("transcript:"):
                        transcript_id = parent.replace("transcript:", "")
                    else:
                        transcript_id = attrs.get("transcript_id")

                    if transcript_id is None:
                        continue

                    seqid, start, end, strand = [fields[i] for i in [0, 3, 4, 6]]
                    transcript_entry = transcripts[transcript_id]
                    transcript_entry["exons"].append((int(start), int(end)))

                    if "seqid" not in transcript_entry["info"]:
                        transcript_entry["info"]["seqid"] = seqid
                        transcript_entry["info"]["strand"] = strand
            except (IndexError, ValueError):
                continue

    return transcripts


def get_splice_junctions(
    transcripts: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Find coordinates (start/end) of splice junctions and whether
      they are donor or acceptor regions.
    Args:
        transcripts (Dict[str, Dict[str, Any]]): A dictionary where keys are transcript IDs.
        Each value is another dictionary with two keys:
            'info': {'seqid': str, 'strand': str} - Chromosome/contig and strand.
            'exons': List[Tuple[int, int]] - A list of (start, end) tuples
                     for each exon belonging to the transcript. Coordinates are integers.

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
        exon_index: int,
    ) -> None:
        junctions.append(
            {
                "id": f"{transcript_id_val}_{junction_type_val}_{exon_index}",
                "seqid": seqid_val,
                "coord": coord_val,
                "strand": strand_val,
                "type": junction_type_val,
            }
        )

    junctions = []
    for transcript_id, transcript_data in transcripts.items():
        transcript_info = transcript_data["info"]
        exons = transcript_data["exons"]
        exons_sorted = sorted(exons, key=lambda exon: exon[0])
        seqid = transcript_info["seqid"]
        strand = transcript_info["strand"]

        if len(exons_sorted) < 2:
            continue

        for i, exon_coords in enumerate(exons_sorted):
            exon_start, exon_end = exon_coords

            if i > 0:  # Acceptor sites
                if strand == "+":
                    acceptor_coord = exon_start
                elif strand == "-":
                    acceptor_coord = exon_end

                _add_junction_entry(
                    transcript_id, seqid, strand, acceptor_coord, "acceptor", i
                )
            if i < len(exons_sorted) - 1:  # Donor sites
                if strand == "+":
                    donor_coord = exon_end
                elif strand == "-":
                    donor_coord = exon_start
                _add_junction_entry(
                    transcript_id, seqid, strand, donor_coord, "donor", i
                )

    return junctions


def get_intron_coords(
    transcripts: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """Find the intron coordinates between all exons in the transcripts.

    Args:
        transcripts (Dict[str, Dict[str, Any]]): A dictionary where keys are transcript IDs.
        Each value is another dictionary with two keys:
            'info': {'seqid': str, 'strand': str} - Chromosome/contig and strand.
            'exons': List[Tuple[int, int]] - A list of (start, end) tuples
                     for each exon belonging to the transcript. Coordinates are integers.

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
        exons = sorted(transcripts[transcript_id]["exons"], key=lambda x: x[0])
        if len(exons) > 1:  # need 2 exons to define the intron between them
            for i in range(len(exons) - 1):
                exon1 = exons[i]
                exon2 = exons[i + 1]
                intron_start = exon1[1] + 1  # start of intron is 1 + end of exon
                intron_end = exon2[0] - 1  # end of intron is the start of next exon - 1
                if intron_start <= intron_end:
                    introns.append((intron_start, intron_end))
        transcripts[transcript_id]["introns"] = introns
    return transcripts


def get_window_coords(
    strand: int, junc_type: str, coord: int, n_exon: int, n_intron: int
) -> Tuple[Optional[int], Optional[int]]:
    """Calculates window for a sequence

    Args:
        junc_type (str): Donor or acceptor
        coord (int): coordinate of splice junction
        n_exon (int): Number of bases into exonic region
        n_intron (int): Number of bases into intronic region

    Returns:
        Tuple[Optional[int], Optional[int]]: Start/end coords of the window
    """
    win_start, win_end = None, None

    if strand == "+":
        if junc_type == "donor":
            win_start = coord - n_exon + 2  # +2 solves an off-by-1 error
            win_end = coord + n_intron + 1
        elif junc_type == "acceptor":
            win_start = coord - n_intron
            win_end = coord + n_exon - 1

    elif strand == "-":
        if junc_type == "donor":
            win_start = coord - n_intron - 1
            win_end = coord + n_exon - 2
        elif junc_type == "acceptor":
            win_start = coord - n_exon + 3
            win_end = coord + n_intron + 2

    return win_start, win_end


def extract_sequence(
    fasta: pysam.FastaFile, seq_id: str, win_start: int, win_end: int
) -> Optional[str]:
    """Obtain a sequence from an indexed FASTA file, given window coordinates.

    Args:
        fasta (pysam.FastaFile): pysam FASTA object
        seq_id (str): sequence ID for the desired sequence
        win_start (int): start coordinate (1-based inclusive)
        win_end (int): end coordinate (1-based inclusive)

    Returns:
        Optional[str]: Sequence of the desired window, or None if invalid
    """
    if win_start > win_end:
        return None
    if win_start < 1:
        win_start = 1
        if win_start > win_end:
            return None

    if seq_id not in fasta.references:
        return None

    seq_len = fasta.get_reference_length(seq_id)
    # Convert to 0-based coordinates for pysam
    win_start_0based = win_start - 1
    win_end_0based = min(win_end, seq_len)

    if win_start_0based >= win_end_0based:
        return None

    seq = fasta.fetch(seq_id, win_start_0based, win_end_0based)
    return seq if seq else None


def write_seq_to_file(
    seq: str,
    seqid: str,
    junc_type: str,
    strand: str,
    win_start: int,
    win_end: int,
    f: TextIO,
) -> None:
    """Write sequences to a file in FASTA format.

    Args:
        seq (str): sequence extracted from FASTA
        seqid (str): sequence/chromosome ID
        junc_type (str): type of junction (donor/acceptor/intron)
        strand (str): strand information
        win_start (int): window start coordinate
        win_end (int): window end coordinate
        f (TextIO): Output file handle
    """
    if seq:
        f.write(f">{seqid}_{junc_type}_{strand}_{win_start}_{win_end}\n")
        f.write(f"{seq}\n")


if __name__ == "__main__":
    main()
