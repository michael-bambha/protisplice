"""
File: splice_sites_protists.py
Author: Michael Bambha
Contact: bambha.m@northeastern.edu
Description: A python script for an ML workflow predicting alternative splice sites in protists.
"""
import argparse
from collections import defaultdict
import re
from typing import Dict, Any


def get_cli_args():
    """
    Generate command line args.
    """
    pass


def parse_transcript_id(gtf_str: str) -> Dict[str, str]:
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
                transcript_id = attrs["transcript_id"]
                seqid, start, end, strand = [fields[i] for i in [0, 3, 4, 6]]
                transcript_entry = transcripts[transcript_id]
                transcript_entry['exons'].append((int(start), int(end)))
                if not transcript_entry['info']:
                    transcript_entry['info']['seqid'] = seqid
                    transcript_entry['info']['strand'] = strand
    return transcripts


def get_splice_junctions(transcripts):
    for transcript in transcripts:
        exons = sorted()
    pass

print(group_exons_by_transcript("transcripts.gtf"))
