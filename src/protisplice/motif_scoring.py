"""
File: motif_scoring.py
Description: Generate PFM/PWM/PPM from aligned DNA sequences
"""

from typing import List, Dict
from Bio import motifs
from Bio.Seq import Seq
import pandas as pd
import numpy as np


def generate_ppm(seqs: List[Seq]) -> pd.DataFrame:
    """Generate a position probability matrix from a list
    of aligned sequences

    Args:
        seqs (List[Seq]): List of aligned sequences (Bio.Seq objects).

    Returns:
        pd.DataFrame: Pandas dataframe containing the probability for each sequence
    """
    motif = motifs.create(seqs)
    counts_df = pd.DataFrame.from_dict(motif.counts, orient="index")
    ppm = counts_df.div(counts_df.sum(axis=0), axis=1)
    return ppm


def generate_pfm(seqs: List[Seq]) -> pd.DataFrame:
    """Generate a position frequency matrix from a list
    of aligned sequences

    Args:
        seqs (List[Seq]): List of aligned sequences (Bio.Seq objects).

    Returns:
        pd.DataFrame: Pandas dataframe containing the probability of each nucleotide
        at each position
    """
    motif = motifs.create(seqs)
    pfm = pd.DataFrame.from_dict(motif.counts, orient="index")
    return pfm.astype("int64")


def generate_pwm(
    seqs: List[Seq],
    background_freq: Dict[str, float] = None,
) -> pd.DataFrame:
    """Generate a position weight matrix from a list of aligned sequences.

    Args:
        seqs (List[Seq]): List of aligned sequences (Bio.Seq objects).
        background_freq (Dict[str, float], optional): Background frequency of the genome/chromosome
        of origin of the sequences. If None is passed in, equal probabilities of 0.25 will be used
        for all nucleotides.
        epsilon (float, optional): Pseudocount to prevent divide by 0 errors. Defaults to 1e-10.

    Returns:
        pd.DataFrame: Pandas dataframe containing the log-odds for each nucleotide at each position.
    """
    motif = motifs.create(seqs)
    counts_df = pd.DataFrame.from_dict(motif.counts, orient="index")

    # add pseudocounts only if a column sum is zero
    col_sums = counts_df.sum(axis=0)
    col_sums = col_sums.replace(0, np.nan)  # mark zeros as NaN
    ppm = counts_df.div(col_sums, axis=1)

    # fill any NaN columns (where sum was 0) with uniform distribution
    ppm = ppm.fillna(1.0 / len(counts_df))

    if background_freq is None:
        background_freq = {base: 0.25 for base in "ACGT"}

    pwm = ppm.div(pd.Series(background_freq), axis=0).map(np.log2)

    return pwm
