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
    counts_df = pd.DataFrame(motif.counts)
    ppm = counts_df.div(counts_df.sum(axis=0), axis=1)
    return ppm


def generate_pfm(seqs: List[Seq]) -> pd.DataFrame:
    """Generate a position frequency matrix from a list
    of aligned sequences

    Args:
        seqs (List[Seq]): List of aligned sequences (Bio.Seq objects).

    Returns:
        pd.DataFrame: Pandas dataframe containing the probability for each sequence
    """
    motif = motifs.create(seqs)
    pfm = pd.DataFrame(motif.counts)
    return pfm


def generate_pwm(
    seqs: List[Seq],
    background_freq: Dict[str, float] = None,
    epsilon: float = 1e-10,
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
    if background_freq is None:
        background_freq = {'A': 0.25, 'G': 0.25, 'C': 0.25, 'T': 0.25}
    bg = pd.Series(background_freq)
    ppm = generate_ppm(seqs)
    if not set(bg.index).issubset(set(ppm.index)):
        ppm = ppm.T
    ppm = ppm.loc[bg.index]
    pwm = (ppm + epsilon).divide(bg, axis=0).apply(np.log2)
    return pwm
