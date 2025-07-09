"""
File: splice_site_model.py
Description: Machine learning model for splice site prediction using extracted sequences
"""

from typing import Dict, Tuple, List
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from scipy.sparse import csr_matrix


class FeatureExtractor:

    def __init__(self):
        pass

    def extract_all(self, seq: str) -> Dict:
        features = {
            **self.get_base_content(seq),
            **self.get_kmers(seq, k=3),
            **self.one_hot_encode(seq)
        }
        return features

    @staticmethod
    def get_kmers(seq: str, k: int) -> List:
        return [seq[i:i + k] for i in range(len(seq) - k + 1)]

    @staticmethod
    def one_hot_encode(seq: str) -> np.ndarray:
        base_map = {
            "A": [1, 0, 0, 0],
            "G": [0, 1, 0, 0],
            "T": [0, 0, 1, 0],
            "C": [0, 0, 0, 1],
            "N": [0, 0, 0, 0],
        }
        encoded_seq = []

        for base in seq:
            encoded_seq.append(base_map[base])

        return np.array(encoded_seq)

    @staticmethod
    def get_base_content(seq: str) -> Dict:
        return {
            "A": (seq.count("A") / len(seq)),
            "G": (seq.count("G") / len(seq)),
            "T": (seq.count("T") / len(seq)),
            "C": (seq.count("C") / len(seq)),
        }


class SpliceModel:

    def __init__(self, features):

        pass

    def train_model(self, features: List[str], labels: List[str]):
        X_train, X_val, y_train, y_val = train_test_split(features, labels, test_size=0.2, random_state=100)


