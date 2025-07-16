"""
File: splice_site_model.py
Description: Machine learning model for splice site prediction using extracted sequences
"""
# pylint: disable=invalid-name
from typing import Dict, Tuple, List
from collections import Counter
from itertools import product
import numpy as np
from Bio import SeqIO
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score


class FeatureExtractor:
    """Finds features for DNA sequences like kmers and base content"""

    def __init__(self, kmer_k: int = 3, flatten_one_hot: bool = True):
        self.kmer_k = kmer_k
        self.flatten_one_hot = flatten_one_hot

    @staticmethod
    def load_data(
        pos_fasta: str, neg_fasta: str
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Load data and create labels

        Args:
            pos_fasta (str): Path to file containing true splice sites
            neg_fasta (str): Path to file containing negative splice sites

        Returns:
            Tuple[np.ndarray, np.ndarray)]: Tuple of sequences, labels
        """
        pos_seqs = [str(record.seq) for record in SeqIO.parse(pos_fasta, "fasta")]
        neg_seqs = [str(record.seq) for record in SeqIO.parse(neg_fasta, "fasta")]
        sequences = np.array(pos_seqs + neg_seqs)
        labels = np.array([1] * len(pos_seqs) + [0] * len(neg_seqs))
        return sequences, labels

    def extract_all(self, seq: str) -> np.ndarray:
        """Combines all feature extraction methods for simpler workflow

        Args:
            seq (str): seq to extract

        Returns:
            np.ndarray: Numpy array of features from the sequence
        """
        base_content = list(self.get_base_content(seq).values())
        kmer_counts = list(self.count_kmers(seq, self.kmer_k).values())
        one_hot = self.one_hot_encode(seq, flatten=self.flatten_one_hot)
        return np.concatenate([base_content, kmer_counts, one_hot])

    def batch_extract_all(self, sequences: List[str]) -> np.ndarray:
        """Run all feature extraction methods on a list of sequences

        Args:
            sequences (List[str]): List of DNA sequences

        Returns:
            np.ndarray: Numpy array of features from each sequence
        """
        return np.array([self.extract_all(seq) for seq in sequences])

    @staticmethod
    def count_kmers(seq: str, k: int) -> Dict[str, int]:
        """Count kmers of a given k in a DNA sequence

        Args:
            seq (str): Sequence to count
            k (int): length k for kmer

        Returns:
            Dict[str, int]: Count dictionary of kmers from the input seq
        """
        all_kmers = [''.join(p) for p in product("ACGT", repeat=k)]
        kmers = [seq[i: i + k] for i in range(len(seq) - k + 1)]
        total = len(kmers)
        counts = Counter(kmers)
        return {kmer: counts.get(kmer, 0) / total for kmer in all_kmers}

    @staticmethod
    def one_hot_encode(seq: str, flatten: bool = True) -> np.ndarray:
        """One-hot encode a DNA sequence
        Args:
            seq (str): DNA sequence
            flatten (bool, optional): Flattens array with np.flatten() if true. Defaults to True.

        Returns:
            np.ndarray: One-hot encoded sequence
        """
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
        arr = np.array(encoded_seq)
        return arr.flatten() if flatten else arr

    @staticmethod
    def get_base_content(seq: str) -> Dict:
        """Count base percentages in the sequence

        Args:
            seq (str): DNA sequence to analyze

        Returns:
            Dict: Count dictionary of base: % composition
        """
        return {
            "A": (seq.count("A") / len(seq)),
            "G": (seq.count("G") / len(seq)),
            "T": (seq.count("T") / len(seq)),
            "C": (seq.count("C") / len(seq)),
        }


class SpliceModel:
    """Scikit-Learn Wrapper Class for ML Workflow"""

    def __init__(
        self,
        test_size: int,
        random_state: int,
        model_cls=LogisticRegression,
        model_params: Dict = None,
    ):
        self.test_size = test_size
        self.random_state = random_state
        self.model_params = model_params
        self.model_cls = model_cls
        self.model = None

    def train_model(self, features: np.ndarray, labels: List[int]):
        """Sklearn wrapper to train any sklearn model

        Args:
            features (np.ndarray): np.array of data
            labels (List[int]): 1/0 for true/false splice sites
        """
        X_train, X_val, y_train, y_val = train_test_split(
            features, labels, test_size=self.test_size, random_state=self.random_state
        )
        model = self.model_cls(**(self.model_params or {}))
        self.model = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', model)
        ])
        self.model.fit(X_train, y_train)
        return X_train, X_val, y_train, y_val

    def predict(self, X_test: np.ndarray) -> List[int]:
        """Sklearn wrapper for model prediction

        Args:
            X_test (np.ndarray): feature test set (no labels)

        Raises:
            ValueError: If model is not trained / does not exist.

        Returns:
            List[int]: List of classification predictions (0/1)
        """
        if not self.model:
            raise ValueError("Model not trained - run train_model first.")
        return self.model.predict(X_test)

    def predict_proba(self, X_test: np.ndarray) -> List[int]:
        """Sklearn wrapper for predict_proba
        Args:
            X_test (np.ndarray): Feature test set

        Raises:
            ValueError: If model is not trained, does not exist
            ValueError: If model does not support predict_proba

        Returns:
            List[int]: List of probabilities
        """
        if not self.model:
            raise ValueError("Model not trained - run train_model first.")
        if not hasattr(self.model, "predict_proba"):
            raise ValueError("This model does not support predict_proba.")
        return self.model.predict_proba(X_test)

    def evaluate(
        self, features: np.ndarray, labels: List[str], metrics: List[str] = ["accuracy"]
    ) -> Dict[str, float]:
        """Sklearn wrapper for evaluating a model given a feature set, labels, and metrics.
        Supported metrics are accuracy, f1, and roc_auc.

        Args:
            features (np.ndarray): Features to evaluate
            labels (List[str]): Labels of features passed in
            metrics (List[str], optional): accuracy, f1, or roc_auc. Defaults to ["accuracy"].

        Raises:
            ValueError: If model is not trained or does not exist
            ValueError: If roc_auc is called on a model that does not support predict_proba

        Returns:
            Dict[str, float]: Dictionary of metric: score
        """
        if self.model is None:
            raise ValueError("Model has not been trained. Call train_model() first.")

        preds = self.model.predict(features)
        results = {}
        for metric in metrics:
            if metric == "accuracy":
                results["accuracy"] = accuracy_score(labels, preds)
            elif metric == "f1":
                results["f1"] = f1_score(labels, preds)
            elif metric == "roc_auc":
                if hasattr(self.model, "predict_proba"):
                    probs = self.model.predict_proba(features)[:, 1]
                    results["roc_auc"] = roc_auc_score(labels, probs)
            else:
                raise ValueError(f"Unsupported metric: {metric}. ")

        return results
