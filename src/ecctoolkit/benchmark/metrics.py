"""
Metrics data structures for benchmark evaluation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class TypeMetrics:
    """Metrics for a specific eccDNA type (U/M/C)."""

    ecc_type: str  # "UeccDNA", "MeccDNA", "CeccDNA"
    truth_count: int = 0
    detected_count: int = 0
    true_positive: int = 0
    false_positive: int = 0
    false_negative: int = 0
    misclassified_as: dict[str, int] = field(default_factory=dict)
    misclassified_from: dict[str, int] = field(default_factory=dict)

    @property
    def precision(self) -> float:
        """Precision = TP / (TP + FP)"""
        if self.detected_count == 0:
            return 0.0
        return self.true_positive / self.detected_count

    @property
    def recall(self) -> float:
        """Recall = TP / (TP + FN) = TP / truth_count"""
        if self.truth_count == 0:
            return 0.0
        return self.true_positive / self.truth_count

    @property
    def f1_score(self) -> float:
        """F1 = 2 * (precision * recall) / (precision + recall)"""
        p, r = self.precision, self.recall
        if p + r == 0:
            return 0.0
        return 2 * p * r / (p + r)

    @property
    def false_positive_rate(self) -> float:
        """FPR = FP / detected_count"""
        if self.detected_count == 0:
            return 0.0
        return self.false_positive / self.detected_count


@dataclass
class StateMetrics:
    """Metrics for detection state (Confirmed/Inferred)."""

    state: str  # "Confirmed", "Inferred"
    total_count: int = 0
    matched_count: int = 0
    unmatched_count: int = 0

    @property
    def precision(self) -> float:
        if self.total_count == 0:
            return 0.0
        return self.matched_count / self.total_count

    @property
    def false_positive_rate(self) -> float:
        if self.total_count == 0:
            return 0.0
        return self.unmatched_count / self.total_count


@dataclass
class LengthBinMetrics:
    """Metrics for a length bin."""

    bin_label: str  # e.g., "<500bp", "500-1000bp"
    bin_min: int
    bin_max: int
    truth_count: int = 0
    detected_count: int = 0
    true_positive: int = 0

    @property
    def recall(self) -> float:
        if self.truth_count == 0:
            return 0.0
        return self.true_positive / self.truth_count


@dataclass
class ConfusionMatrix:
    """Confusion matrix for eccDNA classification."""

    # Rows: predicted, Columns: actual
    # matrix[predicted][actual] = count
    matrix: dict[str, dict[str, int]] = field(default_factory=dict)
    types: list[str] = field(default_factory=lambda: ["UeccDNA", "MeccDNA", "CeccDNA"])

    def __post_init__(self):
        for t in self.types:
            if t not in self.matrix:
                self.matrix[t] = {tt: 0 for tt in self.types}

    def add(self, predicted: str, actual: str, count: int = 1):
        """Add count to matrix[predicted][actual]."""
        if predicted not in self.matrix:
            self.matrix[predicted] = {t: 0 for t in self.types}
        if actual not in self.matrix[predicted]:
            self.matrix[predicted][actual] = 0
        self.matrix[predicted][actual] += count

    def to_dataframe(self):
        """Convert to pandas DataFrame."""
        import pandas as pd

        return pd.DataFrame(self.matrix).T.reindex(
            index=self.types, columns=self.types, fill_value=0
        )


@dataclass
class BenchmarkMetrics:
    """Complete benchmark metrics."""

    # Overall stats
    total_truth: int = 0
    total_detected: int = 0
    total_true_positive: int = 0
    total_false_positive: int = 0
    total_false_negative: int = 0
    total_misclassified: int = 0

    # By type
    by_type: dict[str, TypeMetrics] = field(default_factory=dict)

    # By state
    by_state: dict[str, StateMetrics] = field(default_factory=dict)

    # By length
    by_length: dict[str, LengthBinMetrics] = field(default_factory=dict)

    # Confusion matrix
    confusion_matrix: ConfusionMatrix = field(default_factory=ConfusionMatrix)

    # Matched/unmatched lists for detailed analysis
    matched_truth_ids: set[str] = field(default_factory=set)
    matched_detected_ids: set[str] = field(default_factory=set)
    unmatched_truth_ids: set[str] = field(default_factory=set)
    unmatched_detected_ids: set[str] = field(default_factory=set)

    @property
    def overall_precision(self) -> float:
        if self.total_detected == 0:
            return 0.0
        return self.total_true_positive / self.total_detected

    @property
    def overall_recall(self) -> float:
        if self.total_truth == 0:
            return 0.0
        return self.total_true_positive / self.total_truth

    @property
    def overall_f1(self) -> float:
        p, r = self.overall_precision, self.overall_recall
        if p + r == 0:
            return 0.0
        return 2 * p * r / (p + r)

    def init_types(self, types: Optional[list[str]] = None):
        """Initialize type metrics."""
        if types is None:
            types = ["UeccDNA", "MeccDNA", "CeccDNA"]
        for t in types:
            if t not in self.by_type:
                self.by_type[t] = TypeMetrics(ecc_type=t)

    def init_states(self, states: Optional[list[str]] = None):
        """Initialize state metrics."""
        if states is None:
            states = ["Confirmed", "Inferred"]
        for s in states:
            if s not in self.by_state:
                self.by_state[s] = StateMetrics(state=s)

    def init_length_bins(
        self,
        bins: Optional[list[tuple[str, int, int]]] = None,
    ):
        """Initialize length bin metrics."""
        if bins is None:
            bins = [
                ("<500bp", 0, 500),
                ("500bp-1kb", 500, 1000),
                ("1kb-5kb", 1000, 5000),
                ("5kb-50kb", 5000, 50000),
                (">50kb", 50000, float("inf")),
            ]
        for label, min_len, max_len in bins:
            if label not in self.by_length:
                self.by_length[label] = LengthBinMetrics(
                    bin_label=label, bin_min=min_len, bin_max=max_len
                )

    def summary_dict(self) -> dict:
        """Return summary as dictionary."""
        return {
            "total_truth": self.total_truth,
            "total_detected": self.total_detected,
            "true_positive": self.total_true_positive,
            "false_positive": self.total_false_positive,
            "false_negative": self.total_false_negative,
            "misclassified": self.total_misclassified,
            "precision": self.overall_precision,
            "recall": self.overall_recall,
            "f1_score": self.overall_f1,
        }
