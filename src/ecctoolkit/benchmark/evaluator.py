"""
Core evaluation logic for benchmark analysis.
"""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from itertools import permutations
from pathlib import Path
from typing import Optional

import pandas as pd

from ecctoolkit.benchmark.metrics import (
    BenchmarkMetrics,
    ConfusionMatrix,
    LengthBinMetrics,
    StateMetrics,
    TypeMetrics,
)

logger = logging.getLogger(__name__)


class BenchmarkEvaluator:
    """Evaluator for comparing detection results against ground truth."""

    def __init__(
        self,
        overlap_threshold: float = 0.8,
        length_bins: Optional[list[tuple[str, int, int]]] = None,
    ):
        """
        Initialize evaluator.

        Args:
            overlap_threshold: Minimum overlap ratio for matching (default: 0.8)
            length_bins: List of (label, min_len, max_len) tuples
        """
        self.overlap_threshold = overlap_threshold
        self.length_bins = length_bins or [
            ("<500bp", 0, 500),
            ("500bp-1kb", 500, 1000),
            ("1kb-5kb", 1000, 5000),
            ("5kb-50kb", 5000, 50000),
            (">50kb", 50000, float("inf")),
        ]

        # Data storage
        self.truth_u: Optional[pd.DataFrame] = None
        self.truth_m: Optional[pd.DataFrame] = None
        self.truth_c: Optional[pd.DataFrame] = None
        self.detected: Optional[pd.DataFrame] = None

        # Indices for fast lookup
        self._truth_index: dict[str, list[tuple]] = defaultdict(list)
        # Separate index for chimeric truth: name -> (segments, length, fusion_type)
        self._chimeric_truth_index: dict[str, tuple[list[tuple[str, int, int]], int, str]] = {}
        self._last_matched_truth: Optional[set[str]] = None
        self._last_unmatched_detected: Optional[set[str]] = None

    def load_truth(
        self,
        unique_bed: Optional[str] = None,
        multi_bed: Optional[str] = None,
        chimeric_bed: Optional[str] = None,
        truth_dir: Optional[str] = None,
    ):
        """
        Load ground truth BED files.

        Args:
            unique_bed: Path to unique.bed file
            multi_bed: Path to multi.bed file
            chimeric_bed: Path to chimeric.bed file
            truth_dir: Directory containing *.unique.bed, *.multi.bed, *.chimeric.bed
        """
        if truth_dir:
            truth_path = Path(truth_dir)
            # Find BED files by pattern
            unique_files = list(truth_path.glob("*.unique.bed"))
            multi_files = list(truth_path.glob("*.multi.bed"))
            chimeric_files = list(truth_path.glob("*.chimeric.bed"))

            if unique_files:
                unique_bed = str(unique_files[0])
            if multi_files:
                multi_bed = str(multi_files[0])
            if chimeric_files:
                chimeric_bed = str(chimeric_files[0])

        bed_cols = ["chrom", "start", "end", "name", "length", "strand", "type", "fragments"]

        def _load_bed(path: str, label: str) -> pd.DataFrame:
            logger.info(f"Loading {label} truth from: {path}")
            df = pd.read_csv(path, sep="\t", comment="#", names=bed_cols, header=None)
            if df.empty:
                logger.warning(f"{label} truth BED is empty: {path}")
                return df
            if df.iloc[0]["chrom"] == "chrom":
                df = df.iloc[1:].reset_index(drop=True)
            if df.empty:
                logger.warning(f"{label} truth BED has header only: {path}")
                return df
            df["length"] = df["length"].astype(int)
            df["start"] = df["start"].astype(int)
            df["end"] = df["end"].astype(int)
            return df

        if unique_bed:
            self.truth_u = _load_bed(unique_bed, "UeccDNA")

        if multi_bed:
            self.truth_m = _load_bed(multi_bed, "MeccDNA")

        if chimeric_bed:
            self.truth_c = _load_bed(chimeric_bed, "CeccDNA")

        # Build index
        self._build_truth_index()

    def _build_truth_index(self):
        """Build spatial index for fast overlap queries."""
        self._truth_index.clear()
        self._chimeric_truth_index.clear()

        for df, etype in [
            (self.truth_u, "U"),
            (self.truth_m, "M"),
            (self.truth_c, "C"),
        ]:
            if df is None:
                continue
            for _, row in df.iterrows():
                if etype == "C":
                    # Build chimeric-specific index from fragments column
                    segments = self._parse_truth_segments(str(row["fragments"]))
                    if segments:
                        fusion_type = self._classify_chimeric_segments(segments)
                        self._chimeric_truth_index[row["name"]] = (
                            segments, row["length"], fusion_type
                        )
                    # Also add primary region to general index for backward compat
                self._truth_index[row["chrom"]].append(
                    (row["start"], row["end"], row["name"], etype, row["length"])
                )

    def load_detected(self, result_path: str):
        """
        Load detection results (CircleSeeker merged_output.csv).

        Args:
            result_path: Path to merged_output.csv
        """
        logger.info(f"Loading detection results from: {result_path}")
        self.detected = pd.read_csv(result_path)

    @staticmethod
    def _parse_regions(region_str: str) -> list[tuple[str, int, int]]:
        """Parse Regions column to list of (chrom, start, end).

        Handles formats:
        - chr1:100-200
        - chr1:100-200(+)  (with strand)
        - Multiple segments separated by ; or |
        """
        regions = []
        region_str = str(region_str)

        # Handle both ; and | separators
        for sep in [";", "|"]:
            if sep in region_str:
                parts = region_str.split(sep)
                break
        else:
            parts = [region_str]

        for r in parts:
            r = r.strip()
            if ":" in r and "-" in r:
                # Use regex to handle optional strand suffix
                m = re.match(r"([^:]+):(\d+)-(\d+)", r)
                if m:
                    regions.append((m.group(1), int(m.group(2)), int(m.group(3))))
        return regions

    @staticmethod
    def _overlap_ratio(start1: int, end1: int, start2: int, end2: int) -> float:
        """Calculate overlap ratio between two intervals."""
        overlap = max(0, min(end1, end2) - max(start1, start2))
        len1 = end1 - start1
        len2 = end2 - start2
        if len1 == 0 or len2 == 0:
            return 0.0
        return overlap / min(len1, len2)

    @staticmethod
    def _parse_truth_segments(fragments_str: str) -> list[tuple[str, int, int]]:
        """Parse truth fragments column (e.g., 'Chr2:100-200;Chr1:300-400')."""
        segments = []
        if not fragments_str or fragments_str in (".", "nan"):
            return segments
        for seg in fragments_str.split(";"):
            seg = seg.strip()
            # Match format: chr:start-end (no strand info in truth)
            m = re.match(r"([^:]+):(\d+)-(\d+)", seg)
            if m:
                segments.append((m.group(1), int(m.group(2)), int(m.group(3))))
        return segments

    @staticmethod
    def _classify_chimeric_segments(segments: list[tuple[str, int, int]]) -> str:
        """Classify chimeric segments as inter- or intra-chromosomal."""
        if len(segments) < 2:
            return "unknown"
        chroms = set(s[0] for s in segments)
        return "inter-chromosomal" if len(chroms) > 1 else "intra-chromosomal"

    def _match_chimeric_pair(
        self,
        truth_segments: list[tuple[str, int, int]],
        detected_segments: list[tuple[str, int, int]],
    ) -> tuple[bool, int]:
        """
        Check if two CeccDNA match at segment level.

        Tries all permutations of detected segments to find the best pairing
        with truth segments. Each paired segment must have same chromosome
        and overlap >= threshold.

        Args:
            truth_segments: [(chrom, start, end), ...]
            detected_segments: [(chrom, start, end), ...]

        Returns:
            (is_match, matched_segment_count)
        """
        n_truth = len(truth_segments)
        n_det = len(detected_segments)

        if n_truth == 0 or n_det == 0:
            return False, 0

        # Allow segment count difference of at most 1
        if abs(n_truth - n_det) > 1:
            return False, 0

        # Try all permutations of the smaller list against the larger
        # to find the best pairing
        best_matched = 0

        if n_det <= n_truth:
            # Try permutations of detected against truth (in order)
            for perm in permutations(range(n_det)):
                matched = 0
                for t_idx in range(min(n_truth, n_det)):
                    d_idx = perm[t_idx]
                    t_chr, t_start, t_end = truth_segments[t_idx]
                    d_chr, d_start, d_end = detected_segments[d_idx]
                    if t_chr == d_chr:
                        ratio = self._overlap_ratio(t_start, t_end, d_start, d_end)
                        if ratio >= self.overlap_threshold:
                            matched += 1
                best_matched = max(best_matched, matched)
        else:
            # More detected than truth - try permutations of truth against detected
            for perm in permutations(range(n_truth)):
                matched = 0
                for d_idx in range(min(n_truth, n_det)):
                    t_idx = perm[d_idx]
                    t_chr, t_start, t_end = truth_segments[t_idx]
                    d_chr, d_start, d_end = detected_segments[d_idx]
                    if t_chr == d_chr:
                        ratio = self._overlap_ratio(t_start, t_end, d_start, d_end)
                        if ratio >= self.overlap_threshold:
                            matched += 1
                best_matched = max(best_matched, matched)

        # Full match: all segments in the smaller set are matched
        min_segs = min(n_truth, n_det)
        is_match = best_matched == min_segs and best_matched > 0

        return is_match, best_matched

    def _find_chimeric_truth_match(
        self, detected_segments: list[tuple[str, int, int]]
    ) -> Optional[tuple[str, int, str, int, bool]]:
        """
        Find matching chimeric truth entry for detected segments.

        Returns:
            (truth_name, truth_length, fusion_type, matched_segments, is_partial)
            or None
        """
        best_match = None
        best_matched_count = 0

        for t_name, (t_segments, t_length, fusion_type) in self._chimeric_truth_index.items():
            is_match, matched_count = self._match_chimeric_pair(t_segments, detected_segments)

            if is_match:
                if matched_count > best_matched_count:
                    best_matched_count = matched_count
                    best_match = (t_name, t_length, fusion_type, matched_count, False)
            elif matched_count > 0 and best_match is None:
                # Partial match (some segments matched but not all)
                best_match = (t_name, t_length, fusion_type, matched_count, True)

        return best_match

    def _find_truth_match(
        self, regions: list[tuple[str, int, int]]
    ) -> Optional[tuple[str, str, int]]:
        """
        Find matching truth entry for given regions.

        Returns:
            (truth_name, truth_type, truth_length) or None
        """
        for chrom, start, end in regions:
            for t_start, t_end, t_name, t_type, t_length in self._truth_index.get(chrom, []):
                ratio = self._overlap_ratio(start, end, t_start, t_end)
                if ratio >= self.overlap_threshold:
                    return (t_name, t_type, t_length)
        return None

    def _get_length_bin(self, length: int) -> str:
        """Get length bin label for a given length."""
        for label, min_len, max_len in self.length_bins:
            if min_len <= length < max_len:
                return label
        return self.length_bins[-1][0]  # Last bin

    def evaluate(self) -> BenchmarkMetrics:
        """
        Run evaluation and return metrics.

        Returns:
            BenchmarkMetrics with all computed metrics
        """
        if self.detected is None:
            raise ValueError("No detection results loaded. Call load_detected() first.")

        metrics = BenchmarkMetrics()
        metrics.init_types()
        metrics.init_states()
        metrics.init_length_bins(self.length_bins)

        # Count truth
        truth_counts = {"UeccDNA": 0, "MeccDNA": 0, "CeccDNA": 0}
        truth_by_length: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))

        for df, etype in [
            (self.truth_u, "UeccDNA"),
            (self.truth_m, "MeccDNA"),
            (self.truth_c, "CeccDNA"),
        ]:
            if df is not None:
                truth_counts[etype] = len(df)
                metrics.by_type[etype].truth_count = len(df)
                metrics.total_truth += len(df)

                # Count by length bin
                for _, row in df.iterrows():
                    bin_label = self._get_length_bin(row["length"])
                    truth_by_length[etype][bin_label] += 1
                    if bin_label in metrics.by_length:
                        metrics.by_length[bin_label].truth_count += 1

        # Count CeccDNA inter/intra chromosomal truth and total segments
        cecc_metrics = metrics.by_type["CeccDNA"]
        for t_name, (t_segments, t_length, fusion_type) in self._chimeric_truth_index.items():
            cecc_metrics.total_truth_segments += len(t_segments)
            if fusion_type == "inter-chromosomal":
                cecc_metrics.inter_chr_truth += 1
            elif fusion_type == "intra-chromosomal":
                cecc_metrics.intra_chr_truth += 1

        # Track matched truth IDs
        matched_truth: set[str] = set()

        # Process each detection
        for _, row in self.detected.iterrows():
            ecc_id = row["eccDNA_id"]
            ecc_type = row["eccDNA_type"]
            state = row["State"]
            length = row["Length"]

            # Update detected counts
            metrics.total_detected += 1
            if ecc_type in metrics.by_type:
                metrics.by_type[ecc_type].detected_count += 1
            if state in metrics.by_state:
                metrics.by_state[state].total_count += 1

            # Find truth match
            regions = self._parse_regions(row["Regions"])

            # Use dedicated chimeric matching for CeccDNA
            if ecc_type == "CeccDNA" and self._chimeric_truth_index:
                chimeric_match = self._find_chimeric_truth_match(regions)
                if chimeric_match:
                    truth_name, truth_length, fusion_type, matched_segs, is_partial = chimeric_match

                    # Classify detected CeccDNA fusion type
                    det_fusion = self._classify_chimeric_segments(regions)
                    if det_fusion == "inter-chromosomal":
                        cecc_metrics.inter_chr_detected += 1
                    elif det_fusion == "intra-chromosomal":
                        cecc_metrics.intra_chr_detected += 1

                    if is_partial:
                        # Partial match: some segments matched but not all
                        cecc_metrics.partial_match_count += 1
                        cecc_metrics.total_matched_segments += matched_segs
                        # Treat partial match as false positive
                        metrics.total_false_positive += 1
                        metrics.by_type[ecc_type].false_positive += 1
                        metrics.unmatched_detected_ids.add(ecc_id)
                        if state in metrics.by_state:
                            metrics.by_state[state].unmatched_count += 1
                        logger.debug(
                            f"CeccDNA {ecc_id} partial match to {truth_name}: "
                            f"{matched_segs} segments matched"
                        )
                    else:
                        # Full match
                        matched_truth.add(truth_name)
                        metrics.matched_detected_ids.add(ecc_id)
                        cecc_metrics.total_matched_segments += matched_segs

                        # Correct type classification (detected as CeccDNA, truth is CeccDNA)
                        metrics.total_true_positive += 1
                        metrics.by_type[ecc_type].true_positive += 1
                        if state in metrics.by_state:
                            metrics.by_state[state].matched_count += 1

                        # Update inter/intra chr TP
                        if fusion_type == "inter-chromosomal":
                            cecc_metrics.inter_chr_tp += 1
                        elif fusion_type == "intra-chromosomal":
                            cecc_metrics.intra_chr_tp += 1

                        # Update length bin
                        bin_label = self._get_length_bin(truth_length)
                        if bin_label in metrics.by_length:
                            metrics.by_length[bin_label].true_positive += 1
                            metrics.by_length[bin_label].detected_count += 1

                        # Update confusion matrix diagonal
                        # (will be added in bulk later, skip here)
                else:
                    # No chimeric match found; try general index as fallback
                    # (handles case where CeccDNA is misclassified as U/M)
                    general_match = self._find_truth_match(regions)
                    if general_match:
                        truth_name, truth_type, truth_length = general_match
                        type_map = {"U": "UeccDNA", "M": "MeccDNA", "C": "CeccDNA"}
                        actual_type = type_map.get(truth_type, truth_type)

                        if actual_type != "CeccDNA":
                            # Detected as CeccDNA but truth is U or M
                            matched_truth.add(truth_name)
                            metrics.matched_detected_ids.add(ecc_id)
                            metrics.total_misclassified += 1
                            metrics.by_type[ecc_type].false_positive += 1

                            if ecc_type in metrics.by_type:
                                if actual_type not in metrics.by_type[ecc_type].misclassified_from:
                                    metrics.by_type[ecc_type].misclassified_from[actual_type] = 0
                                metrics.by_type[ecc_type].misclassified_from[actual_type] += 1

                            if actual_type in metrics.by_type:
                                if ecc_type not in metrics.by_type[actual_type].misclassified_as:
                                    metrics.by_type[actual_type].misclassified_as[ecc_type] = 0
                                metrics.by_type[actual_type].misclassified_as[ecc_type] += 1

                            metrics.confusion_matrix.add(ecc_type, actual_type)

                            if state in metrics.by_state:
                                metrics.by_state[state].matched_count += 1
                        else:
                            # Chimeric matching failed but general index matched CeccDNA
                            # This shouldn't normally happen; treat as FP
                            metrics.total_false_positive += 1
                            metrics.by_type[ecc_type].false_positive += 1
                            metrics.unmatched_detected_ids.add(ecc_id)
                            if state in metrics.by_state:
                                metrics.by_state[state].unmatched_count += 1
                    else:
                        # No match at all
                        metrics.total_false_positive += 1
                        metrics.by_type[ecc_type].false_positive += 1
                        metrics.unmatched_detected_ids.add(ecc_id)
                        if state in metrics.by_state:
                            metrics.by_state[state].unmatched_count += 1
            else:
                # Non-chimeric matching (UeccDNA, MeccDNA)
                match = self._find_truth_match(regions)

                if match:
                    truth_name, truth_type, truth_length = match
                    matched_truth.add(truth_name)
                    metrics.matched_detected_ids.add(ecc_id)

                    # Map short type to full type
                    type_map = {"U": "UeccDNA", "M": "MeccDNA", "C": "CeccDNA"}
                    actual_type = type_map.get(truth_type, truth_type)

                    if actual_type == ecc_type:
                        # Correct classification
                        metrics.total_true_positive += 1
                        metrics.by_type[ecc_type].true_positive += 1
                        if state in metrics.by_state:
                            metrics.by_state[state].matched_count += 1

                        # Update length bin
                        bin_label = self._get_length_bin(truth_length)
                        if bin_label in metrics.by_length:
                            metrics.by_length[bin_label].true_positive += 1
                            metrics.by_length[bin_label].detected_count += 1
                    else:
                        # Misclassification
                        metrics.total_misclassified += 1
                        metrics.by_type[ecc_type].false_positive += 1

                        # Track misclassification direction
                        if ecc_type in metrics.by_type:
                            if actual_type not in metrics.by_type[ecc_type].misclassified_from:
                                metrics.by_type[ecc_type].misclassified_from[actual_type] = 0
                            metrics.by_type[ecc_type].misclassified_from[actual_type] += 1

                        if actual_type in metrics.by_type:
                            if ecc_type not in metrics.by_type[actual_type].misclassified_as:
                                metrics.by_type[actual_type].misclassified_as[ecc_type] = 0
                            metrics.by_type[actual_type].misclassified_as[ecc_type] += 1

                        # Update confusion matrix
                        metrics.confusion_matrix.add(ecc_type, actual_type)

                        if state in metrics.by_state:
                            metrics.by_state[state].matched_count += 1
                else:
                    # False positive (no match)
                    metrics.total_false_positive += 1
                    metrics.by_type[ecc_type].false_positive += 1
                    metrics.unmatched_detected_ids.add(ecc_id)

                    if state in metrics.by_state:
                        metrics.by_state[state].unmatched_count += 1

        # Calculate false negatives
        all_truth_ids: set[str] = set()
        for df in [self.truth_u, self.truth_m, self.truth_c]:
            if df is not None:
                all_truth_ids.update(df["name"].tolist())

        metrics.matched_truth_ids = matched_truth
        metrics.unmatched_truth_ids = all_truth_ids - matched_truth
        metrics.total_false_negative = len(metrics.unmatched_truth_ids)

        # Calculate false negatives by type
        for df, etype in [
            (self.truth_u, "UeccDNA"),
            (self.truth_m, "MeccDNA"),
            (self.truth_c, "CeccDNA"),
        ]:
            if df is not None:
                truth_ids = set(df["name"].tolist())
                missed = truth_ids - matched_truth
                metrics.by_type[etype].false_negative = len(missed)

        # Add diagonal to confusion matrix (correct classifications)
        for etype in ["UeccDNA", "MeccDNA", "CeccDNA"]:
            metrics.confusion_matrix.add(etype, etype, metrics.by_type[etype].true_positive)

        self._last_matched_truth = matched_truth
        self._last_unmatched_detected = metrics.unmatched_detected_ids

        logger.info(f"Evaluation complete: {metrics.total_true_positive} TP, "
                    f"{metrics.total_false_positive} FP, {metrics.total_false_negative} FN")

        return metrics

    def get_unmatched_truth_details(self) -> pd.DataFrame:
        """Get details of unmatched (missed) truth entries."""
        matched = self._last_matched_truth
        if matched is None:
            raise ValueError("No evaluation results available. Call evaluate() first.")
        rows = []
        for df, etype in [
            (self.truth_u, "UeccDNA"),
            (self.truth_m, "MeccDNA"),
            (self.truth_c, "CeccDNA"),
        ]:
            if df is None:
                continue
            for _, row in df.iterrows():
                if row["name"] not in matched:
                    rows.append({
                        "name": row["name"],
                        "type": etype,
                        "chrom": row["chrom"],
                        "start": row["start"],
                        "end": row["end"],
                        "length": row["length"],
                    })
        return pd.DataFrame(rows)

    def get_false_positive_details(self) -> pd.DataFrame:
        """Get details of false positive detections."""
        if self.detected is None:
            return pd.DataFrame()

        fp_ids = self._last_unmatched_detected
        if fp_ids is None:
            raise ValueError("No evaluation results available. Call evaluate() first.")
        return self.detected[self.detected["eccDNA_id"].isin(fp_ids)].copy()
