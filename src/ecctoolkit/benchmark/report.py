"""
Report generation for benchmark results.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

from ecctoolkit.benchmark.metrics import BenchmarkMetrics

logger = logging.getLogger(__name__)


class BenchmarkReporter:
    """Generate benchmark reports in various formats."""

    def __init__(self, metrics: BenchmarkMetrics):
        """
        Initialize reporter with metrics.

        Args:
            metrics: BenchmarkMetrics object from evaluator
        """
        self.metrics = metrics

    def generate_text_report(self) -> str:
        """Generate plain text report."""
        lines = []
        m = self.metrics

        lines.append("=" * 70)
        lines.append("eccDNA Detection Benchmark Report")
        lines.append("=" * 70)
        lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("")

        # Overall metrics
        lines.append("-" * 70)
        lines.append("1. OVERALL METRICS")
        lines.append("-" * 70)
        lines.append(f"Total Ground Truth:    {m.total_truth:,}")
        lines.append(f"Total Detected:        {m.total_detected:,}")
        lines.append(f"True Positives:        {m.total_true_positive:,}")
        lines.append(f"False Positives:       {m.total_false_positive:,}")
        lines.append(f"False Negatives:       {m.total_false_negative:,}")
        lines.append(f"Misclassified:         {m.total_misclassified:,}")
        lines.append("")
        lines.append(f"Precision:             {m.overall_precision:.2%}")
        lines.append(f"Recall:                {m.overall_recall:.2%}")
        lines.append(f"F1 Score:              {m.overall_f1:.2%}")
        lines.append("")

        # By type
        lines.append("-" * 70)
        lines.append("2. METRICS BY TYPE")
        lines.append("-" * 70)
        lines.append(f"{'Type':<12} {'Truth':>8} {'Detected':>10} {'TP':>8} {'FP':>8} "
                     f"{'FN':>8} {'Prec':>8} {'Recall':>8} {'F1':>8}")
        lines.append("-" * 70)

        for etype in ["UeccDNA", "MeccDNA", "CeccDNA"]:
            if etype in m.by_type:
                t = m.by_type[etype]
                lines.append(
                    f"{etype:<12} {t.truth_count:>8} {t.detected_count:>10} "
                    f"{t.true_positive:>8} {t.false_positive:>8} {t.false_negative:>8} "
                    f"{t.precision:>7.1%} {t.recall:>7.1%} {t.f1_score:>7.1%}"
                )
        lines.append("")

        # By state
        lines.append("-" * 70)
        lines.append("3. METRICS BY DETECTION SOURCE")
        lines.append("-" * 70)
        lines.append(f"{'State':<12} {'Total':>10} {'Matched':>10} {'Unmatched':>10} "
                     f"{'Precision':>12} {'FP Rate':>10}")
        lines.append("-" * 70)

        for state in ["Confirmed", "Inferred"]:
            if state in m.by_state:
                s = m.by_state[state]
                lines.append(
                    f"{state:<12} {s.total_count:>10} {s.matched_count:>10} "
                    f"{s.unmatched_count:>10} {s.precision:>11.1%} "
                    f"{s.false_positive_rate:>9.1%}"
                )
        lines.append("")

        # By length
        lines.append("-" * 70)
        lines.append("4. RECALL BY LENGTH")
        lines.append("-" * 70)
        lines.append(f"{'Length Bin':<15} {'Truth':>10} {'Detected':>10} {'TP':>10} {'Recall':>10}")
        lines.append("-" * 70)

        for bin_label in ["<500bp", "500bp-1kb", "1kb-5kb", "5kb-50kb", ">50kb"]:
            if bin_label in m.by_length:
                b = m.by_length[bin_label]
                lines.append(
                    f"{bin_label:<15} {b.truth_count:>10} {b.detected_count:>10} "
                    f"{b.true_positive:>10} {b.recall:>9.1%}"
                )
        lines.append("")

        # Confusion matrix
        lines.append("-" * 70)
        lines.append("5. CONFUSION MATRIX (Predicted vs Actual)")
        lines.append("-" * 70)
        lines.append(f"{'Predicted↓ Actual→':<20} {'UeccDNA':>10} {'MeccDNA':>10} {'CeccDNA':>10}")
        lines.append("-" * 70)

        for pred in ["UeccDNA", "MeccDNA", "CeccDNA"]:
            row = m.confusion_matrix.matrix.get(pred, {})
            lines.append(
                f"{pred:<20} {row.get('UeccDNA', 0):>10} "
                f"{row.get('MeccDNA', 0):>10} {row.get('CeccDNA', 0):>10}"
            )
        lines.append("")

        # Misclassification details
        lines.append("-" * 70)
        lines.append("6. MISCLASSIFICATION DETAILS")
        lines.append("-" * 70)

        for etype in ["UeccDNA", "MeccDNA", "CeccDNA"]:
            if etype in m.by_type:
                t = m.by_type[etype]
                if t.misclassified_from:
                    lines.append(f"{etype} <- misclassified from:")
                    for src, cnt in t.misclassified_from.items():
                        lines.append(f"    {src}: {cnt}")
                if t.misclassified_as:
                    lines.append(f"{etype} -> misclassified as:")
                    for dst, cnt in t.misclassified_as.items():
                        lines.append(f"    {dst}: {cnt}")
        lines.append("")

        # CeccDNA-specific metrics
        if "CeccDNA" in m.by_type:
            c = m.by_type["CeccDNA"]
            if c.truth_count > 0 or c.detected_count > 0:
                lines.append("-" * 70)
                lines.append("7. CeccDNA CHIMERIC MATCHING DETAILS")
                lines.append("-" * 70)
                lines.append(f"Segment Accuracy:      {c.segment_accuracy:.2%} "
                             f"({c.total_matched_segments}/{c.total_truth_segments} segments)")
                lines.append(f"Partial Matches:       {c.partial_match_count}")
                lines.append("")
                lines.append(f"{'Subtype':<20} {'Truth':>8} {'Detected':>10} {'TP':>8} "
                             f"{'Prec':>8} {'Recall':>8}")
                lines.append("-" * 70)
                lines.append(
                    f"{'Inter-chromosomal':<20} {c.inter_chr_truth:>8} "
                    f"{c.inter_chr_detected:>10} {c.inter_chr_tp:>8} "
                    f"{c.inter_chr_precision:>7.1%} {c.inter_chr_recall:>7.1%}"
                )
                lines.append(
                    f"{'Intra-chromosomal':<20} {c.intra_chr_truth:>8} "
                    f"{c.intra_chr_detected:>10} {c.intra_chr_tp:>8} "
                    f"{c.intra_chr_precision:>7.1%} {c.intra_chr_recall:>7.1%}"
                )
                lines.append("")

        lines.append("=" * 70)
        lines.append("End of Report")
        lines.append("=" * 70)

        return "\n".join(lines)

    def generate_json_report(self) -> dict:
        """Generate JSON-serializable report."""
        m = self.metrics

        report = {
            "generated_at": datetime.now().isoformat(),
            "overall": {
                "total_truth": m.total_truth,
                "total_detected": m.total_detected,
                "true_positive": m.total_true_positive,
                "false_positive": m.total_false_positive,
                "false_negative": m.total_false_negative,
                "misclassified": m.total_misclassified,
                "precision": m.overall_precision,
                "recall": m.overall_recall,
                "f1_score": m.overall_f1,
            },
            "by_type": {},
            "by_state": {},
            "by_length": {},
            "confusion_matrix": m.confusion_matrix.matrix,
        }

        for etype, t in m.by_type.items():
            type_report = {
                "truth_count": t.truth_count,
                "detected_count": t.detected_count,
                "true_positive": t.true_positive,
                "false_positive": t.false_positive,
                "false_negative": t.false_negative,
                "precision": t.precision,
                "recall": t.recall,
                "f1_score": t.f1_score,
                "misclassified_from": t.misclassified_from,
                "misclassified_as": t.misclassified_as,
            }
            if etype == "CeccDNA":
                type_report["chimeric_details"] = {
                    "segment_accuracy": t.segment_accuracy,
                    "total_matched_segments": t.total_matched_segments,
                    "total_truth_segments": t.total_truth_segments,
                    "partial_match_count": t.partial_match_count,
                    "inter_chromosomal": {
                        "truth": t.inter_chr_truth,
                        "detected": t.inter_chr_detected,
                        "true_positive": t.inter_chr_tp,
                        "precision": t.inter_chr_precision,
                        "recall": t.inter_chr_recall,
                    },
                    "intra_chromosomal": {
                        "truth": t.intra_chr_truth,
                        "detected": t.intra_chr_detected,
                        "true_positive": t.intra_chr_tp,
                        "precision": t.intra_chr_precision,
                        "recall": t.intra_chr_recall,
                    },
                }
            report["by_type"][etype] = type_report

        for state, s in m.by_state.items():
            report["by_state"][state] = {
                "total_count": s.total_count,
                "matched_count": s.matched_count,
                "unmatched_count": s.unmatched_count,
                "precision": s.precision,
                "false_positive_rate": s.false_positive_rate,
            }

        for bin_label, b in m.by_length.items():
            report["by_length"][bin_label] = {
                "truth_count": b.truth_count,
                "detected_count": b.detected_count,
                "true_positive": b.true_positive,
                "recall": b.recall,
            }

        return report

    def save_text_report(self, output_path: str):
        """Save text report to file."""
        report = self.generate_text_report()
        Path(output_path).write_text(report)
        logger.info(f"Text report saved to: {output_path}")

    def save_json_report(self, output_path: str):
        """Save JSON report to file."""
        report = self.generate_json_report()
        with open(output_path, "w") as f:
            json.dump(report, f, indent=2)
        logger.info(f"JSON report saved to: {output_path}")

    def save_csv_summary(self, output_path: str):
        """Save CSV summary with key metrics."""
        import pandas as pd

        m = self.metrics
        rows = []

        # Overall row
        rows.append({
            "Category": "Overall",
            "Type": "All",
            "Truth": m.total_truth,
            "Detected": m.total_detected,
            "TP": m.total_true_positive,
            "FP": m.total_false_positive,
            "FN": m.total_false_negative,
            "Precision": m.overall_precision,
            "Recall": m.overall_recall,
            "F1": m.overall_f1,
        })

        # By type rows
        for etype in ["UeccDNA", "MeccDNA", "CeccDNA"]:
            if etype in m.by_type:
                t = m.by_type[etype]
                rows.append({
                    "Category": "By Type",
                    "Type": etype,
                    "Truth": t.truth_count,
                    "Detected": t.detected_count,
                    "TP": t.true_positive,
                    "FP": t.false_positive,
                    "FN": t.false_negative,
                    "Precision": t.precision,
                    "Recall": t.recall,
                    "F1": t.f1_score,
                })

        # By state rows
        for state in ["Confirmed", "Inferred"]:
            if state in m.by_state:
                s = m.by_state[state]
                rows.append({
                    "Category": "By State",
                    "Type": state,
                    "Truth": None,
                    "Detected": s.total_count,
                    "TP": s.matched_count,
                    "FP": s.unmatched_count,
                    "FN": None,
                    "Precision": s.precision,
                    "Recall": None,
                    "F1": None,
                })

        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False)
        logger.info(f"CSV summary saved to: {output_path}")

    def print_summary(self):
        """Print summary to console."""
        m = self.metrics

        print("\n" + "=" * 60)
        print("BENCHMARK SUMMARY")
        print("=" * 60)
        print(f"Overall Precision: {m.overall_precision:.2%}")
        print(f"Overall Recall:    {m.overall_recall:.2%}")
        print(f"Overall F1:        {m.overall_f1:.2%}")
        print()

        print("By Type:")
        print(f"  {'Type':<12} {'Precision':>12} {'Recall':>12} {'F1':>12}")
        for etype in ["UeccDNA", "MeccDNA", "CeccDNA"]:
            if etype in m.by_type:
                t = m.by_type[etype]
                print(f"  {etype:<12} {t.precision:>11.1%} {t.recall:>11.1%} {t.f1_score:>11.1%}")
        print()

        print("By Detection Source:")
        print(f"  {'State':<12} {'Precision':>12} {'FP Rate':>12}")
        for state in ["Confirmed", "Inferred"]:
            if state in m.by_state:
                s = m.by_state[state]
                print(f"  {state:<12} {s.precision:>11.1%} {s.false_positive_rate:>11.1%}")
        print("=" * 60 + "\n")
