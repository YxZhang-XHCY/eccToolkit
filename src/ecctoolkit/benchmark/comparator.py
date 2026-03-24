"""
Multi-tool benchmark comparison for eccDNA detection.

Evaluates multiple detection tools against the same ground truth
and generates comparison reports.
"""

from __future__ import annotations

import logging
import os
from collections import OrderedDict
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd

from ecctoolkit.benchmark.evaluator import BenchmarkEvaluator
from ecctoolkit.benchmark.metrics import BenchmarkMetrics
from ecctoolkit.benchmark.parsers import parse_tool_output

logger = logging.getLogger(__name__)


class BenchmarkComparator:
    """Compare multiple eccDNA detection tools against the same ground truth."""

    def __init__(self, overlap_threshold: float = 0.8):
        self.overlap_threshold = overlap_threshold
        self._truth_args: Optional[dict] = None
        self._tools: OrderedDict[str, pd.DataFrame] = OrderedDict()
        self._results: OrderedDict[str, BenchmarkMetrics] = OrderedDict()

    def load_truth(
        self,
        unique_bed: Optional[str] = None,
        multi_bed: Optional[str] = None,
        chimeric_bed: Optional[str] = None,
        truth_dir: Optional[str] = None,
    ):
        """Load ground truth BED files.

        Args:
            unique_bed: Path to unique.bed file
            multi_bed: Path to multi.bed file
            chimeric_bed: Path to chimeric.bed file
            truth_dir: Directory containing *.unique.bed, *.multi.bed, *.chimeric.bed
        """
        self._truth_args = {
            "unique_bed": unique_bed,
            "multi_bed": multi_bed,
            "chimeric_bed": chimeric_bed,
            "truth_dir": truth_dir,
        }
        # Validate by loading once
        evaluator = BenchmarkEvaluator(overlap_threshold=self.overlap_threshold)
        evaluator.load_truth(**self._truth_args)
        logger.info("Ground truth loaded successfully")

    def add_tool(self, name: str, result_file: str, fmt: str = "circleseeker"):
        """Add a tool's detection results.

        Args:
            name: Tool display name (e.g., "CircleSeeker", "CReSIL")
            result_file: Path to the tool's output file
            fmt: Format name (circleseeker, circlemap, cresil, eccfinder, bed, csv)
        """
        logger.info(f"Adding tool: {name} (format: {fmt}, file: {result_file})")
        detected_df = parse_tool_output(result_file, fmt)
        self._tools[name] = detected_df
        logger.info(f"Tool '{name}': {len(detected_df)} detections loaded")

    def evaluate_all(self) -> OrderedDict[str, BenchmarkMetrics]:
        """Run benchmark for all added tools against the same truth.

        Returns:
            OrderedDict mapping tool name to BenchmarkMetrics
        """
        if self._truth_args is None:
            raise ValueError("No ground truth loaded. Call load_truth() first.")
        if not self._tools:
            raise ValueError("No tools added. Call add_tool() first.")

        self._results.clear()

        for tool_name, detected_df in self._tools.items():
            logger.info(f"Evaluating tool: {tool_name}")
            evaluator = BenchmarkEvaluator(overlap_threshold=self.overlap_threshold)
            evaluator.load_truth(**self._truth_args)
            evaluator.detected = detected_df
            metrics = evaluator.evaluate()
            self._results[tool_name] = metrics
            logger.info(
                f"  {tool_name}: P={metrics.overall_precision:.3f} "
                f"R={metrics.overall_recall:.3f} F1={metrics.overall_f1:.3f}"
            )

        return self._results

    def generate_comparison_report(self, output_dir: str):
        """Generate multi-tool comparison outputs.

        Files created:
            comparison_summary.csv: One row per tool per type
            comparison_by_size.csv: Performance by length bin
            comparison_matrix.csv: Wide format for easy plotting
            comparison_report.txt: Human-readable summary
        """
        if not self._results:
            raise ValueError("No results available. Call evaluate_all() first.")

        os.makedirs(output_dir, exist_ok=True)

        self._save_comparison_summary(output_dir)
        self._save_comparison_by_size(output_dir)
        self._save_comparison_matrix(output_dir)
        self._save_comparison_text(output_dir)

        logger.info(f"Comparison reports saved to: {output_dir}")

    def _save_comparison_summary(self, output_dir: str):
        """Save comparison_summary.csv: tool, type, metrics."""
        rows = []
        for tool_name, metrics in self._results.items():
            # Overall row
            rows.append({
                "tool": tool_name,
                "type": "Overall",
                "truth_count": metrics.total_truth,
                "detected_count": metrics.total_detected,
                "tp": metrics.total_true_positive,
                "fp": metrics.total_false_positive,
                "fn": metrics.total_false_negative,
                "precision": round(metrics.overall_precision, 4),
                "recall": round(metrics.overall_recall, 4),
                "f1": round(metrics.overall_f1, 4),
            })
            # By type
            for etype in ["UeccDNA", "MeccDNA", "CeccDNA"]:
                if etype in metrics.by_type:
                    t = metrics.by_type[etype]
                    # Map type names for readability
                    type_label = {"UeccDNA": "Unique", "MeccDNA": "Multi", "CeccDNA": "Chimeric"}
                    rows.append({
                        "tool": tool_name,
                        "type": type_label.get(etype, etype),
                        "truth_count": t.truth_count,
                        "detected_count": t.detected_count,
                        "tp": t.true_positive,
                        "fp": t.false_positive,
                        "fn": t.false_negative,
                        "precision": round(t.precision, 4),
                        "recall": round(t.recall, 4),
                        "f1": round(t.f1_score, 4),
                    })

        df = pd.DataFrame(rows)
        path = os.path.join(output_dir, "comparison_summary.csv")
        df.to_csv(path, index=False)
        logger.info(f"Saved: {path}")

    def _save_comparison_by_size(self, output_dir: str):
        """Save comparison_by_size.csv: tool, size_bin, metrics."""
        rows = []
        for tool_name, metrics in self._results.items():
            for bin_label, b in metrics.by_length.items():
                rows.append({
                    "tool": tool_name,
                    "size_bin": bin_label,
                    "truth_count": b.truth_count,
                    "detected_count": b.detected_count,
                    "tp": b.true_positive,
                    "recall": round(b.recall, 4),
                })

        df = pd.DataFrame(rows)
        path = os.path.join(output_dir, "comparison_by_size.csv")
        df.to_csv(path, index=False)
        logger.info(f"Saved: {path}")

    def _save_comparison_matrix(self, output_dir: str):
        """Save comparison_matrix.csv: wide format for plotting.

        Columns: tool, overall_precision, overall_recall, overall_f1,
                 unique_precision, unique_recall, unique_f1, ...
        """
        rows = []
        for tool_name, metrics in self._results.items():
            row = {
                "tool": tool_name,
                "overall_precision": round(metrics.overall_precision, 4),
                "overall_recall": round(metrics.overall_recall, 4),
                "overall_f1": round(metrics.overall_f1, 4),
            }
            type_map = {"UeccDNA": "unique", "MeccDNA": "multi", "CeccDNA": "chimeric"}
            for etype, prefix in type_map.items():
                if etype in metrics.by_type:
                    t = metrics.by_type[etype]
                    row[f"{prefix}_precision"] = round(t.precision, 4)
                    row[f"{prefix}_recall"] = round(t.recall, 4)
                    row[f"{prefix}_f1"] = round(t.f1_score, 4)
                else:
                    row[f"{prefix}_precision"] = None
                    row[f"{prefix}_recall"] = None
                    row[f"{prefix}_f1"] = None
            rows.append(row)

        df = pd.DataFrame(rows)
        path = os.path.join(output_dir, "comparison_matrix.csv")
        df.to_csv(path, index=False)
        logger.info(f"Saved: {path}")

    def _save_comparison_text(self, output_dir: str):
        """Save comparison_report.txt: human-readable summary."""
        lines = []
        lines.append("=" * 70)
        lines.append("Multi-Tool eccDNA Detection Benchmark Comparison")
        lines.append("=" * 70)
        lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"Overlap threshold: {self.overlap_threshold}")
        lines.append(f"Tools compared: {len(self._results)}")
        lines.append("")

        # Overall comparison table
        lines.append("-" * 70)
        lines.append("1. OVERALL COMPARISON")
        lines.append("-" * 70)
        lines.append(
            f"{'Tool':<20} {'Detected':>10} {'TP':>8} {'FP':>8} {'FN':>8} "
            f"{'Prec':>8} {'Recall':>8} {'F1':>8}"
        )
        lines.append("-" * 70)

        for tool_name, m in self._results.items():
            lines.append(
                f"{tool_name:<20} {m.total_detected:>10} "
                f"{m.total_true_positive:>8} {m.total_false_positive:>8} "
                f"{m.total_false_negative:>8} "
                f"{m.overall_precision:>7.1%} {m.overall_recall:>7.1%} "
                f"{m.overall_f1:>7.1%}"
            )
        lines.append("")

        # By type comparison
        for etype, type_label in [
            ("UeccDNA", "Unique eccDNA (UeccDNA)"),
            ("MeccDNA", "Multi-fragment eccDNA (MeccDNA)"),
            ("CeccDNA", "Chimeric eccDNA (CeccDNA)"),
        ]:
            lines.append("-" * 70)
            lines.append(f"2. {type_label}")
            lines.append("-" * 70)
            lines.append(
                f"{'Tool':<20} {'Truth':>8} {'Det':>8} {'TP':>8} {'FP':>8} "
                f"{'FN':>8} {'Prec':>8} {'Recall':>8} {'F1':>8}"
            )
            lines.append("-" * 70)

            for tool_name, m in self._results.items():
                if etype in m.by_type:
                    t = m.by_type[etype]
                    lines.append(
                        f"{tool_name:<20} {t.truth_count:>8} "
                        f"{t.detected_count:>8} {t.true_positive:>8} "
                        f"{t.false_positive:>8} {t.false_negative:>8} "
                        f"{t.precision:>7.1%} {t.recall:>7.1%} "
                        f"{t.f1_score:>7.1%}"
                    )
            lines.append("")

        # By length bin
        lines.append("-" * 70)
        lines.append("3. RECALL BY LENGTH BIN")
        lines.append("-" * 70)

        # Collect all bin labels
        all_bins = []
        for m in self._results.values():
            for bl in m.by_length:
                if bl not in all_bins:
                    all_bins.append(bl)

        header = f"{'Tool':<20}"
        for bl in all_bins:
            header += f" {bl:>12}"
        lines.append(header)
        lines.append("-" * 70)

        for tool_name, m in self._results.items():
            row = f"{tool_name:<20}"
            for bl in all_bins:
                if bl in m.by_length:
                    row += f" {m.by_length[bl].recall:>11.1%}"
                else:
                    row += f" {'N/A':>12}"
            lines.append(row)
        lines.append("")

        lines.append("=" * 70)
        lines.append("End of Comparison Report")
        lines.append("=" * 70)

        text = "\n".join(lines)
        path = os.path.join(output_dir, "comparison_report.txt")
        Path(path).write_text(text)
        logger.info(f"Saved: {path}")

    def print_summary(self):
        """Print comparison summary to console."""
        if not self._results:
            print("No results available.")
            return

        print()
        print("=" * 60)
        print("MULTI-TOOL BENCHMARK COMPARISON")
        print("=" * 60)

        print(f"\n{'Tool':<20} {'Precision':>12} {'Recall':>12} {'F1':>12}")
        print("-" * 56)
        for tool_name, m in self._results.items():
            print(
                f"{tool_name:<20} {m.overall_precision:>11.1%} "
                f"{m.overall_recall:>11.1%} {m.overall_f1:>11.1%}"
            )

        print("\nBy Type:")
        for etype in ["UeccDNA", "MeccDNA", "CeccDNA"]:
            print(f"\n  {etype}:")
            print(f"  {'Tool':<20} {'Precision':>12} {'Recall':>12} {'F1':>12}")
            print("  " + "-" * 56)
            for tool_name, m in self._results.items():
                if etype in m.by_type:
                    t = m.by_type[etype]
                    print(
                        f"  {tool_name:<20} {t.precision:>11.1%} "
                        f"{t.recall:>11.1%} {t.f1_score:>11.1%}"
                    )

        print("\n" + "=" * 60)
