#!/usr/bin/env python3
"""
Benchmark result visualization for eccToolkit iMeta paper.

Reads benchmark output CSV/JSON files and generates publication-quality figures.

Usage:
    python scripts/benchmark_plot.py <BENCHMARK_DIR> [--output-dir FIGURE_DIR]

Expected directory structure:
    <BENCHMARK_DIR>/
        exp1_standard/benchmark/   -> *_summary.csv, *_report.json
        exp2_coverage/benchmark/   -> *_summary.csv, *_report.json
        exp3_complexity/benchmark/ -> *_summary.csv, *_report.json
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

# ============================================================================
# Configuration
# ============================================================================

# Tool display names and colors
TOOL_COLORS = {
    "CircleSeeker": "#E64B35",
    "CReSIL": "#4DBBD5",
    "Circle-Map": "#00A087",
    "ecc_finder": "#3C5488",
}

# eccDNA type colors
TYPE_COLORS = {
    "UeccDNA": "#2166AC",
    "MeccDNA": "#B2182B",
    "CeccDNA": "#35978F",
}

# Platform display names
PLATFORM_NAMES = {
    "hifi": "PacBio HiFi",
    "ngs": "Illumina NGS",
    "ont": "Oxford Nanopore",
}

# Standard figure settings
DPI = 300
FONT_FAMILY = "Arial"


def setup_matplotlib():
    """Configure matplotlib for publication-quality figures."""
    plt.rcParams.update({
        "font.family": FONT_FAMILY,
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.1,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


# ============================================================================
# Data loading helpers
# ============================================================================

def load_benchmark_csv(csv_path: Path) -> pd.DataFrame:
    """Load a benchmark summary CSV file."""
    return pd.read_csv(csv_path)


def load_benchmark_json(json_path: Path) -> dict:
    """Load a benchmark report JSON file."""
    with open(json_path) as f:
        return json.load(f)


def collect_standard_results(benchmark_dir: Path) -> pd.DataFrame:
    """Collect Experiment 1 results into a unified DataFrame.

    Expected files: cs_hifi_summary.csv, cresil_hifi_summary.csv, etc.
    """
    exp_dir = benchmark_dir / "exp1_standard" / "benchmark"
    if not exp_dir.exists():
        print(f"[WARN] Experiment 1 directory not found: {exp_dir}")
        return pd.DataFrame()

    rows = []
    # Map prefix patterns to (tool, platform)
    prefix_map = {
        "cs_hifi": ("CircleSeeker", "hifi"),
        "cs_ngs": ("CircleSeeker", "ngs"),
        "cs_ont": ("CircleSeeker", "ont"),
        "cresil_hifi": ("CReSIL", "hifi"),
        "circlemap_ngs": ("Circle-Map", "ngs"),
        "eccfinder_hifi": ("ecc_finder", "hifi"),
        "eccfinder_ngs": ("ecc_finder", "ngs"),
    }

    for prefix, (tool, platform) in prefix_map.items():
        csv_path = exp_dir / f"{prefix}_summary.csv"
        json_path = exp_dir / f"{prefix}_report.json"

        if json_path.exists():
            data = load_benchmark_json(json_path)
            # Overall metrics
            rows.append({
                "Tool": tool,
                "Platform": platform,
                "Type": "Overall",
                "Precision": data["overall"]["precision"],
                "Recall": data["overall"]["recall"],
                "F1": data["overall"]["f1_score"],
            })
            # By type
            for etype, tdata in data.get("by_type", {}).items():
                rows.append({
                    "Tool": tool,
                    "Platform": platform,
                    "Type": etype,
                    "Precision": tdata["precision"],
                    "Recall": tdata["recall"],
                    "F1": tdata["f1_score"],
                })
        elif csv_path.exists():
            df = load_benchmark_csv(csv_path)
            for _, row in df.iterrows():
                if row["Category"] in ("Overall", "By Type"):
                    rows.append({
                        "Tool": tool,
                        "Platform": platform,
                        "Type": row["Type"],
                        "Precision": row["Precision"],
                        "Recall": row["Recall"],
                        "F1": row["F1"],
                    })

    return pd.DataFrame(rows)


def collect_coverage_results(benchmark_dir: Path) -> pd.DataFrame:
    """Collect Experiment 2 coverage gradient results."""
    exp_dir = benchmark_dir / "exp2_coverage" / "benchmark"
    if not exp_dir.exists():
        print(f"[WARN] Experiment 2 directory not found: {exp_dir}")
        return pd.DataFrame()

    rows = []
    coverages = [5, 10, 30, 50, 100]
    tool_prefixes = {"cs": "CircleSeeker", "cresil": "CReSIL"}

    for cov in coverages:
        for prefix, tool in tool_prefixes.items():
            json_path = exp_dir / f"{prefix}_{cov}x_report.json"
            if json_path.exists():
                data = load_benchmark_json(json_path)
                # Overall
                rows.append({
                    "Tool": tool,
                    "Coverage": cov,
                    "Type": "Overall",
                    "Precision": data["overall"]["precision"],
                    "Recall": data["overall"]["recall"],
                    "F1": data["overall"]["f1_score"],
                })
                # By type
                for etype, tdata in data.get("by_type", {}).items():
                    rows.append({
                        "Tool": tool,
                        "Coverage": cov,
                        "Type": etype,
                        "Precision": tdata["precision"],
                        "Recall": tdata["recall"],
                        "F1": tdata["f1_score"],
                    })

    return pd.DataFrame(rows)


def collect_complexity_results(benchmark_dir: Path) -> pd.DataFrame:
    """Collect Experiment 3 complexity gradient results."""
    exp_dir = benchmark_dir / "exp3_complexity" / "benchmark"
    if not exp_dir.exists():
        print(f"[WARN] Experiment 3 directory not found: {exp_dir}")
        return pd.DataFrame()

    rows = []
    totals = [1000, 5000, 10000, 50000]

    for total in totals:
        json_path = exp_dir / f"cs_n{total}_report.json"
        if json_path.exists():
            data = load_benchmark_json(json_path)
            rows.append({
                "Tool": "CircleSeeker",
                "N_total": total,
                "Type": "Overall",
                "Precision": data["overall"]["precision"],
                "Recall": data["overall"]["recall"],
                "F1": data["overall"]["f1_score"],
            })
            for etype, tdata in data.get("by_type", {}).items():
                rows.append({
                    "Tool": "CircleSeeker",
                    "N_total": total,
                    "Type": etype,
                    "Precision": tdata["precision"],
                    "Recall": tdata["recall"],
                    "F1": tdata["f1_score"],
                })

    return pd.DataFrame(rows)


def collect_length_results(benchmark_dir: Path) -> pd.DataFrame:
    """Collect length-stratified results from Experiment 1 JSON reports."""
    exp_dir = benchmark_dir / "exp1_standard" / "benchmark"
    if not exp_dir.exists():
        return pd.DataFrame()

    rows = []
    prefix_map = {
        "cs_hifi": ("CircleSeeker", "hifi"),
        "cresil_hifi": ("CReSIL", "hifi"),
        "circlemap_ngs": ("Circle-Map", "ngs"),
        "eccfinder_hifi": ("ecc_finder", "hifi"),
    }

    for prefix, (tool, platform) in prefix_map.items():
        json_path = exp_dir / f"{prefix}_report.json"
        if json_path.exists():
            data = load_benchmark_json(json_path)
            for bin_label, bdata in data.get("by_length", {}).items():
                rows.append({
                    "Tool": tool,
                    "Platform": platform,
                    "Length_Bin": bin_label,
                    "Truth": bdata["truth_count"],
                    "Detected": bdata["detected_count"],
                    "TP": bdata["true_positive"],
                    "Recall": bdata["recall"],
                })

    return pd.DataFrame(rows)


def collect_confusion_matrix(benchmark_dir: Path, prefix: str = "cs_hifi") -> dict:
    """Load confusion matrix from a specific benchmark JSON report."""
    json_path = benchmark_dir / "exp1_standard" / "benchmark" / f"{prefix}_report.json"
    if json_path.exists():
        data = load_benchmark_json(json_path)
        return data.get("confusion_matrix", {})
    return {}


# ============================================================================
# Plot functions
# ============================================================================

def plot_standard_comparison(df: pd.DataFrame, output_dir: Path):
    """Plot 1: Tool comparison grouped bar chart (iMeta Fig 3a candidate).

    One panel per platform. Within each panel: tools on x-axis, P/R/F1 bars
    grouped by eccDNA type.
    """
    if df.empty:
        print("[SKIP] No standard benchmark data available.")
        return

    platforms = df["Platform"].unique()
    metrics = ["Precision", "Recall", "F1"]

    for platform in platforms:
        pdf = df[(df["Platform"] == platform) & (df["Type"] != "Overall")]
        if pdf.empty:
            continue

        tools = pdf["Tool"].unique()
        types = [t for t in ["UeccDNA", "MeccDNA", "CeccDNA"] if t in pdf["Type"].values]

        n_tools = len(tools)
        n_types = len(types)
        if n_tools == 0 or n_types == 0:
            continue

        fig, axes = plt.subplots(1, 3, figsize=(12, 4))

        for ax_idx, metric in enumerate(metrics):
            ax = axes[ax_idx]
            x = np.arange(n_tools)
            width = 0.8 / n_types

            for i, etype in enumerate(types):
                values = []
                for tool in tools:
                    row = pdf[(pdf["Tool"] == tool) & (pdf["Type"] == etype)]
                    val = row[metric].values[0] if len(row) > 0 else 0
                    values.append(val)
                offset = (i - n_types / 2 + 0.5) * width
                bars = ax.bar(x + offset, values, width * 0.9,
                              label=etype, color=TYPE_COLORS.get(etype, "#999999"))
                # Add value labels on bars
                for bar, val in zip(bars, values):
                    if val > 0:
                        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                                f"{val:.1%}", ha="center", va="bottom", fontsize=7)

            ax.set_xlabel("")
            ax.set_ylabel(metric)
            ax.set_title(metric)
            ax.set_xticks(x)
            ax.set_xticklabels(tools, rotation=30, ha="right")
            ax.set_ylim(0, 1.15)
            ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))

            if ax_idx == 0:
                ax.legend(loc="upper left", frameon=False)

        platform_name = PLATFORM_NAMES.get(platform, platform)
        fig.suptitle(f"Tool Comparison - {platform_name}", fontsize=13, y=1.02)
        plt.tight_layout()

        out_path = output_dir / f"fig3a_comparison_{platform}.pdf"
        fig.savefig(out_path, dpi=DPI)
        plt.close(fig)
        print(f"  Saved: {out_path}")


def plot_standard_overall(df: pd.DataFrame, output_dir: Path):
    """Plot: Overall P/R/F1 for each tool across platforms."""
    if df.empty:
        return

    odf = df[df["Type"] == "Overall"]
    if odf.empty:
        return

    platforms = odf["Platform"].unique()
    metrics = ["Precision", "Recall", "F1"]

    for platform in platforms:
        pdf = odf[odf["Platform"] == platform]
        tools = pdf["Tool"].unique()
        n_tools = len(tools)
        if n_tools == 0:
            continue

        fig, ax = plt.subplots(figsize=(6, 4))
        x = np.arange(n_tools)
        width = 0.25

        for i, metric in enumerate(metrics):
            values = [pdf[pdf["Tool"] == t][metric].values[0] if len(pdf[pdf["Tool"] == t]) > 0 else 0
                      for t in tools]
            offset = (i - 1) * width
            bars = ax.bar(x + offset, values, width * 0.9, label=metric,
                          color=["#2166AC", "#B2182B", "#35978F"][i])
            for bar, val in zip(bars, values):
                if val > 0:
                    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                            f"{val:.1%}", ha="center", va="bottom", fontsize=8)

        platform_name = PLATFORM_NAMES.get(platform, platform)
        ax.set_title(f"Overall Performance - {platform_name}")
        ax.set_xticks(x)
        ax.set_xticklabels(tools, rotation=30, ha="right")
        ax.set_ylim(0, 1.15)
        ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))
        ax.legend(frameon=False)
        plt.tight_layout()

        out_path = output_dir / f"fig3a_overall_{platform}.pdf"
        fig.savefig(out_path, dpi=DPI)
        plt.close(fig)
        print(f"  Saved: {out_path}")


def plot_coverage_gradient(df: pd.DataFrame, output_dir: Path):
    """Plot 2: Coverage gradient line plot (iMeta Fig 4a candidate).

    x: coverage, y: metric (Recall/F1)
    Lines for each tool, separate panels for U/M/C types.
    """
    if df.empty:
        print("[SKIP] No coverage gradient data available.")
        return

    types = [t for t in ["UeccDNA", "MeccDNA", "CeccDNA"] if t in df["Type"].values]
    tools = df["Tool"].unique()

    for metric in ["Recall", "F1"]:
        fig, axes = plt.subplots(1, len(types), figsize=(4 * len(types), 4), sharey=True)
        if len(types) == 1:
            axes = [axes]

        for ax_idx, etype in enumerate(types):
            ax = axes[ax_idx]
            tdf = df[df["Type"] == etype]

            for tool in tools:
                tool_df = tdf[tdf["Tool"] == tool].sort_values("Coverage")
                if tool_df.empty:
                    continue
                ax.plot(tool_df["Coverage"], tool_df[metric],
                        marker="o", linewidth=2, markersize=6,
                        label=tool, color=TOOL_COLORS.get(tool, "#999999"))

            ax.set_xlabel("Coverage (x)")
            if ax_idx == 0:
                ax.set_ylabel(metric)
            ax.set_title(etype, color=TYPE_COLORS.get(etype, "black"))
            ax.set_xscale("log")
            ax.set_xticks([5, 10, 30, 50, 100])
            ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
            ax.set_ylim(-0.05, 1.05)
            ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))
            ax.grid(True, alpha=0.3, linestyle="--")

            if ax_idx == len(types) - 1:
                ax.legend(loc="lower right", frameon=False)

        fig.suptitle(f"{metric} vs Coverage", fontsize=13, y=1.02)
        plt.tight_layout()

        out_path = output_dir / f"fig4a_coverage_{metric.lower()}.pdf"
        fig.savefig(out_path, dpi=DPI)
        plt.close(fig)
        print(f"  Saved: {out_path}")


def plot_size_stratified(df: pd.DataFrame, output_dir: Path):
    """Plot 3: Size-stratified performance heatmap (iMeta Fig 3b candidate).

    Heatmap: rows = tools, columns = size bins, color = Recall.
    """
    if df.empty:
        print("[SKIP] No length-stratified data available.")
        return

    bin_order = ["<500bp", "500bp-1kb", "1kb-5kb", "5kb-50kb", ">50kb"]
    tools = df["Tool"].unique()

    # Build recall matrix
    recall_matrix = []
    tool_labels = []
    for tool in tools:
        row = []
        for bin_label in bin_order:
            match = df[(df["Tool"] == tool) & (df["Length_Bin"] == bin_label)]
            val = match["Recall"].values[0] if len(match) > 0 else np.nan
            row.append(val)
        recall_matrix.append(row)
        tool_labels.append(tool)

    recall_matrix = np.array(recall_matrix)

    fig, ax = plt.subplots(figsize=(8, max(3, len(tools) * 0.8 + 1)))
    im = ax.imshow(recall_matrix, cmap="RdYlGn", aspect="auto", vmin=0, vmax=1)

    # Annotate cells
    for i in range(len(tools)):
        for j in range(len(bin_order)):
            val = recall_matrix[i, j]
            if not np.isnan(val):
                text_color = "white" if val < 0.3 or val > 0.85 else "black"
                ax.text(j, i, f"{val:.1%}", ha="center", va="center",
                        fontsize=9, color=text_color)

    ax.set_xticks(range(len(bin_order)))
    ax.set_xticklabels(bin_order)
    ax.set_yticks(range(len(tools)))
    ax.set_yticklabels(tool_labels)
    ax.set_xlabel("eccDNA Length")
    ax.set_title("Recall by eccDNA Size")

    cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.04)
    cbar.set_label("Recall")
    cbar.ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))

    plt.tight_layout()

    out_path = output_dir / "fig3b_size_recall_heatmap.pdf"
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_confusion_matrix(cm_data: dict, output_dir: Path, tool_name: str = "CircleSeeker"):
    """Plot 4: Confusion matrix (Supplementary).

    3x3 matrix showing type classification accuracy.
    """
    if not cm_data:
        print("[SKIP] No confusion matrix data available.")
        return

    types = ["UeccDNA", "MeccDNA", "CeccDNA"]
    matrix = np.zeros((3, 3), dtype=int)

    for i, pred in enumerate(types):
        for j, actual in enumerate(types):
            matrix[i, j] = cm_data.get(pred, {}).get(actual, 0)

    # Normalize rows for display (predicted -> actual distribution)
    row_sums = matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    matrix_norm = matrix / row_sums

    fig, ax = plt.subplots(figsize=(5, 4))
    im = ax.imshow(matrix_norm, cmap="Blues", vmin=0, vmax=1)

    for i in range(3):
        for j in range(3):
            count = matrix[i, j]
            pct = matrix_norm[i, j]
            text_color = "white" if pct > 0.6 else "black"
            ax.text(j, i, f"{count}\n({pct:.1%})", ha="center", va="center",
                    fontsize=10, color=text_color)

    ax.set_xticks(range(3))
    ax.set_xticklabels(types, rotation=30, ha="right")
    ax.set_yticks(range(3))
    ax.set_yticklabels(types)
    ax.set_xlabel("Actual Type")
    ax.set_ylabel("Predicted Type")
    ax.set_title(f"Type Classification - {tool_name}")

    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Proportion")
    plt.tight_layout()

    out_path = output_dir / f"suppl_confusion_{tool_name.lower().replace('-', '_')}.pdf"
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_complexity_gradient(df: pd.DataFrame, output_dir: Path):
    """Plot 5: Complexity gradient (Supplementary).

    Bar chart: N_total on x-axis, P/R/F1 on y-axis, grouped by eccDNA type.
    """
    if df.empty:
        print("[SKIP] No complexity gradient data available.")
        return

    types = [t for t in ["UeccDNA", "MeccDNA", "CeccDNA"] if t in df["Type"].values]
    totals = sorted(df["N_total"].unique())

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    metrics = ["Precision", "Recall", "F1"]

    for ax_idx, metric in enumerate(metrics):
        ax = axes[ax_idx]
        x = np.arange(len(totals))
        n_types = len(types)
        width = 0.8 / n_types

        for i, etype in enumerate(types):
            values = []
            for total in totals:
                row = df[(df["N_total"] == total) & (df["Type"] == etype)]
                val = row[metric].values[0] if len(row) > 0 else 0
                values.append(val)
            offset = (i - n_types / 2 + 0.5) * width
            ax.bar(x + offset, values, width * 0.9,
                   label=etype, color=TYPE_COLORS.get(etype, "#999999"))

        ax.set_xlabel("Number of eccDNA")
        ax.set_title(metric)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{t // 1000}K" if t >= 1000 else str(t) for t in totals])
        ax.set_ylim(0, 1.15)
        ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))

        if ax_idx == 0:
            ax.legend(frameon=False)

    fig.suptitle("Performance vs eccDNA Complexity (CircleSeeker)", fontsize=13, y=1.02)
    plt.tight_layout()

    out_path = output_dir / "suppl_complexity_gradient.pdf"
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate benchmark figures for eccToolkit iMeta paper"
    )
    parser.add_argument("benchmark_dir", type=Path,
                        help="Base benchmark output directory")
    parser.add_argument("--output-dir", "-o", type=Path, default=None,
                        help="Figure output directory (default: <benchmark_dir>/figures)")
    args = parser.parse_args()

    benchmark_dir = args.benchmark_dir
    if not benchmark_dir.exists():
        print(f"ERROR: Benchmark directory not found: {benchmark_dir}")
        sys.exit(1)

    output_dir = args.output_dir or benchmark_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    setup_matplotlib()

    print("=" * 60)
    print("eccToolkit Benchmark Figure Generation")
    print("=" * 60)
    print(f"Input:  {benchmark_dir}")
    print(f"Output: {output_dir}")
    print()

    # --- Figure 3a: Standard tool comparison ---
    print("--- Figure 3a: Tool comparison ---")
    std_df = collect_standard_results(benchmark_dir)
    plot_standard_comparison(std_df, output_dir)
    plot_standard_overall(std_df, output_dir)

    # --- Figure 3b: Size-stratified performance ---
    print("--- Figure 3b: Size-stratified recall ---")
    len_df = collect_length_results(benchmark_dir)
    plot_size_stratified(len_df, output_dir)

    # --- Figure 4a: Coverage gradient ---
    print("--- Figure 4a: Coverage gradient ---")
    cov_df = collect_coverage_results(benchmark_dir)
    plot_coverage_gradient(cov_df, output_dir)

    # --- Supplementary: Confusion matrix ---
    print("--- Supplementary: Confusion matrix ---")
    cm = collect_confusion_matrix(benchmark_dir, prefix="cs_hifi")
    plot_confusion_matrix(cm, output_dir, tool_name="CircleSeeker")

    # --- Supplementary: Complexity gradient ---
    print("--- Supplementary: Complexity gradient ---")
    comp_df = collect_complexity_results(benchmark_dir)
    plot_complexity_gradient(comp_df, output_dir)

    print()
    print(f"All figures saved to: {output_dir}")


if __name__ == "__main__":
    main()
