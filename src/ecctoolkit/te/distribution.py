"""TE percentage distribution analysis."""

from __future__ import annotations

import logging
import os
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def analyze_te_distribution(
    input_file: str,
    output_dir: str,
    bins: int = 5,
    group_by: Optional[str] = None,
) -> None:
    """Analyze TE percentage distribution in quantile bins.

    Bins eccDNA by TE coverage percentage into equal-width bins
    (default 5 bins: 0-20%, 20-40%, 40-60%, 60-80%, 80-100%).

    Args:
        input_file: Per-eccDNA TE annotation CSV (from te-analyze output)
            with at least columns: eccDNA_name, te_pct.
        output_dir: Output directory.
        bins: Number of bins (default 5).
        group_by: Optional column name to group by (e.g., "eccDNA_type",
            "te_class"). When set to a TE-level column like "te_class",
            the distribution is computed per TE class.
    """
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(input_file)
    logger.info(f"Loaded {len(df)} rows from {input_file}")

    # Ensure te_pct column exists
    if "te_pct" not in df.columns:
        if "te_bp" in df.columns and "eccDNA_len" in df.columns:
            df["te_pct"] = (df["te_bp"] / df["eccDNA_len"] * 100).clip(0, 100)
        elif "class_pct" in df.columns:
            # TE class breakdown format
            df["te_pct"] = df["class_pct"].clip(0, 100)
        else:
            raise ValueError(
                f"Cannot find or compute te_pct. Need 'te_pct' or "
                f"'te_bp'+'eccDNA_len' columns. Found: {list(df.columns)}"
            )

    # Build bin edges and labels
    bin_edges = np.linspace(0, 100, bins + 1)
    bin_labels = [
        f"{int(bin_edges[i])}-{int(bin_edges[i+1])}%"
        for i in range(bins)
    ]

    if group_by and group_by in df.columns:
        # Grouped distribution
        results = []
        for group_val, group_df in df.groupby(group_by):
            dist = _compute_distribution(group_df["te_pct"], bin_edges, bin_labels)
            dist[group_by] = group_val
            results.append(dist)
        result_df = pd.concat(results, ignore_index=True)
        # Reorder columns
        cols = [group_by] + [c for c in result_df.columns if c != group_by]
        result_df = result_df[cols]
    else:
        result_df = _compute_distribution(df["te_pct"], bin_edges, bin_labels)

    out_path = os.path.join(output_dir, "te_distribution.csv")
    result_df.to_csv(out_path, index=False)
    logger.info(f"Saved TE distribution to {out_path}")

    # Overall summary
    summary = {
        "total_eccDNA": len(df),
        "mean_te_pct": round(df["te_pct"].mean(), 2),
        "median_te_pct": round(df["te_pct"].median(), 2),
        "std_te_pct": round(df["te_pct"].std(), 2),
        "min_te_pct": round(df["te_pct"].min(), 2),
        "max_te_pct": round(df["te_pct"].max(), 2),
        "n_bins": bins,
    }
    summary_df = pd.DataFrame([summary])
    summary_out = os.path.join(output_dir, "te_distribution_summary.csv")
    summary_df.to_csv(summary_out, index=False)
    logger.info(f"Saved distribution summary to {summary_out}")


def _compute_distribution(
    te_pct_series: pd.Series,
    bin_edges: np.ndarray,
    bin_labels: list[str],
) -> pd.DataFrame:
    """Compute distribution counts and percentages for a series of TE percentages."""
    total = len(te_pct_series)
    counts = []
    for i in range(len(bin_labels)):
        low = bin_edges[i]
        high = bin_edges[i + 1]
        if i == len(bin_labels) - 1:
            # Last bin is inclusive on the right
            n = ((te_pct_series >= low) & (te_pct_series <= high)).sum()
        else:
            n = ((te_pct_series >= low) & (te_pct_series < high)).sum()
        counts.append(n)

    return pd.DataFrame({
        "bin": bin_labels,
        "count": counts,
        "pct": [round(c / total * 100, 2) if total > 0 else 0 for c in counts],
    })
