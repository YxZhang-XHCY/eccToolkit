"""Composite TE composition and motif combination analysis."""

from __future__ import annotations

import logging
import os

import pandas as pd

logger = logging.getLogger(__name__)


def analyze_te_composition(
    input_file: str,
    output_dir: str,
    min_count: int = 2,
) -> None:
    """Analyze composite TE motif combinations and their frequencies.

    For each eccDNA classified as composite-TE (or containing multiple TE
    families), extracts the combination of TE families and counts how often
    each combination occurs.

    Args:
        input_file: Either te_classification.csv (from te-classify) or
            eccdna_te_overlaps.csv / te_family_breakdown.csv (from te-analyze).
            Must contain eccDNA_name and te_family columns.
        output_dir: Output directory.
        min_count: Minimum count to include a combination in the output.
    """
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(input_file)
    logger.info(f"Loaded {len(df)} rows from {input_file}")

    # Determine format and extract per-eccDNA family lists
    if "te_family" in df.columns and "eccDNA_name" in df.columns:
        # Overlaps or family breakdown format: multiple rows per eccDNA
        # Optionally filter to only significant families
        if "family_pct" in df.columns:
            # family_breakdown format: keep families with >= 10% relative coverage
            df = df[df["family_pct"] >= 0.1]

        family_lists = (
            df.groupby("eccDNA_name")["te_family"]
            .apply(lambda x: tuple(sorted(set(x))))
            .reset_index()
        )
        family_lists.columns = ["eccDNA_name", "te_combination"]
    elif "te_families" in df.columns:
        # Per-eccDNA summary format: te_families is comma-separated
        df = df[df["te_families"].notna() & (df["te_families"] != "")]
        family_lists = df[["eccDNA_name", "te_families"]].copy()
        family_lists["te_combination"] = family_lists["te_families"].apply(
            lambda x: tuple(sorted(set(s.strip() for s in str(x).split(",") if s.strip())))
        )
        family_lists = family_lists[["eccDNA_name", "te_combination"]]
    elif "classification" in df.columns and "dominant_te" in df.columns:
        # Classification output: need to load overlaps/breakdown separately
        overlaps_path = os.path.join(os.path.dirname(input_file), "te_family_breakdown.csv")
        if not os.path.exists(overlaps_path):
            overlaps_path = os.path.join(os.path.dirname(input_file), "eccdna_te_overlaps.csv")
        if os.path.exists(overlaps_path):
            return analyze_te_composition(overlaps_path, output_dir, min_count)
        raise ValueError(
            f"Cannot find family breakdown data. Provide overlaps or family breakdown CSV."
        )
    else:
        raise ValueError(
            f"Input file must have 'te_family' or 'te_families' column. "
            f"Found: {list(df.columns)}"
        )

    # Filter to multi-family eccDNA (composite)
    family_lists["n_families"] = family_lists["te_combination"].apply(len)
    composite = family_lists[family_lists["n_families"] >= 2].copy()

    if composite.empty:
        logger.warning("No composite-TE eccDNA found")
        pd.DataFrame(columns=["te_combination", "count", "pct"]).to_csv(
            os.path.join(output_dir, "te_combinations.csv"), index=False
        )
        return

    # Count combination frequencies
    combo_str = composite["te_combination"].apply(lambda x: "+".join(x))
    combo_counts = combo_str.value_counts().reset_index()
    combo_counts.columns = ["te_combination", "count"]
    combo_counts = combo_counts[combo_counts["count"] >= min_count]
    combo_counts["pct"] = (
        combo_counts["count"] / len(composite) * 100
    ).round(2)
    combo_counts["n_families"] = combo_counts["te_combination"].apply(
        lambda x: len(x.split("+"))
    )
    combo_counts = combo_counts.sort_values("count", ascending=False).reset_index(drop=True)

    combo_out = os.path.join(output_dir, "te_combinations.csv")
    combo_counts.to_csv(combo_out, index=False)
    logger.info(f"Saved {len(combo_counts)} TE combinations to {combo_out}")

    # Per-eccDNA combination assignment
    composite["te_combination_str"] = combo_str.values
    per_ecc_out = os.path.join(output_dir, "per_eccdna_te_combination.csv")
    composite[["eccDNA_name", "te_combination_str", "n_families"]].to_csv(
        per_ecc_out, index=False
    )
    logger.info(f"Saved per-eccDNA combinations to {per_ecc_out}")

    # Summary
    logger.info(
        f"Composite TE analysis: {len(composite)} composite eccDNA, "
        f"{len(combo_counts)} unique combinations (min_count={min_count})"
    )
    for _, row in combo_counts.head(10).iterrows():
        logger.info(f"  {row['te_combination']}: {row['count']} ({row['pct']}%)")
