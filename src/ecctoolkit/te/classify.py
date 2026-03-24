"""Classify eccDNA by TE composition type."""

from __future__ import annotations

import logging
import os

import pandas as pd

logger = logging.getLogger(__name__)


def classify_te_composition(
    input_file: str,
    output_dir: str,
    min_te_fraction: float = 0.1,
    dominant_threshold: float = 0.8,
) -> None:
    """Classify eccDNA as single-TE, composite-TE, no-TE, or partial-TE.

    Uses per-eccDNA TE class breakdown (output of te-analyze) to determine
    the TE composition type of each eccDNA.

    Classification rules:
        - no_TE: total TE coverage < min_te_fraction (default 10%)
        - single_TE: one TE family covers >= dominant_threshold (default 80%)
        - composite_TE: 2+ TE families each covering >= min_te_fraction
        - partial_TE: has TE but doesn't fit above categories

    Args:
        input_file: Per-eccDNA TE annotation CSV (from te-analyze output,
            with columns: eccDNA_name, eccDNA_len, te_bp, te_pct, te_families).
            OR the raw overlaps CSV with eccDNA_name, te_family, overlap_bp columns.
        output_dir: Output directory.
        min_te_fraction: Minimum TE fraction (0-1) to count as TE-positive.
        dominant_threshold: Fraction threshold (0-1) for single-TE classification.
    """
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(input_file)
    logger.info(f"Loaded {len(df)} rows from {input_file}")

    # Determine input type: per-eccDNA summary or raw overlaps
    if "overlap_bp" in df.columns and "te_family" in df.columns:
        # Raw overlaps format -- aggregate per eccDNA + family
        per_ecc_summary, family_breakdown = _classify_from_overlaps(
            df, min_te_fraction, dominant_threshold
        )
    elif "te_pct" in df.columns:
        # Per-eccDNA summary format -- need overlaps for family breakdown
        # Try to load overlaps from same directory
        overlaps_path = os.path.join(os.path.dirname(input_file), "eccdna_te_overlaps.csv")
        if os.path.exists(overlaps_path):
            overlaps = pd.read_csv(overlaps_path)
            per_ecc_summary, family_breakdown = _classify_from_overlaps(
                overlaps, min_te_fraction, dominant_threshold
            )
            # Merge eccDNA with no TE from the summary
            no_te_names = set(df["eccDNA_name"]) - set(per_ecc_summary["eccDNA_name"])
            if no_te_names:
                no_te = df[df["eccDNA_name"].isin(no_te_names)][
                    ["eccDNA_name"]
                ].copy()
                no_te["classification"] = "no_TE"
                no_te["dominant_te"] = ""
                no_te["n_families"] = 0
                no_te["te_pct"] = 0.0
                per_ecc_summary = pd.concat(
                    [per_ecc_summary, no_te], ignore_index=True
                )
        else:
            # Fall back to simple classification from te_pct and te_families
            per_ecc_summary, family_breakdown = _classify_from_summary(
                df, min_te_fraction, dominant_threshold
            )
    else:
        raise ValueError(
            f"Input file must have either 'overlap_bp'+'te_family' columns "
            f"(raw overlaps) or 'te_pct'+'te_families' columns (per-eccDNA summary). "
            f"Found columns: {list(df.columns)}"
        )

    # Save classification
    class_out = os.path.join(output_dir, "te_classification.csv")
    per_ecc_summary.to_csv(class_out, index=False)
    logger.info(f"Saved TE classification to {class_out}")

    # Save family breakdown if available
    if family_breakdown is not None and not family_breakdown.empty:
        family_out = os.path.join(output_dir, "te_family_breakdown.csv")
        family_breakdown.to_csv(family_out, index=False)
        logger.info(f"Saved family breakdown to {family_out}")

    # Summary statistics
    counts = per_ecc_summary["classification"].value_counts()
    total = len(per_ecc_summary)
    logger.info(f"Classification results ({total} eccDNA):")
    for cls in ["single_TE", "composite_TE", "partial_TE", "no_TE"]:
        n = counts.get(cls, 0)
        logger.info(f"  {cls}: {n} ({n / total * 100:.1f}%)")

    # Save summary
    summary = counts.reset_index()
    summary.columns = ["classification", "count"]
    summary["pct"] = (summary["count"] / total * 100).round(2)
    summary_out = os.path.join(output_dir, "te_classification_summary.csv")
    summary.to_csv(summary_out, index=False)
    logger.info(f"Saved classification summary to {summary_out}")


def _classify_from_overlaps(
    overlaps: pd.DataFrame,
    min_te_fraction: float,
    dominant_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Classify eccDNA from raw overlap data."""
    # Family-level bp per eccDNA
    family_bp = overlaps.groupby(["eccDNA_name", "te_family"]).agg(
        family_bp=("overlap_bp", "sum"),
        eccDNA_len=("eccDNA_len", "first"),
    ).reset_index()
    family_bp["family_pct"] = family_bp["family_bp"] / family_bp["eccDNA_len"]

    results = []
    for name, group in family_bp.groupby("eccDNA_name"):
        ecc_len = group["eccDNA_len"].iloc[0]
        total_te_bp = group["family_bp"].sum()
        total_te_frac = min(total_te_bp / ecc_len, 1.0)

        if total_te_frac < min_te_fraction:
            results.append({
                "eccDNA_name": name,
                "classification": "no_TE",
                "dominant_te": "",
                "n_families": 0,
                "te_pct": round(total_te_frac * 100, 2),
            })
            continue

        # Compute family fractions relative to total TE coverage
        group = group.copy()
        group["rel_frac"] = group["family_bp"] / total_te_bp

        # Families with significant contribution
        significant = group[group["rel_frac"] >= min_te_fraction]
        n_sig = len(significant)

        dominant = group.loc[group["rel_frac"].idxmax()]
        dominant_frac = dominant["rel_frac"]

        if dominant_frac >= dominant_threshold:
            classification = "single_TE"
        elif n_sig >= 2:
            classification = "composite_TE"
        else:
            classification = "partial_TE"

        results.append({
            "eccDNA_name": name,
            "classification": classification,
            "dominant_te": dominant["te_family"],
            "n_families": n_sig,
            "te_pct": round(total_te_frac * 100, 2),
        })

    result_df = pd.DataFrame(results)
    return result_df, family_bp


def _classify_from_summary(
    df: pd.DataFrame,
    min_te_fraction: float,
    dominant_threshold: float,
) -> tuple[pd.DataFrame, None]:
    """Classify eccDNA from per-eccDNA summary (fallback when overlaps not available)."""
    results = []
    for _, row in df.iterrows():
        name = row["eccDNA_name"]
        te_pct = row.get("te_pct", 0)
        te_families_str = str(row.get("te_families", ""))

        if te_pct < min_te_fraction * 100:
            classification = "no_TE"
            dominant_te = ""
            n_families = 0
        else:
            families = [f.strip() for f in te_families_str.split(",") if f.strip()]
            n_families = len(families)
            if n_families == 1:
                classification = "single_TE"
                dominant_te = families[0]
            elif n_families >= 2:
                classification = "composite_TE"
                dominant_te = families[0]
            else:
                classification = "partial_TE"
                dominant_te = ""

        results.append({
            "eccDNA_name": name,
            "classification": classification,
            "dominant_te": dominant_te,
            "n_families": n_families,
            "te_pct": te_pct,
        })

    return pd.DataFrame(results), None
