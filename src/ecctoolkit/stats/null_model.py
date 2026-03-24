"""Null model analysis for eccDNA inter-replicate overlap.

Tests whether low inter-replicate overlap can be explained by incomplete
sampling from a common pool (sampling artifact) vs biologically private
eccDNA generation (micro-random model).

Methods:
  1. Exponential saturation model: y = N * (1 - exp(-c * x))
  2. Geometric extrapolation from final increments
  3. Lincoln-Petersen capture-recapture between replicates
  4. Expected overlap under random sampling (hypergeometric)
"""

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

logger = logging.getLogger(__name__)


def _exp_saturation(x: np.ndarray, n: float, c: float) -> np.ndarray:
    """Exponential saturation model: y = N * (1 - exp(-c * x))."""
    return n * (1 - np.exp(-c * x))


def _fit_exponential(
    fractions: np.ndarray,
    counts: np.ndarray,
) -> Dict[str, float]:
    """Fit exponential saturation model to estimate within-replicate pool size.

    Returns dict with N_estimated, c_estimated, r_squared, saturation_pct.
    """
    p0 = [counts[-1] * 1.5, 2.0]
    popt, pcov = curve_fit(
        _exp_saturation, fractions, counts, p0=p0, maxfev=10000,
    )
    n_est, c_est = popt
    perr = np.sqrt(np.diag(pcov))

    residuals = counts - _exp_saturation(fractions, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((counts - np.mean(counts)) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return {
        "method": "exponential",
        "N_estimated": float(n_est),
        "N_stderr": float(perr[0]),
        "c_estimated": float(c_est),
        "r_squared": float(r_squared),
        "saturation_pct": float(counts[-1] / n_est * 100) if n_est > 0 else 0.0,
    }


def _fit_geometric(
    fractions: np.ndarray,
    counts: np.ndarray,
) -> Dict[str, float]:
    """Geometric extrapolation from marginal discovery rates.

    Uses the ratio of the last 3 increments to estimate decay rate and
    extrapolate remaining eccDNA.
    """
    increments = np.diff(np.insert(counts, 0, 0))
    ratios = increments[1:] / np.maximum(increments[:-1], 1)
    last_3_ratios = ratios[-3:]
    mean_ratio = float(np.mean(last_3_ratios))
    last_increment = float(increments[-1])

    if 0 < mean_ratio < 1:
        remaining = last_increment * mean_ratio / (1 - mean_ratio)
    else:
        remaining = float("inf")

    n_geo = float(counts[-1]) + remaining
    saturation_pct = float(counts[-1]) / n_geo * 100 if np.isfinite(n_geo) and n_geo > 0 else 0.0

    return {
        "method": "geometric",
        "N_estimated": n_geo,
        "decay_ratio": mean_ratio,
        "last_increment": last_increment,
        "saturation_pct": saturation_pct,
    }


def _lincoln_petersen(
    n_a: int,
    n_b: int,
    n_shared: int,
) -> float:
    """Lincoln-Petersen capture-recapture estimator.

    N = (n_A * n_B) / m, where m is the number of shared elements.
    """
    if n_shared == 0:
        return float("inf")
    return (n_a * n_b) / n_shared


def _expected_overlap_fraction(
    n_a: int,
    n_b: int,
    pool_size: float,
) -> float:
    """Expected fraction of A found in B under random sampling from pool.

    Hypergeometric approximation: E[overlap] = n_B / N.
    """
    if pool_size <= 0 or not np.isfinite(pool_size):
        return 0.0
    return n_b / pool_size


def run_null_model(
    saturation_file: str,
    overlap_file: str,
    output_file: str,
) -> None:
    """
    Run null model analysis for eccDNA inter-replicate overlap.

    Estimates within-replicate pool size from saturation curve data,
    between-replicate pool size from Lincoln-Petersen capture-recapture,
    and compares expected vs observed overlap.

    Args:
        saturation_file: CSV with saturation curve data.
            Required columns: sample, fraction, n_eccdna.
            Multiple samples (replicates) should be present.
        overlap_file: CSV with pairwise overlap data.
            Required columns: sample_a, sample_b, n_a, n_b,
            n_shared, overlap_fraction.
        output_file: Output CSV for analysis results.
    """
    # Load saturation data
    sat_df = pd.read_csv(saturation_file)
    required_sat_cols = ["sample", "fraction", "n_eccdna"]
    missing = [c for c in required_sat_cols if c not in sat_df.columns]
    if missing:
        raise ValueError(
            f"Saturation file missing columns: {missing}. "
            f"Required: {required_sat_cols}"
        )

    # Load overlap data
    ov_df = pd.read_csv(overlap_file)
    required_ov_cols = ["sample_a", "sample_b", "n_a", "n_b", "n_shared", "overlap_fraction"]
    missing = [c for c in required_ov_cols if c not in ov_df.columns]
    if missing:
        raise ValueError(
            f"Overlap file missing columns: {missing}. "
            f"Required: {required_ov_cols}"
        )

    results: List[dict] = []

    # Analysis 1: Within-replicate pool size estimation
    logger.info("=== Analysis 1: Within-replicate pool size (saturation) ===")
    samples = sat_df["sample"].unique()
    within_estimates = {}

    for sample in samples:
        sample_data = sat_df[sat_df["sample"] == sample].sort_values("fraction")
        fractions = sample_data["fraction"].values
        counts = sample_data["n_eccdna"].values.astype(float)

        if len(fractions) < 3:
            logger.warning(f"Skipping {sample}: insufficient data points")
            continue

        # Method A: Exponential
        try:
            exp_result = _fit_exponential(fractions, counts)
            exp_result["sample"] = sample
            exp_result["analysis"] = "within_replicate"
            results.append(exp_result)
            within_estimates[sample] = exp_result["N_estimated"]
            logger.info(
                f"  {sample} (exponential): N = {exp_result['N_estimated']:,.0f}, "
                f"R^2 = {exp_result['r_squared']:.6f}, "
                f"saturation = {exp_result['saturation_pct']:.1f}%"
            )
        except (RuntimeError, ValueError) as e:
            logger.warning(f"  {sample} exponential fit failed: {e}")

        # Method B: Geometric
        geo_result = _fit_geometric(fractions, counts)
        geo_result["sample"] = sample
        geo_result["analysis"] = "within_replicate"
        results.append(geo_result)
        logger.info(
            f"  {sample} (geometric): N = {geo_result['N_estimated']:,.0f}, "
            f"decay = {geo_result['decay_ratio']:.4f}"
        )

    # Analysis 2: Between-replicate pool size (Lincoln-Petersen)
    logger.info("=== Analysis 2: Between-replicate pool size (capture-recapture) ===")
    lp_estimates = []

    for _, row in ov_df.iterrows():
        n_shared = int(row["n_shared"])
        n_a = int(row["n_a"])
        n_b = int(row["n_b"])
        pair_label = f"{row['sample_a']} vs {row['sample_b']}"

        n_lp = _lincoln_petersen(n_a, n_b, n_shared)
        lp_estimates.append(n_lp)

        results.append({
            "analysis": "between_replicate",
            "method": "lincoln_petersen",
            "pair": pair_label,
            "n_a": n_a,
            "n_b": n_b,
            "n_shared": n_shared,
            "N_estimated": n_lp,
        })

        logger.info(
            f"  {pair_label}: shared = {n_shared:,}, N_between = {n_lp:,.0f}"
        )

    if lp_estimates:
        mean_lp = np.mean([x for x in lp_estimates if np.isfinite(x)])
        logger.info(f"  Mean N_between = {mean_lp:,.0f}")

    # Analysis 3: Expected vs observed overlap
    logger.info("=== Analysis 3: Expected vs observed overlap ===")
    mean_within = np.mean(list(within_estimates.values())) if within_estimates else 0

    for _, row in ov_df.iterrows():
        n_b = int(row["n_b"])
        observed = float(row["overlap_fraction"])
        pair_label = f"{row['sample_a']} vs {row['sample_b']}"

        if mean_within > 0:
            expected = _expected_overlap_fraction(
                int(row["n_a"]), n_b, mean_within,
            )
            ratio = expected / observed if observed > 0 else float("inf")

            results.append({
                "analysis": "expected_overlap",
                "method": "hypergeometric_approx",
                "pair": pair_label,
                "N_pool": mean_within,
                "expected_overlap": expected,
                "observed_overlap": observed,
                "ratio_expected_to_observed": ratio,
            })

            logger.info(
                f"  {pair_label}: expected = {expected*100:.1f}%, "
                f"observed = {observed*100:.1f}%, "
                f"ratio = {ratio:.1f}x"
            )

    # Save results
    df_results = pd.DataFrame(results)
    df_results.to_csv(output_file, index=False)
    logger.info(f"Results saved: {output_file}")

    # Summary
    if mean_within > 0 and lp_estimates:
        mean_lp_finite = np.mean([x for x in lp_estimates if np.isfinite(x)])
        ratio = mean_lp_finite / mean_within if mean_within > 0 else float("inf")
        logger.info("=== Summary ===")
        logger.info(f"  Mean N_within (exponential) = {mean_within:,.0f}")
        logger.info(f"  Mean N_between (Lincoln-Petersen) = {mean_lp_finite:,.0f}")
        logger.info(f"  N_between / N_within = {ratio:.1f}x")
        if ratio > 3:
            logger.info(
                "  Conclusion: N_between >> N_within suggests biologically "
                "private eccDNA, NOT a sampling artifact."
            )
        else:
            logger.info(
                "  Conclusion: N_between ~ N_within is consistent with "
                "incomplete sampling from a common pool."
            )
