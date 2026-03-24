"""Quantify eccDNA length distribution patterns."""

import logging
import os

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _compute_length_stats(lengths: pd.Series) -> dict:
    """Compute basic statistics for eccDNA lengths."""
    return {
        "count": len(lengths),
        "mean": lengths.mean(),
        "median": lengths.median(),
        "std": lengths.std(),
        "min": lengths.min(),
        "max": lengths.max(),
        "q25": lengths.quantile(0.25),
        "q75": lengths.quantile(0.75),
    }


def _histogram_binning(
    lengths: pd.Series, bin_size: int, max_length: int
) -> pd.DataFrame:
    """Bin length distribution into histogram."""
    bins = np.arange(0, max_length + bin_size, bin_size)
    counts, bin_edges = np.histogram(lengths.clip(upper=max_length), bins=bins)

    return pd.DataFrame({
        "bin_start": bin_edges[:-1].astype(int),
        "bin_end": bin_edges[1:].astype(int),
        "count": counts,
    })


def _periodicity_analysis(
    hist_counts: np.ndarray, bin_size: int
) -> pd.DataFrame:
    """Perform periodicity analysis using ACF and FFT.

    Returns DataFrame with dominant periodicities.
    """
    if len(hist_counts) < 10:
        logger.warning("Too few bins for periodicity analysis")
        return pd.DataFrame(columns=["period_bp", "fft_magnitude", "acf_value"])

    # Normalize counts
    counts = hist_counts.astype(float)
    counts = counts - counts.mean()

    # ACF (autocorrelation function)
    n = len(counts)
    acf_values = np.correlate(counts, counts, mode="full")
    acf_values = acf_values[n - 1:]  # Keep positive lags only
    if acf_values[0] > 0:
        acf_values = acf_values / acf_values[0]  # Normalize

    # FFT
    fft_result = np.fft.rfft(counts)
    fft_magnitude = np.abs(fft_result)
    freqs = np.fft.rfftfreq(n)

    # Find peaks in FFT (skip DC component at index 0)
    if len(fft_magnitude) > 1:
        mag = fft_magnitude[1:]  # Skip DC
        freq = freqs[1:]

        # Find top peaks
        peak_indices = np.argsort(mag)[::-1]
        results = []
        for idx in peak_indices[:10]:  # Top 10
            f = freq[idx]
            if f > 0:
                period_bins = 1.0 / f
                period_bp = period_bins * bin_size
                acf_lag = int(round(period_bins))
                acf_val = acf_values[acf_lag] if acf_lag < len(acf_values) else np.nan

                results.append({
                    "period_bp": round(period_bp, 1),
                    "period_bins": round(period_bins, 2),
                    "fft_magnitude": round(mag[idx], 4),
                    "fft_frequency": round(f, 6),
                    "acf_value": round(acf_val, 4) if not np.isnan(acf_val) else None,
                })

        return pd.DataFrame(results)

    return pd.DataFrame(columns=["period_bp", "fft_magnitude", "acf_value"])


def run_quantification(
    input_file: str,
    output_dir: str,
    bin_size: int = 100,
    max_length: int = 50000,
) -> None:
    """
    Quantify eccDNA length distribution with periodicity analysis.

    Performs ACF (autocorrelation) and FFT analysis to detect
    periodic patterns in eccDNA length distribution.

    Args:
        input_file: eccDNA CSV file with eLength or length column
        output_dir: Output directory for results
        bin_size: Histogram bin size in bp (default: 100)
        max_length: Maximum length to consider in bp (default: 50000)
    """
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    logger.info(f"Loading eccDNA data from {input_file}")
    sep = "\t" if input_file.endswith((".tsv", ".bed", ".txt")) else ","
    df = pd.read_csv(input_file, sep=sep)

    # Find length column
    length_col = None
    for c in df.columns:
        cl = c.lower().strip()
        if cl in ("elength", "length", "size", "ecc_length"):
            length_col = c
            break

    if length_col is None:
        # Try computing from start/end
        col_map = {}
        for c in df.columns:
            cl = c.lower().strip()
            if cl in ("start", "estart"):
                col_map["start"] = c
            elif cl in ("end", "eend"):
                col_map["end"] = c
        if "start" in col_map and "end" in col_map:
            df["_length"] = df[col_map["end"]] - df[col_map["start"]]
            length_col = "_length"
        else:
            raise ValueError(
                f"Cannot find length column in {input_file}. "
                f"Available columns: {list(df.columns)}"
            )

    lengths = pd.to_numeric(df[length_col], errors="coerce").dropna()
    lengths = lengths[lengths > 0]
    logger.info(f"Loaded {len(lengths)} valid length values")

    # Basic statistics
    stats = _compute_length_stats(lengths)
    logger.info(
        f"Length stats: mean={stats['mean']:.0f}, median={stats['median']:.0f}, "
        f"range=[{stats['min']:.0f}, {stats['max']:.0f}]"
    )

    stats_df = pd.DataFrame([stats])
    stats_out = os.path.join(output_dir, "length_statistics.csv")
    stats_df.to_csv(stats_out, index=False)
    logger.info(f"Saved length statistics to {stats_out}")

    # Histogram
    hist_df = _histogram_binning(lengths, bin_size, max_length)
    hist_out = os.path.join(output_dir, "length_distribution.csv")
    hist_df.to_csv(hist_out, index=False)
    logger.info(f"Saved length distribution ({len(hist_df)} bins) to {hist_out}")

    # Periodicity analysis
    logger.info("Running periodicity analysis (ACF + FFT)")
    period_df = _periodicity_analysis(hist_df["count"].values, bin_size)

    if len(period_df) > 0:
        period_out = os.path.join(output_dir, "periodicity_results.csv")
        period_df.to_csv(period_out, index=False)
        logger.info(f"Saved periodicity results to {period_out}")

        # Report top periods
        top = period_df.head(3)
        for _, row in top.iterrows():
            logger.info(
                f"  Period: {row['period_bp']:.0f} bp "
                f"(FFT magnitude: {row['fft_magnitude']:.4f})"
            )
    else:
        logger.info("No significant periodicities detected")

    logger.info("Length distribution quantification completed!")
