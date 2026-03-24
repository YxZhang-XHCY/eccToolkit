"""CNV-controlled correlation analysis of eccDNA genomic distribution.

Tests whether inter-replicate eccDNA density correlations are driven by
copy number variation (CNV) through three analyses:
  1. CNV-normalized density correlation at multiple window sizes
  2. Chromosome exclusion (remove strongest CNV outlier)
  3. Intra-chromosome correlation vs CNV
"""

import logging
import os
from collections import defaultdict
from itertools import combinations
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

logger = logging.getLogger(__name__)

# Default window sizes for multi-scale analysis
DEFAULT_WINDOW_SIZES = [10_000, 50_000, 100_000, 500_000, 1_000_000]


def _format_window_size(ws: int) -> str:
    """Format window size for display (e.g., 10000 -> '10kb')."""
    if ws >= 1_000_000:
        return f"{ws // 1_000_000}Mb"
    return f"{ws // 1_000}kb"


def _load_genome_file(genome_file: str) -> Dict[str, int]:
    """Load genome file (chrom<TAB>size) into dict."""
    chrom_sizes = {}
    with open(genome_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes


def _build_bin_index(
    window_size: int,
    chrom_sizes: Dict[str, int],
    chroms: Optional[List[str]] = None,
) -> Tuple[list, dict]:
    """Build genomic window bins and index.

    Returns:
        bins: list of (chrom, bin_idx) tuples
        idx_map: dict mapping (chrom, bin_idx) to flat index
    """
    if chroms is None:
        chroms = sorted(chrom_sizes.keys(), key=lambda c: (len(c), c))
    else:
        chroms = [c for c in chroms if c in chrom_sizes]

    bins = []
    for chrom in chroms:
        n_bins = int(np.ceil(chrom_sizes[chrom] / window_size))
        for i in range(n_bins):
            bins.append((chrom, i))

    idx_map = {b: i for i, b in enumerate(bins)}
    return bins, idx_map


def _load_eccdna_midpoints(
    filepath: str,
    chrom_sizes: Dict[str, int],
    chr_col: str = "eChr",
    start_col: str = "eStart",
    end_col: str = "eEnd",
) -> List[Tuple[str, int]]:
    """Load eccDNA records and compute midpoints.

    Returns list of (chrom, midpoint) tuples.
    """
    df = pd.read_csv(filepath)

    # Auto-detect column names
    if chr_col not in df.columns:
        for alt in ["chr", "Chr", "chrom", "chromosome"]:
            if alt in df.columns:
                chr_col = alt
                break
    if start_col not in df.columns:
        for alt in ["start", "Start", "chromStart"]:
            if alt in df.columns:
                start_col = alt
                break
    if end_col not in df.columns:
        for alt in ["end", "End", "chromEnd"]:
            if alt in df.columns:
                end_col = alt
                break

    df = df[df[chr_col].isin(chrom_sizes)].copy()
    df[start_col] = pd.to_numeric(df[start_col], errors="coerce")
    df[end_col] = pd.to_numeric(df[end_col], errors="coerce")
    df = df.dropna(subset=[start_col, end_col])

    midpoints = list(
        zip(df[chr_col], ((df[start_col] + df[end_col]) / 2).astype(int))
    )
    return midpoints


def _bin_midpoints(
    midpoints: List[Tuple[str, int]],
    window_size: int,
    idx_map: dict,
    total_bins: int,
) -> np.ndarray:
    """Count eccDNA midpoints per genomic window."""
    counts = np.zeros(total_bins, dtype=np.float64)
    for chrom, mid in midpoints:
        bin_idx = mid // window_size
        key = (chrom, bin_idx)
        if key in idx_map:
            counts[idx_map[key]] += 1
    return counts


def _load_cnv_bed(cnv_file: str, chrom_sizes: Dict[str, int]) -> list:
    """Load CNV data from BED file (chrom, start, end, copy_ratio)."""
    segments = []
    with open(cnv_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            chrom = fields[0]
            if chrom not in chrom_sizes:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
                ratio = float(fields[3])
                segments.append((chrom, start, end, ratio))
            except (ValueError, IndexError):
                continue
    return segments


def _cnv_to_window_values(
    segments: list,
    window_size: int,
    bins: list,
    idx_map: dict,
    chrom_sizes: Dict[str, int],
) -> np.ndarray:
    """Compute weighted-average CNV ratio per genomic window."""
    total_bins = len(bins)
    weighted_sum = np.zeros(total_bins, dtype=np.float64)
    overlap_bp = np.zeros(total_bins, dtype=np.float64)

    for chrom, seg_start, seg_end, ratio in segments:
        first_bin = seg_start // window_size
        last_bin = (seg_end - 1) // window_size
        for b in range(first_bin, last_bin + 1):
            key = (chrom, b)
            if key not in idx_map:
                continue
            flat_idx = idx_map[key]
            bin_start = b * window_size
            bin_end = min((b + 1) * window_size, chrom_sizes.get(chrom, 0))
            ov_start = max(seg_start, bin_start)
            ov_end = min(seg_end, bin_end)
            ov = max(0, ov_end - ov_start)
            weighted_sum[flat_idx] += ratio * ov
            overlap_bp[flat_idx] += ov

    cnv_values = np.ones(total_bins, dtype=np.float64)
    covered = overlap_bp > 0
    cnv_values[covered] = weighted_sum[covered] / overlap_bp[covered]
    return cnv_values


def _compute_chromosome_mean_cnv(
    segments: list,
    chrom_order: List[str],
) -> Dict[str, float]:
    """Compute per-chromosome weighted-average CNV ratio."""
    chrom_data: Dict[str, Dict[str, float]] = defaultdict(
        lambda: {"weighted_sum": 0.0, "total_bp": 0}
    )
    for chrom, start, end, ratio in segments:
        span = end - start
        chrom_data[chrom]["weighted_sum"] += ratio * span
        chrom_data[chrom]["total_bp"] += span

    result = {}
    for chrom in chrom_order:
        d = chrom_data.get(chrom)
        if d and d["total_bp"] > 0:
            result[chrom] = d["weighted_sum"] / d["total_bp"]
        else:
            result[chrom] = 1.0
    return result


def _safe_pearsonr(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """Pearson r that returns NaN for constant arrays."""
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan, 1.0
    return pearsonr(x, y)


def _safe_spearmanr(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """Spearman rho that returns NaN for constant arrays."""
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan, 1.0
    return spearmanr(x, y)


def run_cnv_correlation(
    eccdna_files: List[str],
    cnv_file: str,
    output_dir: str,
    genome_file: str,
    window_sizes: Optional[List[int]] = None,
    sample_names: Optional[List[str]] = None,
) -> None:
    """
    Run CNV-controlled correlation analysis of eccDNA distribution.

    Computes raw and CNV-normalized inter-replicate correlations at
    multiple genomic window sizes, tests chromosome exclusion effects,
    and analyzes intra-chromosome correlation vs CNV.

    Args:
        eccdna_files: List of eccDNA CSV files (one per replicate)
        cnv_file: CNV BED file (chrom, start, end, copy_ratio)
        output_dir: Output directory
        genome_file: Genome file (chrom<TAB>size)
        window_sizes: Window sizes in bp (default: 10kb to 1Mb)
        sample_names: Sample names (default: derived from filenames)
    """
    if len(eccdna_files) < 2:
        raise ValueError("At least 2 eccDNA files are required")

    os.makedirs(output_dir, exist_ok=True)

    if window_sizes is None:
        window_sizes = DEFAULT_WINDOW_SIZES

    if sample_names is None:
        from pathlib import Path
        sample_names = [Path(f).stem for f in eccdna_files]

    chrom_sizes = _load_genome_file(genome_file)
    chrom_order = sorted(chrom_sizes.keys(), key=lambda c: (len(c), c))

    logger.info("=== CNV-Controlled Correlation Analysis ===")
    logger.info(f"Samples: {sample_names}")
    logger.info(f"CNV file: {cnv_file}")
    logger.info(f"Window sizes: {[_format_window_size(w) for w in window_sizes]}")

    # Load eccDNA midpoints
    midpoints = {}
    for filepath, name in zip(eccdna_files, sample_names):
        mps = _load_eccdna_midpoints(filepath, chrom_sizes)
        midpoints[name] = mps
        logger.info(f"  {name}: {len(mps):,} eccDNA")

    # Load CNV data
    cnv_segments = _load_cnv_bed(cnv_file, chrom_sizes)
    cnv_chrom_mean = _compute_chromosome_mean_cnv(cnv_segments, chrom_order)
    logger.info(f"  CNV: {len(cnv_segments)} segments loaded")

    # Generate sample pairs
    pairs = list(combinations(sample_names, 2))
    all_results: List[dict] = []

    CNV_FLOOR = 0.1

    # Analysis 1: Multi-scale correlation (raw and CNV-normalized)
    logger.info("Analysis 1: Multi-scale CNV-normalized correlation...")
    for ws in window_sizes:
        ws_label = _format_window_size(ws)
        bins, idx_map = _build_bin_index(ws, chrom_sizes)
        total_bins = len(bins)

        # Bin eccDNA counts
        binned_raw = {}
        for name in sample_names:
            binned_raw[name] = _bin_midpoints(
                midpoints[name], ws, idx_map, total_bins,
            )

        # Compute CNV per window
        cnv_window = _cnv_to_window_values(
            cnv_segments, ws, bins, idx_map, chrom_sizes,
        )

        # CNV-normalize
        binned_norm = {}
        cnv_floored = np.maximum(cnv_window, CNV_FLOOR)
        for name in sample_names:
            binned_norm[name] = binned_raw[name] / cnv_floored

        for s1, s2 in pairs:
            pr_raw, _ = _safe_pearsonr(binned_raw[s1], binned_raw[s2])
            sr_raw, _ = _safe_spearmanr(binned_raw[s1], binned_raw[s2])
            pr_norm, _ = _safe_pearsonr(binned_norm[s1], binned_norm[s2])
            sr_norm, _ = _safe_spearmanr(binned_norm[s1], binned_norm[s2])

            all_results.append({
                "analysis": "multiscale",
                "window_size": ws_label,
                "window_size_bp": ws,
                "pair": f"{s1} vs {s2}",
                "chr_set": "all",
                "pearson_raw": round(pr_raw, 6),
                "spearman_raw": round(sr_raw, 6),
                "pearson_cnv_norm": round(pr_norm, 6),
                "spearman_cnv_norm": round(sr_norm, 6),
            })

            logger.info(
                f"  {ws_label} {s1} vs {s2}: "
                f"raw r={pr_raw:.4f}, norm r={pr_norm:.4f}"
            )

    # Analysis 2: Chromosome exclusion
    logger.info("Analysis 2: Chromosome exclusion...")

    # Find chromosome with highest CNV deviation
    cnv_deviations = {
        c: abs(cnv_chrom_mean.get(c, 1.0) - 1.0) for c in chrom_order
    }
    if cnv_deviations:
        exclude_chrom = max(cnv_deviations, key=cnv_deviations.get)
        logger.info(
            f"  Excluding {exclude_chrom} "
            f"(CNV = {cnv_chrom_mean.get(exclude_chrom, 1.0):.3f})"
        )

        remaining_chroms = [c for c in chrom_order if c != exclude_chrom]

        for ws in window_sizes:
            ws_label = _format_window_size(ws)
            bins_excl, idx_map_excl = _build_bin_index(
                ws, chrom_sizes, chroms=remaining_chroms,
            )
            total_bins_excl = len(bins_excl)

            for s1, s2 in pairs:
                v1 = _bin_midpoints(
                    midpoints[s1], ws, idx_map_excl, total_bins_excl,
                )
                v2 = _bin_midpoints(
                    midpoints[s2], ws, idx_map_excl, total_bins_excl,
                )

                pr, _ = _safe_pearsonr(v1, v2)
                sr, _ = _safe_spearmanr(v1, v2)

                all_results.append({
                    "analysis": "chr_excluded",
                    "window_size": ws_label,
                    "window_size_bp": ws,
                    "pair": f"{s1} vs {s2}",
                    "chr_set": f"excl_{exclude_chrom}",
                    "pearson_raw": round(pr, 6),
                    "spearman_raw": round(sr, 6),
                    "pearson_cnv_norm": np.nan,
                    "spearman_cnv_norm": np.nan,
                })

    # Analysis 3: Intra-chromosome correlation at 50kb
    logger.info("Analysis 3: Intra-chromosome correlation at 50kb...")
    INTRA_WS = 50_000
    intra_results: List[dict] = []

    for chrom in chrom_order:
        bins_chr, idx_map_chr = _build_bin_index(
            INTRA_WS, chrom_sizes, chroms=[chrom],
        )
        total_bins_chr = len(bins_chr)

        if total_bins_chr < 20:
            continue

        pair_pearsons = []
        for s1, s2 in pairs:
            v1 = _bin_midpoints(midpoints[s1], INTRA_WS, idx_map_chr, total_bins_chr)
            v2 = _bin_midpoints(midpoints[s2], INTRA_WS, idx_map_chr, total_bins_chr)
            pr, _ = _safe_pearsonr(v1, v2)
            pair_pearsons.append(pr)

            all_results.append({
                "analysis": "intra_chromosome",
                "window_size": "50kb",
                "window_size_bp": INTRA_WS,
                "pair": f"{s1} vs {s2}",
                "chr_set": chrom,
                "pearson_raw": round(pr, 6),
                "spearman_raw": np.nan,
                "pearson_cnv_norm": np.nan,
                "spearman_cnv_norm": np.nan,
            })

        mean_pr = float(np.nanmean(pair_pearsons))
        chr_cnv = cnv_chrom_mean.get(chrom, 1.0)

        intra_results.append({
            "chrom": chrom,
            "n_windows": total_bins_chr,
            "mean_pearson_r": round(mean_pr, 6),
            "chr_cnv_ratio": round(chr_cnv, 4),
        })

        logger.debug(f"  {chrom}: CNV={chr_cnv:.3f}, mean r={mean_pr:.4f}")

    # Save results
    df_all = pd.DataFrame(all_results)
    output_csv = os.path.join(output_dir, "cnv_controlled_correlation.csv")
    df_all.to_csv(output_csv, index=False)
    logger.info(f"Correlation results saved: {output_csv}")

    if intra_results:
        df_intra = pd.DataFrame(intra_results)
        intra_csv = os.path.join(output_dir, "cnv_controlled_intra_chr.csv")
        df_intra.to_csv(intra_csv, index=False)
        logger.info(f"Intra-chromosome results saved: {intra_csv}")

    # Summary
    logger.info("=== Summary ===")
    df_ms = df_all[df_all["analysis"] == "multiscale"]
    for ws in window_sizes:
        ws_sub = df_ms[df_ms["window_size_bp"] == ws]
        if ws_sub.empty:
            continue
        mean_raw = ws_sub["pearson_raw"].mean()
        mean_norm = ws_sub["pearson_cnv_norm"].mean()
        logger.info(
            f"  {_format_window_size(ws)}: "
            f"raw r={mean_raw:.4f}, CNV-norm r={mean_norm:.4f}, "
            f"delta={mean_norm-mean_raw:+.4f}"
        )

    logger.info("=== CNV-Controlled Correlation Analysis Complete ===")
