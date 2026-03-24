"""Window-based feature enrichment analysis for eccDNA hotspots."""

from __future__ import annotations

import logging
import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)


def _load_genome_sizes(genome_file: str) -> Dict[str, int]:
    """Load chromosome sizes from genome file (chrom<tab>size)."""
    sizes = {}
    with open(genome_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 2:
                try:
                    sizes[fields[0]] = int(fields[1])
                except ValueError:
                    continue
    return sizes


def _build_windows(
    chrom_sizes: Dict[str, int], window_size: int
) -> List[Tuple[str, int, int]]:
    """Build genomic windows of specified size."""
    chrom_order = sorted(
        chrom_sizes.keys(), key=lambda c: (len(c), c)
    )
    windows = []
    for chrom in chrom_order:
        size = chrom_sizes[chrom]
        for start in range(0, size, window_size):
            end = min(start + window_size, size)
            windows.append((chrom, start, end))
    return windows


def _load_bed_regions(filepath: str) -> List[Tuple[str, int, int]]:
    """Load regions from BED-like format (handles plain BED, narrowPeak, etc.)."""
    regions = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            # Standard BED: chrom, start, end
            if fields[0].startswith("chr"):
                chrom, start, end = fields[0], fields[1], fields[2]
            # UCSC table format: bin, chrom, start, end
            elif len(fields) >= 4 and fields[1].startswith("chr"):
                chrom, start, end = fields[1], fields[2], fields[3]
            else:
                continue
            try:
                regions.append((chrom, int(start), int(end)))
            except ValueError:
                continue
    return regions


def _count_eccdna_per_window(
    eccdna_file: str,
    windows: List[Tuple[str, int, int]],
    window_size: int,
    chrom_sizes: Dict[str, int],
) -> Dict[Tuple[str, int, int], int]:
    """Count eccDNA regions per window using midpoint assignment."""
    ext = Path(eccdna_file).suffix.lower()
    if ext == ".bed":
        df = pd.read_csv(eccdna_file, sep="\t", header=None, comment="#")
        df = df.iloc[:, :3]
        df.columns = ["chrom", "start", "end"]
    else:
        df = pd.read_csv(eccdna_file)
        chr_col = next(
            (c for c in ["eChr", "chr", "chrom"] if c in df.columns), None
        )
        start_col = next(
            (c for c in ["eStart", "start", "chromStart"] if c in df.columns), None
        )
        end_col = next(
            (c for c in ["eEnd", "end", "chromEnd"] if c in df.columns), None
        )
        if chr_col is None or start_col is None or end_col is None:
            raise ValueError(f"Cannot find chr/start/end columns in {eccdna_file}")
        df = df[[chr_col, start_col, end_col]].copy()
        df.columns = ["chrom", "start", "end"]

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna().astype({"start": int, "end": int})

    # Create window set for fast lookup
    window_set = set(windows)

    counts: Dict[Tuple[str, int, int], int] = {w: 0 for w in windows}
    for _, row in df.iterrows():
        chrom = str(row["chrom"])
        mid = (row["start"] + row["end"]) // 2
        w_idx = mid // window_size
        w_start = w_idx * window_size
        w_end = min((w_idx + 1) * window_size, chrom_sizes.get(chrom, 0))
        w = (chrom, w_start, w_end)
        if w in counts:
            counts[w] += 1

    return counts


def _coverage_per_window(
    regions: List[Tuple[str, int, int]],
    windows: List[Tuple[str, int, int]],
) -> Dict[Tuple[str, int, int], float]:
    """Compute fraction of each window covered by regions."""
    chrom_regions: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for chrom, start, end in regions:
        chrom_regions[chrom].append((start, end))
    for chrom in chrom_regions:
        chrom_regions[chrom].sort()

    coverage = {}
    for chrom, wstart, wend in windows:
        covered_bp = 0
        for rstart, rend in chrom_regions.get(chrom, []):
            if rstart >= wend:
                break
            if rend > wstart:
                ov_start = max(rstart, wstart)
                ov_end = min(rend, wend)
                covered_bp += ov_end - ov_start
        coverage[(chrom, wstart, wend)] = covered_bp / (wend - wstart) if wend > wstart else 0.0

    return coverage


def _identify_hotspots(
    eccdna_counts: Dict[Tuple[str, int, int], int],
    threshold: float,
) -> set:
    """Identify hotspot windows based on fold-above-median threshold."""
    values = np.array(list(eccdna_counts.values()), dtype=float)
    nonzero = values[values > 0]
    if len(nonzero) == 0:
        return set()
    median_val = np.median(nonzero)
    if median_val <= 0:
        return set()
    hotspots = set()
    for w, count in eccdna_counts.items():
        if count >= threshold * median_val:
            hotspots.add(w)
    return hotspots


def run_feature_enrichment(
    eccdna_file: str,
    feature_files: Dict[str, str],
    output_dir: str,
    genome_file: str,
    window_size: int = 100000,
    hotspot_file: Optional[str] = None,
    hotspot_threshold: float = 3.0,
) -> None:
    """Analyze genomic feature enrichment in eccDNA hotspot windows.

    For each genomic window, computes feature coverage/count, then compares
    hotspot vs non-hotspot windows using Mann-Whitney U test.

    Args:
        eccdna_file: eccDNA regions CSV/BED file.
        feature_files: Mapping of feature name to BED file path.
        output_dir: Output directory.
        genome_file: Genome sizes file (chrom<tab>size).
        window_size: Genomic window size in bp (default: 100kb).
        hotspot_file: Optional pre-computed hotspot regions BED file.
            If not provided, hotspots are identified from eccdna_file.
        hotspot_threshold: Fold above median for hotspot identification
            (only used when hotspot_file is not provided).
    """
    os.makedirs(output_dir, exist_ok=True)

    # Load genome and build windows
    chrom_sizes = _load_genome_sizes(genome_file)
    windows = _build_windows(chrom_sizes, window_size)
    logger.info(f"Built {len(windows)} windows of {window_size // 1000}kb")

    # Determine hotspot vs non-hotspot windows
    if hotspot_file:
        # Load pre-computed hotspot regions
        hotspot_regions = _load_bed_regions(hotspot_file)
        hotspot_windows = set()
        for chrom, start, end in hotspot_regions:
            # Map to windows
            w_start_idx = start // window_size
            w_end_idx = (end - 1) // window_size
            for w_idx in range(w_start_idx, w_end_idx + 1):
                ws = w_idx * window_size
                we = min((w_idx + 1) * window_size, chrom_sizes.get(chrom, 0))
                w = (chrom, ws, we)
                if w in set(windows):
                    hotspot_windows.add(w)
        logger.info(f"Loaded {len(hotspot_windows)} hotspot windows from {hotspot_file}")
    else:
        # Identify hotspots from eccDNA density
        eccdna_counts = _count_eccdna_per_window(
            eccdna_file, windows, window_size, chrom_sizes
        )
        hotspot_windows = _identify_hotspots(eccdna_counts, hotspot_threshold)
        logger.info(
            f"Identified {len(hotspot_windows)} hotspot windows "
            f"(>= {hotspot_threshold}x median)"
        )

    non_hotspot_windows = set(windows) - hotspot_windows
    if len(hotspot_windows) < 5 or len(non_hotspot_windows) < 5:
        logger.error(
            f"Too few windows for comparison: "
            f"{len(hotspot_windows)} hotspot, {len(non_hotspot_windows)} non-hotspot"
        )
        return

    # Load features and compute enrichment
    all_results = []
    feature_matrix_data = {"chrom": [], "start": [], "end": [], "is_hotspot": []}

    for w in windows:
        feature_matrix_data["chrom"].append(w[0])
        feature_matrix_data["start"].append(w[1])
        feature_matrix_data["end"].append(w[2])
        feature_matrix_data["is_hotspot"].append(1 if w in hotspot_windows else 0)

    for feat_name, feat_file in sorted(feature_files.items()):
        logger.info(f"Processing feature: {feat_name}")

        if not os.path.exists(feat_file):
            logger.warning(f"  Feature file not found: {feat_file}")
            continue

        regions = _load_bed_regions(feat_file)
        if len(regions) == 0:
            logger.warning(f"  No regions loaded from {feat_file}")
            continue
        logger.info(f"  Loaded {len(regions)} regions")

        cov = _coverage_per_window(regions, windows)

        # Store in matrix
        feature_matrix_data[feat_name] = [cov.get(w, 0.0) for w in windows]

        # Split by hotspot status
        hot_vals = np.array([cov.get(w, 0.0) for w in windows if w in hotspot_windows])
        nonhot_vals = np.array(
            [cov.get(w, 0.0) for w in windows if w in non_hotspot_windows]
        )

        mean_hot = np.mean(hot_vals)
        mean_nonhot = np.mean(nonhot_vals)
        fold_enrichment = mean_hot / mean_nonhot if mean_nonhot > 0 else float("inf")

        try:
            stat, pval = mannwhitneyu(hot_vals, nonhot_vals, alternative="two-sided")
        except ValueError:
            stat, pval = 0, 1.0

        all_results.append({
            "feature": feat_name,
            "n_hotspot": len(hot_vals),
            "n_non_hotspot": len(nonhot_vals),
            "mean_hotspot": round(mean_hot, 6),
            "mean_non_hotspot": round(mean_nonhot, 6),
            "fold_enrichment": round(fold_enrichment, 4),
            "mann_whitney_U": stat,
            "p_value": pval,
        })

    # Apply BH FDR correction
    if all_results:
        p_values = [r["p_value"] for r in all_results]
        _, q_values, _, _ = multipletests(p_values, method="fdr_bh")
        for r, q in zip(all_results, q_values):
            r["q_value"] = round(q, 6)

    # Save enrichment results
    df_enrich = pd.DataFrame(all_results)
    if not df_enrich.empty:
        df_enrich = df_enrich.sort_values("p_value")
    enrich_out = os.path.join(output_dir, "feature_enrichment.csv")
    df_enrich.to_csv(enrich_out, index=False)
    logger.info(f"Saved enrichment results to {enrich_out}")

    # Save feature matrix
    df_matrix = pd.DataFrame(feature_matrix_data)
    matrix_out = os.path.join(output_dir, "feature_matrix.csv")
    df_matrix.to_csv(matrix_out, index=False)
    logger.info(f"Saved feature matrix to {matrix_out}")

    # Log summary
    for r in all_results:
        sig = "***" if r.get("q_value", 1) < 0.001 else (
            "**" if r.get("q_value", 1) < 0.01 else (
                "*" if r.get("q_value", 1) < 0.05 else "ns"
            )
        )
        logger.info(
            f"  {r['feature']}: fold={r['fold_enrichment']:.2f}, "
            f"p={r['p_value']:.2e}, q={r.get('q_value', 'N/A')} ({sig})"
        )
