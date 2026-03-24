"""Hotspot significance testing using genome-wide permutation.

Tests whether eccDNA are non-uniformly distributed by comparing observed
window counts against a null distribution from randomly shuffled eccDNA
positions within valid (non-excluded) genomic regions.

Algorithm:
  1. Load eccDNA midpoints and count per non-overlapping window.
  2. Build exclusion regions from gap/centromere BED files (optional).
  3. Build valid chromosome segments for uniform random sampling.
  4. Per-chromosome permutation: redistribute midpoints uniformly within
     valid regions, preserving per-chromosome totals.
  5. Repeat N times to build null distribution per window.
  6. Empirical p-value = (exceed_count + 1) / (N_perm + 1) (conservative).
  7. Benjamini-Hochberg FDR correction.
"""

import os
import logging
from typing import List, Optional, Dict, Tuple
from multiprocessing import Pool
from functools import partial

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_exclusion_regions(
    exclude_files: List[str],
    valid_chroms: set,
) -> Dict[str, List[Tuple[int, int]]]:
    """Load and merge exclusion regions from BED-like files.

    Supports plain BED (chrom, start, end) and UCSC table formats where the
    chrom is in column 0 or column 1. For cytoBand files, only 'acen' bands
    are included.

    Args:
        exclude_files: List of file paths (BED, gap.txt, cytoBand.txt).
        valid_chroms: Set of chromosome names to keep.

    Returns:
        Dict mapping chrom -> sorted, merged list of (start, end) intervals.
    """
    from collections import defaultdict
    raw = defaultdict(list)

    for fpath in exclude_files:
        if not os.path.exists(fpath):
            logger.warning(f"Exclusion file not found, skipping: {fpath}")
            continue

        basename = os.path.basename(fpath).lower()
        is_cytoband = "cytoband" in basename or "cyto" in basename

        with open(fpath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.split("\t")

                # cytoBand: only keep 'acen' entries
                if is_cytoband:
                    if len(parts) >= 5 and parts[4] == "acen":
                        chrom = parts[0]
                        if chrom in valid_chroms:
                            try:
                                raw[chrom].append((int(parts[1]), int(parts[2])))
                            except ValueError:
                                pass
                    continue

                # Generic BED or UCSC table
                if len(parts) < 3:
                    continue
                if parts[0].startswith("chr") and parts[0] in valid_chroms:
                    try:
                        raw[parts[0]].append((int(parts[1]), int(parts[2])))
                    except ValueError:
                        pass
                elif len(parts) >= 4 and parts[1].startswith("chr") and parts[1] in valid_chroms:
                    try:
                        raw[parts[1]].append((int(parts[2]), int(parts[3])))
                    except ValueError:
                        pass

    # Merge overlapping regions per chromosome
    merged = {}
    for chrom in valid_chroms:
        regions = sorted(raw.get(chrom, []))
        if not regions:
            merged[chrom] = []
            continue
        m = [regions[0]]
        for s, e in regions[1:]:
            if s <= m[-1][1]:
                m[-1] = (m[-1][0], max(m[-1][1], e))
            else:
                m.append((s, e))
        merged[chrom] = m
    return merged


def _build_valid_segments(
    chrom_sizes: Dict[str, int],
    exclude_regions: Dict[str, List[Tuple[int, int]]],
) -> Dict[str, Tuple[List[Tuple[int, int]], np.ndarray, float]]:
    """Build valid (non-excluded) segments per chromosome for random sampling.

    Returns:
        Dict: chrom -> (segments, cumulative_lengths, total_valid_length)
            segments: list of (start, end) valid intervals
            cumulative_lengths: np.array for inverse CDF sampling
            total_valid_length: total bp available for sampling
    """
    valid = {}
    for chrom, csize in chrom_sizes.items():
        excl = exclude_regions.get(chrom, [])
        segments = []
        pos = 0
        for es, ee in excl:
            if es > pos:
                segments.append((pos, es))
            pos = max(pos, ee)
        if csize > pos:
            segments.append((pos, csize))

        if not segments:
            segments = [(0, csize)]

        lengths = np.array([e - s for s, e in segments], dtype=np.float64)
        cum = np.cumsum(lengths)
        total = cum[-1] if len(cum) > 0 else 0
        valid[chrom] = (segments, cum, total)
    return valid


def _sample_uniform_positions(
    valid_info: Tuple[List[Tuple[int, int]], np.ndarray, float],
    n: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample n positions uniformly from valid segments using inverse CDF.

    Args:
        valid_info: (segments, cumulative_lengths, total_valid_length).
        n: Number of positions to sample.
        rng: numpy random Generator.

    Returns:
        Array of sampled genomic positions.
    """
    segments, cum, total = valid_info
    if total == 0 or n == 0:
        return np.array([], dtype=np.int64)

    u = rng.uniform(0, total, size=n)
    seg_idx = np.searchsorted(cum, u, side="right")
    seg_idx = np.clip(seg_idx, 0, len(segments) - 1)

    seg_starts = np.array([segments[i][0] for i in seg_idx])
    offsets = u - np.where(seg_idx > 0, cum[seg_idx - 1], 0)
    return (seg_starts + offsets).astype(np.int64)


def _count_windows(midpoints: np.ndarray, chrom_size: int, window_size: int) -> np.ndarray:
    """Count midpoints per non-overlapping window."""
    n_windows = (chrom_size + window_size - 1) // window_size
    if len(midpoints) == 0:
        return np.zeros(n_windows, dtype=np.int32)
    bins = midpoints // window_size
    counts = np.bincount(bins, minlength=n_windows)[:n_windows]
    return counts.astype(np.int32)


def _benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction. Returns q-values."""
    n = len(p_values)
    if n == 0:
        return np.array([])
    p_arr = np.asarray(p_values, dtype=np.float64)
    order = np.argsort(p_arr)
    q = np.empty(n, dtype=np.float64)
    cummin = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        rank = i + 1
        val = min(p_arr[idx] * n / rank, 1.0)
        cummin = min(cummin, val)
        q[idx] = cummin
    return q


def _permute_chromosome(args):
    """Worker function for parallelized per-chromosome permutation.

    Args:
        args: Tuple of (chrom, n_ecc, obs_counts, valid_info, chrom_size,
               window_size, n_perm, seed).

    Returns:
        Tuple: (chrom, exceed_counts) where exceed_counts[i] = number of
               permutations where random count >= observed count for window i.
    """
    (chrom, n_ecc, obs_counts, valid_info, chrom_size,
     window_size, n_perm, seed) = args

    rng = np.random.default_rng(seed)
    n_win = len(obs_counts)
    exceed = np.zeros(n_win, dtype=np.int32)

    for _ in range(n_perm):
        rand_pos = _sample_uniform_positions(valid_info, n_ecc, rng)
        rand_counts = _count_windows(rand_pos, chrom_size, window_size)
        exceed += (rand_counts >= obs_counts).astype(np.int32)

    return chrom, exceed


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_hotspot_permtest(
    input_file: str,
    fai_file: str,
    output_dir: str,
    window_size: int = 100_000,
    n_perm: int = 1000,
    n_cores: int = 4,
    fdr_threshold: float = 0.05,
    fold_threshold: float = 3.0,
    exclude_files: Optional[List[str]] = None,
    chrom_col: Optional[str] = None,
    start_col: Optional[str] = None,
    end_col: Optional[str] = None,
    chrom_filter: Optional[List[str]] = None,
    seed: int = 42,
) -> str:
    """Test hotspot significance using genome-wide permutation with FDR.

    For each non-overlapping window, tests whether the observed eccDNA count
    is significantly higher than expected under a uniform null distribution.
    Uses per-chromosome shuffling to preserve chromosome-level totals.

    Args:
        input_file: eccDNA CSV or BED file.
        fai_file: Reference genome .fai index for chromosome sizes.
        output_dir: Output directory.
        window_size: Window size in bp (default: 100000).
        n_perm: Number of permutations (default: 1000).
        n_cores: Number of CPU cores for parallelization (default: 4).
        fdr_threshold: FDR q-value threshold for significance (default: 0.05).
        fold_threshold: Minimum fold-above-median to be called hotspot (default: 3.0).
        exclude_files: List of BED/gap/cytoBand files defining exclusion regions.
        chrom_col: Column name for chromosome (auto-detected if None).
        start_col: Column name for start (auto-detected if None).
        end_col: Column name for end (auto-detected if None).
        chrom_filter: List of chromosomes to include (None = all from fai).
        seed: Random seed for reproducibility.

    Returns:
        Path to output CSV with columns: chrom, start, end, observed_count,
        empirical_pvalue, q_value, is_hotspot.
    """
    from ecctoolkit.hotspot.matrix import _load_chrom_sizes, _load_eccdna_midpoints

    os.makedirs(output_dir, exist_ok=True)

    # Load chromosome sizes
    chrom_sizes = _load_chrom_sizes(fai_file)
    if chrom_filter is not None:
        chrom_sizes = {c: s for c, s in chrom_sizes.items() if c in chrom_filter}
    chroms = sorted(chrom_sizes.keys(),
                    key=lambda c: (not c.startswith("chr"),
                                   c.replace("chr", "").zfill(5)))
    valid_chroms = set(chroms)

    # Build exclusion regions
    if exclude_files:
        exclude_regions = _load_exclusion_regions(exclude_files, valid_chroms)
    else:
        exclude_regions = {c: [] for c in chroms}

    valid_segments = _build_valid_segments(chrom_sizes, exclude_regions)

    for chrom in chroms[:3]:
        seg, cum, total = valid_segments[chrom]
        pct = total / chrom_sizes[chrom] * 100
        logger.info(
            f"  {chrom}: {len(seg)} valid segments, "
            f"{total / 1e6:.1f} Mb valid ({pct:.1f}%)"
        )

    # Load eccDNA midpoints
    df = _load_eccdna_midpoints(input_file, chrom_col, start_col, end_col,
                                valid_chroms)
    logger.info(f"Loaded {len(df)} eccDNA records")

    midpoints_by_chrom = {}
    for chrom in chroms:
        chrom_df = df[df["chrom"] == chrom]
        midpoints_by_chrom[chrom] = np.sort(chrom_df["midpoint"].values)

    ws_label = (f"{window_size // 1000}kb" if window_size < 1_000_000
                else f"{window_size // 1_000_000}Mb")
    logger.info(
        f"Permutation test: window={ws_label}, n_perm={n_perm}, "
        f"FDR<{fdr_threshold}, fold>={fold_threshold}"
    )

    # Step 1: Observed counts
    obs_by_chrom = {}
    for chrom in chroms:
        csize = chrom_sizes[chrom]
        obs_by_chrom[chrom] = _count_windows(
            midpoints_by_chrom[chrom], csize, window_size
        )

    # Step 2: Permutation test (parallelized across chromosomes)
    rng_master = np.random.default_rng(seed)
    worker_args = []
    for chrom in chroms:
        n_ecc = len(midpoints_by_chrom[chrom])
        chrom_seed = int(rng_master.integers(0, 2**31))
        worker_args.append((
            chrom, n_ecc, obs_by_chrom[chrom],
            valid_segments[chrom], chrom_sizes[chrom],
            window_size, n_perm, chrom_seed,
        ))

    if n_cores > 1 and len(chroms) > 1:
        logger.info(f"Running permutation with {n_cores} processes...")
        with Pool(processes=min(n_cores, len(chroms))) as pool:
            results = pool.map(_permute_chromosome, worker_args)
    else:
        logger.info("Running permutation (single process)...")
        results = [_permute_chromosome(a) for a in worker_args]

    exceed_by_chrom = dict(results)

    # Step 3: Empirical p-values (conservative: +1 / +1)
    all_pvals = []
    window_index = []  # (chrom, window_idx, obs_count)
    for chrom in chroms:
        obs = obs_by_chrom[chrom]
        exc = exceed_by_chrom[chrom]
        pv = (exc + 1).astype(np.float64) / (n_perm + 1)
        for wi in range(len(obs)):
            all_pvals.append(pv[wi])
            window_index.append((chrom, wi, int(obs[wi])))

    all_pvals = np.array(all_pvals)

    # Step 4: FDR correction (genome-wide)
    q_values = _benjamini_hochberg(all_pvals)

    # Step 5: Determine significance with fold-threshold pre-filter
    all_obs = np.concatenate([obs_by_chrom[c] for c in chroms])
    non_empty = all_obs[all_obs > 0]
    median_val = float(np.median(non_empty)) if len(non_empty) > 0 else 0
    count_threshold = fold_threshold * median_val

    logger.info(
        f"Median (non-empty windows) = {median_val:.1f}, "
        f"count threshold ({fold_threshold}x) = {count_threshold:.1f}"
    )

    # Build output records
    records = []
    n_sig = 0
    for gi, (chrom, wi, obs_count) in enumerate(window_index):
        wstart = wi * window_size
        wend = min(wstart + window_size, chrom_sizes[chrom])
        is_hot = (q_values[gi] < fdr_threshold and obs_count >= count_threshold)
        if is_hot:
            n_sig += 1
        records.append({
            "chrom": chrom,
            "start": wstart,
            "end": wend,
            "observed_count": obs_count,
            "empirical_pvalue": float(all_pvals[gi]),
            "q_value": float(q_values[gi]),
            "is_hotspot": is_hot,
        })

    logger.info(f"Significant hotspot windows (q<{fdr_threshold}): {n_sig}")

    # Save output
    df_out = pd.DataFrame(records)
    out_path = os.path.join(output_dir, f"permtest_{ws_label}.csv")
    df_out.to_csv(out_path, index=False)
    logger.info(f"Saved: {out_path}")

    return out_path
