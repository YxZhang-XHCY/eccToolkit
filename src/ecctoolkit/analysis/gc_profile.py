"""GC content analysis for eccDNA regions and breakpoint flanking sequences.

Provides three levels of GC analysis:
  1. Per-eccDNA GC content distribution
  2. Breakpoint-flanking GC profile (metagene-style, ±N bp around junctions)
  3. Comparison with genomic background (random region sampling)
"""

import logging
import os
import random
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_fai(genome_file: str) -> Dict[str, int]:
    """Load chromosome sizes from .fai index."""
    fai_path = genome_file + ".fai"
    sizes: Dict[str, int] = {}
    if not os.path.exists(fai_path):
        raise FileNotFoundError(
            f"FAI index not found: {fai_path}. "
            "Run 'samtools faidx <genome.fa>' first."
        )
    with open(fai_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def _parse_input(input_file: str) -> pd.DataFrame:
    """Parse eccDNA input (CSV/TSV/BED) into DataFrame with chr, start, end.

    Supports:
      - CSV/TSV with chr/chrom, start, end columns
      - BED format (tab-separated, no header)
      - seqname format (chr1:100-200)
    """
    # Try CSV/TSV first
    sep = "\t" if input_file.endswith((".bed", ".tsv")) else ","
    try:
        df = pd.read_csv(input_file, sep=sep, comment="#")
    except Exception:
        df = pd.read_csv(input_file, sep="\t", comment="#", header=None)

    # Normalize column names
    col_map = {}
    for c in df.columns:
        cl = str(c).lower().strip()
        if cl in ("chr", "chrom", "chromosome", "#chrom", "#chr"):
            col_map[c] = "chr"
        elif cl == "start":
            col_map[c] = "start"
        elif cl == "end":
            col_map[c] = "end"
        elif cl in ("name", "eccdna_name", "eccdna_id", "id"):
            col_map[c] = "name"
        elif cl in ("seqname",):
            col_map[c] = "seqname"
    df = df.rename(columns=col_map)

    # Handle BED with no header
    if "chr" not in df.columns and df.shape[1] >= 3:
        df.columns = ["chr", "start", "end"] + [
            f"col{i}" for i in range(3, df.shape[1])
        ]

    # Handle seqname format
    if "chr" not in df.columns and "seqname" in df.columns:
        import re
        pattern = re.compile(r"^(.+?):(\d+)-(\d+)")
        parsed = df["seqname"].str.extract(pattern)
        df["chr"] = parsed[0]
        df["start"] = parsed[1].astype(int)
        df["end"] = parsed[2].astype(int)

    if "chr" not in df.columns or "start" not in df.columns:
        raise ValueError(
            "Cannot find chr/start/end columns. "
            "Expected CSV/TSV with chr,start,end or BED format."
        )

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["length"] = df["end"] - df["start"]

    if "name" not in df.columns:
        df["name"] = [f"ecc_{i}" for i in range(len(df))]

    return df[["name", "chr", "start", "end", "length"]].copy()


def _fetch_sequence(genome_fh, chrom: str, start: int, end: int) -> str:
    """Fetch sequence from indexed FASTA using pysam."""
    try:
        return genome_fh.fetch(chrom, start, end).upper()
    except (KeyError, ValueError):
        return ""


def _gc_content(seq: str) -> float:
    """Calculate GC content of a sequence."""
    if not seq:
        return np.nan
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    valid = gc + seq.count("A") + seq.count("T")
    return gc / valid if valid > 0 else np.nan


def _gc_per_position(seq: str) -> np.ndarray:
    """Return per-base GC indicator (1=G/C, 0=A/T, nan=other)."""
    arr = np.full(len(seq), np.nan)
    for i, base in enumerate(seq.upper()):
        if base in ("G", "C"):
            arr[i] = 1.0
        elif base in ("A", "T"):
            arr[i] = 0.0
    return arr


# ---------------------------------------------------------------------------
# Core analysis functions
# ---------------------------------------------------------------------------

def _analyze_eccdna_gc(
    df: pd.DataFrame,
    genome_fh,
) -> pd.DataFrame:
    """Compute overall GC content for each eccDNA region.

    Returns DataFrame with: name, chr, start, end, length, gc_content.
    """
    gc_values = []
    for _, row in df.iterrows():
        seq = _fetch_sequence(genome_fh, row["chr"], row["start"], row["end"])
        gc_values.append(_gc_content(seq))

    result = df.copy()
    result["gc_content"] = gc_values
    return result


def _analyze_breakpoint_gc_profile(
    df: pd.DataFrame,
    genome_fh,
    chrom_sizes: Dict[str, int],
    flank_size: int = 150,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute metagene-style GC profile around breakpoints.

    For each eccDNA, extracts ±flank_size bp around start (5' junction)
    and end (3' junction), computes per-position GC fraction.

    Returns:
        start_profile: shape (2 * flank_size,) - mean GC at each position around start
        end_profile: shape (2 * flank_size,) - mean GC at each position around end
    """
    window = 2 * flank_size
    start_matrix = []
    end_matrix = []

    for _, row in df.iterrows():
        chrom = row["chr"]
        ecc_start = int(row["start"])
        ecc_end = int(row["end"])
        chrom_len = chrom_sizes.get(chrom, float("inf"))

        # Start junction: [start - flank_size, start + flank_size)
        s_left = max(0, ecc_start - flank_size)
        s_right = min(chrom_len, ecc_start + flank_size)
        seq_start = _fetch_sequence(genome_fh, chrom, s_left, s_right)

        if len(seq_start) == window:
            start_matrix.append(_gc_per_position(seq_start))

        # End junction: [end - flank_size, end + flank_size)
        e_left = max(0, ecc_end - flank_size)
        e_right = min(chrom_len, ecc_end + flank_size)
        seq_end = _fetch_sequence(genome_fh, chrom, e_left, e_right)

        if len(seq_end) == window:
            end_matrix.append(_gc_per_position(seq_end))

    if start_matrix:
        start_profile = np.nanmean(np.array(start_matrix), axis=0)
    else:
        start_profile = np.full(window, np.nan)

    if end_matrix:
        end_profile = np.nanmean(np.array(end_matrix), axis=0)
    else:
        end_profile = np.full(window, np.nan)

    return start_profile, end_profile


def _sample_background_gc(
    chrom_sizes: Dict[str, int],
    genome_fh,
    lengths: List[int],
    n_samples: int = 10000,
    seed: int = 42,
) -> List[float]:
    """Sample random genomic regions and compute their GC content.

    Samples regions matching the length distribution of input eccDNA
    for fair comparison.
    """
    rng = random.Random(seed)
    chroms = [c for c in chrom_sizes if not c.startswith(("chrM", "chrUn", "chr.*_"))]
    if not chroms:
        chroms = list(chrom_sizes.keys())

    # Weight chromosomes by size
    total_size = sum(chrom_sizes[c] for c in chroms)
    weights = [chrom_sizes[c] / total_size for c in chroms]

    gc_values = []
    attempts = 0
    max_attempts = n_samples * 5

    while len(gc_values) < n_samples and attempts < max_attempts:
        attempts += 1
        chrom = rng.choices(chroms, weights=weights, k=1)[0]
        length = rng.choice(lengths)
        chrom_len = chrom_sizes[chrom]

        if length >= chrom_len:
            continue

        start = rng.randint(0, chrom_len - length)
        end = start + length
        seq = _fetch_sequence(genome_fh, chrom, start, end)
        gc = _gc_content(seq)
        if not np.isnan(gc):
            gc_values.append(gc)

    return gc_values


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_gc_profile(
    input_file: str,
    genome_file: str,
    output_dir: str,
    flank_size: int = 150,
    n_background: int = 10000,
    seed: int = 42,
    smooth_window: int = 10,
) -> None:
    """Run GC content analysis for eccDNA regions.

    Performs three analyses:
      1. Per-eccDNA GC content distribution
      2. Breakpoint-flanking GC profile (metagene-style, ±flank_size bp)
      3. Comparison with random genomic background

    Args:
        input_file: eccDNA regions file (CSV/TSV/BED with chr, start, end).
        genome_file: Reference genome FASTA (must have .fai index).
        output_dir: Output directory for results.
        flank_size: Flanking size (bp) around breakpoints (default: 150).
        n_background: Number of random background regions to sample (default: 10000).
        seed: Random seed for background sampling (default: 42).
        smooth_window: Smoothing window size for profile (default: 10).
    """
    import pysam

    os.makedirs(output_dir, exist_ok=True)

    # Load data
    logger.info(f"Loading eccDNA regions from {input_file}")
    df = _parse_input(input_file)
    logger.info(f"Loaded {len(df)} eccDNA regions")

    logger.info(f"Loading reference genome: {genome_file}")
    chrom_sizes = _load_fai(genome_file)
    genome_fh = pysam.FastaFile(genome_file)

    # -----------------------------------------------------------------------
    # Analysis 1: Per-eccDNA GC content
    # -----------------------------------------------------------------------
    logger.info("Analysis 1: Computing per-eccDNA GC content")
    gc_df = _analyze_eccdna_gc(df, genome_fh)
    gc_out = os.path.join(output_dir, "eccdna_gc_content.csv")
    gc_df.to_csv(gc_out, index=False)
    logger.info(f"  Mean GC: {gc_df['gc_content'].mean():.4f}")
    logger.info(f"  Median GC: {gc_df['gc_content'].median():.4f}")
    logger.info(f"  Saved to {gc_out}")

    # -----------------------------------------------------------------------
    # Analysis 2: Breakpoint flanking GC profile
    # -----------------------------------------------------------------------
    logger.info(f"Analysis 2: Computing breakpoint GC profile (±{flank_size}bp)")
    start_profile, end_profile = _analyze_breakpoint_gc_profile(
        df, genome_fh, chrom_sizes, flank_size
    )

    # Build position axis: relative to breakpoint
    positions = np.arange(-flank_size, flank_size)

    # Smooth profiles
    def _smooth(arr, w):
        if w <= 1 or len(arr) < w:
            return arr
        kernel = np.ones(w) / w
        return np.convolve(arr, kernel, mode="same")

    start_smooth = _smooth(start_profile, smooth_window)
    end_smooth = _smooth(end_profile, smooth_window)

    profile_df = pd.DataFrame({
        "position": positions,
        "start_junction_gc": start_profile,
        "start_junction_gc_smooth": start_smooth,
        "end_junction_gc": end_profile,
        "end_junction_gc_smooth": end_smooth,
    })
    profile_out = os.path.join(output_dir, "breakpoint_gc_profile.csv")
    profile_df.to_csv(profile_out, index=False)
    logger.info(f"  Saved breakpoint GC profile to {profile_out}")

    # Combined profile (average of start and end junctions)
    combined_profile = np.nanmean(
        np.array([start_profile, end_profile]), axis=0
    )
    combined_smooth = _smooth(combined_profile, smooth_window)

    # Inside vs outside GC comparison
    inside_gc = np.nanmean(combined_profile[flank_size:])
    outside_gc = np.nanmean(combined_profile[:flank_size])
    logger.info(f"  Inside eccDNA mean GC: {inside_gc:.4f}")
    logger.info(f"  Outside eccDNA mean GC: {outside_gc:.4f}")

    # -----------------------------------------------------------------------
    # Analysis 3: Background comparison
    # -----------------------------------------------------------------------
    logger.info(f"Analysis 3: Sampling {n_background} background regions")
    lengths = df["length"].tolist()
    bg_gc = _sample_background_gc(
        chrom_sizes, genome_fh, lengths, n_background, seed
    )
    logger.info(f"  Sampled {len(bg_gc)} background regions")
    logger.info(f"  Background mean GC: {np.mean(bg_gc):.4f}")

    bg_df = pd.DataFrame({"gc_content": bg_gc})
    bg_out = os.path.join(output_dir, "background_gc_content.csv")
    bg_df.to_csv(bg_out, index=False)
    logger.info(f"  Saved background GC to {bg_out}")

    # -----------------------------------------------------------------------
    # Summary statistics
    # -----------------------------------------------------------------------
    summary = {
        "n_eccdna": len(df),
        "eccdna_gc_mean": gc_df["gc_content"].mean(),
        "eccdna_gc_median": gc_df["gc_content"].median(),
        "eccdna_gc_std": gc_df["gc_content"].std(),
        "background_gc_mean": np.mean(bg_gc),
        "background_gc_median": np.median(bg_gc),
        "background_gc_std": np.std(bg_gc),
        "breakpoint_inside_gc": inside_gc,
        "breakpoint_outside_gc": outside_gc,
        "flank_size": flank_size,
        "n_background": len(bg_gc),
    }

    # Mann-Whitney U test: eccDNA GC vs background GC
    from scipy.stats import mannwhitneyu
    ecc_gc_valid = gc_df["gc_content"].dropna().values
    if len(ecc_gc_valid) > 0 and len(bg_gc) > 0:
        stat, pval = mannwhitneyu(ecc_gc_valid, bg_gc, alternative="two-sided")
        summary["mannwhitney_statistic"] = stat
        summary["mannwhitney_pvalue"] = pval
        logger.info(f"  eccDNA vs background GC: U={stat:.0f}, p={pval:.2e}")

    summary_df = pd.DataFrame([summary])
    summary_out = os.path.join(output_dir, "gc_analysis_summary.csv")
    summary_df.to_csv(summary_out, index=False)
    logger.info(f"  Saved summary to {summary_out}")

    genome_fh.close()
    logger.info("GC profile analysis completed!")
