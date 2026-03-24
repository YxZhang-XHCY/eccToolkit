"""Generate multi-scale window count matrices for hotspot analysis.

Counts eccDNA midpoints in non-overlapping genomic windows at multiple
resolutions (default: 10kb, 50kb, 100kb). Uses vectorized numpy bincount
for fast counting without external dependencies like bedtools.
"""

import os
import logging
from typing import List, Optional, Dict, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _load_chrom_sizes(fai_file: str) -> Dict[str, int]:
    """Load chromosome sizes from a FASTA index (.fai) file.

    Args:
        fai_file: Path to .fai file (tab-separated: name, length, ...).

    Returns:
        Dict mapping chromosome name to size.
    """
    chrom_sizes = {}
    with open(fai_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes


def _load_eccdna_midpoints(
    input_file: str,
    chrom_col: Optional[str] = None,
    start_col: Optional[str] = None,
    end_col: Optional[str] = None,
    valid_chroms: Optional[set] = None,
) -> pd.DataFrame:
    """Load eccDNA records and compute midpoints.

    Supports CSV (with header) and BED (no header, tab-separated) formats.
    Auto-detects column names for common eccDNA CSV formats.

    Args:
        input_file: Path to eccDNA CSV or BED file.
        chrom_col: Column name for chromosome (auto-detected if None).
        start_col: Column name for start position (auto-detected if None).
        end_col: Column name for end position (auto-detected if None).
        valid_chroms: Set of chromosome names to keep (None = keep all).

    Returns:
        DataFrame with columns: chrom, start, end, midpoint.
    """
    ext = os.path.splitext(input_file)[1].lower()

    if ext == ".bed":
        df = pd.read_csv(input_file, sep="\t", header=None,
                         usecols=[0, 1, 2], names=["chrom", "start", "end"],
                         dtype={"chrom": str})
    else:
        df = pd.read_csv(input_file, low_memory=False)
        # Auto-detect column names
        col_map = {}
        for col in df.columns:
            cl = col.lower()
            if cl in ("chr", "chrom", "chromosome"):
                col_map["chrom"] = col
            elif cl in ("start", "chromstart"):
                col_map["start"] = col
            elif cl in ("end", "chromend"):
                col_map["end"] = col

        # Allow caller overrides
        if chrom_col:
            col_map["chrom"] = chrom_col
        if start_col:
            col_map["start"] = start_col
        if end_col:
            col_map["end"] = end_col

        if not all(k in col_map for k in ("chrom", "start", "end")):
            raise ValueError(
                f"Cannot detect chrom/start/end columns in {input_file}. "
                f"Available columns: {list(df.columns)}. "
                f"Specify chrom_col, start_col, end_col explicitly."
            )

        df = df.rename(columns={
            col_map["chrom"]: "chrom",
            col_map["start"]: "start",
            col_map["end"]: "end",
        })

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)

    if valid_chroms is not None:
        df = df[df["chrom"].isin(valid_chroms)]

    df["midpoint"] = ((df["start"] + df["end"]) // 2).astype(np.int64)
    return df[["chrom", "start", "end", "midpoint"]].reset_index(drop=True)


def _count_midpoints_in_windows(
    midpoints: np.ndarray,
    chrom_size: int,
    window_size: int,
) -> np.ndarray:
    """Count midpoints per non-overlapping window using integer division.

    Args:
        midpoints: Sorted array of midpoint positions.
        chrom_size: Chromosome size in bp.
        window_size: Window size in bp.

    Returns:
        Array of counts, one per window.
    """
    n_windows = (chrom_size + window_size - 1) // window_size
    if len(midpoints) == 0:
        return np.zeros(n_windows, dtype=np.int32)
    bins = midpoints // window_size
    counts = np.bincount(bins, minlength=n_windows)[:n_windows]
    return counts.astype(np.int32)


def generate_multiscale_matrix(
    input_file: str,
    fai_file: str,
    output_dir: str,
    window_sizes: Optional[List[int]] = None,
    chrom_col: Optional[str] = None,
    start_col: Optional[str] = None,
    end_col: Optional[str] = None,
    chrom_filter: Optional[List[str]] = None,
) -> Dict[int, str]:
    """Generate multi-scale window count matrices.

    Creates eccDNA density matrices at different resolutions for downstream
    hotspot detection. For each window size, counts eccDNA midpoints in
    non-overlapping genomic windows using vectorized numpy operations.

    Args:
        input_file: eccDNA CSV or BED file (must have chrom/start/end columns).
        fai_file: Reference genome .fai index file for chromosome sizes.
        output_dir: Output directory for CSV files.
        window_sizes: List of window sizes in bp (default: [10000, 50000, 100000]).
        chrom_col: Column name for chromosome (auto-detected if None).
        start_col: Column name for start (auto-detected if None).
        end_col: Column name for end (auto-detected if None).
        chrom_filter: List of chromosomes to include (None = use all from fai).

    Returns:
        Dict mapping window_size -> output CSV file path.
    """
    if window_sizes is None:
        window_sizes = [10_000, 50_000, 100_000]

    os.makedirs(output_dir, exist_ok=True)

    # Load chromosome sizes
    chrom_sizes = _load_chrom_sizes(fai_file)
    if chrom_filter is not None:
        chrom_sizes = {c: s for c, s in chrom_sizes.items() if c in chrom_filter}
    chroms = sorted(chrom_sizes.keys(),
                    key=lambda c: (not c.startswith("chr"),
                                   c.replace("chr", "").zfill(5)))
    valid_chroms = set(chroms)

    logger.info(f"Loaded {len(chroms)} chromosomes from {fai_file}")

    # Load eccDNA midpoints
    df = _load_eccdna_midpoints(input_file, chrom_col, start_col, end_col,
                                valid_chroms)
    logger.info(f"Loaded {len(df)} eccDNA records from {input_file}")

    # Group midpoints by chromosome (sorted arrays)
    midpoints_by_chrom: Dict[str, np.ndarray] = {}
    for chrom in chroms:
        chrom_df = df[df["chrom"] == chrom]
        midpoints_by_chrom[chrom] = np.sort(chrom_df["midpoint"].values)

    # Generate per-window-size matrices
    output_files = {}
    all_records = []

    for ws in window_sizes:
        ws_label = f"{ws // 1000}kb" if ws < 1_000_000 else f"{ws // 1_000_000}Mb"
        logger.info(f"Counting at {ws_label} resolution...")

        records = []
        for chrom in chroms:
            csize = chrom_sizes[chrom]
            mids = midpoints_by_chrom.get(chrom, np.array([], dtype=np.int64))
            counts = _count_midpoints_in_windows(mids, csize, ws)

            for wi, cnt in enumerate(counts):
                wstart = wi * ws
                wend = min(wstart + ws, csize)
                records.append({
                    "chrom": chrom,
                    "start": wstart,
                    "end": wend,
                    "count": int(cnt),
                })

        df_out = pd.DataFrame(records)
        out_path = os.path.join(output_dir, f"window_counts_{ws_label}.csv")
        df_out.to_csv(out_path, index=False)
        output_files[ws] = out_path

        n_nonzero = (df_out["count"] > 0).sum()
        total_count = df_out["count"].sum()
        logger.info(
            f"  {ws_label}: {len(df_out)} windows, "
            f"{n_nonzero} non-empty, total count = {total_count}"
        )

        # Add window_size column for combined output
        df_out_combined = df_out.copy()
        df_out_combined.insert(0, "window_size", ws_label)
        all_records.append(df_out_combined)

    # Save combined matrix
    df_combined = pd.concat(all_records, ignore_index=True)
    combined_path = os.path.join(output_dir, "window_counts_combined.csv")
    df_combined.to_csv(combined_path, index=False)
    output_files["combined"] = combined_path
    logger.info(f"Combined matrix: {combined_path} ({len(df_combined)} rows)")

    return output_files
