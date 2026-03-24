"""Pairwise reciprocal overlap analysis for eccDNA intervals."""

from __future__ import annotations

import bisect
import csv
import logging
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _detect_format(filepath: str) -> str:
    """Detect input file format based on extension and content."""
    ext = Path(filepath).suffix.lower()
    if ext == ".bed":
        return "bed"
    if ext in (".csv", ".tsv"):
        return "csv"
    # Peek at first line to guess
    with open(filepath) as f:
        first_line = f.readline().strip()
    if "\t" in first_line and not first_line.startswith("#"):
        fields = first_line.split("\t")
        if len(fields) >= 3 and fields[0].startswith("chr"):
            return "bed"
    return "csv"


def _load_intervals(
    filepath: str,
    fmt: str = "auto",
) -> List[Tuple[str, int, int]]:
    """Load genomic intervals from CSV or BED file.

    CSV files are expected to have columns: chr/eChr, start/eStart, end/eEnd.
    BED files use the first three tab-separated columns.

    Returns:
        List of (chrom, start, end) tuples.
    """
    if fmt == "auto":
        fmt = _detect_format(filepath)

    intervals = []

    if fmt == "bed":
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                fields = line.split("\t")
                if len(fields) < 3:
                    continue
                try:
                    intervals.append((fields[0], int(fields[1]), int(fields[2])))
                except ValueError:
                    continue
    else:
        df = pd.read_csv(filepath)
        # Try standard eccDNA column names, then generic
        chr_col = next(
            (c for c in ["eChr", "chr", "chrom", "Chr"] if c in df.columns), None
        )
        start_col = next(
            (c for c in ["eStart", "start", "chromStart", "Start"] if c in df.columns),
            None,
        )
        end_col = next(
            (c for c in ["eEnd", "end", "chromEnd", "End"] if c in df.columns), None
        )
        if chr_col is None or start_col is None or end_col is None:
            raise ValueError(
                f"Cannot find chr/start/end columns in {filepath}. "
                f"Available columns: {list(df.columns)}"
            )
        for _, row in df.iterrows():
            try:
                intervals.append((str(row[chr_col]), int(row[start_col]), int(row[end_col])))
            except (ValueError, TypeError):
                continue

    logger.info(f"Loaded {len(intervals)} intervals from {Path(filepath).name}")
    return intervals


def _group_and_sort(
    intervals: List[Tuple[str, int, int]],
) -> Tuple[Dict[str, List[Tuple[int, int]]], Dict[str, List[int]]]:
    """Group intervals by chromosome and sort by start position.

    Returns:
        intervals_by_chr: chrom -> [(start, end), ...] sorted by start
        starts_by_chr: chrom -> [start, ...] aligned with intervals_by_chr
    """
    by_chr: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for chrom, start, end in intervals:
        by_chr[chrom].append((start, end))

    starts_by_chr: Dict[str, List[int]] = {}
    for chrom, ivs in by_chr.items():
        ivs.sort()
        starts_by_chr[chrom] = [s for s, _ in ivs]

    return dict(by_chr), starts_by_chr


def _recip_overlap(
    a_start: int, a_end: int, b_start: int, b_end: int, threshold: float
) -> bool:
    """Check if two intervals meet reciprocal overlap threshold."""
    ov = min(a_end, b_end) - max(a_start, b_start)
    if ov <= 0:
        return False
    len_a = a_end - a_start
    len_b = b_end - b_start
    if len_a <= 0 or len_b <= 0:
        return False
    return ov >= threshold * len_a and ov >= threshold * len_b


def _count_overlap(
    a_by_chr: Dict[str, List[Tuple[int, int]]],
    b_by_chr: Dict[str, List[Tuple[int, int]]],
    b_starts_by_chr: Dict[str, List[int]],
    recip: float,
) -> int:
    """Count how many intervals in A have reciprocal overlap with any in B.

    Uses binary search and length constraints for efficient pruning.
    """
    inv = 1.0 / recip
    count = 0

    for chrom, a_intervals in a_by_chr.items():
        b_intervals = b_by_chr.get(chrom)
        if not b_intervals:
            continue
        b_starts = b_starts_by_chr[chrom]

        for a_start, a_end in a_intervals:
            len_a = a_end - a_start
            if len_a <= 0:
                continue

            min_len_b = len_a * recip
            max_len_b = len_a * inv

            # Prune search range using start coordinate constraints
            left = bisect.bisect_left(b_starts, a_start - int(math.ceil(max_len_b)))
            right = bisect.bisect_left(b_starts, a_end)

            found = False
            for k in range(left, right):
                b_start, b_end = b_intervals[k]
                if b_end <= a_start:
                    continue
                len_b = b_end - b_start
                if len_b <= 0 or len_b < min_len_b or len_b > max_len_b:
                    continue
                if _recip_overlap(a_start, a_end, b_start, b_end, recip):
                    found = True
                    break
            if found:
                count += 1

    return count


def compute_reciprocal_overlap(
    input_files: List[str],
    output_file: str,
    sample_names: Optional[List[str]] = None,
    min_reciprocal: float = 0.5,
    input_format: str = "auto",
    output_matrix: bool = False,
) -> None:
    """Compute pairwise reciprocal overlap between eccDNA interval sets.

    For each pair of samples (A, B), counts how many intervals in A have
    at least one interval in B with reciprocal overlap >= threshold.

    Outputs a long-format CSV with one row per pair. When output_matrix is
    True, also writes an NxN overlap percentage matrix CSV alongside the
    main output (with "_matrix" suffix).

    Args:
        input_files: List of input file paths (CSV or BED).
        output_file: Output CSV path for the pairwise overlap table.
        sample_names: Optional sample labels (default: derived from filenames).
        min_reciprocal: Reciprocal overlap threshold (default: 0.5 = 50%).
        input_format: Input format: "csv", "bed", or "auto" (detect).
        output_matrix: Also write NxN percentage matrix CSV (default: False).
    """
    if sample_names is None:
        sample_names = [Path(f).stem for f in input_files]

    if len(sample_names) != len(input_files):
        raise ValueError("Number of sample names must match number of input files")

    os.makedirs(Path(output_file).parent, exist_ok=True)

    # Load all samples
    data: Dict[str, Dict[str, List[Tuple[int, int]]]] = {}
    starts: Dict[str, Dict[str, List[int]]] = {}
    sizes: Dict[str, int] = {}

    for name, filepath in zip(sample_names, input_files):
        intervals = _load_intervals(filepath, fmt=input_format)
        by_chr, starts_by_chr = _group_and_sort(intervals)
        data[name] = by_chr
        starts[name] = starts_by_chr
        sizes[name] = len(intervals)
        logger.info(f"  {name}: {len(intervals):,} intervals")

    n = len(sample_names)

    # Compute pairwise overlap
    overlap_pct = np.zeros((n, n), dtype=np.float64)
    rows = []
    for i in range(n):
        overlap_pct[i, i] = 100.0
        for j in range(n):
            if i == j:
                continue
            a, b = sample_names[i], sample_names[j]
            cnt = _count_overlap(data[a], data[b], starts[b], recip=min_reciprocal)
            frac = cnt / sizes[a] if sizes[a] > 0 else 0.0
            overlap_pct[i, j] = frac * 100.0
            rows.append({
                "sample_A": a,
                "sample_B": b,
                "n_A": sizes[a],
                "n_B": sizes[b],
                "overlap_count": cnt,
                "overlap_fraction": round(frac, 6),
            })
            logger.info(
                f"  {a} vs {b}: {cnt:,}/{sizes[a]:,} = {frac:.4f}"
            )

    # Save long-format output
    df = pd.DataFrame(rows)
    df.to_csv(output_file, index=False)
    logger.info(f"Saved pairwise overlap to {output_file}")

    # Optionally save NxN percentage matrix
    if output_matrix:
        df_matrix = pd.DataFrame(
            overlap_pct, index=sample_names, columns=sample_names
        )
        base, ext = os.path.splitext(output_file)
        matrix_file = f"{base}_matrix{ext}"
        df_matrix.to_csv(matrix_file, index_label="sample")
        logger.info(f"Saved overlap matrix to {matrix_file}")
