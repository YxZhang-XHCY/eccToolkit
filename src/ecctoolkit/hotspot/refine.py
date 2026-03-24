"""Refine hotspot boundaries at high resolution.

Takes coarse hotspot regions (e.g., 100kb windows) and re-counts eccDNA
at finer resolution (e.g., 1kb sub-windows) within each hotspot to trim
boundaries where signal density drops below a threshold.
"""

import os
import logging
from typing import Optional, List

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def refine_hotspot_boundaries(
    input_file: str,
    eccdna_file: str,
    output_file: str,
    fine_window: int = 1000,
    density_fraction: float = 0.1,
    chrom_col: Optional[str] = None,
    start_col: Optional[str] = None,
    end_col: Optional[str] = None,
) -> str:
    """Refine hotspot boundaries at high resolution.

    For each coarse hotspot region, divides it into fine sub-windows and
    counts eccDNA midpoints. Then trims from both ends any sub-windows
    whose count is below a fraction of the peak sub-window count within
    that hotspot.

    Args:
        input_file: Hotspot CSV with columns: chrom, start, end
            (and optionally other columns that will be preserved).
        eccdna_file: eccDNA CSV or BED file with chrom/start/end.
        output_file: Output CSV path for refined boundaries.
        fine_window: Sub-window size in bp for refinement (default: 1000).
        density_fraction: Fraction of peak sub-window count below which
            boundary sub-windows are trimmed (default: 0.1, i.e., 10%).
        chrom_col: Column name for chromosome in eccdna_file (auto-detected).
        start_col: Column name for start in eccdna_file (auto-detected).
        end_col: Column name for end in eccdna_file (auto-detected).

    Returns:
        Path to output CSV with refined boundaries.
    """
    from ecctoolkit.hotspot.matrix import _load_eccdna_midpoints

    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)

    # Load hotspot regions
    df_hot = pd.read_csv(input_file)

    # Detect hotspot chrom/start/end columns
    hot_chrom_col = None
    hot_start_col = None
    hot_end_col = None
    for col in df_hot.columns:
        cl = col.lower()
        if cl in ("chr", "chrom", "chromosome"):
            hot_chrom_col = col
        elif cl in ("start", "chromstart", "window_start"):
            hot_start_col = col
        elif cl in ("end", "chromend", "window_end"):
            hot_end_col = col

    if not all([hot_chrom_col, hot_start_col, hot_end_col]):
        raise ValueError(
            f"Cannot detect chrom/start/end columns in hotspot file. "
            f"Available: {list(df_hot.columns)}"
        )

    # Load eccDNA midpoints
    df_ecc = _load_eccdna_midpoints(eccdna_file, chrom_col, start_col, end_col)
    logger.info(
        f"Loaded {len(df_ecc)} eccDNA records and {len(df_hot)} hotspot regions"
    )

    # Index eccDNA midpoints by chromosome for fast lookup
    midpoints_by_chrom = {}
    for chrom, group in df_ecc.groupby("chrom"):
        midpoints_by_chrom[chrom] = np.sort(group["midpoint"].values)

    # Refine each hotspot
    refined_records = []
    n_trimmed = 0

    for idx, row in df_hot.iterrows():
        chrom = row[hot_chrom_col]
        h_start = int(row[hot_start_col])
        h_end = int(row[hot_end_col])

        mids = midpoints_by_chrom.get(chrom, np.array([], dtype=np.int64))

        # Select midpoints within this hotspot region
        mask = (mids >= h_start) & (mids < h_end)
        region_mids = mids[mask]

        if len(region_mids) == 0:
            # No eccDNA in this hotspot, keep original boundaries
            record = row.to_dict()
            record["refined_start"] = h_start
            record["refined_end"] = h_end
            record["refined"] = False
            refined_records.append(record)
            continue

        # Divide hotspot into fine sub-windows
        n_sub = max(1, (h_end - h_start + fine_window - 1) // fine_window)
        sub_counts = np.zeros(n_sub, dtype=np.int32)

        # Count midpoints in each sub-window
        sub_bins = (region_mids - h_start) // fine_window
        sub_bins = np.clip(sub_bins, 0, n_sub - 1)
        for b in sub_bins:
            sub_counts[b] += 1

        peak_count = sub_counts.max()
        threshold = density_fraction * peak_count

        # Trim from left: find first sub-window above threshold
        left_idx = 0
        while left_idx < n_sub and sub_counts[left_idx] < threshold:
            left_idx += 1

        # Trim from right: find last sub-window above threshold
        right_idx = n_sub - 1
        while right_idx >= left_idx and sub_counts[right_idx] < threshold:
            right_idx -= 1

        refined_start = h_start + left_idx * fine_window
        refined_end = min(h_start + (right_idx + 1) * fine_window, h_end)

        was_trimmed = (refined_start != h_start or refined_end != h_end)
        if was_trimmed:
            n_trimmed += 1

        record = row.to_dict()
        record["refined_start"] = refined_start
        record["refined_end"] = refined_end
        record["refined"] = was_trimmed

        refined_records.append(record)

    df_refined = pd.DataFrame(refined_records)
    df_refined.to_csv(output_file, index=False)

    logger.info(
        f"Refined {n_trimmed}/{len(df_hot)} hotspot boundaries. "
        f"Output: {output_file}"
    )

    return output_file
