"""Combine FLED and CircleMap results."""

import logging
import os
from pathlib import Path
from typing import List

import pandas as pd

logger = logging.getLogger(__name__)


def _find_fled_files(fled_dir: str) -> List[Path]:
    """Find FLED eccDNA output files in directory."""
    p = Path(fled_dir)
    patterns = ["*.csv", "*.tsv", "*.txt", "*.bed"]
    files = []
    for pat in patterns:
        files.extend(p.glob(pat))
    if not files:
        raise FileNotFoundError(f"No FLED result files found in {fled_dir}")
    return sorted(files)


def _find_circlemap_files(circlemap_dir: str) -> List[Path]:
    """Find CircleMap output files (BED format) in directory."""
    p = Path(circlemap_dir)
    patterns = ["*.bed", "*.tsv", "*.txt"]
    files = []
    for pat in patterns:
        files.extend(p.glob(pat))
    if not files:
        raise FileNotFoundError(f"No CircleMap result files found in {circlemap_dir}")
    return sorted(files)


def _load_fled(filepath: Path) -> pd.DataFrame:
    """Load and standardize a FLED result file."""
    sep = "\t" if filepath.suffix in (".tsv", ".bed", ".txt") else ","
    df = pd.read_csv(filepath, sep=sep, comment="#")

    col_map = {}
    for c in df.columns:
        cl = str(c).lower().strip()
        if cl in ("chr", "chrom", "chromosome", "#chrom", "#chr"):
            col_map[c] = "chr"
        elif cl in ("start", "estart"):
            col_map[c] = "start"
        elif cl in ("end", "eend"):
            col_map[c] = "end"
        elif cl in ("length", "elength"):
            col_map[c] = "length"
    df = df.rename(columns=col_map)

    # BED with no header
    if "chr" not in df.columns and df.shape[1] >= 3:
        df.columns = ["chr", "start", "end"] + [
            f"col{i}" for i in range(3, df.shape[1])
        ]

    if "chr" not in df.columns or "start" not in df.columns or "end" not in df.columns:
        raise ValueError(f"Cannot parse FLED file {filepath}: missing chr/start/end columns")

    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["chr", "start", "end"])

    if "length" not in df.columns:
        df["length"] = df["end"] - df["start"]

    df["tool"] = "FLED"
    df["sample"] = filepath.stem
    return df[["chr", "start", "end", "length", "tool", "sample"]].copy()


def _load_circlemap(filepath: Path) -> pd.DataFrame:
    """Load and standardize a CircleMap result file."""
    df = pd.read_csv(filepath, sep="\t", comment="#", header=None)

    if df.shape[1] >= 3:
        df.columns = ["chr", "start", "end"] + [
            f"col{i}" for i in range(3, df.shape[1])
        ]
    else:
        raise ValueError(f"CircleMap file {filepath} has fewer than 3 columns")

    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["chr", "start", "end"])

    df["length"] = df["end"] - df["start"]
    df["tool"] = "CircleMap"
    df["sample"] = filepath.stem
    return df[["chr", "start", "end", "length", "tool", "sample"]].copy()


def _reciprocal_overlap(row1: pd.Series, row2: pd.Series, threshold: float = 0.8) -> bool:
    """Check if two intervals have reciprocal overlap >= threshold."""
    if row1["chr"] != row2["chr"]:
        return False

    s1, e1 = int(row1["start"]), int(row1["end"])
    s2, e2 = int(row2["start"]), int(row2["end"])

    overlap_start = max(s1, s2)
    overlap_end = min(e1, e2)
    overlap_len = max(0, overlap_end - overlap_start)

    len1 = e1 - s1
    len2 = e2 - s2

    if len1 == 0 or len2 == 0:
        return False

    frac1 = overlap_len / len1
    frac2 = overlap_len / len2

    return frac1 >= threshold and frac2 >= threshold


def _deduplicate_by_overlap(df: pd.DataFrame, threshold: float = 0.8) -> pd.DataFrame:
    """Deduplicate entries with reciprocal overlap >= threshold.

    For overlapping entries from different tools, keeps both but marks them as merged.
    For overlapping entries from the same tool, keeps the first.
    """
    if len(df) == 0:
        return df

    df = df.sort_values(["chr", "start", "end"]).reset_index(drop=True)
    keep = [True] * len(df)
    merged_tools = list(df["tool"])
    merged_samples = list(df["sample"])

    for i in range(len(df)):
        if not keep[i]:
            continue
        for j in range(i + 1, len(df)):
            if not keep[j]:
                continue
            # Early exit: if different chromosome or too far apart
            if df.at[i, "chr"] != df.at[j, "chr"]:
                break
            if int(df.at[j, "start"]) > int(df.at[i, "end"]):
                break

            if _reciprocal_overlap(df.iloc[i], df.iloc[j], threshold):
                # Merge: keep the first, annotate with both tools
                if merged_tools[i] != df.at[j, "tool"]:
                    merged_tools[i] = "FLED+CircleMap"
                if merged_samples[i] != df.at[j, "sample"]:
                    merged_samples[i] = f"{merged_samples[i]};{df.at[j, 'sample']}"
                keep[j] = False

    result = df[keep].copy()
    result["tool"] = [merged_tools[i] for i in range(len(df)) if keep[i]]
    result["sample"] = [merged_samples[i] for i in range(len(df)) if keep[i]]
    return result.reset_index(drop=True)


def combine_fled_circlemap(
    fled_dir: str,
    circlemap_dir: str,
    output_dir: str,
) -> None:
    """
    Combine FLED and CircleMap results into unified format.

    Loads results from both tools, standardizes columns, merges,
    and deduplicates by reciprocal overlap >= 0.8.

    Args:
        fled_dir: Directory containing FLED result files
        circlemap_dir: Directory containing CircleMap result files (BED)
        output_dir: Output directory
    """
    os.makedirs(output_dir, exist_ok=True)

    # Load FLED results
    logger.info(f"Searching for FLED results in {fled_dir}")
    fled_files = _find_fled_files(fled_dir)
    logger.info(f"Found {len(fled_files)} FLED files")

    fled_dfs = []
    for f in fled_files:
        try:
            fled_dfs.append(_load_fled(f))
            logger.info(f"  Loaded {f.name}: {len(fled_dfs[-1])} records")
        except Exception as e:
            logger.warning(f"  Skipping {f.name}: {e}")

    # Load CircleMap results
    logger.info(f"Searching for CircleMap results in {circlemap_dir}")
    cm_files = _find_circlemap_files(circlemap_dir)
    logger.info(f"Found {len(cm_files)} CircleMap files")

    cm_dfs = []
    for f in cm_files:
        try:
            cm_dfs.append(_load_circlemap(f))
            logger.info(f"  Loaded {f.name}: {len(cm_dfs[-1])} records")
        except Exception as e:
            logger.warning(f"  Skipping {f.name}: {e}")

    # Combine
    all_dfs = fled_dfs + cm_dfs
    if not all_dfs:
        raise ValueError("No valid result files loaded from either tool")

    combined = pd.concat(all_dfs, ignore_index=True)
    logger.info(f"Combined total: {len(combined)} records before deduplication")

    # Deduplicate
    deduped = _deduplicate_by_overlap(combined, threshold=0.8)
    logger.info(f"After deduplication: {len(deduped)} records")

    # Tool source summary
    tool_counts = deduped["tool"].value_counts()
    for tool, count in tool_counts.items():
        logger.info(f"  {tool}: {count} records")

    # Save
    out_file = os.path.join(output_dir, "combined_eccdna.csv")
    deduped.to_csv(out_file, index=False)
    logger.info(f"Saved combined results to {out_file}")
