"""File I/O utilities for eccToolkit."""

import logging
from pathlib import Path
from typing import List, Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)

# Standard eccDNA column names
ECCDNA_REQUIRED_COLS = ["eChr", "eStart", "eEnd"]
ECCDNA_OPTIONAL_COLS = ["eName", "eLength"]
BED_COLS = ["chrom", "start", "end", "name", "score", "strand"]


def load_eccdna_csv(
    filepath: Union[str, Path],
    required_cols: Optional[List[str]] = None,
    calculate_length: bool = True,
) -> pd.DataFrame:
    """
    Load eccDNA data from CSV file.

    Args:
        filepath: Path to CSV file
        required_cols: List of required columns (default: eChr, eStart, eEnd)
        calculate_length: Calculate eLength if missing

    Returns:
        DataFrame with eccDNA data

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If required columns are missing
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    if required_cols is None:
        required_cols = ECCDNA_REQUIRED_COLS

    df = pd.read_csv(filepath)

    # Check required columns
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Calculate eLength if missing
    if calculate_length and "eLength" not in df.columns:
        if "eStart" in df.columns and "eEnd" in df.columns:
            df["eLength"] = df["eEnd"] - df["eStart"]

    logger.info(f"Loaded {len(df)} records from {filepath.name}")
    return df


def load_bed(
    filepath: Union[str, Path],
    names: Optional[List[str]] = None,
    sep: str = "\t",
    header: Optional[int] = None,
) -> pd.DataFrame:
    """
    Load BED format file.

    Args:
        filepath: Path to BED file
        names: Column names (default: chrom, start, end, ...)
        sep: Field separator
        header: Header row number (default: None)

    Returns:
        DataFrame with BED data
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    # Detect number of columns
    with open(filepath, "r") as f:
        first_line = f.readline().strip()
        n_cols = len(first_line.split(sep))

    # Use default names based on column count
    if names is None:
        names = BED_COLS[:n_cols]

    df = pd.read_csv(filepath, sep=sep, header=header, names=names)
    logger.info(f"Loaded {len(df)} records from {filepath.name}")
    return df


def save_bed(
    df: pd.DataFrame,
    filepath: Union[str, Path],
    cols: Optional[List[str]] = None,
    header: bool = False,
) -> None:
    """
    Save DataFrame to BED format.

    Args:
        df: DataFrame to save
        filepath: Output file path
        cols: Columns to write (default: first 3 columns)
        header: Write header row
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    if cols is None:
        # Try standard eccDNA columns first
        if all(c in df.columns for c in ["eChr", "eStart", "eEnd"]):
            cols = ["eChr", "eStart", "eEnd"]
            if "eName" in df.columns:
                cols.append("eName")
        else:
            cols = df.columns[:3].tolist()

    df[cols].to_csv(filepath, sep="\t", index=False, header=header)
    logger.info(f"Saved {len(df)} records to {filepath.name}")


def create_output_dirs(
    base_dir: Union[str, Path],
    subdirs: Optional[List[str]] = None,
    prefix: Optional[str] = None,
) -> dict:
    """
    Create output directory structure.

    Args:
        base_dir: Base output directory
        subdirs: List of subdirectory names
        prefix: Optional prefix for directory names

    Returns:
        Dictionary mapping subdir names to Path objects
    """
    base_dir = Path(base_dir)

    if subdirs is None:
        subdirs = ["results", "figures", "plot_data", "temp"]

    dirs = {}
    for subdir in subdirs:
        if prefix:
            dir_name = f"{prefix}_{subdir}"
        else:
            dir_name = subdir
        dir_path = base_dir / dir_name
        dir_path.mkdir(parents=True, exist_ok=True)
        dirs[subdir] = dir_path

    return dirs


def merge_csv_files(
    input_dir: Union[str, Path],
    pattern: str = "*.csv",
    sample_col: str = "sample",
) -> pd.DataFrame:
    """
    Merge multiple CSV files and add sample column.

    Args:
        input_dir: Directory containing CSV files
        pattern: Glob pattern for files
        sample_col: Name of sample column to add

    Returns:
        Merged DataFrame
    """
    input_dir = Path(input_dir)
    files = sorted(input_dir.glob(pattern))

    if not files:
        raise FileNotFoundError(f"No files matching {pattern} in {input_dir}")

    dfs = []
    for f in files:
        df = pd.read_csv(f)
        df[sample_col] = f.stem
        dfs.append(df)

    merged = pd.concat(dfs, ignore_index=True)
    logger.info(f"Merged {len(files)} files, {len(merged)} total records")
    return merged
