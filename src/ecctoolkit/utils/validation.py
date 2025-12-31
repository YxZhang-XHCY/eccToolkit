"""Data validation utilities for eccToolkit."""

import logging
import shutil
from typing import List, Optional

import pandas as pd

logger = logging.getLogger(__name__)


def validate_eccdna_columns(
    df: pd.DataFrame,
    required_cols: Optional[List[str]] = None,
) -> bool:
    """
    Validate that DataFrame has required eccDNA columns.

    Args:
        df: DataFrame to validate
        required_cols: List of required column names

    Returns:
        True if valid, raises ValueError otherwise
    """
    if required_cols is None:
        required_cols = ["eChr", "eStart", "eEnd"]

    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    return True


def filter_chromosomes(
    df: pd.DataFrame,
    chr_col: str = "eChr",
    keep_sex: bool = False,
    keep_mito: bool = False,
    keep_alt: bool = False,
) -> pd.DataFrame:
    """
    Filter DataFrame by chromosome.

    Args:
        df: DataFrame with chromosome column
        chr_col: Name of chromosome column
        keep_sex: Keep sex chromosomes (chrX, chrY)
        keep_mito: Keep mitochondrial chromosome
        keep_alt: Keep alternative/unplaced contigs (containing "_")

    Returns:
        Filtered DataFrame
    """
    original_len = len(df)
    filtered = df.copy()

    # Filter mitochondrial
    if not keep_mito:
        filtered = filtered[~filtered[chr_col].str.contains("chrM", na=False)]

    # Filter sex chromosomes
    if not keep_sex:
        filtered = filtered[~filtered[chr_col].str.contains("chrX|chrY", na=False)]

    # Filter alternative/unplaced contigs
    if not keep_alt:
        filtered = filtered[~filtered[chr_col].str.contains("_", na=False)]

    removed = original_len - len(filtered)
    if removed > 0:
        logger.info(f"Filtered {removed} records by chromosome")

    return filtered


def check_dependencies(tools: List[str]) -> dict:
    """
    Check if required command-line tools are installed.

    Args:
        tools: List of tool names to check

    Returns:
        Dictionary mapping tool names to availability status
    """
    status = {}
    for tool in tools:
        path = shutil.which(tool)
        status[tool] = path is not None
        if path:
            logger.debug(f"{tool}: found at {path}")
        else:
            logger.warning(f"{tool}: not found in PATH")

    return status


def require_dependencies(tools: List[str]) -> None:
    """
    Require that all specified tools are available.

    Args:
        tools: List of required tool names

    Raises:
        RuntimeError: If any tool is missing
    """
    status = check_dependencies(tools)
    missing = [tool for tool, available in status.items() if not available]
    if missing:
        raise RuntimeError(
            f"Missing required tools: {missing}. "
            "Please install them and ensure they are in your PATH."
        )


def validate_file_exists(filepath: str, description: str = "File") -> None:
    """
    Validate that a file exists.

    Args:
        filepath: Path to check
        description: Description for error message

    Raises:
        FileNotFoundError: If file doesn't exist
    """
    from pathlib import Path

    if not Path(filepath).exists():
        raise FileNotFoundError(f"{description} not found: {filepath}")
