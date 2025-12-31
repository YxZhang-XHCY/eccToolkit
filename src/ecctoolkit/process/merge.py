"""Merge CSV files."""

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def merge_files(
    input_dir: str,
    output_file: str,
    pattern: str = "*.csv",
) -> None:
    """Merge multiple CSV files with sample column."""
    input_path = Path(input_dir)
    files = sorted(input_path.glob(pattern))

    if not files:
        raise FileNotFoundError(f"No files matching {pattern} in {input_dir}")

    dfs = []
    for f in files:
        df = pd.read_csv(f)
        df["sample"] = f.stem
        dfs.append(df)
        logger.info(f"Loaded {len(df)} records from {f.name}")

    merged = pd.concat(dfs, ignore_index=True)
    merged.to_csv(output_file, index=False)
    logger.info(f"Merged {len(files)} files, {len(merged)} total records -> {output_file}")
