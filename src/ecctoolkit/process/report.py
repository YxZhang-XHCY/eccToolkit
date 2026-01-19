"""Generate summary reports."""

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def generate_reports(input_file: str, output_prefix: str) -> None:
    """Generate summary reports from combined data."""
    df = pd.read_csv(input_file)
    output_path = Path(output_prefix)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Sample counts
    if "sample" in df.columns:
        sample_counts = df.groupby("sample").size().reset_index(name="count")
        sample_counts.to_csv(f"{output_prefix}_sample_counts.csv", index=False)
        logger.info(f"Generated sample counts report")

    # Length distribution
    if "eLength" in df.columns:
        length_stats = df["eLength"].describe()
        length_stats.to_csv(f"{output_prefix}_length_stats.csv")
        logger.info(f"Generated length statistics report")

    # Summary
    summary = {
        "total_records": len(df),
        "columns": list(df.columns),
    }
    if "sample" in df.columns:
        summary["n_samples"] = df["sample"].nunique()

    logger.info(f"Reports generated with prefix: {output_prefix}")
