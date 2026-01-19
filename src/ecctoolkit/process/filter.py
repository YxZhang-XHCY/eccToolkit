"""Filter eccDNA data by annotation percentage."""

import logging

import pandas as pd

logger = logging.getLogger(__name__)


def filter_by_annotation(
    input_file: str,
    output_file: str,
    min_percent: float = 80,
) -> None:
    """Filter by minimum annotation percentage."""
    df = pd.read_csv(input_file)

    if "anno_Percent" not in df.columns:
        raise ValueError("Column 'anno_Percent' not found in input file")

    original_count = len(df)
    filtered = df[df["anno_Percent"] >= min_percent]
    filtered.to_csv(output_file, index=False)

    logger.info(
        f"Filtered {original_count} -> {len(filtered)} records "
        f"(anno_Percent >= {min_percent})"
    )
