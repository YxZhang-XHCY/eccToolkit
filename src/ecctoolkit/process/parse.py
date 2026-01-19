"""Parse eccDNA CSV files."""

import logging
import re

import pandas as pd

logger = logging.getLogger(__name__)


def parse_eccdna(input_file: str, output_file: str) -> None:
    """Parse eccDNA CSV and extract coordinates from seqname."""
    df = pd.read_csv(input_file)

    # Extract chr, start, end from seqname if needed
    if "seqname" in df.columns and "eChr" not in df.columns:
        pattern = r"([^:]+):(\d+)-(\d+)"

        def extract_coords(seqname):
            match = re.match(pattern, str(seqname))
            if match:
                return match.groups()
            return None, None, None

        coords = df["seqname"].apply(extract_coords)
        df["eChr"] = [c[0] for c in coords]
        df["eStart"] = [int(c[1]) if c[1] else None for c in coords]
        df["eEnd"] = [int(c[2]) if c[2] else None for c in coords]

    # Calculate length if missing
    if "eLength" not in df.columns and "eStart" in df.columns and "eEnd" in df.columns:
        df["eLength"] = df["eEnd"] - df["eStart"]

    df.to_csv(output_file, index=False)
    logger.info(f"Parsed {len(df)} records -> {output_file}")
