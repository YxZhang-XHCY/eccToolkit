"""Data format conversion utilities."""

import logging
from pathlib import Path
from typing import List

import pandas as pd

logger = logging.getLogger(__name__)


def convert_gff_to_csv(
    input_dir: str,
    output_file: str,
) -> None:
    """
    Batch convert GFF files to merged CSV.

    Parses GFF attributes and combines all files into single CSV.

    Args:
        input_dir: Directory containing GFF files
        output_file: Output CSV file
    """
    input_path = Path(input_dir)
    gff_files = list(input_path.glob("*.gff")) + list(input_path.glob("*.gff3"))

    if not gff_files:
        raise FileNotFoundError(f"No GFF files found in {input_dir}")

    logger.info(f"Found {len(gff_files)} GFF files")

    all_records: List[dict] = []
    for gff_file in gff_files:
        # TODO: Implement GFF parsing
        logger.info(f"Processing {gff_file.name}")

    # TODO: Migrate from phase1_parse_and_stats.py
    raise NotImplementedError("This module is pending migration")


def convert_cnv_to_bed(
    input_file: str,
    output_dir: str,
) -> None:
    """
    Convert CNV TSV to per-cell-line BED files.

    Converts 1-based coordinates to 0-based BED format.

    Args:
        input_file: CNV TSV file (1-based coordinates)
        output_dir: Output directory for BED files
    """
    logger.info("CNV to BED conversion - placeholder")
    logger.info(f"Input: {input_file}")
    # TODO: Migrate from cnv_to_bed.py
    raise NotImplementedError("This module is pending migration")
