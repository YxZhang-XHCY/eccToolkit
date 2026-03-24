"""
Format parsers for various eccDNA detection tool outputs.

Each parser reads a tool-specific output format and converts it to a
standardized pandas DataFrame compatible with BenchmarkEvaluator.

Required output columns:
    eccDNA_id, Regions, eccDNA_type, State, Length
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def parse_circleseeker(filepath: str) -> pd.DataFrame:
    """Parse CircleSeeker merged_output.csv.

    Expected columns: eccDNA_id, Regions, eccDNA_type, State, Length
    This is the native format — pass through with minimal validation.
    """
    logger.info(f"Parsing CircleSeeker output: {filepath}")
    df = pd.read_csv(filepath)

    required = {"eccDNA_id", "Regions", "eccDNA_type", "State", "Length"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"CircleSeeker file missing columns: {missing}. "
            f"Found: {list(df.columns)}"
        )
    return df


def parse_circlemap(filepath: str) -> pd.DataFrame:
    """Parse Circle-Map Realign output (BED-like format).

    Circle-Map outputs BED with columns:
        chrom, start, end, discordants, soft_clipped, ...
    All entries treated as UeccDNA (single-region).
    """
    logger.info(f"Parsing Circle-Map output: {filepath}")

    # Try reading with header first, fall back to no header
    df = pd.read_csv(filepath, sep="\t", header=None)
    if df.shape[1] < 3:
        raise ValueError(
            f"Circle-Map file has fewer than 3 columns: {filepath}"
        )

    chrom = df.iloc[:, 0]
    start = df.iloc[:, 1].astype(int)
    end = df.iloc[:, 2].astype(int)

    rows = []
    for i in range(len(df)):
        c, s, e = chrom.iloc[i], start.iloc[i], end.iloc[i]
        length = e - s
        region_str = f"{c}:{s}-{e}"
        rows.append({
            "eccDNA_id": f"circlemap_{i+1}",
            "Regions": region_str,
            "eccDNA_type": "UeccDNA",
            "State": "Confirmed",
            "Length": length,
        })

    result = pd.DataFrame(rows)
    logger.info(f"Parsed {len(result)} entries from Circle-Map output")
    return result


def parse_cresil(filepath: str) -> pd.DataFrame:
    """Parse CReSIL output.

    CReSIL outputs TSV/CSV with columns including:
        merge_region (e.g. "chr2:47235815-47237959_-"), merge_len, num_region, ...

    Entries with num_region > 1 are classified as MeccDNA (multi-fragment),
    otherwise as UeccDNA.
    """
    logger.info(f"Parsing CReSIL output: {filepath}")

    # Auto-detect separator
    with open(filepath) as f:
        first_line = f.readline()
    sep = "\t" if "\t" in first_line else ","

    df = pd.read_csv(filepath, sep=sep)

    # Find the relevant columns
    region_col = None
    len_col = None
    num_region_col = None

    for col in df.columns:
        col_lower = col.lower().strip()
        if "merge_region" in col_lower or col_lower == "region":
            region_col = col
        elif "merge_len" in col_lower or col_lower == "length":
            len_col = col
        elif "num_region" in col_lower:
            num_region_col = col

    if region_col is None:
        raise ValueError(
            f"CReSIL file missing region column. Found: {list(df.columns)}"
        )

    rows = []
    for i, row in df.iterrows():
        raw_region = str(row[region_col])
        # Parse CReSIL region format: "chr2:47235815-47237959_-" or multiple separated by ;
        regions = _parse_cresil_regions(raw_region)
        if not regions:
            continue

        region_str = ";".join(f"{c}:{s}-{e}" for c, s, e in regions)
        length = int(row[len_col]) if len_col and pd.notna(row[len_col]) else sum(e - s for _, s, e in regions)

        num_regions = 1
        if num_region_col and pd.notna(row[num_region_col]):
            num_regions = int(row[num_region_col])

        ecc_type = "MeccDNA" if num_regions > 1 else "UeccDNA"

        rows.append({
            "eccDNA_id": f"cresil_{i+1}",
            "Regions": region_str,
            "eccDNA_type": ecc_type,
            "State": "Confirmed",
            "Length": length,
        })

    result = pd.DataFrame(rows)
    logger.info(f"Parsed {len(result)} entries from CReSIL output")
    return result


def _parse_cresil_regions(region_str: str) -> list[tuple[str, int, int]]:
    """Parse CReSIL region format.

    Handles formats like:
        "chr2:47235815-47237959_-"
        "chr1:100-200_+;chr2:300-400_-"
    """
    regions = []
    # Split by ; for multiple regions
    parts = region_str.split(";")
    for part in parts:
        part = part.strip()
        # Remove strand suffix (_+, _-, etc.)
        part = re.sub(r"_[+-]$", "", part)
        if ":" in part and "-" in part:
            try:
                chrom, pos = part.rsplit(":", 1)
                start_str, end_str = pos.split("-", 1)
                regions.append((chrom, int(start_str), int(end_str)))
            except (ValueError, IndexError):
                pass
    return regions


def parse_eccfinder(filepath: str) -> pd.DataFrame:
    """Parse ecc_finder output (BED-like format).

    ecc_finder outputs BED format:
        chrom, start, end, name, score, strand
    All entries treated as UeccDNA.
    """
    logger.info(f"Parsing ecc_finder output: {filepath}")

    df = pd.read_csv(filepath, sep="\t", header=None)
    if df.shape[1] < 3:
        raise ValueError(
            f"ecc_finder file has fewer than 3 columns: {filepath}"
        )

    chrom = df.iloc[:, 0]
    start = df.iloc[:, 1].astype(int)
    end = df.iloc[:, 2].astype(int)

    rows = []
    for i in range(len(df)):
        c, s, e = chrom.iloc[i], start.iloc[i], end.iloc[i]
        length = e - s
        region_str = f"{c}:{s}-{e}"
        rows.append({
            "eccDNA_id": f"eccfinder_{i+1}",
            "Regions": region_str,
            "eccDNA_type": "UeccDNA",
            "State": "Confirmed",
            "Length": length,
        })

    result = pd.DataFrame(rows)
    logger.info(f"Parsed {len(result)} entries from ecc_finder output")
    return result


def parse_generic_bed(filepath: str) -> pd.DataFrame:
    """Parse generic BED file (chrom, start, end, ...).

    All entries treated as UeccDNA.
    """
    logger.info(f"Parsing generic BED: {filepath}")

    df = pd.read_csv(filepath, sep="\t", header=None, comment="#")
    if df.shape[1] < 3:
        raise ValueError(f"BED file has fewer than 3 columns: {filepath}")

    # Skip header row if present
    if str(df.iloc[0, 0]).lower() in ("chrom", "chr", "#chrom"):
        df = df.iloc[1:].reset_index(drop=True)

    chrom = df.iloc[:, 0]
    start = df.iloc[:, 1].astype(int)
    end = df.iloc[:, 2].astype(int)

    rows = []
    for i in range(len(df)):
        c, s, e = chrom.iloc[i], start.iloc[i], end.iloc[i]
        length = e - s
        region_str = f"{c}:{s}-{e}"
        rows.append({
            "eccDNA_id": f"bed_{i+1}",
            "Regions": region_str,
            "eccDNA_type": "UeccDNA",
            "State": "Confirmed",
            "Length": length,
        })

    result = pd.DataFrame(rows)
    logger.info(f"Parsed {len(result)} entries from BED file")
    return result


def parse_generic_csv(filepath: str) -> pd.DataFrame:
    """Parse generic CSV with chr/chrom, start, end columns.

    All entries treated as UeccDNA.
    """
    logger.info(f"Parsing generic CSV: {filepath}")

    df = pd.read_csv(filepath)

    # Find chrom/start/end columns (case-insensitive)
    col_map = {}
    for col in df.columns:
        cl = col.lower().strip()
        if cl in ("chr", "chrom", "chromosome", "#chrom"):
            col_map["chrom"] = col
        elif cl in ("start", "chromstart", "begin"):
            col_map["start"] = col
        elif cl in ("end", "chromend", "stop"):
            col_map["end"] = col

    if "chrom" not in col_map or "start" not in col_map or "end" not in col_map:
        raise ValueError(
            f"CSV file missing required columns (chr/chrom, start, end). "
            f"Found: {list(df.columns)}"
        )

    chrom = df[col_map["chrom"]]
    start = df[col_map["start"]].astype(int)
    end = df[col_map["end"]].astype(int)

    rows = []
    for i in range(len(df)):
        c, s, e = chrom.iloc[i], start.iloc[i], end.iloc[i]
        length = e - s
        region_str = f"{c}:{s}-{e}"
        rows.append({
            "eccDNA_id": f"csv_{i+1}",
            "Regions": region_str,
            "eccDNA_type": "UeccDNA",
            "State": "Confirmed",
            "Length": length,
        })

    result = pd.DataFrame(rows)
    logger.info(f"Parsed {len(result)} entries from CSV file")
    return result


# Format name to parser function mapping
FORMAT_PARSERS = {
    "circleseeker": parse_circleseeker,
    "circlemap": parse_circlemap,
    "cresil": parse_cresil,
    "eccfinder": parse_eccfinder,
    "bed": parse_generic_bed,
    "csv": parse_generic_csv,
}

SUPPORTED_FORMATS = list(FORMAT_PARSERS.keys())


def parse_tool_output(filepath: str, fmt: str) -> pd.DataFrame:
    """Parse a tool's output file using the specified format.

    Args:
        filepath: Path to the tool's output file
        fmt: Format name (one of SUPPORTED_FORMATS)

    Returns:
        Standardized DataFrame with columns:
        eccDNA_id, Regions, eccDNA_type, State, Length
    """
    fmt = fmt.lower().strip()
    if fmt not in FORMAT_PARSERS:
        raise ValueError(
            f"Unsupported format: '{fmt}'. "
            f"Supported formats: {SUPPORTED_FORMATS}"
        )
    return FORMAT_PARSERS[fmt](filepath)
