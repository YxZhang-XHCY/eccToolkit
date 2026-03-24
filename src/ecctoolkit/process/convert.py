"""Data format conversion utilities."""

import logging
import os
from pathlib import Path
from typing import Dict, List

import pandas as pd

logger = logging.getLogger(__name__)


def _parse_gff_attributes(attr_string: str) -> Dict[str, str]:
    """Parse GFF3 attribute string (key=value;key=value;...) into dict."""
    attrs = {}
    if not attr_string or pd.isna(attr_string) or attr_string == ".":
        return attrs
    for item in str(attr_string).split(";"):
        item = item.strip()
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs


def _parse_gff_file(filepath: Path) -> List[dict]:
    """Parse a single GFF3 file into list of records."""
    records = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue

            attrs = _parse_gff_attributes(parts[8])

            record = {
                "seqid": parts[0],
                "source": parts[1],
                "type": parts[2],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "score": parts[5] if parts[5] != "." else None,
                "strand": parts[6] if parts[6] != "." else None,
                "phase": parts[7] if parts[7] != "." else None,
                "attributes_raw": parts[8],
            }
            # Extract common attributes as separate columns
            for key in ("ID", "Name", "Parent", "gene_id", "gene_name",
                        "transcript_id", "gene_biotype"):
                if key in attrs:
                    record[key] = attrs[key]

            record["source_file"] = filepath.name
            records.append(record)

    return records


def convert_gff_to_csv(
    input_dir: str,
    output_file: str,
) -> None:
    """
    Batch convert GFF files to merged CSV.

    Parses GFF3 attributes and combines all files into single CSV.

    Args:
        input_dir: Directory containing GFF files
        output_file: Output CSV file
    """
    input_path = Path(input_dir)
    gff_files = sorted(
        list(input_path.glob("*.gff")) + list(input_path.glob("*.gff3"))
    )

    if not gff_files:
        raise FileNotFoundError(f"No GFF files found in {input_dir}")

    logger.info(f"Found {len(gff_files)} GFF files")

    all_records: List[dict] = []
    for gff_file in gff_files:
        logger.info(f"Processing {gff_file.name}")
        records = _parse_gff_file(gff_file)
        all_records.extend(records)
        logger.info(f"  Parsed {len(records)} records")

    if not all_records:
        logger.warning("No records parsed from GFF files")
        return

    df = pd.DataFrame(all_records)
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    df.to_csv(output_file, index=False)
    logger.info(f"Saved {len(df)} records to {output_file}")

    # Summary by type
    type_counts = df["type"].value_counts()
    for feat_type, count in type_counts.head(10).items():
        logger.info(f"  {feat_type}: {count}")


def convert_cnv_to_bed(
    input_file: str,
    output_dir: str,
) -> None:
    """
    Convert CNV TSV to per-cell-line BED files.

    Converts 1-based coordinates to 0-based BED format.
    If a cell_line column exists, splits into separate BED files.

    Args:
        input_file: CNV TSV file (1-based coordinates)
        output_dir: Output directory for BED files
    """
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(input_file, sep="\t")
    logger.info(f"Loaded {len(df)} CNV records from {input_file}")

    # Normalize column names
    col_map = {}
    for c in df.columns:
        cl = c.lower().strip()
        if cl in ("chr", "chrom", "chromosome", "#chrom", "#chr"):
            col_map[c] = "chr"
        elif cl == "start":
            col_map[c] = "start"
        elif cl == "end":
            col_map[c] = "end"
        elif cl in ("cell_line", "cellline", "sample", "cell"):
            col_map[c] = "cell_line"
        elif cl in ("name", "id", "cnv_id"):
            col_map[c] = "name"
        elif cl in ("cn", "copy_number", "copynumber"):
            col_map[c] = "copy_number"
    df = df.rename(columns=col_map)

    if "chr" not in df.columns or "start" not in df.columns or "end" not in df.columns:
        raise ValueError(
            f"Cannot find chr/start/end columns in {input_file}. "
            f"Available columns: {list(df.columns)}"
        )

    # Convert 1-based to 0-based (BED format)
    df["start"] = df["start"].astype(int) - 1
    df["end"] = df["end"].astype(int)

    # Build BED columns
    bed_cols = ["chr", "start", "end"]
    if "name" in df.columns:
        bed_cols.append("name")
    if "copy_number" in df.columns:
        bed_cols.append("copy_number")

    if "cell_line" in df.columns:
        # Split by cell_line
        groups = df.groupby("cell_line")
        logger.info(f"Splitting into {len(groups)} cell-line BED files")
        for cell_line, group_df in groups:
            safe_name = str(cell_line).replace("/", "_").replace(" ", "_")
            out_file = os.path.join(output_dir, f"{safe_name}.bed")
            group_df[bed_cols].to_csv(out_file, sep="\t", index=False, header=False)
            logger.info(f"  {cell_line}: {len(group_df)} records -> {out_file}")
    else:
        # Single output file
        out_file = os.path.join(output_dir, "cnv.bed")
        df[bed_cols].to_csv(out_file, sep="\t", index=False, header=False)
        logger.info(f"Saved {len(df)} records to {out_file}")
