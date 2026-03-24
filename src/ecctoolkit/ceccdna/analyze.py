"""CeccDNA (chimeric eccDNA) analysis: extraction, classification, and summary statistics."""

import logging
import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _parse_location(loc_str: str) -> List[Dict]:
    """Parse location string like 'chr17:1000-2000(-);chr9:5000-6000(+)' into segment dicts."""
    segments = []
    if pd.isna(loc_str) or loc_str == ".":
        return segments
    for seg in loc_str.split(";"):
        seg = seg.strip()
        m = re.match(r"([^:]+):(\d+)-(\d+)\(([+-])\)", seg)
        if m:
            start = int(m.group(2))
            end = int(m.group(3))
            segments.append({
                "chr": m.group(1),
                "start": start,
                "end": end,
                "strand": m.group(4),
                "seg_length": end - start,
            })
    return segments


def _classify_fusion(segments: List[Dict]) -> str:
    """Classify a CeccDNA as inter-chromosomal or intra-chromosomal."""
    if len(segments) < 2:
        return "unknown"
    chrs = set(s["chr"] for s in segments)
    return "inter-chromosomal" if len(chrs) > 1 else "intra-chromosomal"


def _compute_intra_distance(segments: List[Dict]) -> float:
    """For intra-chromosomal CeccDNAs, compute max distance between segments."""
    if len(segments) < 2:
        return np.nan
    chrs = set(s["chr"] for s in segments)
    if len(chrs) > 1:
        return np.nan
    positions = []
    for s in segments:
        positions.extend([s["start"], s["end"]])
    return max(positions) - min(positions)


def _get_chromosomes_involved(segments: List[Dict]) -> List[str]:
    """Get sorted list of chromosomes involved."""
    return sorted(set(s["chr"] for s in segments))


def _process_single_sample(
    input_file: str,
    sample_name: str,
    location_col: str,
    type_col: str,
    type_value: str,
) -> Tuple[List[Dict], Dict]:
    """Process a single sample CSV and return (cecc_records, sample_summary)."""
    df = pd.read_csv(input_file)
    total = len(df)

    # Filter for CeccDNA
    if type_col in df.columns:
        cecc_df = df[df[type_col] == type_value].copy()
    else:
        logger.warning(
            f"Column '{type_col}' not found in {input_file}. "
            f"Available columns: {list(df.columns)}. Using all rows."
        )
        cecc_df = df.copy()

    cecc_n = len(cecc_df)
    logger.info(f"{sample_name}: {total} total eccDNA, {cecc_n} CeccDNA")

    # Parse each CeccDNA
    cecc_records = []
    for _, row in cecc_df.iterrows():
        loc_str = row.get(location_col, "")
        segments = _parse_location(str(loc_str))
        fusion_type = _classify_fusion(segments)
        intra_dist = _compute_intra_distance(segments)
        chrs_involved = _get_chromosomes_involved(segments)

        rec = {
            "sample": sample_name,
            "location": loc_str,
            "segment_count": len(segments),
            "total_length": sum(s["seg_length"] for s in segments),
            "fusion_type": fusion_type,
            "intra_distance": intra_dist,
            "n_chromosomes": len(chrs_involved),
            "chromosomes": ",".join(chrs_involved),
        }
        # Carry over useful columns from original data
        for col in ["eccDNA_id", "length", "read_count", "confidence"]:
            if col in row.index:
                rec[col] = row[col]

        cecc_records.append(rec)

    # Compute sample-level summary
    if cecc_n > 0:
        lengths = cecc_df["length"].values if "length" in cecc_df.columns else np.array([r["total_length"] for r in cecc_records])
        seg_counts = np.array([r["segment_count"] for r in cecc_records])
        fusion_types = [r["fusion_type"] for r in cecc_records]
        inter_n = sum(1 for f in fusion_types if f == "inter-chromosomal")
        intra_n = sum(1 for f in fusion_types if f == "intra-chromosomal")
        intra_dists = [r["intra_distance"] for r in cecc_records if not np.isnan(r["intra_distance"])]
    else:
        lengths = np.array([])
        seg_counts = np.array([])
        inter_n = intra_n = 0
        intra_dists = []

    summary = {
        "sample": sample_name,
        "total_eccDNA": total,
        "CeccDNA_n": cecc_n,
        "CeccDNA_pct": cecc_n / total * 100 if total > 0 else 0,
        "CeccDNA_mean_length": np.mean(lengths) if len(lengths) > 0 else np.nan,
        "CeccDNA_median_length": np.median(lengths) if len(lengths) > 0 else np.nan,
        "CeccDNA_min_length": np.min(lengths) if len(lengths) > 0 else np.nan,
        "CeccDNA_max_length": np.max(lengths) if len(lengths) > 0 else np.nan,
        "mean_segment_count": np.mean(seg_counts) if len(seg_counts) > 0 else np.nan,
        "seg2_count": int(np.sum(seg_counts == 2)) if len(seg_counts) > 0 else 0,
        "seg3_count": int(np.sum(seg_counts == 3)) if len(seg_counts) > 0 else 0,
        "seg4plus_count": int(np.sum(seg_counts >= 4)) if len(seg_counts) > 0 else 0,
        "inter_chromosomal_n": inter_n,
        "intra_chromosomal_n": intra_n,
        "inter_chromosomal_pct": inter_n / cecc_n * 100 if cecc_n > 0 else 0,
        "intra_chromosomal_pct": intra_n / cecc_n * 100 if cecc_n > 0 else 0,
        "intra_median_distance": np.median(intra_dists) if intra_dists else np.nan,
        "intra_mean_distance": np.mean(intra_dists) if intra_dists else np.nan,
    }

    return cecc_records, summary


def run_ceccdna_analysis(
    input_files: List[str],
    output_dir: str,
    sample_names: Optional[List[str]] = None,
    location_col: str = "location",
    type_col: str = "type",
    type_value: str = "Cecc",
) -> None:
    """Run CeccDNA analysis on one or more eccDNA CSV files.

    For each input file, filters for CeccDNA records (type == type_value),
    parses the location string to extract segment coordinates, classifies
    each CeccDNA as inter- or intra-chromosomal, and computes per-sample
    and cross-sample summary statistics.

    Args:
        input_files: List of paths to eccDNA CSV files.
        output_dir: Directory to write output files.
        sample_names: Optional sample names corresponding to input_files.
            If None, derived from input file basenames.
        location_col: Column name containing the location string.
        type_col: Column name containing the eccDNA type.
        type_value: Value in type_col that identifies CeccDNA records.

    Output files:
        ceccdna_records.csv: One row per CeccDNA with parsed segment info.
        ceccdna_sample_summary.csv: Per-sample statistics.
        ceccdna_segment_stats.csv: Segment count distribution per sample.
    """
    os.makedirs(output_dir, exist_ok=True)

    if sample_names is None:
        sample_names = [
            os.path.splitext(os.path.basename(f))[0] for f in input_files
        ]

    if len(sample_names) != len(input_files):
        raise ValueError(
            f"Number of sample names ({len(sample_names)}) does not match "
            f"number of input files ({len(input_files)})"
        )

    all_records = []
    all_summaries = []

    for input_file, sample_name in zip(input_files, sample_names):
        logger.info(f"Processing {sample_name}: {input_file}")
        records, summary = _process_single_sample(
            input_file, sample_name, location_col, type_col, type_value
        )
        all_records.extend(records)
        all_summaries.append(summary)

    # Save CeccDNA records
    records_path = os.path.join(output_dir, "ceccdna_records.csv")
    if all_records:
        pd.DataFrame(all_records).to_csv(records_path, index=False)
    else:
        pd.DataFrame().to_csv(records_path, index=False)
    logger.info(f"Wrote {len(all_records)} CeccDNA records to {records_path}")

    # Save sample summary
    summary_path = os.path.join(output_dir, "ceccdna_sample_summary.csv")
    pd.DataFrame(all_summaries).to_csv(summary_path, index=False)
    logger.info(f"Wrote sample summary to {summary_path}")

    # Save segment count distribution
    seg_stats = []
    for s in all_summaries:
        seg_stats.append({
            "sample": s["sample"],
            "CeccDNA_n": s["CeccDNA_n"],
            "seg2_count": s["seg2_count"],
            "seg3_count": s["seg3_count"],
            "seg4plus_count": s["seg4plus_count"],
            "seg2_pct": s["seg2_count"] / s["CeccDNA_n"] * 100 if s["CeccDNA_n"] > 0 else 0,
            "seg3_pct": s["seg3_count"] / s["CeccDNA_n"] * 100 if s["CeccDNA_n"] > 0 else 0,
            "seg4plus_pct": s["seg4plus_count"] / s["CeccDNA_n"] * 100 if s["CeccDNA_n"] > 0 else 0,
        })
    seg_stats_path = os.path.join(output_dir, "ceccdna_segment_stats.csv")
    pd.DataFrame(seg_stats).to_csv(seg_stats_path, index=False)
    logger.info(f"Wrote segment stats to {seg_stats_path}")
