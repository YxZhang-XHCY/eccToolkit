"""TE composition analysis for eccDNA using RepeatMasker/GFF/BED annotation."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# TE annotation loaders
# ---------------------------------------------------------------------------

def _detect_te_format(filepath: str) -> str:
    """Auto-detect TE annotation format from file content."""
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            # UCSC rmsk.txt: 17 tab-separated fields, first field is bin (int)
            if len(fields) >= 17:
                try:
                    int(fields[0])  # bin
                    int(fields[1])  # swScore
                    return "rmsk"
                except ValueError:
                    pass
            # GFF/GTF: 9 tab-separated fields with attributes
            if len(fields) >= 9 and fields[2] in (
                "repeat_region", "transposable_element", "dispersed_repeat",
                "match", "similarity", "gene", "exon",
            ):
                return "gff"
            # RepeatMasker .out: space-delimited, starts with score
            space_fields = line.split()
            if len(space_fields) >= 15:
                try:
                    int(space_fields[0])  # SW score
                    float(space_fields[1])  # perc_div
                    return "rmsk_out"
                except ValueError:
                    pass
            # BED format: 4+ tab fields, second and third are int
            if len(fields) >= 4:
                try:
                    int(fields[1])
                    int(fields[2])
                    return "bed"
                except ValueError:
                    pass
            break
    return "bed"


def _load_te_rmsk(filepath: str) -> pd.DataFrame:
    """Load UCSC rmsk.txt format.

    Columns: bin, swScore, milliDiv, milliDel, milliIns, genoName,
             genoStart, genoEnd, genoLeft, strand, repName, repClass,
             repFamily, repStart, repEnd, repLeft, id
    """
    col_names = [
        "bin", "swScore", "milliDiv", "milliDel", "milliIns",
        "chrom", "start", "end", "genoLeft",
        "strand", "te_name", "te_class", "te_family",
        "repStart", "repEnd", "repLeft", "id",
    ]
    df = pd.read_csv(filepath, sep="\t", header=None, comment="#",
                     names=col_names, low_memory=False)
    df["pct_div"] = df["milliDiv"] / 10.0
    return df[["chrom", "start", "end", "strand", "te_name",
               "te_class", "te_family", "pct_div"]].copy()


def _load_te_rmsk_out(filepath: str) -> pd.DataFrame:
    """Load RepeatMasker .out format (space-delimited)."""
    records = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("SW") or line.startswith("score"):
                continue
            if line.endswith("*"):
                line = line[:-1].strip()
            fields = line.split()
            if len(fields) < 15:
                continue
            try:
                pct_div = float(fields[1])
                chrom = fields[4]
                start = int(fields[5]) - 1  # 1-based to 0-based
                end = int(fields[6])
                strand = "-" if fields[8] == "C" else "+"
                te_name = fields[9]
                class_family = fields[10]
                if "/" in class_family:
                    te_class, te_family = class_family.split("/", 1)
                else:
                    te_class = class_family
                    te_family = class_family
                records.append({
                    "chrom": chrom, "start": start, "end": end,
                    "strand": strand, "te_name": te_name,
                    "te_class": te_class, "te_family": te_family,
                    "pct_div": pct_div,
                })
            except (ValueError, IndexError):
                continue
    return pd.DataFrame(records)


def _load_te_gff(filepath: str) -> pd.DataFrame:
    """Load GFF/GFF3 format TE annotation."""
    records = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            chrom = fields[0]
            start = int(fields[3]) - 1  # GFF is 1-based
            end = int(fields[4])
            strand = fields[6]
            attrs = fields[8]
            te_name = _parse_gff_attr(attrs, "Name") or _parse_gff_attr(attrs, "ID") or ""
            te_class = _parse_gff_attr(attrs, "class") or _parse_gff_attr(attrs, "repeat_class") or fields[2]
            te_family = _parse_gff_attr(attrs, "family") or _parse_gff_attr(attrs, "repeat_family") or te_class
            pct_div_str = _parse_gff_attr(attrs, "pct_div") or _parse_gff_attr(attrs, "divergence")
            pct_div = float(pct_div_str) if pct_div_str else np.nan
            records.append({
                "chrom": chrom, "start": start, "end": end,
                "strand": strand, "te_name": te_name,
                "te_class": te_class, "te_family": te_family,
                "pct_div": pct_div,
            })
    return pd.DataFrame(records)


def _parse_gff_attr(attrs: str, key: str) -> Optional[str]:
    """Extract attribute value from GFF attributes string."""
    # GFF3: key=value;key2=value2
    for part in attrs.split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            if k.strip() == key:
                return v.strip()
        # GTF: key "value"
        elif part.startswith(f'{key} "'):
            return part.split('"')[1]
    return None


def _load_te_bed(filepath: str) -> pd.DataFrame:
    """Load BED format TE annotation.

    Supports two layouts:
      - Standard BED6+: chr, start, end, te_name, score, strand, te_class, te_family
      - Minimal BED4:   chr, start, end, te_name
    """
    df = pd.read_csv(filepath, sep="\t", header=None, comment="#")
    n_cols = df.shape[1]

    # Standard BED: col0=chr, col1=start, col2=end, col3=name, col4=score, col5=strand
    # TE-extended:  col6=te_class, col7=te_family
    cols = {0: "chrom", 1: "start", 2: "end"}
    if n_cols >= 4:
        cols[3] = "te_name"
    if n_cols >= 6:
        cols[5] = "strand"
    if n_cols >= 8:
        # Standard BED6 + te_class + te_family in cols 6-7
        cols[6] = "te_class"
        cols[7] = "te_family"
    elif n_cols >= 5:
        # Check if col4 looks like a score (numeric) or a TE class name
        sample_val = str(df.iloc[0, 4]) if len(df) > 0 else ""
        try:
            float(sample_val)
            is_score = True
        except ValueError:
            is_score = False

        if is_score and n_cols >= 7:
            # BED6 + te_class in col6
            cols[6] = "te_class"
        elif not is_score:
            # Non-standard: col4=te_class, col5=te_family (no score/strand)
            cols[4] = "te_class"
            if n_cols >= 6:
                cols[5] = "te_family"

    df = df.rename(columns=cols)
    for col in ["te_name", "te_class", "te_family", "strand"]:
        if col not in df.columns:
            df[col] = "Unknown" if col != "strand" else "+"
    df["pct_div"] = np.nan
    return df[["chrom", "start", "end", "strand", "te_name",
               "te_class", "te_family", "pct_div"]].copy()


def load_te_annotation(filepath: str, te_format: str = "auto") -> pd.DataFrame:
    """Load TE annotation from various formats.

    Args:
        filepath: Path to TE annotation file.
        te_format: One of "rmsk", "rmsk_out", "gff", "bed", "auto".

    Returns:
        DataFrame with columns: chrom, start, end, strand, te_name,
        te_class, te_family, pct_div.
    """
    if te_format == "auto":
        te_format = _detect_te_format(filepath)
        logger.info(f"Auto-detected TE format: {te_format}")

    loaders = {
        "rmsk": _load_te_rmsk,
        "rmsk_out": _load_te_rmsk_out,
        "gff": _load_te_gff,
        "bed": _load_te_bed,
    }
    if te_format not in loaders:
        raise ValueError(f"Unknown TE format: {te_format}. Use one of {list(loaders.keys())}")
    df = loaders[te_format](filepath)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    logger.info(f"Loaded {len(df)} TE annotations from {filepath}")
    return df


# ---------------------------------------------------------------------------
# eccDNA loader
# ---------------------------------------------------------------------------

def _load_eccdna(filepath: str) -> pd.DataFrame:
    """Load eccDNA regions from CSV or BED.

    Returns DataFrame with columns: chrom, start, end, name, type (if available).
    """
    ext = Path(filepath).suffix.lower()
    if ext == ".bed":
        df = pd.read_csv(filepath, sep="\t", header=None, comment="#")
        df = df.iloc[:, :min(df.shape[1], 6)]
        col_map = {0: "chrom", 1: "start", 2: "end"}
        if df.shape[1] >= 4:
            col_map[3] = "name"
        if df.shape[1] >= 5:
            col_map[4] = "score"
        if df.shape[1] >= 6:
            col_map[5] = "strand"
        df = df.rename(columns=col_map)
    else:
        df = pd.read_csv(filepath)
        # Normalize column names
        col_remap = {}
        for target, candidates in [
            ("chrom", ["eChr", "chr", "chrom", "chromosome"]),
            ("start", ["eStart", "start", "chromStart"]),
            ("end", ["eEnd", "end", "chromEnd"]),
            ("name", ["eccDNA_name", "name", "eccDNA_id", "id"]),
            ("type", ["type", "eccDNA_type", "ecc_type"]),
        ]:
            for c in candidates:
                if c in df.columns:
                    col_remap[c] = target
                    break
        df = df.rename(columns=col_remap)

    if "chrom" not in df.columns or "start" not in df.columns or "end" not in df.columns:
        raise ValueError(f"Cannot find chrom/start/end columns in {filepath}")

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    if "name" not in df.columns:
        df["name"] = [f"eccDNA_{i}" for i in range(len(df))]

    return df


# ---------------------------------------------------------------------------
# Interval overlap
# ---------------------------------------------------------------------------

def _interval_overlap(s1: int, e1: int, s2: int, e2: int) -> int:
    """Compute overlap length between two intervals."""
    return max(0, min(e1, e2) - max(s1, s2))


def _intersect_eccdna_te(
    eccdna_df: pd.DataFrame, te_df: pd.DataFrame
) -> pd.DataFrame:
    """Intersect eccDNA regions with TE annotations.

    Uses a chromosome-based sweep for efficiency.

    Returns a DataFrame with one row per eccDNA-TE overlap, containing:
    eccDNA columns + te columns + overlap_bp.
    """
    te_by_chrom = {}
    for chrom, group in te_df.groupby("chrom"):
        sorted_group = group.sort_values("start")
        te_by_chrom[chrom] = (
            sorted_group["start"].values,
            sorted_group["end"].values,
            sorted_group["te_name"].values,
            sorted_group["te_class"].values,
            sorted_group["te_family"].values,
            sorted_group["pct_div"].values,
        )

    records = []
    for _, ecc in eccdna_df.iterrows():
        chrom = ecc["chrom"]
        ecc_start = ecc["start"]
        ecc_end = ecc["end"]
        ecc_name = ecc["name"]
        ecc_type = ecc.get("type", "")
        ecc_len = ecc_end - ecc_start

        if chrom not in te_by_chrom:
            continue

        te_starts, te_ends, te_names, te_classes, te_families, te_divs = te_by_chrom[chrom]

        # Binary search for first TE that could overlap
        idx = np.searchsorted(te_ends, ecc_start, side="right")
        for i in range(idx, len(te_starts)):
            if te_starts[i] >= ecc_end:
                break
            ov = _interval_overlap(ecc_start, ecc_end, te_starts[i], te_ends[i])
            if ov > 0:
                records.append({
                    "eccDNA_name": ecc_name,
                    "eccDNA_chr": chrom,
                    "eccDNA_start": ecc_start,
                    "eccDNA_end": ecc_end,
                    "eccDNA_len": ecc_len,
                    "eccDNA_type": ecc_type,
                    "te_name": te_names[i],
                    "te_class": str(te_classes[i]),
                    "te_family": str(te_families[i]),
                    "te_pct_div": te_divs[i],
                    "overlap_bp": ov,
                })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def run_te_analysis(
    input_file: str,
    output_dir: str,
    te_annotation: str = "",
    te_format: str = "auto",
    threads: int = 8,
) -> None:
    """Analyze TE composition of eccDNA using TE annotation.

    Intersects eccDNA regions with TE annotations to compute per-eccDNA
    TE coverage, class breakdown, and family breakdown. Optionally compares
    TE content between eccDNA types (e.g. Mecc vs Uecc).

    Args:
        input_file: eccDNA regions file (CSV or BED with chr, start, end).
        output_dir: Output directory for results.
        te_annotation: TE annotation file (RepeatMasker .out, UCSC rmsk.txt,
            GFF, or BED). If empty, input_file is treated as pre-computed
            TE annotation CSV.
        te_format: Format of TE annotation file. One of "rmsk", "rmsk_out",
            "gff", "bed", "auto".
        threads: Number of threads (reserved for future use).
    """
    os.makedirs(output_dir, exist_ok=True)

    if not te_annotation:
        raise ValueError("te_annotation is required: path to RepeatMasker/GFF/BED file")

    # Load data
    logger.info(f"Loading eccDNA regions from {input_file}")
    eccdna_df = _load_eccdna(input_file)
    logger.info(f"Loaded {len(eccdna_df)} eccDNA regions")

    logger.info(f"Loading TE annotations from {te_annotation}")
    te_df = load_te_annotation(te_annotation, te_format)

    # Intersect
    logger.info("Intersecting eccDNA with TE annotations...")
    overlaps = _intersect_eccdna_te(eccdna_df, te_df)
    logger.info(f"Found {len(overlaps)} eccDNA-TE overlaps")

    # --- Per-eccDNA TE annotation ---
    if overlaps.empty:
        logger.warning("No overlaps found between eccDNA and TE annotations")
        per_ecc = eccdna_df[["name", "chrom", "start", "end"]].copy()
        per_ecc = per_ecc.rename(columns={"name": "eccDNA_name", "chrom": "eccDNA_chr",
                                           "start": "eccDNA_start", "end": "eccDNA_end"})
        per_ecc["eccDNA_len"] = per_ecc["eccDNA_end"] - per_ecc["eccDNA_start"]
        per_ecc["te_bp"] = 0
        per_ecc["te_pct"] = 0.0
        per_ecc["n_te_hits"] = 0
        per_ecc["te_classes"] = ""
        per_ecc["te_families"] = ""
        per_ecc.to_csv(os.path.join(output_dir, "per_eccdna_te_annotation.csv"), index=False)
        pd.DataFrame().to_csv(os.path.join(output_dir, "te_summary.csv"), index=False)
        return

    # Save raw overlaps
    overlaps.to_csv(os.path.join(output_dir, "eccdna_te_overlaps.csv"), index=False)

    # Per-eccDNA summary
    per_ecc = overlaps.groupby("eccDNA_name").agg(
        eccDNA_chr=("eccDNA_chr", "first"),
        eccDNA_start=("eccDNA_start", "first"),
        eccDNA_end=("eccDNA_end", "first"),
        eccDNA_len=("eccDNA_len", "first"),
        eccDNA_type=("eccDNA_type", "first"),
        te_bp=("overlap_bp", "sum"),
        n_te_hits=("overlap_bp", "count"),
        te_classes=("te_class", lambda x: ",".join(sorted(set(x)))),
        te_families=("te_family", lambda x: ",".join(sorted(set(x)))),
    ).reset_index()

    # Add eccDNA with no TE overlap
    no_overlap_names = set(eccdna_df["name"]) - set(per_ecc["eccDNA_name"])
    if no_overlap_names:
        no_ov = eccdna_df[eccdna_df["name"].isin(no_overlap_names)].copy()
        no_ov_records = []
        for _, row in no_ov.iterrows():
            no_ov_records.append({
                "eccDNA_name": row["name"],
                "eccDNA_chr": row["chrom"],
                "eccDNA_start": row["start"],
                "eccDNA_end": row["end"],
                "eccDNA_len": row["end"] - row["start"],
                "eccDNA_type": row.get("type", ""),
                "te_bp": 0, "n_te_hits": 0,
                "te_classes": "", "te_families": "",
            })
        per_ecc = pd.concat([per_ecc, pd.DataFrame(no_ov_records)], ignore_index=True)

    # Cap te_bp at eccDNA length (overlapping TEs)
    per_ecc["te_bp"] = per_ecc[["te_bp", "eccDNA_len"]].min(axis=1)
    per_ecc["te_pct"] = (per_ecc["te_bp"] / per_ecc["eccDNA_len"] * 100).round(2)

    per_ecc_out = os.path.join(output_dir, "per_eccdna_te_annotation.csv")
    per_ecc.to_csv(per_ecc_out, index=False)
    logger.info(f"Saved per-eccDNA TE annotation to {per_ecc_out}")

    # --- Per-eccDNA TE class breakdown ---
    class_bp = overlaps.groupby(["eccDNA_name", "te_class"])["overlap_bp"].sum().reset_index()
    class_bp = class_bp.rename(columns={"overlap_bp": "class_bp"})
    # Merge eccDNA length
    class_bp = class_bp.merge(
        per_ecc[["eccDNA_name", "eccDNA_len"]], on="eccDNA_name", how="left"
    )
    class_bp["class_pct"] = (class_bp["class_bp"] / class_bp["eccDNA_len"] * 100).round(2)
    class_bp_out = os.path.join(output_dir, "per_eccdna_te_class_breakdown.csv")
    class_bp.to_csv(class_bp_out, index=False)
    logger.info(f"Saved per-eccDNA TE class breakdown to {class_bp_out}")

    # --- Summary statistics ---
    summary_records = []
    total_ecc = len(per_ecc)
    te_positive = (per_ecc["te_bp"] > 0).sum()
    summary_records.append({
        "metric": "total_eccDNA", "value": total_ecc,
    })
    summary_records.append({
        "metric": "te_positive_eccDNA", "value": te_positive,
    })
    summary_records.append({
        "metric": "te_positive_pct", "value": round(te_positive / total_ecc * 100, 2),
    })
    summary_records.append({
        "metric": "mean_te_pct", "value": round(per_ecc["te_pct"].mean(), 2),
    })
    summary_records.append({
        "metric": "median_te_pct", "value": round(per_ecc["te_pct"].median(), 2),
    })

    # TE class summary across all eccDNA
    class_summary = overlaps.groupby("te_class")["overlap_bp"].sum()
    total_te_bp = class_summary.sum()
    for cls, bp in class_summary.items():
        summary_records.append({
            "metric": f"te_class_{cls}_bp", "value": bp,
        })
        summary_records.append({
            "metric": f"te_class_{cls}_pct",
            "value": round(bp / total_te_bp * 100, 2) if total_te_bp > 0 else 0,
        })

    # By eccDNA type if available
    if "eccDNA_type" in per_ecc.columns and per_ecc["eccDNA_type"].notna().any():
        types = per_ecc["eccDNA_type"].dropna().unique()
        types = [t for t in types if t != ""]
        if len(types) > 1:
            for etype in sorted(types):
                subset = per_ecc[per_ecc["eccDNA_type"] == etype]
                summary_records.append({
                    "metric": f"type_{etype}_count", "value": len(subset),
                })
                summary_records.append({
                    "metric": f"type_{etype}_mean_te_pct",
                    "value": round(subset["te_pct"].mean(), 2),
                })

    summary_df = pd.DataFrame(summary_records)
    summary_out = os.path.join(output_dir, "te_summary.csv")
    summary_df.to_csv(summary_out, index=False)
    logger.info(f"Saved TE summary to {summary_out}")

    logger.info(
        f"TE analysis complete: {te_positive}/{total_ecc} eccDNA "
        f"({te_positive / total_ecc * 100:.1f}%) overlap with TEs, "
        f"mean TE coverage {per_ecc['te_pct'].mean():.1f}%"
    )
