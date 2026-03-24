"""Multi-tool detection overlap comparison for eccDNA."""

from __future__ import annotations

import json
import logging
import os
import re
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


def _detect_format(filepath: str) -> str:
    """Detect input file format based on extension and content."""
    ext = Path(filepath).suffix.lower()
    if ext == ".bed":
        return "bed"
    if ext in (".csv",):
        return "csv"
    if ext in (".tsv", ".txt"):
        return "tsv"
    return "csv"


def _load_single_intervals(
    filepath: str, fmt: str
) -> List[Tuple[str, int, int]]:
    """Load single-segment intervals from various formats."""
    intervals = []

    if fmt == "bed":
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                fields = line.split("\t")
                if len(fields) < 3:
                    continue
                try:
                    chrom = fields[0]
                    if not chrom.startswith("chr"):
                        # Try UCSC table format (bin, chrom, start, end)
                        if len(fields) >= 4 and fields[1].startswith("chr"):
                            chrom = fields[1]
                            intervals.append((chrom, int(fields[2]), int(fields[3])))
                            continue
                        continue
                    intervals.append((chrom, int(fields[1]), int(fields[2])))
                except ValueError:
                    continue
    elif fmt == "tsv":
        df = pd.read_csv(filepath, sep="\t")
        chr_col = next(
            (c for c in ["chr", "chrom", "Chr", "eChr"] if c in df.columns), None
        )
        start_col = next(
            (c for c in ["start", "chromStart", "Start", "eStart"] if c in df.columns),
            None,
        )
        end_col = next(
            (c for c in ["end", "chromEnd", "End", "eEnd"] if c in df.columns), None
        )
        if chr_col and start_col and end_col:
            for _, row in df.iterrows():
                try:
                    intervals.append(
                        (str(row[chr_col]), int(row[start_col]), int(row[end_col]))
                    )
                except (ValueError, TypeError):
                    continue
    else:
        df = pd.read_csv(filepath)
        chr_col = next(
            (c for c in ["eChr", "chr", "chrom", "Chr"] if c in df.columns), None
        )
        start_col = next(
            (c for c in ["eStart", "start", "chromStart", "Start"] if c in df.columns),
            None,
        )
        end_col = next(
            (c for c in ["eEnd", "end", "chromEnd", "End"] if c in df.columns), None
        )
        if chr_col and start_col and end_col:
            for _, row in df.iterrows():
                try:
                    intervals.append(
                        (str(row[chr_col]), int(row[start_col]), int(row[end_col]))
                    )
                except (ValueError, TypeError):
                    continue

    return intervals


def _load_chimeric_segments(
    filepath: str, fmt: str
) -> List[List[Tuple[str, int, int]]]:
    """Load chimeric (multi-segment) entries.

    Each entry is a list of (chrom, start, end) segments.
    Attempts to parse segment information from location/region columns.
    """
    entries = []
    pat = re.compile(r"(chr[\dXYMmt]+):(\d+)-(\d+)")

    if fmt == "csv":
        df = pd.read_csv(filepath)
    elif fmt in ("tsv", "txt"):
        df = pd.read_csv(filepath, sep="\t")
    else:
        return entries

    # Look for a location-like column with segment info
    loc_col = next(
        (c for c in ["location", "merge_region", "segments", "region"] if c in df.columns),
        None,
    )
    if loc_col is None:
        return entries

    for _, row in df.iterrows():
        loc_str = str(row.get(loc_col, ""))
        segs = [
            (m.group(1), int(m.group(2)), int(m.group(3)))
            for m in pat.finditer(loc_str)
        ]
        if len(segs) >= 2:
            entries.append(segs)

    return entries


class _IntervalIndex:
    """Chromosome-partitioned sorted interval index with binary search."""

    def __init__(self, intervals: List[Tuple[str, int, int]]):
        self.data: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        for chrom, start, end in intervals:
            self.data[chrom].append((start, end))
        self.starts: Dict[str, List[int]] = {}
        for chrom in self.data:
            self.data[chrom].sort()
            self.starts[chrom] = [iv[0] for iv in self.data[chrom]]

    def has_reciprocal_overlap(
        self, chrom: str, qs: int, qe: int, threshold: float
    ) -> bool:
        """Check if any indexed interval has reciprocal overlap >= threshold."""
        if chrom not in self.data:
            return False
        starts = self.starts[chrom]
        ivs = self.data[chrom]
        qlen = qe - qs
        if qlen <= 0:
            return False

        slack = (1.0 - threshold) * qlen
        lo = bisect_left(starts, qs - qlen)
        hi = bisect_right(starts, qs + int(slack) + 1)

        for i in range(lo, min(hi, len(ivs))):
            ts, te = ivs[i]
            tlen = te - ts
            if tlen <= 0:
                continue
            ov = min(qe, te) - max(qs, ts)
            if ov > 0 and ov / qlen >= threshold and ov / tlen >= threshold:
                return True
        return False


def _reciprocal_overlap(
    s1: int, e1: int, s2: int, e2: int, threshold: float
) -> bool:
    """Check reciprocal overlap between two intervals."""
    ov = min(e1, e2) - max(s1, s2)
    if ov <= 0:
        return False
    l1, l2 = e1 - s1, e2 - s2
    return l1 > 0 and l2 > 0 and ov / l1 >= threshold and ov / l2 >= threshold


def _segments_match(
    segs_a: List[Tuple[str, int, int]],
    segs_b: List[Tuple[str, int, int]],
    threshold: float,
) -> bool:
    """Check if two chimeric entries match via 1-to-1 segment matching."""
    if len(segs_a) != len(segs_b):
        return False
    used = set()
    for ca, sa, ea in segs_a:
        found = False
        for j, (cb, sb, eb) in enumerate(segs_b):
            if j not in used and ca == cb and _reciprocal_overlap(sa, ea, sb, eb, threshold):
                used.add(j)
                found = True
                break
        if not found:
            return False
    return True


def _compute_pairwise_single(
    tools: Dict[str, List[Tuple[str, int, int]]],
    threshold: float,
) -> pd.DataFrame:
    """Compute pairwise overlap matrix for single-segment intervals."""
    tool_names = sorted(tools.keys())
    indices = {name: _IntervalIndex(ivs) for name, ivs in tools.items()}

    rows = []
    for a_name in tool_names:
        for b_name in tool_names:
            if a_name == b_name:
                continue
            a_intervals = tools[a_name]
            b_idx = indices[b_name]
            count = sum(
                1
                for chrom, start, end in a_intervals
                if b_idx.has_reciprocal_overlap(chrom, start, end, threshold)
            )
            n_a = len(a_intervals)
            rows.append({
                "tool_A": a_name,
                "tool_B": b_name,
                "n_A": n_a,
                "n_B": len(tools[b_name]),
                "overlap_count": count,
                "overlap_fraction": round(count / n_a, 6) if n_a > 0 else 0.0,
            })

    return pd.DataFrame(rows)


def _compute_pairwise_chimeric(
    tools: Dict[str, List[List[Tuple[str, int, int]]]],
    threshold: float,
) -> pd.DataFrame:
    """Compute pairwise overlap for chimeric (multi-segment) entries."""
    tool_names = sorted(tools.keys())
    rows = []

    for a_name in tool_names:
        for b_name in tool_names:
            if a_name == b_name:
                continue
            a_entries = tools[a_name]
            b_entries = tools[b_name]

            # Group by segment count for efficiency
            b_by_n: Dict[int, List[Tuple[int, List[Tuple[str, int, int]]]]] = defaultdict(list)
            for idx, segs in enumerate(b_entries):
                b_by_n[len(segs)].append((idx, segs))

            matched_a = 0
            matched_b = set()

            for segs_a in a_entries:
                n_segs = len(segs_a)
                found = False
                for j, segs_b in b_by_n.get(n_segs, []):
                    if j in matched_b:
                        continue
                    if _segments_match(segs_a, segs_b, threshold):
                        matched_a += 1
                        matched_b.add(j)
                        found = True
                        break

            n_a = len(a_entries)
            rows.append({
                "tool_A": a_name,
                "tool_B": b_name,
                "n_A": n_a,
                "n_B": len(b_entries),
                "overlap_count": matched_a,
                "overlap_fraction": round(matched_a / n_a, 6) if n_a > 0 else 0.0,
            })

    return pd.DataFrame(rows)


def _classify_nway(
    tools: Dict[str, List[Tuple[str, int, int]]],
    threshold: float,
) -> pd.DataFrame:
    """Classify each interval by which tools detect it (N-way classification).

    For each tool's intervals, checks which other tools also detect the same
    region (using reciprocal overlap). This produces a detailed classification
    suitable for Venn/UpSet diagram visualization.

    Returns:
        DataFrame with columns: tool, chrom, start, end, matched_tools, n_matched.
        matched_tools is a semicolon-separated list of other tools that match.
    """
    tool_names = sorted(tools.keys())
    indices = {name: _IntervalIndex(ivs) for name, ivs in tools.items()}

    rows = []
    for src_name in tool_names:
        other_names = [t for t in tool_names if t != src_name]
        for chrom, start, end in tools[src_name]:
            matched = [
                t for t in other_names
                if indices[t].has_reciprocal_overlap(chrom, start, end, threshold)
            ]
            rows.append({
                "tool": src_name,
                "chrom": chrom,
                "start": start,
                "end": end,
                "matched_tools": ";".join(matched) if matched else "",
                "n_matched": len(matched),
            })

    return pd.DataFrame(rows)


def _summarize_nway(df_classified: pd.DataFrame, tool_names: List[str]) -> pd.DataFrame:
    """Summarize N-way classification into category counts per tool.

    For each source tool, counts intervals that are unique to it, shared with
    exactly one other tool, shared with two, etc. Also produces combination labels.

    Returns:
        DataFrame with columns: tool, category, count.
    """
    rows = []
    for tool in tool_names:
        tool_df = df_classified[df_classified["tool"] == tool]
        # Group by matched_tools combination
        grouped = tool_df.groupby("matched_tools").size().reset_index(name="count")
        for _, row in grouped.iterrows():
            matched_str = row["matched_tools"]
            if matched_str == "":
                category = f"{tool}_only"
            else:
                category = f"{tool}+{matched_str.replace(';', '+')}"
            rows.append({
                "tool": tool,
                "category": category,
                "count": row["count"],
            })
    return pd.DataFrame(rows)


def compare_detection_tools(
    tool_files: Dict[str, str],
    output_dir: str,
    min_reciprocal: float = 0.9,
    input_formats: Optional[Dict[str, str]] = None,
    classify_nway: bool = True,
) -> None:
    """Compare eccDNA detection results across multiple tools.

    Computes pairwise reciprocal overlap between all tool pairs for both
    single-segment and chimeric (multi-segment) entries. When classify_nway
    is True, also classifies each interval by which combination of tools
    detects it (for Venn/UpSet diagrams).

    Args:
        tool_files: Mapping of tool name to file path.
        output_dir: Output directory for results.
        min_reciprocal: Reciprocal overlap threshold (default: 0.9 = 90%).
        input_formats: Optional per-tool format override (default: auto-detect).
        classify_nway: Produce N-way classification per interval (default: True).
    """
    os.makedirs(output_dir, exist_ok=True)

    if input_formats is None:
        input_formats = {}

    # Load intervals for each tool
    single_tools: Dict[str, List[Tuple[str, int, int]]] = {}
    chimeric_tools: Dict[str, List[List[Tuple[str, int, int]]]] = {}

    for tool_name, filepath in tool_files.items():
        fmt = input_formats.get(tool_name, _detect_format(filepath))
        logger.info(f"Loading {tool_name} from {filepath} (format: {fmt})...")

        single = _load_single_intervals(filepath, fmt)
        chimeric = _load_chimeric_segments(filepath, fmt)

        single_tools[tool_name] = single
        chimeric_tools[tool_name] = chimeric

        logger.info(
            f"  {tool_name}: {len(single)} single-segment, "
            f"{len(chimeric)} chimeric entries"
        )

    # Compute single-segment pairwise overlap
    logger.info("Computing single-segment pairwise overlap...")
    df_single = _compute_pairwise_single(single_tools, min_reciprocal)
    single_out = os.path.join(output_dir, "single_segment_overlap.csv")
    df_single.to_csv(single_out, index=False)
    logger.info(f"Saved single-segment overlap to {single_out}")

    # N-way classification for single-segment intervals
    if classify_nway and len(single_tools) >= 2:
        logger.info("Computing N-way interval classification...")
        df_classified = _classify_nway(single_tools, min_reciprocal)
        classified_out = os.path.join(output_dir, "nway_classification.csv")
        df_classified.to_csv(classified_out, index=False)
        logger.info(f"Saved N-way classification to {classified_out}")

        tool_names = sorted(single_tools.keys())
        df_nway_summary = _summarize_nway(df_classified, tool_names)
        nway_summary_out = os.path.join(output_dir, "nway_summary.csv")
        df_nway_summary.to_csv(nway_summary_out, index=False)
        logger.info(f"Saved N-way summary to {nway_summary_out}")

    # Compute chimeric pairwise overlap (only if any tool has chimeric entries)
    has_chimeric = any(len(v) > 0 for v in chimeric_tools.values())
    if has_chimeric:
        logger.info("Computing chimeric pairwise overlap...")
        df_chimeric = _compute_pairwise_chimeric(chimeric_tools, min_reciprocal)
        chimeric_out = os.path.join(output_dir, "chimeric_overlap.csv")
        df_chimeric.to_csv(chimeric_out, index=False)
        logger.info(f"Saved chimeric overlap to {chimeric_out}")

    # Save summary
    summary = {
        "min_reciprocal": min_reciprocal,
        "tools": {
            name: {
                "file": filepath,
                "n_single": len(single_tools[name]),
                "n_chimeric": len(chimeric_tools[name]),
            }
            for name, filepath in tool_files.items()
        },
    }
    summary_out = os.path.join(output_dir, "overlap_summary.json")
    with open(summary_out, "w") as f:
        json.dump(summary, f, indent=2)
    logger.info(f"Saved summary to {summary_out}")
