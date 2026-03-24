"""Genomic linkage visualization for Circos-style plots."""

import logging
import os
import random
from itertools import combinations
from typing import List, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


def _random_rgb(rng: random.Random) -> Tuple[int, int, int]:
    """Generate a random RGB color tuple."""
    return (rng.randint(50, 230), rng.randint(50, 230), rng.randint(50, 230))


def _format_rgb(rgb: Tuple[int, int, int]) -> str:
    """Format RGB tuple as Circos color string."""
    return f"{rgb[0]},{rgb[1]},{rgb[2]}"


def generate_genolink(
    input_file: str,
    output_file: str,
    chromosomes: Optional[List[str]] = None,
    color_mode: str = "random",
) -> None:
    """
    Generate genomic link data for Circos-style visualization.

    Groups by eName and creates pairwise links between eccDNA positions
    with RGB colors for visualization.

    Args:
        input_file: eccDNA CSV file with eName column
        output_file: Output TSV file in Circos link format
        chromosomes: List of chromosomes to include (default: all)
        color_mode: Color assignment mode - "random" or "highlight-max"
    """
    logger.info(f"Loading eccDNA data from {input_file}")
    sep = "\t" if input_file.endswith((".tsv", ".bed", ".txt")) else ","
    df = pd.read_csv(input_file, sep=sep)

    # Normalize column names
    col_map = {}
    for c in df.columns:
        cl = c.lower().strip()
        if cl in ("chr", "chrom", "chromosome", "#chrom", "#chr"):
            col_map[c] = "chr"
        elif cl in ("start", "estart"):
            col_map[c] = "start"
        elif cl in ("end", "eend"):
            col_map[c] = "end"
        elif cl in ("ename", "eccdna_name", "eccdna_id", "name", "id"):
            col_map[c] = "eName"
    df = df.rename(columns=col_map)

    if "eName" not in df.columns:
        raise ValueError(
            f"Cannot find eName column in {input_file}. "
            f"Available columns: {list(df.columns)}"
        )

    for required in ("chr", "start", "end"):
        if required not in df.columns:
            raise ValueError(f"Missing required column: {required}")

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chr", "start", "end", "eName"])

    # Filter by chromosomes if specified
    if chromosomes:
        df = df[df["chr"].isin(chromosomes)]
        logger.info(f"Filtered to {len(df)} records on chromosomes: {chromosomes}")

    # Group by eName and generate pairwise links
    grouped = df.groupby("eName")
    rng = random.Random(42)

    # Pre-compute group sizes for highlight-max mode
    group_sizes = {}
    links_by_group = {}
    for ename, group in grouped:
        if len(group) < 2:
            continue
        entries = group[["chr", "start", "end"]].values.tolist()
        pairs = list(combinations(range(len(entries)), 2))
        links_by_group[ename] = [(entries[i], entries[j]) for i, j in pairs]
        group_sizes[ename] = len(pairs)

    if not links_by_group:
        logger.warning("No multi-segment groups found; no links to generate")
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        with open(output_file, "w") as fh:
            pass  # Empty file
        return

    logger.info(f"Found {len(links_by_group)} groups with >=2 segments")
    total_links = sum(group_sizes.values())
    logger.info(f"Total pairwise links: {total_links}")

    # Assign colors
    if color_mode == "highlight-max":
        max_group = max(group_sizes, key=group_sizes.get)
        highlight_color = (255, 0, 0)  # Red for max group
        gray_color = (180, 180, 180)
        colors = {}
        for ename in links_by_group:
            colors[ename] = highlight_color if ename == max_group else gray_color
        logger.info(f"Highlight-max: group '{max_group}' with {group_sizes[max_group]} links")
    elif color_mode == "random":
        colors = {ename: _random_rgb(rng) for ename in links_by_group}
    else:
        raise ValueError(f"Unknown color_mode: {color_mode}. Use 'random' or 'highlight-max'.")

    # Generate Circos link format output
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    lines = []
    for ename, link_pairs in links_by_group.items():
        rgb = colors[ename]
        color_str = f"color={_format_rgb(rgb)}"
        for entry1, entry2 in link_pairs:
            chr1, s1, e1 = entry1
            chr2, s2, e2 = entry2
            lines.append(
                f"{chr1}\t{int(s1)}\t{int(e1)}\t{chr2}\t{int(s2)}\t{int(e2)}\t{color_str}"
            )

    with open(output_file, "w") as fh:
        fh.write("\n".join(lines))
        if lines:
            fh.write("\n")

    logger.info(f"Saved {len(lines)} links to {output_file}")
