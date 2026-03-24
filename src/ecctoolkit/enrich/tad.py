"""TAD boundary enrichment analysis using permutation test."""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _load_eccdna_as_bed(filepaths: List[str]) -> pd.DataFrame:
    """Load eccDNA regions from CSV/BED files into BED format DataFrame."""
    all_dfs = []
    for filepath in filepaths:
        ext = Path(filepath).suffix.lower()
        if ext == ".bed":
            df = pd.read_csv(filepath, sep="\t", header=None, comment="#")
            df = df.iloc[:, :3]
            df.columns = ["chrom", "start", "end"]
        else:
            df = pd.read_csv(filepath)
            chr_col = next(
                (c for c in ["eChr", "chr", "chrom"] if c in df.columns), None
            )
            start_col = next(
                (c for c in ["eStart", "start", "chromStart"] if c in df.columns), None
            )
            end_col = next(
                (c for c in ["eEnd", "end", "chromEnd"] if c in df.columns), None
            )
            if chr_col is None or start_col is None or end_col is None:
                raise ValueError(f"Cannot find chr/start/end columns in {filepath}")
            df = df[[chr_col, start_col, end_col]].copy()
            df.columns = ["chrom", "start", "end"]

        df["start"] = pd.to_numeric(df["start"], errors="coerce")
        df["end"] = pd.to_numeric(df["end"], errors="coerce")
        df = df.dropna().astype({"start": int, "end": int})
        all_dfs.append(df)

    return pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()


def _load_tad_boundaries(tad_file: str, flank: int = 50000) -> pd.DataFrame:
    """Load TAD boundaries and create flanking regions.

    TAD boundaries are represented as the regions at TAD edges.
    A flank region around each boundary is created.

    Args:
        tad_file: TAD domains BED file (chrom, start, end).
        flank: Size of boundary flanking region in bp.

    Returns:
        DataFrame with boundary regions.
    """
    df = pd.read_csv(tad_file, sep="\t", header=None, comment="#")
    df = df.iloc[:, :3]
    df.columns = ["chrom", "start", "end"]
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna().astype({"start": int, "end": int})

    # Create boundary regions (left and right edges of TADs)
    boundaries = []
    for _, row in df.iterrows():
        chrom = row["chrom"]
        # Left boundary
        b_start = max(0, row["start"] - flank)
        b_end = row["start"] + flank
        boundaries.append({"chrom": chrom, "start": b_start, "end": b_end})
        # Right boundary
        b_start = max(0, row["end"] - flank)
        b_end = row["end"] + flank
        boundaries.append({"chrom": chrom, "start": b_start, "end": b_end})

    return pd.DataFrame(boundaries).drop_duplicates()


def _compute_min_distances(
    eccdna_df: pd.DataFrame, boundaries_df: pd.DataFrame
) -> np.ndarray:
    """Compute minimum distance from each eccDNA midpoint to nearest boundary.

    Uses bedtools closest for efficiency.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write eccDNA midpoints as BED
        midpoints = eccdna_df.copy()
        midpoints["mid"] = (midpoints["start"] + midpoints["end"]) // 2
        midpoints_bed = os.path.join(tmpdir, "midpoints.bed")
        mid_df = midpoints[["chrom", "mid"]].copy()
        mid_df["end"] = mid_df["mid"] + 1
        mid_df = mid_df[["chrom", "mid", "end"]]
        mid_df.columns = ["chrom", "start", "end"]
        mid_df = mid_df.sort_values(["chrom", "start"])
        mid_df.to_csv(midpoints_bed, sep="\t", index=False, header=False)

        # Write boundary midpoints
        bound_mid = boundaries_df.copy()
        bound_mid["mid"] = (bound_mid["start"] + bound_mid["end"]) // 2
        bound_bed = os.path.join(tmpdir, "boundaries.bed")
        b_df = bound_mid[["chrom", "mid"]].copy()
        b_df["end"] = b_df["mid"] + 1
        b_df = b_df[["chrom", "mid", "end"]]
        b_df.columns = ["chrom", "start", "end"]
        b_df = b_df.sort_values(["chrom", "start"])
        b_df.to_csv(bound_bed, sep="\t", index=False, header=False)

        result = subprocess.run(
            ["bedtools", "closest", "-a", midpoints_bed, "-b", bound_bed, "-d"],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            logger.warning(f"bedtools closest failed: {result.stderr.strip()}")
            return np.array([])

        distances = []
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            fields = line.split("\t")
            try:
                dist = int(fields[-1])
                if dist >= 0:
                    distances.append(dist)
            except (ValueError, IndexError):
                continue

    return np.array(distances)


def run_tad_enrichment(
    eccdna_files: List[str],
    tad_file: str,
    output_dir: str,
    genome_file: Optional[str] = None,
    prefix: Optional[str] = None,
    n_shuffles: int = 1000,
    n_cores: int = 8,
    keep_sex_chr: bool = False,
) -> None:
    """Analyze eccDNA enrichment at TAD boundaries.

    Computes the distance from each eccDNA to the nearest TAD boundary,
    then uses a permutation test to assess whether eccDNA are closer to
    TAD boundaries than expected by chance.

    Args:
        eccdna_files: List of eccDNA CSV/BED files.
        tad_file: TAD boundaries/domains BED file.
        output_dir: Output directory.
        genome_file: Genome sizes file.
        prefix: Output file prefix.
        n_shuffles: Number of permutations.
        n_cores: Number of CPU cores.
        keep_sex_chr: Keep sex chromosomes in analysis.
    """
    os.makedirs(output_dir, exist_ok=True)
    out_prefix = prefix or "tad_enrichment"

    # Load data
    eccdna_df = _load_eccdna_as_bed(eccdna_files)
    if not keep_sex_chr:
        eccdna_df = eccdna_df[~eccdna_df["chrom"].isin(["chrX", "chrY"])]
    logger.info(f"Loaded {len(eccdna_df)} eccDNA regions")

    boundaries_df = _load_tad_boundaries(tad_file)
    if not keep_sex_chr:
        boundaries_df = boundaries_df[~boundaries_df["chrom"].isin(["chrX", "chrY"])]
    logger.info(f"Created {len(boundaries_df)} TAD boundary regions")

    # Observed distances
    obs_distances = _compute_min_distances(eccdna_df, boundaries_df)
    if len(obs_distances) == 0:
        logger.error("No distances computed. Check input files.")
        return
    obs_median = np.median(obs_distances)
    obs_mean = np.mean(obs_distances)
    logger.info(
        f"Observed distances: median={obs_median:.0f}, mean={obs_mean:.0f}"
    )

    if genome_file is None:
        raise ValueError("genome_file is required for TAD enrichment analysis")

    # Permutation test: shuffle eccDNA and recompute distances
    with tempfile.TemporaryDirectory() as tmpdir:
        eccdna_bed = os.path.join(tmpdir, "eccdna.bed")
        eccdna_df[["chrom", "start", "end"]].to_csv(
            eccdna_bed, sep="\t", index=False, header=False
        )

        perm_medians = []
        perm_means = []

        for i in range(n_shuffles):
            shuffle_cmd = [
                "bedtools", "shuffle",
                "-i", eccdna_bed,
                "-g", genome_file,
                "-seed", str(i),
                "-noOverlapping",
            ]
            shuffle_result = subprocess.run(
                shuffle_cmd, capture_output=True, text=True
            )
            if shuffle_result.returncode != 0:
                shuffle_cmd = [c for c in shuffle_cmd if c != "-noOverlapping"]
                shuffle_result = subprocess.run(
                    shuffle_cmd, capture_output=True, text=True
                )
                if shuffle_result.returncode != 0:
                    continue

            shuffled_bed = os.path.join(tmpdir, f"shuffled_{i}.bed")
            with open(shuffled_bed, "w") as f:
                f.write(shuffle_result.stdout)

            # Load shuffled as DataFrame
            shuffled_df = pd.read_csv(
                shuffled_bed, sep="\t", header=None, names=["chrom", "start", "end"]
            )
            perm_distances = _compute_min_distances(shuffled_df, boundaries_df)

            if len(perm_distances) > 0:
                perm_medians.append(np.median(perm_distances))
                perm_means.append(np.mean(perm_distances))

            os.remove(shuffled_bed)

            if (i + 1) % 100 == 0:
                logger.info(f"  Completed {i + 1}/{n_shuffles} permutations")

    perm_medians = np.array(perm_medians)
    perm_means = np.array(perm_means)

    # P-value: fraction of permutations with median distance <= observed
    p_value_median = (np.sum(perm_medians <= obs_median) + 1) / (len(perm_medians) + 1)
    p_value_mean = (np.sum(perm_means <= obs_mean) + 1) / (len(perm_means) + 1)

    fold_median = (
        np.mean(perm_medians) / obs_median if obs_median > 0 else float("inf")
    )
    fold_mean = np.mean(perm_means) / obs_mean if obs_mean > 0 else float("inf")

    results = {
        "metric": ["median_distance", "mean_distance"],
        "observed": [obs_median, obs_mean],
        "mean_permuted": [
            round(np.mean(perm_medians), 2),
            round(np.mean(perm_means), 2),
        ],
        "sd_permuted": [
            round(np.std(perm_medians), 2),
            round(np.std(perm_means), 2),
        ],
        "fold_closer": [round(fold_median, 4), round(fold_mean, 4)],
        "p_value": [round(p_value_median, 6), round(p_value_mean, 6)],
        "n_permutations": [n_shuffles, n_shuffles],
    }

    df_results = pd.DataFrame(results)
    output_file = os.path.join(output_dir, f"{out_prefix}.csv")
    df_results.to_csv(output_file, index=False)
    logger.info(f"Saved TAD enrichment results to {output_file}")

    # Save distance distribution
    dist_file = os.path.join(output_dir, f"{out_prefix}_distances.csv")
    pd.DataFrame({"distance_to_tad_boundary": obs_distances}).to_csv(
        dist_file, index=False
    )
    logger.info(f"Saved distance distribution to {dist_file}")

    logger.info(
        f"TAD boundary enrichment: median distance {obs_median:.0f} bp "
        f"(expected {np.mean(perm_medians):.0f}), "
        f"fold={fold_median:.2f}x closer, p={p_value_median:.4f}"
    )
