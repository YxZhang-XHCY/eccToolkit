"""CNV region enrichment analysis using permutation test."""

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
    """Load eccDNA regions from CSV/BED files and combine into BED format."""
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


def _load_cnv_regions(cnv_file: str) -> pd.DataFrame:
    """Load CNV regions from BED file.

    Expected format: chrom, start, end, cnv_value [, additional columns]
    CNV categories are assigned based on cnv_value:
      - gain: cnv_value > 1.2 (above diploid)
      - loss: cnv_value < 0.8 (below diploid)
      - neutral: 0.8 <= cnv_value <= 1.2
    """
    df = pd.read_csv(cnv_file, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        raise ValueError(
            f"CNV file must have at least 4 columns (chrom, start, end, cnv_value): {cnv_file}"
        )
    df = df.iloc[:, :4]
    df.columns = ["chrom", "start", "end", "cnv_value"]
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["cnv_value"] = pd.to_numeric(df["cnv_value"], errors="coerce")
    df = df.dropna()

    # Classify CNV categories
    df["cnv_category"] = "neutral"
    df.loc[df["cnv_value"] > 1.2, "cnv_category"] = "gain"
    df.loc[df["cnv_value"] < 0.8, "cnv_category"] = "loss"

    return df


def _count_intersections_by_category(
    eccdna_bed: str, cnv_df: pd.DataFrame, tmpdir: str
) -> dict:
    """Count eccDNA regions falling in each CNV category."""
    counts = {}
    for category in ["gain", "loss", "neutral"]:
        cat_df = cnv_df[cnv_df["cnv_category"] == category]
        if cat_df.empty:
            counts[category] = 0
            continue
        cat_bed = os.path.join(tmpdir, f"cnv_{category}.bed")
        cat_df[["chrom", "start", "end"]].to_csv(
            cat_bed, sep="\t", index=False, header=False
        )
        result = subprocess.run(
            ["bedtools", "intersect", "-a", eccdna_bed, "-b", cat_bed, "-u", "-wa"],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0 or not result.stdout.strip():
            counts[category] = 0
        else:
            counts[category] = len(result.stdout.strip().split("\n"))
    return counts


def _shuffle_and_count_by_category(
    eccdna_bed: str,
    cnv_df: pd.DataFrame,
    genome_file: str,
    exclusion_file: Optional[str],
    tmpdir: str,
    seed: int,
) -> dict:
    """Shuffle eccDNA and count intersections per CNV category."""
    shuffle_cmd = [
        "bedtools", "shuffle",
        "-i", eccdna_bed,
        "-g", genome_file,
        "-seed", str(seed),
        "-noOverlapping",
    ]
    if exclusion_file:
        shuffle_cmd.extend(["-excl", exclusion_file])

    shuffle_result = subprocess.run(shuffle_cmd, capture_output=True, text=True)
    if shuffle_result.returncode != 0:
        shuffle_cmd = [c for c in shuffle_cmd if c != "-noOverlapping"]
        shuffle_result = subprocess.run(shuffle_cmd, capture_output=True, text=True)
        if shuffle_result.returncode != 0:
            return {"gain": 0, "loss": 0, "neutral": 0}

    # Write shuffled to temp file
    shuffled_bed = os.path.join(tmpdir, f"shuffled_{seed}.bed")
    with open(shuffled_bed, "w") as f:
        f.write(shuffle_result.stdout)

    counts = _count_intersections_by_category(shuffled_bed, cnv_df, tmpdir)

    # Clean up
    if os.path.exists(shuffled_bed):
        os.remove(shuffled_bed)

    return counts


def run_cnv_enrichment(
    eccdna_files: List[str],
    cnv_file: str,
    output_dir: str,
    genome_file: Optional[str] = None,
    prefix: Optional[str] = None,
    n_shuffles: int = 1000,
    n_cores: int = 8,
    keep_sex_chr: bool = False,
) -> None:
    """Analyze eccDNA enrichment in CNV gain/loss/neutral regions.

    Uses bedtools shuffle + intersect for permutation test.

    Args:
        eccdna_files: List of eccDNA CSV/BED files.
        cnv_file: CNV regions BED file (chrom, start, end, cnv_value).
        output_dir: Output directory.
        genome_file: Genome sizes file (auto-detected if not provided).
        prefix: Output file prefix.
        n_shuffles: Number of permutations.
        n_cores: Number of CPU cores.
        keep_sex_chr: Keep sex chromosomes in analysis.
    """
    os.makedirs(output_dir, exist_ok=True)
    out_prefix = prefix or "cnv_enrichment"

    # Load data
    eccdna_df = _load_eccdna_as_bed(eccdna_files)
    if not keep_sex_chr:
        eccdna_df = eccdna_df[~eccdna_df["chrom"].isin(["chrX", "chrY"])]
    logger.info(f"Loaded {len(eccdna_df)} eccDNA regions")

    cnv_df = _load_cnv_regions(cnv_file)
    if not keep_sex_chr:
        cnv_df = cnv_df[~cnv_df["chrom"].isin(["chrX", "chrY"])]
    for cat in ["gain", "loss", "neutral"]:
        n = (cnv_df["cnv_category"] == cat).sum()
        logger.info(f"  CNV {cat}: {n} regions")

    if genome_file is None:
        raise ValueError("genome_file is required for CNV enrichment analysis")

    with tempfile.TemporaryDirectory() as tmpdir:
        eccdna_bed = os.path.join(tmpdir, "eccdna.bed")
        eccdna_df[["chrom", "start", "end"]].to_csv(
            eccdna_bed, sep="\t", index=False, header=False
        )

        # Observed counts
        observed = _count_intersections_by_category(eccdna_bed, cnv_df, tmpdir)
        logger.info(f"Observed counts: {observed}")

        # Permutation test
        permuted_counts = {"gain": [], "loss": [], "neutral": []}
        for i in range(n_shuffles):
            perm = _shuffle_and_count_by_category(
                eccdna_bed, cnv_df, genome_file, None, tmpdir, seed=i
            )
            for cat in permuted_counts:
                permuted_counts[cat].append(perm[cat])
            if (i + 1) % 100 == 0:
                logger.info(f"  Completed {i + 1}/{n_shuffles} permutations")

        # Compute statistics
        results = []
        for category in ["gain", "loss", "neutral"]:
            perm_arr = np.array(permuted_counts[category], dtype=float)
            obs = observed[category]
            mean_perm = np.mean(perm_arr)
            fold = obs / mean_perm if mean_perm > 0 else float("inf")
            p_value = (np.sum(perm_arr >= obs) + 1) / (n_shuffles + 1)

            results.append({
                "cnv_category": category,
                "observed": obs,
                "mean_permuted": round(mean_perm, 2),
                "sd_permuted": round(np.std(perm_arr), 2),
                "fold_enrichment": round(fold, 4),
                "p_value": round(p_value, 6),
                "n_permutations": n_shuffles,
            })
            logger.info(
                f"  {category}: obs={obs}, expected={mean_perm:.1f}, "
                f"fold={fold:.2f}, p={p_value:.4f}"
            )

    # Save results
    df_results = pd.DataFrame(results)
    output_file = os.path.join(output_dir, f"{out_prefix}.csv")
    df_results.to_csv(output_file, index=False)
    logger.info(f"Saved CNV enrichment results to {output_file}")
