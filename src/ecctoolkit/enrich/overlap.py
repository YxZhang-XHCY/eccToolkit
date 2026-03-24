"""General permutation test for genomic region overlap enrichment."""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _load_query_as_bed(input_file: str) -> pd.DataFrame:
    """Load query regions from CSV or BED and return a 3-column BED DataFrame."""
    ext = Path(input_file).suffix.lower()

    if ext == ".bed":
        df = pd.read_csv(input_file, sep="\t", header=None, comment="#")
        df = df.iloc[:, :3]
        df.columns = ["chrom", "start", "end"]
    else:
        df = pd.read_csv(input_file)
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
        if chr_col is None or start_col is None or end_col is None:
            raise ValueError(
                f"Cannot find chr/start/end columns in {input_file}. "
                f"Available: {list(df.columns)}"
            )
        df = df[[chr_col, start_col, end_col]].copy()
        df.columns = ["chrom", "start", "end"]

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna().astype({"start": int, "end": int})
    return df


def _write_temp_bed(df: pd.DataFrame, tmpdir: str, name: str) -> str:
    """Write DataFrame to a temporary BED file."""
    path = os.path.join(tmpdir, f"{name}.bed")
    df.to_csv(path, sep="\t", index=False, header=False)
    return path


def _count_overlaps_bedtools(query_bed: str, annotation_bed: str) -> int:
    """Count query regions overlapping at least one annotation region."""
    result = subprocess.run(
        ["bedtools", "intersect", "-a", query_bed, "-b", annotation_bed, "-u", "-wa"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        logger.warning(f"bedtools intersect failed: {result.stderr.strip()}")
        return 0
    return len(result.stdout.strip().split("\n")) if result.stdout.strip() else 0


def _shuffle_and_count(
    query_bed: str,
    annotation_bed: str,
    genome_file: str,
    exclusion_file: Optional[str],
    seed: int,
) -> int:
    """Shuffle query regions and count overlaps with annotation."""
    shuffle_cmd = [
        "bedtools", "shuffle",
        "-i", query_bed,
        "-g", genome_file,
        "-seed", str(seed),
        "-noOverlapping",
    ]
    if exclusion_file:
        shuffle_cmd.extend(["-excl", exclusion_file])

    shuffle_result = subprocess.run(
        shuffle_cmd, capture_output=True, text=True
    )
    if shuffle_result.returncode != 0:
        # Fall back without -noOverlapping if it fails
        shuffle_cmd = [c for c in shuffle_cmd if c != "-noOverlapping"]
        shuffle_result = subprocess.run(
            shuffle_cmd, capture_output=True, text=True
        )
        if shuffle_result.returncode != 0:
            return 0

    intersect_cmd = [
        "bedtools", "intersect",
        "-a", "stdin", "-b", annotation_bed,
        "-u", "-wa",
    ]
    intersect_result = subprocess.run(
        intersect_cmd,
        input=shuffle_result.stdout,
        capture_output=True,
        text=True,
    )
    if intersect_result.returncode != 0:
        return 0

    return (
        len(intersect_result.stdout.strip().split("\n"))
        if intersect_result.stdout.strip()
        else 0
    )


def _run_one_annotation(
    query_bed: str,
    annotation_file: str,
    genome_file: str,
    n_shuffles: int,
    exclusion_file: Optional[str],
) -> dict:
    """Run permutation test for one annotation file."""
    ann_name = Path(annotation_file).stem

    # Observed overlap
    observed = _count_overlaps_bedtools(query_bed, annotation_file)

    # Permuted overlaps
    permuted = []
    for i in range(n_shuffles):
        cnt = _shuffle_and_count(
            query_bed, annotation_file, genome_file, exclusion_file, seed=i
        )
        permuted.append(cnt)

    permuted = np.array(permuted, dtype=float)
    mean_perm = np.mean(permuted)
    fold_enrichment = observed / mean_perm if mean_perm > 0 else float("inf")
    p_value = (np.sum(permuted >= observed) + 1) / (n_shuffles + 1)

    return {
        "annotation": ann_name,
        "observed": observed,
        "mean_permuted": round(mean_perm, 2),
        "sd_permuted": round(np.std(permuted), 2),
        "fold_enrichment": round(fold_enrichment, 4),
        "p_value": round(p_value, 6),
        "n_permutations": n_shuffles,
    }


def run_overlap_test(
    input_file: str,
    annotation_files: List[str],
    output_dir: str,
    genome_file: str,
    n_shuffles: int = 1000,
    n_cores: int = 8,
    exclusion_file: Optional[str] = None,
) -> None:
    """Run permutation test for genomic region overlap enrichment.

    Tests whether query regions (eccDNA) overlap with each annotation
    more than expected by chance using bedtools shuffle.

    Args:
        input_file: Query regions BED/CSV file.
        annotation_files: List of annotation BED files.
        output_dir: Output directory.
        genome_file: Genome sizes file (chrom<tab>size).
        n_shuffles: Number of shuffle permutations.
        n_cores: Number of CPU cores for parallel execution.
        exclusion_file: Optional BED file of regions to exclude from shuffling.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Load and prepare query regions
    query_df = _load_query_as_bed(input_file)
    logger.info(f"Loaded {len(query_df)} query regions from {input_file}")

    with tempfile.TemporaryDirectory() as tmpdir:
        query_bed = _write_temp_bed(query_df, tmpdir, "query")

        # Run permutation test for each annotation
        results = []
        if n_cores > 1 and len(annotation_files) > 1:
            with ProcessPoolExecutor(max_workers=min(n_cores, len(annotation_files))) as executor:
                futures = {
                    executor.submit(
                        _run_one_annotation,
                        query_bed,
                        ann_file,
                        genome_file,
                        n_shuffles,
                        exclusion_file,
                    ): ann_file
                    for ann_file in annotation_files
                }
                for future in as_completed(futures):
                    ann_file = futures[future]
                    try:
                        result = future.result()
                        results.append(result)
                        logger.info(
                            f"  {result['annotation']}: observed={result['observed']}, "
                            f"fold={result['fold_enrichment']:.2f}, "
                            f"p={result['p_value']:.4f}"
                        )
                    except Exception as e:
                        logger.error(f"Error processing {ann_file}: {e}")
        else:
            for ann_file in annotation_files:
                result = _run_one_annotation(
                    query_bed, ann_file, genome_file, n_shuffles, exclusion_file
                )
                results.append(result)
                logger.info(
                    f"  {result['annotation']}: observed={result['observed']}, "
                    f"fold={result['fold_enrichment']:.2f}, "
                    f"p={result['p_value']:.4f}"
                )

    # Apply BH FDR correction
    if results:
        from statsmodels.stats.multitest import multipletests

        p_values = [r["p_value"] for r in results]
        _, q_values, _, _ = multipletests(p_values, method="fdr_bh")
        for r, q in zip(results, q_values):
            r["q_value"] = round(q, 6)

    # Save results
    df_results = pd.DataFrame(results)
    if not df_results.empty:
        df_results = df_results.sort_values("p_value")
    output_file = os.path.join(output_dir, "overlap_enrichment.csv")
    df_results.to_csv(output_file, index=False)
    logger.info(f"Saved enrichment results to {output_file}")
