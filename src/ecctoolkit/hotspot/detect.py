"""Hotspot detection with replicate and group classification.

Integrates matrix counting and permutation testing to identify significant
eccDNA hotspots. When multiple samples are provided, classifies hotspots by:

  - Replicate reproducibility (within a group):
      core:      significant in ALL replicates
      recurrent: significant in >= 2 replicates but not all
      sporadic:  significant in only 1 replicate

  - Sample specificity (across groups):
      shared:         hotspot in >= 2 groups
      group_specific: hotspot in only 1 group
"""

import os
import logging
from typing import List, Optional, Dict, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _run_single_sample_detection(
    input_file: str,
    fai_file: str,
    output_dir: str,
    window_size: int,
    n_perm: int,
    n_cores: int,
    fdr_threshold: float,
    fold_threshold: float,
    exclude_files: Optional[List[str]],
    chrom_col: Optional[str],
    start_col: Optional[str],
    end_col: Optional[str],
    chrom_filter: Optional[List[str]],
    seed: int,
) -> pd.DataFrame:
    """Run permutation test for a single sample file.

    Returns DataFrame with columns: chrom, start, end, observed_count,
    empirical_pvalue, q_value, is_hotspot.
    """
    from ecctoolkit.hotspot.permtest import run_hotspot_permtest

    out_path = run_hotspot_permtest(
        input_file=input_file,
        fai_file=fai_file,
        output_dir=output_dir,
        window_size=window_size,
        n_perm=n_perm,
        n_cores=n_cores,
        fdr_threshold=fdr_threshold,
        fold_threshold=fold_threshold,
        exclude_files=exclude_files,
        chrom_col=chrom_col,
        start_col=start_col,
        end_col=end_col,
        chrom_filter=chrom_filter,
        seed=seed,
    )
    return pd.read_csv(out_path)


def run_sharp_detection(
    input_files: List[str],
    fai_file: str,
    output_dir: str,
    sample_names: Optional[List[str]] = None,
    group_labels: Optional[List[str]] = None,
    window_size: int = 100_000,
    n_perm: int = 1000,
    n_cores: int = 4,
    fdr_threshold: float = 0.05,
    fold_threshold: float = 3.0,
    exclude_files: Optional[List[str]] = None,
    chrom_col: Optional[str] = None,
    start_col: Optional[str] = None,
    end_col: Optional[str] = None,
    chrom_filter: Optional[List[str]] = None,
    seed: int = 42,
) -> str:
    """Detect hotspots with replicate and group classification.

    For a single sample, runs permutation test and outputs significant windows.
    For multiple samples in the same group, classifies by replicate reproducibility.
    For multiple groups, additionally classifies by sample specificity.

    Args:
        input_files: List of eccDNA CSV/BED files (one per sample).
            For legacy single-file usage, pass a list with one element.
        fai_file: Reference genome .fai index file.
        output_dir: Output directory.
        sample_names: Names for each sample (defaults to filenames).
        group_labels: Group label per sample (e.g., ["HeLa", "HeLa", "HeLa",
            "U87MG", "U87MG", "U87MG"]). If None, all samples are in one group.
        window_size: Window size in bp.
        n_perm: Number of permutations.
        n_cores: CPU cores for parallel permutation.
        fdr_threshold: FDR threshold for significance.
        fold_threshold: Minimum fold-above-median for hotspot calling.
        exclude_files: Exclusion region files (gap/centromere BED).
        chrom_col: Column name for chromosome.
        start_col: Column name for start position.
        end_col: Column name for end position.
        chrom_filter: Chromosomes to include.
        seed: Random seed.

    Returns:
        Path to the main output CSV (hotspot_results.csv).
    """
    os.makedirs(output_dir, exist_ok=True)

    # Handle legacy single-file interface
    if isinstance(input_files, str):
        input_files = [input_files]

    n_samples = len(input_files)

    # Default sample names from filenames
    if sample_names is None:
        sample_names = [
            os.path.splitext(os.path.basename(f))[0] for f in input_files
        ]

    # Default: all samples in one group
    if group_labels is None:
        group_labels = ["default"] * n_samples

    if len(sample_names) != n_samples or len(group_labels) != n_samples:
        raise ValueError(
            "sample_names and group_labels must have same length as input_files"
        )

    # Build group structure: group_name -> [sample_indices]
    groups: Dict[str, List[int]] = {}
    for i, gl in enumerate(group_labels):
        groups.setdefault(gl, []).append(i)

    ws_label = (f"{window_size // 1000}kb" if window_size < 1_000_000
                else f"{window_size // 1_000_000}Mb")

    # Run permutation test per sample
    logger.info(f"Running hotspot detection for {n_samples} sample(s)...")
    sample_results: Dict[str, pd.DataFrame] = {}
    sample_hotsets: Dict[str, set] = {}

    for i, (fpath, sname) in enumerate(zip(input_files, sample_names)):
        logger.info(f"  [{i+1}/{n_samples}] {sname}: {fpath}")
        sample_outdir = os.path.join(output_dir, "per_sample", sname)

        df_res = _run_single_sample_detection(
            input_file=fpath,
            fai_file=fai_file,
            output_dir=sample_outdir,
            window_size=window_size,
            n_perm=n_perm,
            n_cores=n_cores,
            fdr_threshold=fdr_threshold,
            fold_threshold=fold_threshold,
            exclude_files=exclude_files,
            chrom_col=chrom_col,
            start_col=start_col,
            end_col=end_col,
            chrom_filter=chrom_filter,
            seed=seed + i,
        )
        sample_results[sname] = df_res

        # Build set of hotspot window keys
        hot_df = df_res[df_res["is_hotspot"]]
        hotset = set(zip(hot_df["chrom"], hot_df["start"], hot_df["end"]))
        sample_hotsets[sname] = hotset
        logger.info(f"    {sname}: {len(hotset)} hotspot windows")

    # Single sample - just output the permtest result
    if n_samples == 1:
        sname = sample_names[0]
        df_out = sample_results[sname].copy()
        df_out.insert(0, "sample", sname)
        out_path = os.path.join(output_dir, "hotspot_results.csv")
        df_out.to_csv(out_path, index=False)
        logger.info(f"Single sample result saved: {out_path}")
        return out_path

    # Multi-sample: classify by replicate reproducibility within each group
    all_records = []

    # Collect core hotspots per group for cross-group specificity analysis
    group_core_sets: Dict[str, set] = {}

    for group_name, member_indices in groups.items():
        member_names = [sample_names[i] for i in member_indices]
        member_hotsets = [sample_hotsets[s] for s in member_names]
        n_reps = len(member_names)

        # Union of all hotspot windows in this group
        union_hot = set()
        for hs in member_hotsets:
            union_hot |= hs

        # Core = intersection of all replicates
        core = set.intersection(*member_hotsets) if member_hotsets else set()
        group_core_sets[group_name] = core

        logger.info(
            f"  Group '{group_name}': {n_reps} replicates, "
            f"union={len(union_hot)}, core={len(core)}"
        )

        # Classify each window in the union
        for w in union_hot:
            chrom, wstart, wend = w
            n_in = sum(1 for hs in member_hotsets if w in hs)

            if n_in == n_reps:
                rep_category = "core"
            elif n_in >= 2:
                rep_category = "recurrent"
            else:
                rep_category = "sporadic"

            record = {
                "chrom": chrom,
                "start": wstart,
                "end": wend,
                "group": group_name,
                "replicate_category": rep_category,
                "n_replicates_hit": n_in,
                "n_replicates_total": n_reps,
            }

            # Add per-replicate observed counts
            for si, sname in enumerate(member_names):
                df_s = sample_results[sname]
                match = df_s[(df_s["chrom"] == chrom) &
                             (df_s["start"] == wstart)]
                obs = int(match["observed_count"].iloc[0]) if len(match) > 0 else 0
                record[f"rep{si+1}_count"] = obs

            all_records.append(record)

    # Cross-group specificity (only if multiple groups)
    unique_groups = list(groups.keys())
    if len(unique_groups) > 1:
        # For each record, check how many other groups also have it as core
        for record in all_records:
            w = (record["chrom"], record["start"], record["end"])
            current_group = record["group"]

            other_groups_with_core = [
                g for g in unique_groups
                if g != current_group and w in group_core_sets.get(g, set())
            ]

            if other_groups_with_core:
                record["specificity"] = "shared"
                record["shared_with"] = ";".join(other_groups_with_core)
            else:
                record["specificity"] = "group_specific"
                record["shared_with"] = ""
    else:
        for record in all_records:
            record["specificity"] = ""
            record["shared_with"] = ""

    # Sort and save
    df_out = pd.DataFrame(all_records)
    if not df_out.empty:
        df_out = df_out.sort_values(["group", "chrom", "start"]).reset_index(drop=True)

    out_path = os.path.join(output_dir, "hotspot_results.csv")
    df_out.to_csv(out_path, index=False)
    logger.info(f"Hotspot results saved: {out_path} ({len(df_out)} rows)")

    # Summary
    summary_records = []
    for group_name in unique_groups:
        gdf = df_out[df_out["group"] == group_name]
        summary_records.append({
            "group": group_name,
            "total_hotspots": len(gdf),
            "core": (gdf["replicate_category"] == "core").sum(),
            "recurrent": (gdf["replicate_category"] == "recurrent").sum(),
            "sporadic": (gdf["replicate_category"] == "sporadic").sum(),
        })
    df_summary = pd.DataFrame(summary_records)
    summary_path = os.path.join(output_dir, "hotspot_summary.csv")
    df_summary.to_csv(summary_path, index=False)
    logger.info(f"Summary saved: {summary_path}")

    return out_path
