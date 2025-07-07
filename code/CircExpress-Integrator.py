#!/usr/bin/env python3
"""
ECC-DNA enrichment/depletion vs. RNA-seq DEGs – per-sample association analysis
==============================================================================
* English-only figures
* FC-stratified gradient curves with summary heatmaps
* Modular & CLI-driven via argparse
* Enhanced with gene lists output and comprehensive logging
* Dependencies: pandas, numpy, matplotlib, seaborn, scipy, statsmodels
"""

import argparse
import os
import warnings
import logging
from pathlib import Path
from datetime import datetime
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ─────────────────────────────────────────────────────────────────────────────
# 0. Logging setup
# ─────────────────────────────────────────────────────────────────────────────
def setup_logging(outdir: Path, verbose: bool = False):
    """Setup logging configuration."""
    log_level = logging.DEBUG if verbose else logging.INFO
    log_file = outdir / f"analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    
    # Configure root logger
    logger = logging.getLogger()
    logger.setLevel(log_level)
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger

# ─────────────────────────────────────────────────────────────────────────────
# I. I/O helpers
# ─────────────────────────────────────────────────────────────────────────────
def load_eccdna_data(path: str) -> pd.DataFrame:
    """Load eccDNA enrichment/depletion results (CSV)."""
    logger = logging.getLogger(__name__)
    logger.info(f"Loading eccDNA data from: {path}")
    
    df = pd.read_csv(path)
    required = {"Sample", "gene", "direction", "significant"}
    if not required.issubset(df.columns):
        missing = ", ".join(required - set(df.columns))
        raise ValueError(f"eccDNA file is missing required columns: {missing}")
    
    logger.info(f"Loaded {len(df)} eccDNA records from {df['Sample'].nunique()} samples")
    return df


def load_degs(path: str) -> pd.DataFrame:
    """Load DEG table produced by e.g. GEPIA2 (tab-separated)."""
    logger = logging.getLogger(__name__)
    logger.info(f"Loading DEG data from: {path}")
    
    cols = ["gene", "ensembl_id", "tumor_exp", "normal_exp", "log2fc", "qvalue"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols)
    
    logger.info(f"Loaded {len(df)} DEGs")
    logger.info(f"DEG range: log2FC from {df['log2fc'].min():.2f} to {df['log2fc'].max():.2f}")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# II. Pre-compute DEG sets for FC thresholds
# ─────────────────────────────────────────────────────────────────────────────
def build_deg_sets(degs: pd.DataFrame, thresholds=(1, 2, 4, 6)):
    """Return dict fc -> (up_set, down_set)."""
    logger = logging.getLogger(__name__)
    deg_sets = {}
    
    for fc in thresholds:
        up = set(degs[degs["log2fc"] > fc]["gene"])
        down = set(degs[degs["log2fc"] < -fc]["gene"])
        deg_sets[fc] = (up, down)
        logger.debug(f"FC={fc}: {len(up)} up-regulated, {len(down)} down-regulated genes")
    
    return deg_sets


# ─────────────────────────────────────────────────────────────────────────────
# III. Stats helpers
# ─────────────────────────────────────────────────────────────────────────────
def fisher_test(overlap: int, set1_size: int, set2_size: int, total: int):
    """One-tailed Fisher exact test (alternative='greater')."""
    table = [
        [overlap, set1_size - overlap],
        [set2_size - overlap, total - set1_size - set2_size + overlap],
    ]
    or_, p = fisher_exact(table, alternative="greater")
    return or_, p


# ─────────────────────────────────────────────────────────────────────────────
# IV. Overlap analysis per FC threshold (Enhanced with gene lists)
# ─────────────────────────────────────────────────────────────────────────────
def analyze_overlap(sample_ecc: pd.DataFrame, deg_sets: dict, fc: float):
    """Return stats dict for a single FC, including gene lists."""
    ecc_enriched = set(
        sample_ecc[
            (sample_ecc.direction == "enrichment") & sample_ecc.significant
        ]["gene"]
    )
    ecc_depleted = set(
        sample_ecc[
            (sample_ecc.direction == "depletion") & sample_ecc.significant
        ]["gene"]
    )
    deg_up, deg_down = deg_sets[fc]

    # Calculate overlaps
    enriched_up = ecc_enriched & deg_up
    enriched_down = ecc_enriched & deg_down
    depleted_up = ecc_depleted & deg_up
    depleted_down = ecc_depleted & deg_down

    stats = {
        "n_ecc_enriched": len(ecc_enriched),
        "n_ecc_depleted": len(ecc_depleted),
        "n_deg_up": len(deg_up),
        "n_deg_down": len(deg_down),
        "n_enriched_up": len(enriched_up),
        "n_enriched_down": len(enriched_down),
        "n_depleted_up": len(depleted_up),
        "n_depleted_down": len(depleted_down),
        # Add gene lists
        "genes_enriched_up": sorted(list(enriched_up)),
        "genes_enriched_down": sorted(list(enriched_down)),
        "genes_depleted_up": sorted(list(depleted_up)),
        "genes_depleted_down": sorted(list(depleted_down))
    }

    total_genes = len(set(sample_ecc["gene"]).union(*deg_sets[fc]))

    # Fisher – enriched∩up and depleted∩down only
    if stats["n_enriched_up"] > 0:
        or_eu, p_eu = fisher_test(
            stats["n_enriched_up"],
            stats["n_ecc_enriched"],
            stats["n_deg_up"],
            total_genes,
        )
    else:
        or_eu, p_eu = 0, 1

    if stats["n_depleted_down"] > 0:
        or_dd, p_dd = fisher_test(
            stats["n_depleted_down"],
            stats["n_ecc_depleted"],
            stats["n_deg_down"],
            total_genes,
        )
    else:
        or_dd, p_dd = 0, 1

    stats["OR_enriched_up"] = or_eu
    stats["p_enriched_up"] = p_eu
    stats["OR_depleted_down"] = or_dd
    stats["p_depleted_down"] = p_dd

    return stats


# ─────────────────────────────────────────────────────────────────────────────
# V. Gradient analysis per sample
# ─────────────────────────────────────────────────────────────────────────────
def gradient_analysis(
    sample: str,
    eccdna_df: pd.DataFrame,
    deg_sets: dict,
    fc_list=(1, 2, 4, 6),
) -> pd.DataFrame:
    logger = logging.getLogger(__name__)
    logger.info(f"Performing gradient analysis for sample: {sample}")
    
    sample_ecc = eccdna_df[eccdna_df["Sample"] == sample]
    rows = []
    
    for fc in fc_list:
        stats = analyze_overlap(sample_ecc, deg_sets, fc)
        stats["sample"] = sample
        stats["fc_threshold"] = fc
        rows.append(stats)
    
    df = pd.DataFrame(rows)

    # FDR correction for the two p-value columns (per sample)
    mask = df[["p_enriched_up", "p_depleted_down"]].values.flatten() < 1
    if mask.any():
        pvals = df[["p_enriched_up", "p_depleted_down"]].values.flatten()
        valid_idx = np.where(~np.isnan(pvals) & (pvals < 1))[0]
        if len(valid_idx):
            qvals = multipletests(pvals[valid_idx], method="fdr_bh")[1]
            # Map back
            flat_q = np.full_like(pvals, np.nan, dtype=float)
            flat_q[valid_idx] = qvals
            df["q_enriched_up"] = flat_q[0 :: 2]
            df["q_depleted_down"] = flat_q[1 :: 2]
    
    return df


# ─────────────────────────────────────────────────────────────────────────────
# VI. Plotting helpers
# ─────────────────────────────────────────────────────────────────────────────
def plot_gradient(sample: str, grad_df: pd.DataFrame, outdir: Path):
    """Original gradient curves for individual sample."""
    fc = grad_df["fc_threshold"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(
        f"{sample}: eccDNA–DEG association across FC thresholds", fontsize=14
    )

    # (A) Overlap counts
    ax = axes[0, 0]
    ax.plot(fc, grad_df["n_enriched_up"], "o-", label="Enriched ∩ Up", linewidth=2)
    ax.plot(fc, grad_df["n_depleted_down"], "s-", label="Depleted ∩ Down", linewidth=2)
    ax.set_xlabel("FC threshold")
    ax.set_ylabel("# overlapping genes")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (B) Odds ratios
    ax = axes[0, 1]
    ax.plot(fc, grad_df["OR_enriched_up"], "o-", label="Enriched-Up", linewidth=2)
    ax.plot(fc, grad_df["OR_depleted_down"], "s-", label="Depleted-Down", linewidth=2)
    ax.axhline(1, ls="--", c="k", alpha=0.5)
    ax.set_yscale("log")
    ax.set_xlabel("FC threshold")
    ax.set_ylabel("Odds ratio (log scale)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (C) –log10 p
    ax = axes[1, 0]
    # clip p to avoid log10(0)
    p_eu = np.clip(grad_df["p_enriched_up"], 1e-300, 1.0)
    p_dd = np.clip(grad_df["p_depleted_down"], 1e-300, 1.0)
    ax.plot(fc, -np.log10(p_eu), "o-", label="Enriched-Up", linewidth=2)
    ax.plot(fc, -np.log10(p_dd), "s-", label="Depleted-Down", linewidth=2)
    ax.axhline(-np.log10(0.05), ls="--", c="k", alpha=0.5, label="p=0.05")
    ax.set_xlabel("FC threshold")
    ax.set_ylabel("–log10 p-value")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (D) DEG counts
    ax = axes[1, 1]
    ax.plot(fc, grad_df["n_deg_up"], "o-", label="Up-DEGs", linewidth=2)
    ax.plot(fc, grad_df["n_deg_down"], "s-", label="Down-DEGs", linewidth=2)
    ax.set_xlabel("FC threshold")
    ax.set_ylabel("# DEGs")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outdir / f"{sample}_gradient_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()


def plot_summary_heatmap(all_grad: pd.DataFrame, outdir: Path):
    """Create summary heatmaps for all samples."""
    logger = logging.getLogger(__name__)
    logger.info("Creating summary heatmaps")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Summary: eccDNA-DEG associations across all samples", fontsize=16)
    
    # Prepare data for heatmaps
    samples = sorted(all_grad['sample'].unique())
    fc_thresholds = sorted(all_grad['fc_threshold'].unique())
    
    # (A) Odds ratio heatmap - Enriched-Up
    ax = axes[0, 0]
    or_eu_matrix = all_grad.pivot(index='sample', columns='fc_threshold', values='OR_enriched_up')
    sns.heatmap(or_eu_matrix, annot=True, fmt='.2f', cmap='RdYlBu_r', ax=ax, 
                cbar_kws={'label': 'Odds Ratio'})
    ax.set_title('Odds Ratio: Enriched ∩ Up-regulated')
    ax.set_xlabel('FC threshold')
    ax.set_ylabel('Sample')
    
    # (B) Odds ratio heatmap - Depleted-Down
    ax = axes[0, 1]
    or_dd_matrix = all_grad.pivot(index='sample', columns='fc_threshold', values='OR_depleted_down')
    sns.heatmap(or_dd_matrix, annot=True, fmt='.2f', cmap='RdYlBu_r', ax=ax,
                cbar_kws={'label': 'Odds Ratio'})
    ax.set_title('Odds Ratio: Depleted ∩ Down-regulated')
    ax.set_xlabel('FC threshold')
    ax.set_ylabel('Sample')
    
    # (C) -log10(p) heatmap - Enriched-Up
    ax = axes[1, 0]
    p_eu_matrix = all_grad.pivot(index='sample', columns='fc_threshold', values='p_enriched_up')
    p_eu_matrix = -np.log10(np.clip(p_eu_matrix, 1e-300, 1))
    sns.heatmap(p_eu_matrix, annot=True, fmt='.2f', cmap='Reds', ax=ax,
                cbar_kws={'label': '-log10(p-value)'})
    ax.set_title('-log10(p): Enriched ∩ Up-regulated')
    ax.set_xlabel('FC threshold')
    ax.set_ylabel('Sample')
    
    # (D) -log10(p) heatmap - Depleted-Down
    ax = axes[1, 1]
    p_dd_matrix = all_grad.pivot(index='sample', columns='fc_threshold', values='p_depleted_down')
    p_dd_matrix = -np.log10(np.clip(p_dd_matrix, 1e-300, 1))
    sns.heatmap(p_dd_matrix, annot=True, fmt='.2f', cmap='Reds', ax=ax,
                cbar_kws={'label': '-log10(p-value)'})
    ax.set_title('-log10(p): Depleted ∩ Down-regulated')
    ax.set_xlabel('FC threshold')
    ax.set_ylabel('Sample')
    
    plt.tight_layout()
    plt.savefig(outdir / "summary_heatmaps.png", dpi=300, bbox_inches='tight')
    plt.close()


def plot_overlap_summary(all_grad: pd.DataFrame, outdir: Path):
    """Create bar plots summarizing overlap counts across samples."""
    logger = logging.getLogger(__name__)
    logger.info("Creating overlap summary plots")
    
    # Filter for FC=2 as example
    fc2_data = all_grad[all_grad['fc_threshold'] == 2].sort_values('sample')
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Gene overlap counts at FC threshold = 2', fontsize=16)
    
    # Enriched-Up overlaps
    ax1.bar(fc2_data['sample'], fc2_data['n_enriched_up'], color='darkred', alpha=0.7)
    ax1.set_xlabel('Sample')
    ax1.set_ylabel('# genes')
    ax1.set_title('Enriched ∩ Up-regulated')
    ax1.tick_params(axis='x', rotation=45)
    
    # Depleted-Down overlaps
    ax2.bar(fc2_data['sample'], fc2_data['n_depleted_down'], color='darkblue', alpha=0.7)
    ax2.set_xlabel('Sample')
    ax2.set_ylabel('# genes')
    ax2.set_title('Depleted ∩ Down-regulated')
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(outdir / "overlap_summary_fc2.png", dpi=300, bbox_inches='tight')
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# VII. Gene list export
# ─────────────────────────────────────────────────────────────────────────────
def export_gene_lists(grad_df: pd.DataFrame, sample: str, outdir: Path):
    """Export gene lists to separate files for each FC threshold."""
    logger = logging.getLogger(__name__)
    logger.info(f"Exporting gene lists for sample: {sample}")
    
    gene_dir = outdir / "gene_lists"
    gene_dir.mkdir(exist_ok=True)
    
    for _, row in grad_df.iterrows():
        fc = row['fc_threshold']
        
        # Create a dictionary with all gene lists
        gene_lists = {
            'enriched_up': row['genes_enriched_up'],
            'enriched_down': row['genes_enriched_down'],
            'depleted_up': row['genes_depleted_up'],
            'depleted_down': row['genes_depleted_down']
        }
        
        # Save as JSON for easy reading
        json_file = gene_dir / f"{sample}_fc{fc}_gene_lists.json"
        with open(json_file, 'w') as f:
            json.dump(gene_lists, f, indent=2)
        
        # Also save key lists as plain text
        if row['genes_enriched_up']:
            txt_file = gene_dir / f"{sample}_fc{fc}_enriched_up.txt"
            with open(txt_file, 'w') as f:
                f.write('\n'.join(row['genes_enriched_up']))
        
        if row['genes_depleted_down']:
            txt_file = gene_dir / f"{sample}_fc{fc}_depleted_down.txt"
            with open(txt_file, 'w') as f:
                f.write('\n'.join(row['genes_depleted_down']))


# ─────────────────────────────────────────────────────────────────────────────
# VIII. Main workflow
# ─────────────────────────────────────────────────────────────────────────────
def run(args):
    # Setup output directory and logging
    output_root = Path(args.outdir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    
    logger = setup_logging(output_root, args.verbose)
    logger.info("="*80)
    logger.info("Starting eccDNA-DEG association analysis")
    logger.info(f"Output directory: {output_root}")
    logger.info("="*80)
    
    # Load data
    ecc = load_eccdna_data(args.eccdna)
    deg = load_degs(args.degs)

    # Parse FC thresholds
    fc_list = [float(x) for x in args.fc.split(",")]
    logger.info(f"FC thresholds: {fc_list}")
    
    # Build DEG sets
    deg_sets = build_deg_sets(deg, thresholds=fc_list)

    # Get samples
    samples = sorted(ecc["Sample"].unique())
    logger.info(f"Processing {len(samples)} samples: {', '.join(samples)}")

    summary_tables = []

    # Process each sample
    for s in samples:
        logger.info(f"Processing sample: {s}")
        sample_dir = output_root / f"sample_{s}"
        sample_dir.mkdir(exist_ok=True)

        # Gradient analysis
        grad_df = gradient_analysis(s, ecc, deg_sets, fc_list)
        
        # Save detailed results
        grad_df.to_csv(sample_dir / f"{s}_gradient_results.csv", index=False)
        
        # Save summary (without gene lists for readability)
        summary_df = grad_df.drop(columns=[col for col in grad_df.columns if col.startswith('genes_')])
        summary_df.to_csv(sample_dir / f"{s}_gradient_summary.csv", index=False)
        
        # Generate plots
        plot_gradient(s, grad_df, sample_dir)
        
        # Export gene lists
        export_gene_lists(grad_df, s, sample_dir)
        
        summary_tables.append(summary_df)
        logger.info(f"Completed processing for sample: {s}")

    # Combine all results
    logger.info("Creating combined analysis results")
    all_grad = pd.concat(summary_tables, ignore_index=True)
    all_grad.to_csv(output_root / "all_samples_gradient_analysis.csv", index=False)
    
    # Generate summary visualizations
    plot_summary_heatmap(all_grad, output_root)
    plot_overlap_summary(all_grad, output_root)
    
    # Summary statistics
    logger.info("="*80)
    logger.info("Analysis Summary:")
    logger.info(f"Total samples analyzed: {len(samples)}")
    logger.info(f"FC thresholds used: {fc_list}")
    
    # Report significant associations
    sig_associations = all_grad[
        (all_grad['p_enriched_up'] < 0.05) | (all_grad['p_depleted_down'] < 0.05)
    ]
    logger.info(f"Significant associations (p<0.05): {len(sig_associations)}")
    
    for _, row in sig_associations.iterrows():
        if row['p_enriched_up'] < 0.05:
            logger.info(f"  {row['sample']} @ FC={row['fc_threshold']}: "
                       f"Enriched-Up (OR={row['OR_enriched_up']:.2f}, p={row['p_enriched_up']:.3e})")
        if row['p_depleted_down'] < 0.05:
            logger.info(f"  {row['sample']} @ FC={row['fc_threshold']}: "
                       f"Depleted-Down (OR={row['OR_depleted_down']:.2f}, p={row['p_depleted_down']:.3e})")
    
    logger.info("="*80)
    logger.info(f"Analysis complete! Results saved to: {output_root}")


# ─────────────────────────────────────────────────────────────────────────────
# IX. CLI
# ─────────────────────────────────────────────────────────────────────────────
def cli():
    p = argparse.ArgumentParser(
        description="eccDNA & DEG association analysis with enhanced features",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("--eccdna", required=True, help="Merged eccDNA CSV file")
    p.add_argument("--degs", required=True, help="DEG table (TSV)")
    p.add_argument(
        "--outdir", 
        default="eccdna_deg_results", 
        help="Output directory"
    )
    p.add_argument(
        "--fc",
        default="1,2,4,6",
        help="Comma-separated log2FC thresholds",
    )
    p.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    args = p.parse_args()
    run(args)


if __name__ == "__main__":
    cli()
