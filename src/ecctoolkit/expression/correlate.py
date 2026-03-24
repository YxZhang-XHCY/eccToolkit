"""eccDNA-DEG correlation analysis using Fisher's exact test."""

import logging
import os
from typing import List, Optional

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

logger = logging.getLogger(__name__)

DEFAULT_FC_THRESHOLDS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]


def _load_eccdna_genes(eccdna_file: str, enrichment_threshold: float) -> set:
    """Load eccDNA-enriched gene set from CSV.

    Expects columns: gene_name (or gene), and enrichment_score (or eccDNA_count).
    """
    df = pd.read_csv(eccdna_file)
    col_map = {}
    for c in df.columns:
        cl = c.lower().strip()
        if cl in ("gene_name", "gene", "gene_symbol", "symbol"):
            col_map[c] = "gene_name"
        elif cl in ("enrichment_score", "enrichment", "eccdna_count", "count", "score"):
            col_map[c] = "enrichment_score"
    df = df.rename(columns=col_map)

    if "gene_name" not in df.columns:
        raise ValueError(
            f"Cannot find gene name column in {eccdna_file}. "
            "Expected: gene_name, gene, gene_symbol, or symbol."
        )

    if "enrichment_score" in df.columns:
        enriched = df.loc[
            df["enrichment_score"] >= enrichment_threshold, "gene_name"
        ]
    else:
        logger.warning(
            "No enrichment_score column found; treating all genes as enriched"
        )
        enriched = df["gene_name"]

    gene_set = set(enriched.dropna().astype(str).str.strip())
    logger.info(f"Loaded {len(gene_set)} eccDNA-enriched genes (threshold={enrichment_threshold})")
    return gene_set


def _load_deg_results(deg_file: str) -> pd.DataFrame:
    """Load DEG results from TSV/CSV.

    Expects columns: gene_name (or gene), log2FC (or logFC), p_value/padj/FDR.
    """
    sep = "\t" if deg_file.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(deg_file, sep=sep)

    col_map = {}
    for c in df.columns:
        cl = c.lower().strip()
        if cl in ("gene_name", "gene", "gene_symbol", "symbol"):
            col_map[c] = "gene_name"
        elif cl in ("log2fc", "logfc", "log2foldchange"):
            col_map[c] = "log2FC"
        elif cl in ("padj", "fdr", "adj.p.val", "p_adjusted", "adjusted_pvalue"):
            col_map[c] = "padj"
        elif cl in ("p_value", "pvalue", "p.value"):
            col_map[c] = "p_value"
    df = df.rename(columns=col_map)

    if "gene_name" not in df.columns:
        # Try using index (row names) as gene names
        if df.index.dtype == object:
            df["gene_name"] = df.index
        else:
            raise ValueError(
                f"Cannot find gene name column in {deg_file}. "
                "Expected: gene_name, gene, gene_symbol, or symbol."
            )

    if "log2FC" not in df.columns:
        raise ValueError(f"Cannot find log2FC column in {deg_file}.")

    if "padj" not in df.columns and "p_value" in df.columns:
        logger.warning("No padj column found; using p_value instead")
        df["padj"] = df["p_value"]

    if "padj" not in df.columns:
        raise ValueError(f"Cannot find p-value column in {deg_file}.")

    df["gene_name"] = df["gene_name"].astype(str).str.strip()
    return df


def _fisher_test_at_threshold(
    eccdna_genes: set,
    deg_df: pd.DataFrame,
    fc_threshold: float,
    pval_threshold: float,
) -> dict:
    """Run Fisher's exact test at a given FC threshold."""
    all_genes = set(deg_df["gene_name"].dropna().astype(str))
    deg_genes = set(
        deg_df.loc[
            (deg_df["log2FC"].abs() >= fc_threshold)
            & (deg_df["padj"] <= pval_threshold),
            "gene_name",
        ]
    )

    # Only consider genes present in both datasets
    universe = all_genes
    eccdna_in_universe = eccdna_genes & universe

    a = len(eccdna_in_universe & deg_genes)      # eccDNA-enriched AND DEG
    b = len(eccdna_in_universe - deg_genes)       # eccDNA-enriched AND non-DEG
    c = len(deg_genes - eccdna_in_universe)       # non-enriched AND DEG
    d = len(universe - eccdna_in_universe - deg_genes)  # non-enriched AND non-DEG

    table = np.array([[a, b], [c, d]])
    odds_ratio, pvalue = fisher_exact(table, alternative="two-sided")

    return {
        "fc_threshold": fc_threshold,
        "pval_threshold": pval_threshold,
        "n_eccdna_genes": len(eccdna_in_universe),
        "n_deg_genes": len(deg_genes),
        "n_overlap": a,
        "n_eccdna_only": b,
        "n_deg_only": c,
        "n_neither": d,
        "odds_ratio": odds_ratio,
        "fisher_pvalue": pvalue,
        "universe_size": len(universe),
    }


def run_expression_correlation(
    eccdna_file: str,
    deg_file: str,
    output_dir: str,
    mode: str = "gradient",
    fc_thresholds: Optional[List[float]] = None,
    pval_threshold: float = 0.05,
    enrichment_threshold: float = 1.0,
) -> None:
    """
    Analyze correlation between eccDNA enrichment and DEGs using Fisher test.

    Args:
        eccdna_file: eccDNA enrichment CSV file
        deg_file: DEG results TSV file
        output_dir: Output directory
        mode: Analysis mode - "single" (single threshold) or "gradient" (gradient FC)
        fc_thresholds: List of FC thresholds for gradient mode
        pval_threshold: Adjusted p-value threshold for DEG calling
        enrichment_threshold: Minimum enrichment score for eccDNA genes
    """
    os.makedirs(output_dir, exist_ok=True)

    logger.info(f"Loading eccDNA data from {eccdna_file}")
    eccdna_genes = _load_eccdna_genes(eccdna_file, enrichment_threshold)

    logger.info(f"Loading DEG data from {deg_file}")
    deg_df = _load_deg_results(deg_file)
    logger.info(f"Loaded {len(deg_df)} genes from DEG results")

    if fc_thresholds is None:
        fc_thresholds = DEFAULT_FC_THRESHOLDS

    if mode == "single":
        thresholds = [fc_thresholds[0]] if fc_thresholds else [1.0]
    elif mode == "gradient":
        thresholds = sorted(fc_thresholds)
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'single' or 'gradient'.")

    results = []
    for fc in thresholds:
        result = _fisher_test_at_threshold(eccdna_genes, deg_df, fc, pval_threshold)
        results.append(result)
        logger.info(
            f"  FC>={fc}: overlap={result['n_overlap']}, "
            f"OR={result['odds_ratio']:.3f}, p={result['fisher_pvalue']:.2e}"
        )

    results_df = pd.DataFrame(results)
    out_file = os.path.join(output_dir, "correlation_results.csv")
    results_df.to_csv(out_file, index=False)
    logger.info(f"Saved correlation results to {out_file}")

    # Summary
    if mode == "gradient" and len(results) > 1:
        sig_results = results_df[results_df["fisher_pvalue"] < 0.05]
        logger.info(
            f"Gradient summary: {len(sig_results)}/{len(results)} thresholds "
            f"show significant correlation (p<0.05)"
        )
