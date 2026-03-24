"""DEG enrichment analysis at multiple FC thresholds."""

import logging
import os
from typing import List, Optional

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

logger = logging.getLogger(__name__)

DEFAULT_FC_THRESHOLDS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]


def run_deg_enrichment(
    input_file: str,
    deg_file: str,
    output_file: str,
    fc_thresholds: Optional[List[float]] = None,
    pval_threshold: float = 0.05,
) -> None:
    """
    Analyze DEG enrichment at multiple fold-change thresholds.

    For each FC threshold, computes how many eccDNA-associated genes are DEGs
    compared to background, using Fisher's exact test.

    Args:
        input_file: Gene list CSV (eccDNA-associated genes, column: gene_name or gene)
        deg_file: DEG results TSV (columns: gene_name, log2FC, padj/FDR)
        output_file: Output CSV file
        fc_thresholds: List of log2FC thresholds to test
        pval_threshold: Adjusted p-value threshold for DEG calling
    """
    if fc_thresholds is None:
        fc_thresholds = DEFAULT_FC_THRESHOLDS

    # Load eccDNA-associated gene list
    gene_df = pd.read_csv(input_file)
    col_map = {}
    for c in gene_df.columns:
        cl = c.lower().strip()
        if cl in ("gene_name", "gene", "gene_symbol", "symbol"):
            col_map[c] = "gene_name"
    gene_df = gene_df.rename(columns=col_map)

    if "gene_name" not in gene_df.columns:
        raise ValueError(
            f"Cannot find gene name column in {input_file}. "
            "Expected: gene_name, gene, gene_symbol, or symbol."
        )

    eccdna_genes = set(gene_df["gene_name"].dropna().astype(str).str.strip())
    logger.info(f"Loaded {len(eccdna_genes)} eccDNA-associated genes from {input_file}")

    # Load DEG results
    sep = "\t" if deg_file.endswith((".tsv", ".txt")) else ","
    deg_df = pd.read_csv(deg_file, sep=sep)

    col_map = {}
    for c in deg_df.columns:
        cl = c.lower().strip()
        if cl in ("gene_name", "gene", "gene_symbol", "symbol"):
            col_map[c] = "gene_name"
        elif cl in ("log2fc", "logfc", "log2foldchange"):
            col_map[c] = "log2FC"
        elif cl in ("padj", "fdr", "adj.p.val", "p_adjusted", "adjusted_pvalue"):
            col_map[c] = "padj"
        elif cl in ("p_value", "pvalue", "p.value"):
            col_map[c] = "p_value"
    deg_df = deg_df.rename(columns=col_map)

    if "gene_name" not in deg_df.columns:
        if deg_df.index.dtype == object:
            deg_df["gene_name"] = deg_df.index
        else:
            raise ValueError(f"Cannot find gene name column in {deg_file}.")

    if "log2FC" not in deg_df.columns:
        raise ValueError(f"Cannot find log2FC column in {deg_file}.")

    if "padj" not in deg_df.columns and "p_value" in deg_df.columns:
        logger.warning("No padj column found; using p_value instead")
        deg_df["padj"] = deg_df["p_value"]

    if "padj" not in deg_df.columns:
        raise ValueError(f"Cannot find p-value column in {deg_file}.")

    deg_df["gene_name"] = deg_df["gene_name"].astype(str).str.strip()
    all_genes = set(deg_df["gene_name"])
    eccdna_in_universe = eccdna_genes & all_genes
    n_eccdna = len(eccdna_in_universe)
    n_total = len(all_genes)

    logger.info(f"Universe: {n_total} genes, {n_eccdna} eccDNA-associated genes in universe")

    results = []
    for fc in sorted(fc_thresholds):
        deg_genes = set(
            deg_df.loc[
                (deg_df["log2FC"].abs() >= fc) & (deg_df["padj"] <= pval_threshold),
                "gene_name",
            ]
        )
        n_deg = len(deg_genes)

        # 2x2 table
        a = len(eccdna_in_universe & deg_genes)       # eccDNA & DEG
        b = len(eccdna_in_universe - deg_genes)        # eccDNA & non-DEG
        c = len(deg_genes - eccdna_in_universe)        # non-eccDNA & DEG
        d = len(all_genes - eccdna_in_universe - deg_genes)  # non-eccDNA & non-DEG

        # Enrichment ratio
        eccdna_deg_rate = a / n_eccdna if n_eccdna > 0 else 0
        bg_deg_rate = n_deg / n_total if n_total > 0 else 0
        enrichment_ratio = eccdna_deg_rate / bg_deg_rate if bg_deg_rate > 0 else np.inf

        # Fisher's exact test
        table = np.array([[a, b], [c, d]])
        odds_ratio, pvalue = fisher_exact(table, alternative="two-sided")

        results.append({
            "fc_threshold": fc,
            "n_deg_total": n_deg,
            "n_eccdna_deg": a,
            "n_eccdna_non_deg": b,
            "eccdna_deg_rate": round(eccdna_deg_rate, 4),
            "background_deg_rate": round(bg_deg_rate, 4),
            "enrichment_ratio": round(enrichment_ratio, 4),
            "odds_ratio": odds_ratio,
            "fisher_pvalue": pvalue,
        })

        logger.info(
            f"  |log2FC|>={fc}: DEGs={n_deg}, eccDNA_DEGs={a}, "
            f"enrichment={enrichment_ratio:.2f}, p={pvalue:.2e}"
        )

    results_df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    results_df.to_csv(output_file, index=False)
    logger.info(f"Saved enrichment results to {output_file}")
