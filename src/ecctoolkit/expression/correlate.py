"""eccDNA-DEG correlation analysis using Fisher's exact test."""

import logging

logger = logging.getLogger(__name__)


def run_expression_correlation(
    eccdna_file: str,
    deg_file: str,
    output_dir: str,
    mode: str = "gradient",
) -> None:
    """
    Analyze correlation between eccDNA enrichment and DEGs using Fisher test.

    Args:
        eccdna_file: eccDNA enrichment CSV file
        deg_file: DEG results TSV file
        output_dir: Output directory
        mode: Analysis mode - "single" (single threshold) or "gradient" (gradient FC)
    """
    logger.info("Expression correlation analysis - placeholder")
    logger.info(f"eccDNA: {eccdna_file}")
    logger.info(f"DEG: {deg_file}")
    logger.info(f"Mode: {mode}")
    # TODO: Migrate from eccdna_deg_correlation_analysis_*.py
    raise NotImplementedError("This module is pending migration")
