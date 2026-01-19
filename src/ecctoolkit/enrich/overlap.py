"""General permutation test for genomic region overlap."""

import logging
from typing import List

logger = logging.getLogger(__name__)


def run_overlap_test(
    input_file: str,
    annotation_files: List[str],
    output_dir: str,
    genome_file: str,
    n_shuffles: int = 1000,
    n_cores: int = 8,
) -> None:
    """
    General permutation test for genomic region overlap enrichment.

    Tests whether query regions overlap with annotation regions more than
    expected by chance using bedtools shuffle.

    Args:
        input_file: Query regions BED/CSV file
        annotation_files: List of annotation BED files
        output_dir: Output directory
        genome_file: Genome sizes file
        n_shuffles: Number of permutations
        n_cores: Number of CPU cores
    """
    logger.info("General overlap permutation test - placeholder")
    logger.info(f"Query: {input_file}")
    logger.info(f"Annotations: {annotation_files}")
    logger.info(f"Permutations: {n_shuffles}, Cores: {n_cores}")
    # TODO: Migrate from permutation_overlap_test.py
    raise NotImplementedError("This module is pending migration")
