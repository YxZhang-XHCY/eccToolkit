"""Enrichment analysis module - permutation tests for genomic regions."""

from ecctoolkit.enrich.cnv import run_cnv_enrichment
from ecctoolkit.enrich.overlap import run_overlap_test
from ecctoolkit.enrich.tad import run_tad_enrichment

__all__ = [
    "run_cnv_enrichment",
    "run_tad_enrichment",
    "run_overlap_test",
]
