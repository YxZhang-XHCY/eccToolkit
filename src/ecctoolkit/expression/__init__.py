"""Gene expression correlation analysis module."""

from ecctoolkit.expression.correlate import run_expression_correlation
from ecctoolkit.expression.enrich import run_deg_enrichment

__all__ = [
    "run_expression_correlation",
    "run_deg_enrichment",
]
