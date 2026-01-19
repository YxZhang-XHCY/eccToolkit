"""Transposable element (TE) analysis module."""

from ecctoolkit.te.analyze import run_te_analysis
from ecctoolkit.te.classify import classify_te_composition
from ecctoolkit.te.composition import analyze_te_composition
from ecctoolkit.te.distribution import analyze_te_distribution
from ecctoolkit.te.process import process_te_data

__all__ = [
    "run_te_analysis",
    "classify_te_composition",
    "analyze_te_distribution",
    "analyze_te_composition",
    "process_te_data",
]
