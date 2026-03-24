"""
Benchmark module for evaluating eccDNA detection results against ground truth.

This module provides tools to:
- Compare CircleSeeker (or other tool) outputs against simulated ground truth
- Calculate precision, recall, F1 scores
- Analyze performance by eccDNA type (U/M/C) and detection source (Confirmed/Inferred)
- Generate comprehensive benchmark reports
"""

from ecctoolkit.benchmark.comparator import BenchmarkComparator
from ecctoolkit.benchmark.evaluator import BenchmarkEvaluator
from ecctoolkit.benchmark.metrics import BenchmarkMetrics, TypeMetrics
from ecctoolkit.benchmark.parsers import SUPPORTED_FORMATS, parse_tool_output
from ecctoolkit.benchmark.report import BenchmarkReporter

__all__ = [
    "BenchmarkComparator",
    "BenchmarkEvaluator",
    "BenchmarkMetrics",
    "TypeMetrics",
    "BenchmarkReporter",
    "parse_tool_output",
    "SUPPORTED_FORMATS",
]
