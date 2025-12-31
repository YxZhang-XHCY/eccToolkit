"""Hotspot detection and analysis module."""

from ecctoolkit.hotspot.detect import run_sharp_detection
from ecctoolkit.hotspot.matrix import generate_multiscale_matrix
from ecctoolkit.hotspot.permtest import run_hotspot_permtest
from ecctoolkit.hotspot.refine import refine_hotspot_boundaries

__all__ = [
    "generate_multiscale_matrix",
    "run_sharp_detection",
    "run_hotspot_permtest",
    "refine_hotspot_boundaries",
]
