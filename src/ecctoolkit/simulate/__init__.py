"""Simulation module for eccDNA and sequencing reads."""

from ecctoolkit.simulate.reads import run_read_simulation
from ecctoolkit.simulate.region import run_region_simulation

__all__ = [
    "run_region_simulation",
    "run_read_simulation",
]
