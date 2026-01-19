"""Simulation module for eccDNA and sequencing reads.

Provides two main components:
1. Region simulation (sim-region): Generate eccDNA sequences from reference genome
2. Read simulation (readsim): Simulate sequencing reads using external tools (ART, PBSIM2)

The unified `simulate` command combines both components in a single pipeline.
"""

from ecctoolkit.simulate.region import run_region_simulation
from ecctoolkit.simulate.readsim import seqsim, libsim, fqsim
from ecctoolkit.simulate.unified_config import UnifiedSimulateConfig
from ecctoolkit.simulate.pipeline import SimulatePipeline

__all__ = [
    # Region simulation
    "run_region_simulation",
    # Read simulation
    "seqsim",
    "libsim",
    "fqsim",
    # Unified pipeline
    "UnifiedSimulateConfig",
    "SimulatePipeline",
]
