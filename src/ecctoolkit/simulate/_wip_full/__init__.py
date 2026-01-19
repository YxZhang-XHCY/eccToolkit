"""
ecsim full mode - Detailed RCA (Rolling Circle Amplification) simulation.

This module provides a biophysically accurate simulation of eccDNA sequencing,
including RCA amplification, branching, chimera formation, and platform-specific
read generation.
"""

from .reads import run_read_simulation

__all__ = ["run_read_simulation"]
