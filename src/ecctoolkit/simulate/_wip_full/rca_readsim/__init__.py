"""
RCA read simulator.

从 eccDNA 环状模板出发，模拟 RCA 扩增（支链/超支化）和三平台测序。
"""

from .config import SimConfig, get_default_config
from .models import EccDNA, SequencedRead, LinearMolecule, RCAMoleculeGraph
from .io_utils import parse_fasta, write_fastq, write_paired_fastq, build_ecc_db

__version__ = "1.0.0"
__all__ = [
    'SimConfig',
    'get_default_config',
    'EccDNA',
    'SequencedRead',
    'LinearMolecule',
    'RCAMoleculeGraph',
    'parse_fasta',
    'write_fastq',
    'write_paired_fastq',
    'build_ecc_db'
]
