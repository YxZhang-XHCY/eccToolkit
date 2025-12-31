"""
文库生成模块

支持三个平台：
- NGS (Illumina): paired-end短读段
- HiFi (PacBio): 高质量长读段
- ONT (Nanopore): 超长读段
"""

from .ngs import NGSLibraryGenerator, NGSLibraryStats
from .hifi import HiFiLibraryGenerator, HiFiLibraryStats
from .ont import ONTLibraryGenerator, ONTLibraryStats

__all__ = [
    'NGSLibraryGenerator', 'NGSLibraryStats',
    'HiFiLibraryGenerator', 'HiFiLibraryStats',
    'ONTLibraryGenerator', 'ONTLibraryStats'
]
