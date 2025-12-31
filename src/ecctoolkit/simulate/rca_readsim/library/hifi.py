"""
PacBio HiFi文库生成模块

实现§9.4：HiFi文库生成
- DNA打断（fragmentation）模拟
- size selection
- 高质量CCS reads
- truth追踪
"""

from typing import List, Dict, Tuple, Optional
import numpy as np
from dataclasses import dataclass

from ..models import (
    LinearMolecule, SequencedRead, Segment, EccDNA,
    Strand, SegmentType, get_ecc_length
)
from ..config import SimConfig
from ..debranch import LinearMoleculePool
from ..seq_utils import quality_string
from ..error_models import PacBioHiFiErrorModel
from .base import extract_segments, compute_repeat_count, compute_repeat_count_by_source


@dataclass
class HiFiLibraryStats:
    """HiFi文库统计"""
    total_reads: int = 0
    mean_read_length: float = 0.0
    max_read_length: int = 0
    min_read_length: int = 0
    reads_from_trunk: int = 0
    reads_from_branch: int = 0
    reads_with_chimera: int = 0
    molecules_fragmented: int = 0      # 被打断的分子数
    total_fragments_generated: int = 0  # 生成的片段总数

    def summary(self) -> str:
        return (
            f"HiFi: {self.total_reads} reads, "
            f"Length: {self.mean_read_length:.0f}bp "
            f"({self.min_read_length}-{self.max_read_length}), "
            f"Fragmented: {self.molecules_fragmented} molecules -> "
            f"{self.total_fragments_generated} fragments, "
            f"Chimera: {self.reads_with_chimera}"
        )


class HiFiLibraryGenerator:
    """PacBio HiFi文库生成器"""

    def __init__(
        self,
        config: SimConfig,
        rng: Optional[np.random.Generator] = None
    ):
        self.config = config
        self.rng = rng if rng is not None else np.random.default_rng()
        self._read_counter = 0

        # 初始化PacBio HiFi错误模型（非常低的错误率）
        self.error_model = PacBioHiFiErrorModel(
            substitution_rate=0.001,  # 0.1% 替换
            insertion_rate=0.0001,    # 0.01% 插入
            deletion_rate=0.0001,     # 0.01% 删除
            rng=self.rng
        )

    def fragment_molecule(
        self,
        mol: LinearMolecule,
        mol_seq: str
    ) -> List[Tuple[int, int, LinearMolecule]]:
        """
        将长分子随机打断成目标大小的片段

        模拟DNA shearing过程：
        1. 如果分子长度在目标范围内，直接返回
        2. 如果分子太长，随机打断成多个片段
        3. 每个片段的大小服从正态分布

        Args:
            mol: 线性分子
            mol_seq: 分子序列

        Returns:
            [(start, length, fragment_mol), ...] 片段信息列表
        """
        mol_length = len(mol_seq)
        min_len = self.config.constants.hifi_min_length
        max_len = self.config.constants.hifi_max_length
        target_mean = self.config.constants.hifi_target_size_mean
        target_std = self.config.constants.hifi_target_size_std

        # 如果分子长度在目标范围内，直接返回
        if mol_length <= max_len:
            if mol_length >= min_len:
                return [(0, mol_length, mol)]
            else:
                return []  # 太短

        # 长分子需要打断
        fragments = []
        # Randomize start to avoid phase artifacts at molecule boundaries.
        first_frag_size = int(self.rng.normal(target_mean, target_std))
        first_frag_size = max(min_len, min(max_len, first_frag_size))
        start_offset = int(self.rng.integers(0, max(1, first_frag_size)))
        pos = start_offset
        use_first = True

        while pos < mol_length:
            # 采样片段大小
            if use_first:
                frag_size = first_frag_size
                use_first = False
            else:
                frag_size = int(self.rng.normal(target_mean, target_std))
            frag_size = max(min_len, min(max_len, frag_size))  # 限制范围

            # 确保不超过分子剩余长度
            frag_size = min(frag_size, mol_length - pos)

            # 过滤太短的片段
            if frag_size >= min_len:
                # 创建片段的LinearMolecule（提取对应的segments）
                frag_mol = self._create_fragment_molecule(mol, pos, frag_size)
                if frag_mol:
                    fragments.append((pos, frag_size, frag_mol))

            pos += frag_size

        # Add the prefix fragment (0..start_offset) if it is long enough.
        if start_offset >= min_len:
            frag_mol = self._create_fragment_molecule(mol, 0, start_offset)
            if frag_mol:
                fragments.append((0, start_offset, frag_mol))

        return fragments

    def _create_fragment_molecule(
        self,
        mol: LinearMolecule,
        start: int,
        length: int
    ) -> Optional[LinearMolecule]:
        """
        从分子的指定区域创建片段分子

        Args:
            mol: 原始分子
            start: 起始位置
            length: 片段长度

        Returns:
            片段的LinearMolecule
        """
        # 提取对应区域的segments
        frag_segments = self._extract_segments(mol, start, length)

        if not frag_segments:
            return None

        # 检查是否有chimera在片段范围内
        has_chimera = False
        chimera_positions = []
        if mol.has_chimera:
            for chim_pos in mol.chimera_positions:
                if start <= chim_pos < start + length:
                    has_chimera = True
                    chimera_positions.append(chim_pos - start)

        # 计算片段内的活跃锚点
        frag_anchors = [
            a - start for a in mol.active_branch_anchors
            if start <= a < start + length
        ]

        return LinearMolecule(
            molecule_id=f"{mol.molecule_id}_frag_{start}",
            segments=frag_segments,
            source_graph_id=mol.source_graph_id,
            is_from_branch=mol.is_from_branch,
            has_chimera=has_chimera,
            chimera_positions=chimera_positions,
            active_branch_anchors=frag_anchors
        )

    def _extract_segments(
        self,
        mol: LinearMolecule,
        start: int,
        length: int
    ) -> List[Segment]:
        """提取指定范围的segments（委托给共享函数）"""
        return extract_segments(mol, start, length)

    def generate_read(
        self,
        mol: LinearMolecule,
        mol_seq: str,
        frag_start: int = 0,
        frag_length: Optional[int] = None
    ) -> Optional[SequencedRead]:
        """
        从分子（或片段）生成HiFi read

        Args:
            mol: 线性分子（或片段分子）
            mol_seq: 分子序列
            frag_start: 片段在原始分子中的起始位置
            frag_length: 片段长度

        Returns:
            SequencedRead或None
        """
        if frag_length is None:
            frag_length = len(mol_seq)

        self._read_counter += 1
        read_id = f"hifi_read_{self._read_counter}"

        # 提取序列
        read_seq = mol_seq[frag_start:frag_start + frag_length] if frag_start > 0 else mol_seq[:frag_length]

        # 应用PacBio HiFi错误模型（低错误率 + 偶尔indel）
        read_seq, read_qual = self.error_model.apply(read_seq)

        # 使用片段的segments
        segments = list(mol.segments)

        # 收集来源eccDNA
        ecc_ids = list(set(s.ecc_id for s in segments))

        # 检查chimera
        has_chimera = mol.has_chimera
        chimera_bps = mol.chimera_positions.copy() if mol.chimera_positions else []

        # 检查分支锚点
        covered_anchor = len(mol.active_branch_anchors) > 0

        return SequencedRead(
            read_id=read_id,
            platform="HiFi",
            sequence=read_seq,
            quality=read_qual,
            source_molecule_id=mol.molecule_id,
            source_ecc_ids=ecc_ids,
            segments=segments,
            repeat_count_truth=self._compute_repeat_count(segments),
            repeat_count_by_source=self._compute_repeat_count_by_source(segments),
            has_inter_chimera=has_chimera,
            chimera_breakpoints=chimera_bps,
            covered_branch_anchor=covered_anchor,
            truncated=False,
            junction_covered_possible=mol.junction_covered_possible,
            # Background DNA fields
            is_background=mol.is_background,
            background_chrom=mol.background_chrom,
            background_start=mol.background_start,
            background_end=mol.background_end
        )

    def _compute_repeat_count(self, segments: List[Segment]) -> int:
        """计算repeat次数（委托给共享函数）"""
        return compute_repeat_count(segments)

    def _compute_repeat_count_by_source(self, segments: List[Segment]) -> dict:
        """计算各eccDNA来源repeat比例（委托给共享函数）"""
        return compute_repeat_count_by_source(segments)

    def generate(
        self,
        pool: LinearMoleculePool,
        target_reads: int
    ) -> Tuple[List[SequencedRead], HiFiLibraryStats]:
        """
        生成HiFi reads

        流程：
        1. 从分子池抽样分子
        2. 长分子进行打断（fragmentation）
        3. 对每个片段生成read

        Args:
            pool: 线性分子池
            target_reads: 目标read数

        Returns:
            (reads列表, 统计信息)
        """
        stats = HiFiLibraryStats()
        reads = []
        read_lengths = []

        min_len = self.config.constants.hifi_min_length
        max_len = self.config.constants.hifi_max_length

        attempts = 0
        max_attempts = target_reads * 10  # 增加尝试次数

        while stats.total_reads < target_reads and attempts < max_attempts:
            attempts += 1

            # 抽样分子
            mol = pool.sample(1, self.rng)[0]
            mol_seq = pool.get_sequence(mol)
            mol_len = len(mol_seq)

            # 跳过太短的分子
            if mol_len < min_len:
                continue

            # 打断分子（如果需要）
            fragments = self.fragment_molecule(mol, mol_seq)

            if not fragments:
                continue

            # 统计打断情况
            if mol_len > max_len:
                stats.molecules_fragmented += 1
                stats.total_fragments_generated += len(fragments)

            # 对每个片段生成read
            for frag_start, frag_length, frag_mol in fragments:
                if stats.total_reads >= target_reads:
                    break

                # 提取片段序列
                frag_seq = mol_seq[frag_start:frag_start + frag_length]

                # 生成read
                read = self.generate_read(frag_mol, frag_seq, 0, frag_length)

                if read is None:
                    continue

                reads.append(read)
                read_lengths.append(len(read.sequence))

                # 更新统计
                stats.total_reads += 1

                if mol.is_from_branch:
                    stats.reads_from_branch += 1
                else:
                    stats.reads_from_trunk += 1

                if read.has_inter_chimera:
                    stats.reads_with_chimera += 1

        # 检查是否因达到最大尝试次数而退出
        if attempts >= max_attempts and stats.total_reads < target_reads:
            import warnings
            warnings.warn(
                f"HiFi sampling reached max attempts ({max_attempts}). "
                f"Generated {stats.total_reads}/{target_reads} reads. "
                f"Consider checking molecule pool lengths (HiFi needs long molecules)."
            )

        # 计算长度统计
        if read_lengths:
            stats.mean_read_length = np.mean(read_lengths)
            stats.max_read_length = max(read_lengths)
            stats.min_read_length = min(read_lengths)

        return reads, stats
