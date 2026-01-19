"""
分子打断模块

用于 NGS 和 HiFi 平台的分子打断处理：
- NGS: 打断成短片段 (~400bp insert size)
- HiFi: 打断成中等片段 (15-25kb)

打断会自然解除分支结构。
"""

from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass, field
import numpy as np

from .models import (
    EccDNA, LinearMolecule, Segment, RCAMoleculeGraph,
    Strand, SegmentType
)
from .config import SimConfig


@dataclass
class FragmentationStats:
    """打断统计"""
    input_molecules: int = 0
    output_fragments: int = 0
    total_input_length: int = 0
    total_output_length: int = 0
    avg_fragment_length: float = 0.0

    def summary(self) -> str:
        return (
            f"Fragmentation: {self.input_molecules} molecules -> "
            f"{self.output_fragments} fragments, "
            f"avg length: {self.avg_fragment_length:.0f}bp"
        )


class Fragmenter:
    """分子打断器"""

    def __init__(
        self,
        config: SimConfig,
        rng: Optional[np.random.Generator] = None
    ):
        self.config = config
        self.rng = rng if rng is not None else np.random.default_rng()
        self._frag_counter = 0  # 全局唯一 fragment 计数器

    def fragment_for_ngs(
        self,
        molecules: List[LinearMolecule],
        ecc_db: Dict[str, EccDNA],
        insert_mean: float = 400.0,
        insert_std: float = 100.0,
        min_insert: int = 150,
        max_insert: int = 800,
    ) -> Tuple[List[LinearMolecule], FragmentationStats]:
        """
        为 NGS 打断分子

        每个分子被打断成多个短片段，模拟 Illumina 文库制备。
        打断会自然解除分支结构。

        Args:
            molecules: 输入分子列表
            ecc_db: eccDNA 数据库
            insert_mean: insert size 均值
            insert_std: insert size 标准差
            min_insert: 最小 insert size
            max_insert: 最大 insert size

        Returns:
            (片段列表, 统计信息)
        """
        stats = FragmentationStats()
        stats.input_molecules = len(molecules)

        fragments = []
        for mol in molecules:
            mol_len = mol.total_length
            stats.total_input_length += mol_len

            if mol_len < min_insert:
                # 分子太短，跳过
                continue

            # 计算可以产生多少片段
            # 使用泊松过程模拟随机打断点
            n_fragments = max(1, int(mol_len / insert_mean))

            # 生成打断点
            break_points = sorted(self.rng.integers(0, mol_len, size=n_fragments - 1))
            break_points = [0] + list(break_points) + [mol_len]

            for i in range(len(break_points) - 1):
                start = break_points[i]
                end = break_points[i + 1]
                frag_len = end - start

                # 过滤太短或太长的片段
                if frag_len < min_insert or frag_len > max_insert:
                    continue

                # 提取片段的 segments
                frag_segments = mol.extract_subregion(start, frag_len)
                if not frag_segments:
                    continue

                self._frag_counter += 1
                frag = LinearMolecule(
                    molecule_id=f"ngs_frag_{self._frag_counter}",
                    segments=frag_segments,
                    source_graph_id=mol.source_graph_id,
                    is_from_branch=False,
                    repeat_count=mol.repeat_count,
                    source_ecc_length=mol.source_ecc_length,
                    is_background=mol.is_background,
                    background_chrom=mol.background_chrom,
                    background_start=mol.background_start,
                    background_end=mol.background_end,
                    has_chimera=mol.has_chimera,
                )
                fragments.append(frag)
                stats.total_output_length += frag_len

        stats.output_fragments = len(fragments)
        if fragments:
            stats.avg_fragment_length = stats.total_output_length / len(fragments)

        return fragments, stats

    def fragment_for_hifi(
        self,
        molecules: List[LinearMolecule],
        ecc_db: Dict[str, EccDNA],
        target_mean: float = 20000.0,
        target_std: float = 2000.0,
        min_length: int = 10000,
        max_length: int = 30000,
    ) -> Tuple[List[LinearMolecule], FragmentationStats]:
        """
        为 HiFi 打断分子

        每个分子被打断成中等长度片段，模拟 PacBio 文库制备。
        打断会自然解除分支结构。

        Args:
            molecules: 输入分子列表
            ecc_db: eccDNA 数据库
            target_mean: 目标片段长度均值
            target_std: 目标片段长度标准差
            min_length: 最小片段长度
            max_length: 最大片段长度

        Returns:
            (片段列表, 统计信息)
        """
        stats = FragmentationStats()
        stats.input_molecules = len(molecules)

        fragments = []
        for mol in molecules:
            mol_len = mol.total_length
            stats.total_input_length += mol_len

            if mol_len < min_length:
                # 分子太短，整个保留（如果足够长）
                if mol_len >= min_length // 2:
                    # 创建一个完整分子的拷贝
                    self._frag_counter += 1
                    frag = LinearMolecule(
                        molecule_id=f"hifi_frag_{self._frag_counter}",
                        segments=mol.segments.copy(),
                        source_graph_id=mol.source_graph_id,
                        is_from_branch=mol.is_from_branch,
                        repeat_count=mol.repeat_count,
                        source_ecc_length=mol.source_ecc_length,
                        is_background=mol.is_background,
                        background_chrom=mol.background_chrom,
                        background_start=mol.background_start,
                        background_end=mol.background_end,
                        has_chimera=mol.has_chimera,
                    )
                    fragments.append(frag)
                    stats.total_output_length += mol_len
                continue

            # 计算打断点数量
            n_fragments = max(1, int(mol_len / target_mean))

            if n_fragments == 1:
                # 只需要一个片段
                frag_len = min(mol_len, max_length)
                start = self.rng.integers(0, max(1, mol_len - frag_len + 1))
                frag_segments = mol.extract_subregion(start, frag_len)
                if frag_segments:
                    self._frag_counter += 1
                    frag = LinearMolecule(
                        molecule_id=f"hifi_frag_{self._frag_counter}",
                        segments=frag_segments,
                        source_graph_id=mol.source_graph_id,
                        is_from_branch=False,
                        repeat_count=mol.repeat_count,
                        source_ecc_length=mol.source_ecc_length,
                        is_background=mol.is_background,
                        background_chrom=mol.background_chrom,
                        background_start=mol.background_start,
                        background_end=mol.background_end,
                        has_chimera=mol.has_chimera,
                    )
                    fragments.append(frag)
                    stats.total_output_length += frag_len
            else:
                # 多个片段
                break_points = sorted(self.rng.integers(0, mol_len, size=n_fragments - 1))
                break_points = [0] + list(break_points) + [mol_len]

                for i in range(len(break_points) - 1):
                    start = break_points[i]
                    end = break_points[i + 1]
                    frag_len = end - start

                    # 过滤不合适长度的片段
                    if frag_len < min_length or frag_len > max_length:
                        continue

                    frag_segments = mol.extract_subregion(start, frag_len)
                    if not frag_segments:
                        continue

                    self._frag_counter += 1
                    frag = LinearMolecule(
                        molecule_id=f"hifi_frag_{self._frag_counter}",
                        segments=frag_segments,
                        source_graph_id=mol.source_graph_id,
                        is_from_branch=False,
                        repeat_count=mol.repeat_count,
                        source_ecc_length=mol.source_ecc_length,
                        is_background=mol.is_background,
                        background_chrom=mol.background_chrom,
                        background_start=mol.background_start,
                        background_end=mol.background_end,
                        has_chimera=mol.has_chimera,
                    )
                    fragments.append(frag)
                    stats.total_output_length += frag_len

        stats.output_fragments = len(fragments)
        if fragments:
            stats.avg_fragment_length = stats.total_output_length / len(fragments)

        return fragments, stats


def fragment_molecules(
    molecules: List[LinearMolecule],
    ecc_db: Dict[str, EccDNA],
    platform: str,
    config: SimConfig,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[List[LinearMolecule], FragmentationStats]:
    """
    根据平台打断分子

    Args:
        molecules: 输入分子列表
        ecc_db: eccDNA 数据库
        platform: 平台类型 ("NGS", "HiFi")
        config: 配置
        rng: 随机数生成器

    Returns:
        (片段列表, 统计信息)
    """
    fragmenter = Fragmenter(config, rng)

    if platform.upper() == "NGS":
        return fragmenter.fragment_for_ngs(
            molecules,
            ecc_db,
            insert_mean=config.I_mu,
            insert_std=config.constants.insert_size_std,
            min_insert=150,
            max_insert=800,
        )
    elif platform.upper() == "HIFI":
        return fragmenter.fragment_for_hifi(
            molecules,
            ecc_db,
            target_mean=config.constants.hifi_target_size_mean,
            target_std=config.constants.hifi_target_size_std,
            min_length=config.constants.hifi_min_length,
            max_length=config.constants.hifi_max_length,
        )
    else:
        raise ValueError(f"Unknown platform for fragmentation: {platform}")
