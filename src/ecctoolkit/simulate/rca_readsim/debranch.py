"""
去支化/线性化模块

实现§8：从分子图得到可测序的线性dsDNA分子集合
- 去支化模型（边删除模型）
- Nick断裂模型（T7处理后的效应）
- 线性化输出
"""

from typing import List, Dict, Tuple, Optional, Set
import numpy as np
import warnings
from dataclasses import dataclass, field

from .models import (
    EccDNA, RCAMoleculeGraph, LinearMolecule, Segment, BranchNode,
    Strand, SegmentType, ChimeraJunction, ActiveBranchInfo,
    compute_junction_covered_possible
)
from .config import SimConfig


def validate_non_overlapping_junctions(junctions: List[ChimeraJunction]) -> List[ChimeraJunction]:
    """
    验证并修复重叠的chimera junctions

    如果发现重叠，发出警告并移除后一个重叠的junction

    Args:
        junctions: 按position排序的junction列表

    Returns:
        非重叠的junction列表
    """
    if len(junctions) <= 1:
        return junctions

    valid_junctions = [junctions[0]]

    for i in range(1, len(junctions)):
        prev = valid_junctions[-1]
        curr = junctions[i]

        prev_end = prev.position + prev.donor_segment.length
        curr_start = curr.position

        if prev_end > curr_start:
            # 发现重叠
            warnings.warn(
                f"Overlapping chimera junctions detected: "
                f"junction at {prev.position} (end={prev_end}) overlaps with "
                f"junction at {curr_start}. Removing the overlapping junction."
            )
            # 跳过当前junction（保留前一个）
            continue

        valid_junctions.append(curr)

    return valid_junctions


@dataclass
class DebranchStats:
    """去支化统计"""
    total_branches: int = 0
    debranched_count: int = 0
    retained_count: int = 0
    nick_count: int = 0                    # 产生的nick数
    break_count: int = 0                   # 断裂数
    linear_molecules_from_trunk: int = 0
    linear_molecules_from_branch: int = 0
    trunk_fragments_from_breaks: int = 0   # 因断裂产生的主干片段数

    @property
    def debranch_rate(self) -> float:
        if self.total_branches == 0:
            return 0.0
        return self.debranched_count / self.total_branches

    def summary(self) -> str:
        return (
            f"Branches: {self.total_branches} "
            f"(debranched: {self.debranched_count}, retained: {self.retained_count}), "
            f"Rate: {self.debranch_rate:.2%}, "
            f"Nicks: {self.nick_count}, Breaks: {self.break_count}, "
            f"Linear molecules: {self.linear_molecules_from_trunk} trunk "
            f"(+{self.trunk_fragments_from_breaks} fragments) + "
            f"{self.linear_molecules_from_branch} branch"
        )


class Debrancher:
    """去支化处理器"""

    def __init__(self, config: SimConfig, rng: Optional[np.random.Generator] = None):
        self.config = config
        self.rng = rng if rng is not None else np.random.default_rng()

    def debranch_molecule(self, mol: RCAMoleculeGraph) -> None:
        """
        对分子进行去支化处理（原地修改）

        对每个分支，以概率D_eff_actual断开连接。
        如果k_debranch > 0，长分支更难去除：
            D_eff_actual = D_eff * exp(-k_debranch * (L_branch / L_ref))

        Args:
            mol: RCA分子图（会被修改）
        """
        import math

        D_eff_base = self.config.D_eff
        k_debranch = self.config.k_debranch
        L_ref = self.config.L_ref_debranch

        for branch in mol.branches:
            # 计算长度依赖的去支化效率
            branch_length = branch.segment.length if branch.segment else 0

            if k_debranch > 0 and branch_length > 0:
                # 长分支更难去除
                length_factor = math.exp(-k_debranch * (branch_length / L_ref))
                D_eff_actual = D_eff_base * length_factor
            else:
                # 长度无关模式（默认）
                D_eff_actual = D_eff_base

            if self.rng.random() < D_eff_actual:
                branch.debranched = True

    def calculate_break_positions(
        self,
        mol: RCAMoleculeGraph
    ) -> Tuple[List[int], int, int]:
        """
        计算由nick导致的断裂位置

        模型：
        1. 每个去支化的锚点有 p_nick 概率产生nick
        2. 每个nick有 p_break 概率转为双链断裂

        Args:
            mol: RCA分子图

        Returns:
            (断裂位置列表, nick数, 断裂数)
        """
        p_nick = self.config.constants.p_nick_per_debranch
        p_break = self.config.constants.p_nick_to_break

        break_positions = []
        nick_count = 0

        # 只考虑已去支化的锚点
        debranched_anchors = [b.anchor_pos for b in mol.branches if b.debranched]

        for anchor_pos in debranched_anchors:
            # Step 1: 是否产生nick
            if self.rng.random() < p_nick:
                nick_count += 1

                # Step 2: nick是否转为断裂
                if self.rng.random() < p_break:
                    break_positions.append(anchor_pos)

        # 排序并去重（多个锚点可能在同一位置）
        break_positions = sorted(set(break_positions))

        return break_positions, nick_count, len(break_positions)

    def split_trunk_at_breaks(
        self,
        mol: RCAMoleculeGraph,
        trunk_mol: LinearMolecule,
        break_positions: List[int],
        ecc_db: Dict[str, EccDNA]
    ) -> List[LinearMolecule]:
        """
        在断裂位置切分主干分子

        Args:
            mol: 原始RCA分子图
            trunk_mol: 原始主干线性分子
            break_positions: 断裂位置列表
            ecc_db: eccDNA字典

        Returns:
            切分后的线性分子列表
        """
        if not break_positions:
            return [trunk_mol]

        ecc = mol.source_ecc
        fragments = []

        # 添加起始和结束位置
        all_positions = [0] + break_positions + [mol.trunk_length]

        # 收集未去支化的分支信息
        active_branch_list = [
            (b.anchor_pos, ActiveBranchInfo(
                anchor_pos=b.anchor_pos,
                branch_segment=b.segment,
                branch_id=b.branch_id
            ))
            for b in mol.branches if not b.debranched
        ]

        # 收集已去支化的锚点位置（用于解耦的残余结构影响）
        debranched_anchor_list = [b.anchor_pos for b in mol.branches if b.debranched]

        for i in range(len(all_positions) - 1):
            start_pos = all_positions[i]
            end_pos = all_positions[i + 1]
            frag_length = end_pos - start_pos

            if frag_length < 100:  # 忽略太短的片段
                continue

            # 创建片段的segment
            # 需要处理chimera的情况
            frag_segments = self._extract_segments_in_range(
                mol, start_pos, frag_length, ecc_db
            )

            if not frag_segments:
                continue

            # 计算这个片段内的活跃锚点和分支信息
            frag_active_anchors = []
            frag_active_branches = []
            for anchor_pos, branch_info in active_branch_list:
                if start_pos <= anchor_pos < end_pos:
                    # 调整锚点位置为相对于片段的位置
                    adjusted_anchor = anchor_pos - start_pos
                    frag_active_anchors.append(adjusted_anchor)
                    frag_active_branches.append(ActiveBranchInfo(
                        anchor_pos=adjusted_anchor,
                        branch_segment=branch_info.branch_segment,
                        branch_id=branch_info.branch_id
                    ))

            # 计算这个片段内的已去支化锚点
            frag_debranched_anchors = [
                anchor_pos - start_pos
                for anchor_pos in debranched_anchor_list
                if start_pos <= anchor_pos < end_pos
            ]

            # 检查是否有chimera在这个片段内
            has_chimera = False
            chimera_positions = []
            for j in mol.chimera_junctions:
                if start_pos <= j.position < end_pos:
                    has_chimera = True
                    chimera_positions.append(j.position - start_pos)

            frag_mol = LinearMolecule(
                molecule_id=f"{mol.instance_id}_trunk_frag{i}",
                segments=frag_segments,
                source_graph_id=mol.instance_id,
                is_from_branch=False,
                repeat_count=mol.repeat_count,
                source_ecc_length=mol.source_ecc.length,
                has_chimera=has_chimera,
                chimera_positions=chimera_positions,
                active_branch_anchors=frag_active_anchors,
                debranched_anchors=frag_debranched_anchors,
                active_branches=frag_active_branches,
                junction_covered_possible=compute_junction_covered_possible(frag_segments)
            )
            fragments.append(frag_mol)

        return fragments if fragments else [trunk_mol]

    def _extract_segments_in_range(
        self,
        mol: RCAMoleculeGraph,
        start_pos: int,
        length: int,
        ecc_db: Dict[str, EccDNA]
    ) -> List[Segment]:
        """
        从主干的指定范围提取segments

        处理chimera嵌入的情况

        Args:
            mol: RCA分子图
            start_pos: 起始位置
            length: 提取长度
            ecc_db: eccDNA字典

        Returns:
            segment列表
        """
        ecc = mol.source_ecc
        end_pos = start_pos + length
        segments = []

        # 如果没有chimera，直接创建简单segment
        if not mol.chimera_junctions:
            offset = start_pos % ecc.length
            segments.append(Segment(
                ecc_id=ecc.id,
                ecc_offset=offset,
                length=length,
                strand=Strand.FORWARD,
                segment_type=SegmentType.TRUNK,
                parent_offset=start_pos
            ))
            return segments

        # 有chimera，需要处理
        junctions = sorted(mol.chimera_junctions, key=lambda j: j.position)
        # 验证并移除重叠的junctions
        junctions = validate_non_overlapping_junctions(junctions)

        current_pos = start_pos

        for junction in junctions:
            chim_start = junction.position
            chim_end = junction.position + junction.donor_segment.length

            # chimera在范围之前，跳过
            if chim_end <= start_pos:
                continue

            # chimera在范围之后，结束
            if chim_start >= end_pos:
                break

            # 添加chimera之前的主干部分（在范围内的）
            if chim_start > current_pos:
                trunk_start = max(current_pos, start_pos)
                trunk_end = min(chim_start, end_pos)
                if trunk_end > trunk_start:
                    offset = trunk_start % ecc.length
                    segments.append(Segment(
                        ecc_id=ecc.id,
                        ecc_offset=offset,
                        length=trunk_end - trunk_start,
                        strand=Strand.FORWARD,
                        segment_type=SegmentType.TRUNK,
                        parent_offset=trunk_start
                    ))

            # 添加chimera部分（在范围内的）
            if chim_start < end_pos and chim_end > start_pos:
                # 计算chimera在范围内的部分
                chim_local_start = max(0, start_pos - chim_start)
                chim_local_end = min(junction.donor_segment.length, end_pos - chim_start)
                if chim_local_end > chim_local_start:
                    donor_seg = junction.donor_segment
                    # 验证donor ecc存在于数据库中
                    if donor_seg.ecc_id not in ecc_db:
                        # 跳过无效的chimera段，记录警告但不中断
                        continue
                    donor_ecc_len = ecc_db[donor_seg.ecc_id].length
                    # 确保offset在有效范围内
                    safe_offset = (donor_seg.ecc_offset + chim_local_start) % donor_ecc_len
                    # 注意：不截断length到donor_ecc_len，因为chimera可能跨越多个模板周期
                    # get_circular_substr() 会正确处理 length > template_length 的情况
                    actual_length = chim_local_end - chim_local_start
                    segments.append(Segment(
                        ecc_id=donor_seg.ecc_id,
                        ecc_offset=safe_offset,
                        length=actual_length,
                        strand=donor_seg.strand,
                        segment_type=SegmentType.CHIMERA,
                        parent_offset=chim_start + chim_local_start
                    ))

            current_pos = chim_end

        # 添加最后的主干部分
        if current_pos < end_pos:
            trunk_start = max(current_pos, start_pos)
            if trunk_start < end_pos:
                offset = trunk_start % ecc.length
                segments.append(Segment(
                    ecc_id=ecc.id,
                    ecc_offset=offset,
                    length=end_pos - trunk_start,
                    strand=Strand.FORWARD,
                    segment_type=SegmentType.TRUNK,
                    parent_offset=trunk_start
                ))

        return segments

    def linearize_molecule(
        self,
        mol: RCAMoleculeGraph,
        ecc_db: Dict[str, EccDNA],
        break_positions: List[int]
    ) -> List[LinearMolecule]:
        """
        将分子图线性化为可测序分子

        Args:
            mol: RCA分子图
            ecc_db: eccDNA字典
            break_positions: 断裂位置列表

        Returns:
            线性分子列表
        """
        linear_molecules = []
        ecc = mol.source_ecc

        # 1. 创建主干分子（可能被断裂切分）
        trunk_mol = self._create_trunk_molecule(mol, ecc_db)

        # 2. 如果有断裂，切分主干
        if break_positions:
            trunk_fragments = self.split_trunk_at_breaks(
                mol, trunk_mol, break_positions, ecc_db
            )
            linear_molecules.extend(trunk_fragments)
        else:
            linear_molecules.append(trunk_mol)

        # 3. 去支化后的独立分支分子
        for branch in mol.branches:
            if branch.debranched:
                branch_mol = LinearMolecule(
                    molecule_id=f"{mol.instance_id}_br{branch.branch_id}",
                    segments=[branch.segment],
                    source_graph_id=mol.instance_id,
                    is_from_branch=True,
                    repeat_count=mol.repeat_count,
                    source_ecc_length=mol.source_ecc.length,
                    junction_covered_possible=compute_junction_covered_possible([branch.segment])
                )
                linear_molecules.append(branch_mol)

        return linear_molecules

    def _create_trunk_molecule(
        self,
        mol: RCAMoleculeGraph,
        ecc_db: Dict[str, EccDNA]
    ) -> LinearMolecule:
        """
        创建主干线性分子（完整，未考虑断裂）

        Args:
            mol: RCA分子图
            ecc_db: eccDNA字典

        Returns:
            主干线性分子
        """
        ecc = mol.source_ecc

        # 收集未去支化的分支信息（用于假chimera模拟）
        active_anchors = []
        active_branches = []
        debranched_anchors = []
        for b in mol.branches:
            if not b.debranched:
                active_anchors.append(b.anchor_pos)
                active_branches.append(ActiveBranchInfo(
                    anchor_pos=b.anchor_pos,
                    branch_segment=b.segment,
                    branch_id=b.branch_id
                ))
            else:
                debranched_anchors.append(b.anchor_pos)

        # 如果没有chimera，直接创建简单的主干分子
        if not mol.chimera_junctions:
            trunk_segment = Segment(
                ecc_id=ecc.id,
                ecc_offset=0,
                length=mol.trunk_length,
                strand=Strand.FORWARD,
                segment_type=SegmentType.TRUNK
            )
            return LinearMolecule(
                molecule_id=f"{mol.instance_id}_trunk",
                segments=[trunk_segment],
                source_graph_id=mol.instance_id,
                is_from_branch=False,
                repeat_count=mol.repeat_count,
                source_ecc_length=ecc.length,
                has_chimera=False,
                active_branch_anchors=active_anchors,
                debranched_anchors=debranched_anchors,
                active_branches=active_branches,
                junction_covered_possible=compute_junction_covered_possible([trunk_segment])
            )

        # 有chimera，需要构建复合segments
        segments = self._build_chimera_segments(mol, ecc_db)

        chimera_positions = [j.position for j in mol.chimera_junctions]

        return LinearMolecule(
            molecule_id=f"{mol.instance_id}_trunk",
            segments=segments,
            source_graph_id=mol.instance_id,
            is_from_branch=False,
            repeat_count=mol.repeat_count,
            source_ecc_length=ecc.length,
            has_chimera=True,
            chimera_positions=chimera_positions,
            active_branch_anchors=active_anchors,
            debranched_anchors=debranched_anchors,
            active_branches=active_branches,
            junction_covered_possible=compute_junction_covered_possible(segments)
        )

    def _build_chimera_segments(
        self,
        mol: RCAMoleculeGraph,
        ecc_db: Dict[str, EccDNA]
    ) -> List[Segment]:
        """
        构建包含chimera插入的segment列表
        """
        ecc = mol.source_ecc

        # 按位置排序chimera junctions
        junctions = sorted(mol.chimera_junctions, key=lambda j: j.position)
        # 验证并移除重叠的junctions
        junctions = validate_non_overlapping_junctions(junctions)

        segments = []
        current_pos = 0

        for junction in junctions:
            # 添加chimera之前的主干片段
            if junction.position > current_pos:
                pre_len = junction.position - current_pos
                pre_offset = current_pos % ecc.length
                segments.append(Segment(
                    ecc_id=ecc.id,
                    ecc_offset=pre_offset,
                    length=pre_len,
                    strand=Strand.FORWARD,
                    segment_type=SegmentType.TRUNK,
                    parent_offset=current_pos
                ))

            # 添加chimera片段
            segments.append(junction.donor_segment)

            # 更新current_pos为chimera结束位置（修复：之前错误地使用了起始位置）
            current_pos = junction.position + junction.donor_segment.length

        # 添加最后一段主干
        if current_pos < mol.trunk_length:
            remaining = mol.trunk_length - current_pos
            final_offset = current_pos % ecc.length
            segments.append(Segment(
                ecc_id=ecc.id,
                ecc_offset=final_offset,
                length=remaining,
                strand=Strand.FORWARD,
                segment_type=SegmentType.TRUNK,
                parent_offset=current_pos
            ))

        return segments

    def process_all(
        self,
        molecules: List[RCAMoleculeGraph],
        ecc_db: Dict[str, EccDNA]
    ) -> Tuple[List[LinearMolecule], DebranchStats]:
        """
        处理所有分子的去支化和线性化

        Args:
            molecules: RCA分子图列表
            ecc_db: eccDNA字典

        Returns:
            (线性分子列表, 统计信息)
        """
        stats = DebranchStats()
        all_linear = []

        for mol in molecules:
            # 统计分支
            stats.total_branches += len(mol.branches)

            # 去支化
            self.debranch_molecule(mol)

            # 统计去支化结果
            for branch in mol.branches:
                if branch.debranched:
                    stats.debranched_count += 1
                else:
                    stats.retained_count += 1

            # 计算断裂位置
            break_positions, nick_count, break_count = self.calculate_break_positions(mol)
            stats.nick_count += nick_count
            stats.break_count += break_count

            # 线性化（考虑断裂）
            linear_mols = self.linearize_molecule(mol, ecc_db, break_positions)

            for lm in linear_mols:
                if lm.is_from_branch:
                    stats.linear_molecules_from_branch += 1
                else:
                    stats.linear_molecules_from_trunk += 1
                    # 统计因断裂产生的额外片段
                    if '_frag' in lm.molecule_id:
                        stats.trunk_fragments_from_breaks += 1

            all_linear.extend(linear_mols)

        return all_linear, stats


class LinearMoleculePool:
    """
    线性分子池

    管理所有可测序的线性分子，支持按权重抽样
    """

    def __init__(
        self,
        molecules: List[LinearMolecule],
        ecc_db: Dict[str, EccDNA],
        weight_by_length: bool = True
    ):
        self.molecules = molecules
        self.ecc_db = ecc_db
        self.weight_by_length = weight_by_length

        # 预计算权重
        self._compute_weights()

    def _compute_weights(self):
        """计算抽样权重"""
        # 处理空分子池
        if not self.molecules:
            self.weights = np.array([])
            return

        if self.weight_by_length:
            self.weights = np.array([m.total_length for m in self.molecules], dtype=float)
        else:
            self.weights = np.ones(len(self.molecules))

        # 防止除零：如果总权重为0，使用均匀分布
        weight_sum = self.weights.sum()
        if weight_sum > 0:
            self.weights /= weight_sum
        else:
            # 所有分子长度都为0的极端情况，使用均匀分布
            self.weights = np.ones(len(self.molecules)) / len(self.molecules)

    def sample(self, n: int, rng: np.random.Generator) -> List[LinearMolecule]:
        """
        按权重抽样分子

        Args:
            n: 抽样数量
            rng: 随机数生成器

        Returns:
            抽样的分子列表

        Raises:
            ValueError: 如果分子池为空
        """
        if not self.molecules:
            raise ValueError("Cannot sample from empty molecule pool")

        indices = rng.choice(len(self.molecules), size=n, replace=True, p=self.weights)
        return [self.molecules[i] for i in indices]

    def get_sequence(self, mol: LinearMolecule) -> str:
        """获取分子序列"""
        return mol.get_sequence(self.ecc_db)

    @property
    def total_molecules(self) -> int:
        return len(self.molecules)

    @property
    def total_length(self) -> int:
        return sum(m.total_length for m in self.molecules)
