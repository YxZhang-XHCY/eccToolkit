"""
NGS（Illumina）文库生成模块

实现§9.2：NGS文库生成
- insert size采样
- paired-end read生成
- 支链导致的假chimera模拟
- truth追踪
"""

from typing import List, Dict, Tuple, Optional
import numpy as np
from dataclasses import dataclass

from ..models import (
    LinearMolecule, SequencedRead, Segment, EccDNA,
    Strand, SegmentType, ActiveBranchInfo, get_ecc_length
)
from ..config import SimConfig
from ..debranch import LinearMoleculePool
from ..seq_utils import reverse_complement, simulate_quality_decay
from ..error_models import IlluminaErrorModel
from .base import extract_segments, compute_repeat_count, compute_repeat_count_by_source


@dataclass
class NGSLibraryStats:
    """NGS文库统计"""
    total_reads: int = 0
    total_pairs: int = 0
    insert_size_mean: float = 0.0
    insert_size_std: float = 0.0
    reads_from_trunk: int = 0
    reads_from_branch: int = 0
    reads_with_inter_chimera: int = 0      # 跨eccDNA的chimera
    reads_with_branch_chimera: int = 0     # 支链导致的假chimera

    def summary(self) -> str:
        return (
            f"NGS: {self.total_reads} reads ({self.total_pairs} pairs), "
            f"Insert: {self.insert_size_mean:.0f}±{self.insert_size_std:.0f}bp, "
            f"Inter-chimera: {self.reads_with_inter_chimera}, "
            f"Branch-chimera: {self.reads_with_branch_chimera}"
        )


class NGSLibraryGenerator:
    """NGS文库生成器"""

    def __init__(
        self,
        config: SimConfig,
        rng: Optional[np.random.Generator] = None
    ):
        self.config = config
        self.rng = rng if rng is not None else np.random.default_rng()
        self._read_counter = 0
        self._ecc_db: Optional[Dict[str, EccDNA]] = None

        # 初始化Illumina错误模型
        self.error_model = IlluminaErrorModel(
            base_error_rate=0.001,  # 5'端错误率 0.1%
            end_error_rate=0.01,    # 3'端错误率 1%
            rng=self.rng
        )

    def set_ecc_db(self, ecc_db: Dict[str, EccDNA]):
        """设置eccDNA数据库，用于获取分支序列"""
        self._ecc_db = ecc_db

    def sample_insert_size(self) -> int:
        """采样insert size"""
        I_mu = self.config.I_mu
        I_sigma = self.config.I_sigma
        read_length = self.config.constants.read_length
        min_insert = 2 * read_length

        insert = self.rng.normal(I_mu, I_sigma)
        return max(min_insert, int(round(insert)))

    def check_branch_chimera(
        self,
        mol: LinearMolecule,
        insert_start: int,
        insert_size: int
    ) -> Optional[Tuple[ActiveBranchInfo, int]]:
        """
        检查insert是否会在分支锚点处产生假chimera

        Args:
            mol: 线性分子
            insert_start: insert起始位置
            insert_size: insert大小

        Returns:
            如果产生假chimera，返回 (分支信息, 断点在insert中的位置)
            否则返回 None
        """
        if not mol.active_branches:
            return None

        p_chimera = self.config.constants.p_branch_chimera
        window = self.config.constants.branch_chimera_window
        insert_end = insert_start + insert_size

        for branch_info in mol.active_branches:
            anchor = branch_info.anchor_pos

            # 检查锚点是否在insert范围内或窗口内
            # 窗口扩展：insert边界向外扩展window bp
            # 使用<=确保边界条件正确
            in_range = (insert_start - window) <= anchor <= (insert_end + window)

            if in_range:
                # 锚点在扩展范围内，有概率产生chimera
                if self.rng.random() < p_chimera:
                    # 断点位置：如果锚点在insert外，clamp到边界
                    if anchor < insert_start:
                        breakpoint_in_insert = 0  # 在insert起始处断开
                    elif anchor >= insert_end:
                        breakpoint_in_insert = insert_size  # 在insert末端断开
                    else:
                        breakpoint_in_insert = anchor - insert_start
                    return (branch_info, breakpoint_in_insert)

        return None

    def generate_branch_chimera_insert(
        self,
        mol: LinearMolecule,
        mol_seq: str,
        insert_start: int,
        insert_size: int,
        branch_info: ActiveBranchInfo,
        breakpoint: int
    ) -> Tuple[str, List[Segment], int]:
        """
        生成包含分支的假chimera insert序列

        打断可能发生在锚点的两侧，产生两种情况：
        1. 主干前半部分 + 分支后半部分
        2. 分支前半部分 + 主干后半部分

        Args:
            mol: 线性分子
            mol_seq: 分子序列
            insert_start: insert起始
            insert_size: insert大小
            branch_info: 分支信息
            breakpoint: 断点在insert中的位置

        Returns:
            (chimera序列, segments列表, 断点位置)
        """
        if self._ecc_db is None:
            # 没有ecc_db，无法获取分支序列，返回原序列
            orig_seq = mol_seq[insert_start:insert_start + insert_size]
            orig_segments = self._extract_segments(mol, insert_start, insert_size)
            return orig_seq, orig_segments, -1

        # 获取分支序列
        branch_seg = branch_info.branch_segment
        ecc = self._ecc_db.get(branch_seg.ecc_id)
        if ecc is None:
            orig_seq = mol_seq[insert_start:insert_start + insert_size]
            orig_segments = self._extract_segments(mol, insert_start, insert_size)
            return orig_seq, orig_segments, -1

        branch_seq = ecc.get_circular_substr(branch_seg.ecc_offset, branch_seg.length)
        if branch_seg.strand == Strand.REVERSE:
            branch_seq = reverse_complement(branch_seq)

        # 决定chimera的结构（随机选择方向）
        # 50%概率：主干在前，分支在后
        # 50%概率：分支在前，主干在后
        trunk_first = self.rng.random() < 0.5

        if trunk_first:
            # 主干前 + 分支后
            trunk_len = breakpoint
            branch_len = insert_size - breakpoint

            trunk_part = mol_seq[insert_start:insert_start + trunk_len]
            # 从分支序列中取后半部分
            branch_start = min(len(branch_seq) // 2, branch_len)
            branch_part = branch_seq[branch_start:branch_start + branch_len]

            # 如果分支部分不够长，用分支序列填充
            while len(branch_part) < branch_len and len(branch_seq) > 0:
                branch_part += branch_seq[:branch_len - len(branch_part)]

            chimera_seq = trunk_part + branch_part[:branch_len]

            # 构建segments
            segments = []
            # 主干部分
            trunk_segments = self._extract_segments(mol, insert_start, trunk_len)
            for seg in trunk_segments:
                if seg.parent_offset is not None:
                    seg.parent_offset = max(0, seg.parent_offset - insert_start)
            segments.extend(trunk_segments)
            # 分支部分
            if branch_len > 0:
                segments.append(Segment(
                    ecc_id=branch_seg.ecc_id,
                    ecc_offset=(branch_seg.ecc_offset + branch_start) % ecc.length,
                    length=min(branch_len, len(branch_part)),
                    strand=branch_seg.strand,
                    segment_type=SegmentType.BRANCH,
                    parent_offset=trunk_len
                ))
        else:
            # 分支前 + 主干后
            branch_len = breakpoint
            trunk_len = insert_size - breakpoint

            # 从分支序列中取前半部分
            branch_part = branch_seq[:branch_len]
            while len(branch_part) < branch_len and len(branch_seq) > 0:
                branch_part += branch_seq[:branch_len - len(branch_part)]

            trunk_start = insert_start + breakpoint
            trunk_part = mol_seq[trunk_start:trunk_start + trunk_len]

            chimera_seq = branch_part[:branch_len] + trunk_part

            # 构建segments
            segments = []
            # 分支部分
            if branch_len > 0:
                segments.append(Segment(
                    ecc_id=branch_seg.ecc_id,
                    ecc_offset=branch_seg.ecc_offset,
                    length=min(branch_len, len(branch_part)),
                    strand=branch_seg.strand,
                    segment_type=SegmentType.BRANCH,
                    parent_offset=0
                ))
            # 主干部分
            trunk_segments = self._extract_segments(mol, trunk_start, trunk_len)
            for seg in trunk_segments:
                # 修复：确保parent_offset不为负值
                # seg.parent_offset是segment在原分子中的位置，需要转换为在chimera insert中的相对位置
                relative_pos = max(0, seg.parent_offset - trunk_start)
                seg.parent_offset = branch_len + relative_pos
            segments.extend(trunk_segments)

        return chimera_seq[:insert_size], segments, breakpoint

    def generate_paired_reads(
        self,
        mol: LinearMolecule,
        mol_seq: str,
        insert_start: int,
        insert_size: int,
        branch_chimera: Optional[Tuple[str, List[Segment], int]] = None
    ) -> Tuple[SequencedRead, SequencedRead]:
        """
        从insert生成paired-end reads

        Args:
            mol: 线性分子
            mol_seq: 分子序列
            insert_start: insert在分子中的起始位置
            insert_size: insert大小
            branch_chimera: 如果有分支chimera，(序列, segments, 断点)

        Returns:
            (R1 read, R2 read)
        """
        self._read_counter += 1
        read_id_base = f"read_{self._read_counter}"
        read_length = self.config.constants.read_length

        # 确定insert序列和segments
        has_branch_chimera = False
        branch_breakpoint = -1

        if branch_chimera is not None:
            insert_seq, all_segments, branch_breakpoint = branch_chimera
            has_branch_chimera = branch_breakpoint >= 0
        else:
            insert_seq = mol_seq[insert_start:insert_start + insert_size]
            all_segments = None

        # R1: insert前端
        r1_seq = insert_seq[:read_length]
        # 应用Illumina错误模型（替换错误 + 位置相关质量）
        r1_seq, r1_qual = self.error_model.apply(r1_seq)

        # R2: insert后端反向互补
        r2_seq_raw = insert_seq[-read_length:] if len(insert_seq) >= read_length else insert_seq
        r2_seq = reverse_complement(r2_seq_raw)
        # 应用Illumina错误模型
        r2_seq, r2_qual = self.error_model.apply(r2_seq)

        # 计算segments（用于truth追踪）
        if all_segments is not None:
            # 使用chimera的segments
            r1_segments = self._extract_from_segment_list(all_segments, 0, read_length)
            r2_start = insert_size - read_length
            r2_segments = self._extract_from_segment_list(all_segments, max(0, r2_start), read_length)
        else:
            r1_segments = self._extract_segments(mol, insert_start, read_length)
            r2_start = insert_start + insert_size - read_length
            r2_segments = self._extract_segments(mol, max(0, r2_start), read_length)

        # R2的segments需要反向
        for seg in r2_segments:
            seg.strand = Strand.REVERSE if seg.strand == Strand.FORWARD else Strand.FORWARD

        # 计算repeat count
        ecc_ids = list(set(s.ecc_id for s in r1_segments + r2_segments))

        # 检查是否覆盖inter-eccDNA chimera
        has_inter_chimera = False
        chimera_bps = []
        if mol.has_chimera:
            for pos in mol.chimera_positions:
                if insert_start <= pos < insert_start + insert_size:
                    has_inter_chimera = True
                    chimera_bps.append(pos - insert_start)

        # 检查是否覆盖未去支化分支锚点
        covered_anchor = False
        for anchor in mol.active_branch_anchors:
            if insert_start <= anchor < insert_start + insert_size:
                covered_anchor = True
                break

        chimera_breakpoints = list(chimera_bps)
        if has_branch_chimera and branch_breakpoint >= 0:
            if branch_breakpoint not in chimera_breakpoints:
                chimera_breakpoints.append(branch_breakpoint)

        # 创建reads
        r1 = SequencedRead(
            read_id=f"{read_id_base}/1",
            platform="NGS",
            sequence=r1_seq,
            quality=r1_qual,
            source_molecule_id=mol.molecule_id,
            source_ecc_ids=ecc_ids,
            segments=r1_segments,
            repeat_count_truth=self._compute_repeat_count(r1_segments, mol),
            repeat_count_by_source=self._compute_repeat_count_by_source(r1_segments),
            has_inter_chimera=has_inter_chimera,
            chimera_breakpoints=chimera_breakpoints,
            covered_branch_anchor=covered_anchor,
            has_branch_chimera=has_branch_chimera,
            is_paired=True,
            mate_id=f"{read_id_base}/2",
            read_number=1,
            insert_size=insert_size,
            junction_covered_possible=mol.junction_covered_possible,
            # Background DNA fields
            is_background=mol.is_background,
            background_chrom=mol.background_chrom,
            background_start=mol.background_start,
            background_end=mol.background_end
        )

        r2 = SequencedRead(
            read_id=f"{read_id_base}/2",
            platform="NGS",
            sequence=r2_seq,
            quality=r2_qual,
            source_molecule_id=mol.molecule_id,
            source_ecc_ids=ecc_ids,
            segments=r2_segments,
            repeat_count_truth=self._compute_repeat_count(r2_segments, mol),
            repeat_count_by_source=self._compute_repeat_count_by_source(r2_segments),
            has_inter_chimera=has_inter_chimera,
            chimera_breakpoints=chimera_breakpoints,
            covered_branch_anchor=covered_anchor,
            has_branch_chimera=has_branch_chimera,
            is_paired=True,
            mate_id=f"{read_id_base}/1",
            read_number=2,
            insert_size=insert_size,
            junction_covered_possible=mol.junction_covered_possible,
            # Background DNA fields
            is_background=mol.is_background,
            background_chrom=mol.background_chrom,
            background_start=mol.background_start,
            background_end=mol.background_end
        )

        return r1, r2

    def _extract_from_segment_list(
        self,
        segments: List[Segment],
        start: int,
        length: int
    ) -> List[Segment]:
        """从segment列表中提取指定范围"""
        result = []
        current_pos = 0

        for seg in segments:
            # 修复：使用 is not None 而非 falsy 检查，因为 parent_offset=0 是有效值
            seg_start = seg.parent_offset if hasattr(seg, 'parent_offset') and seg.parent_offset is not None else current_pos
            seg_end = seg_start + seg.length

            # 检查是否与目标范围重叠
            if seg_end <= start:
                current_pos = seg_end
                continue
            if seg_start >= start + length:
                break

            # 计算重叠部分
            overlap_start = max(seg_start, start)
            overlap_end = min(seg_end, start + length)
            overlap_len = overlap_end - overlap_start

            if overlap_len > 0:
                local_offset = overlap_start - seg_start
                ecc_len = get_ecc_length(seg.ecc_id)

                new_seg = Segment(
                    ecc_id=seg.ecc_id,
                    ecc_offset=(seg.ecc_offset + local_offset) % ecc_len,
                    length=overlap_len,
                    strand=seg.strand,
                    segment_type=seg.segment_type,
                    parent_offset=overlap_start - start
                )
                result.append(new_seg)

            current_pos = seg_end

        return result

    def _extract_segments(
        self,
        mol: LinearMolecule,
        start: int,
        length: int
    ) -> List[Segment]:
        """从分子中提取指定区域的segments（委托给共享函数）"""
        return extract_segments(mol, start, length)

    def _compute_repeat_count(
        self,
        segments: List[Segment],
        mol: LinearMolecule
    ) -> int:
        """计算read覆盖的repeat次数（委托给共享函数）"""
        return compute_repeat_count(segments)

    def _compute_repeat_count_by_source(
        self,
        segments: List[Segment]
    ) -> dict:
        """计算各eccDNA来源repeat比例（委托给共享函数）"""
        return compute_repeat_count_by_source(segments)

    def generate(
        self,
        pool: LinearMoleculePool,
        target_reads: int
    ) -> Tuple[List[SequencedRead], NGSLibraryStats]:
        """
        生成NGS reads

        Args:
            pool: 线性分子池
            target_reads: 目标read数

        Returns:
            (reads列表, 统计信息)
        """
        # 设置ecc_db
        self._ecc_db = pool.ecc_db

        stats = NGSLibraryStats()
        reads = []
        insert_sizes = []

        target_pairs = target_reads // 2

        # 防止无限循环：最多尝试 target * 10 次
        max_attempts = target_pairs * 10
        attempts = 0
        skipped_too_short = 0

        while stats.total_pairs < target_pairs and attempts < max_attempts:
            attempts += 1

            # 抽样分子
            mol = pool.sample(1, self.rng)[0]
            mol_seq = pool.get_sequence(mol)
            mol_len = len(mol_seq)

            # 采样insert
            insert_size = self.sample_insert_size()

            # 检查分子是否足够长
            if mol_len < insert_size:
                skipped_too_short += 1
                continue

            # 随机选择insert起点
            max_start = mol_len - insert_size
            insert_start = self.rng.integers(0, max_start + 1)

            # 检查是否产生分支chimera
            branch_result = self.check_branch_chimera(mol, insert_start, insert_size)

            if branch_result is not None:
                branch_info, breakpoint = branch_result
                branch_chimera = self.generate_branch_chimera_insert(
                    mol, mol_seq, insert_start, insert_size,
                    branch_info, breakpoint
                )
            else:
                branch_chimera = None

            # 生成paired reads
            r1, r2 = self.generate_paired_reads(
                mol, mol_seq, insert_start, insert_size, branch_chimera
            )
            reads.extend([r1, r2])

            # 更新统计
            stats.total_pairs += 1
            stats.total_reads += 2
            insert_sizes.append(insert_size)

            if mol.is_from_branch:
                stats.reads_from_branch += 2
            else:
                stats.reads_from_trunk += 2

            if r1.has_inter_chimera:
                stats.reads_with_inter_chimera += 2

            if r1.has_branch_chimera:
                stats.reads_with_branch_chimera += 2

        # 检查是否因达到最大尝试次数而退出
        if attempts >= max_attempts and stats.total_pairs < target_pairs:
            import warnings
            warnings.warn(
                f"NGS sampling reached max attempts ({max_attempts}). "
                f"Generated {stats.total_pairs}/{target_pairs} pairs. "
                f"Skipped {skipped_too_short} molecules too short for insert size. "
                f"Consider reducing I_mu or checking molecule pool lengths."
            )

        # 计算insert统计
        if insert_sizes:
            stats.insert_size_mean = np.mean(insert_sizes)
            stats.insert_size_std = np.std(insert_sizes)

        return reads, stats
