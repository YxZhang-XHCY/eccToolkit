"""
Nanopore文库生成模块

实现§9.3：ONT文库生成
- 全长分子测序
- 未去支化分支导致的截断
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
from ..error_models import NanoporeErrorModel
from .base import extract_segments, compute_repeat_count, compute_repeat_count_by_source


@dataclass
class ONTLibraryStats:
    """ONT文库统计"""
    total_reads: int = 0
    mean_read_length: float = 0.0
    max_read_length: int = 0
    truncated_reads: int = 0
    branch_truncated: int = 0      # 分支锚点导致的截断
    pore_clogged: int = 0          # 孔道堵塞导致的截断
    reads_from_trunk: int = 0
    reads_from_branch: int = 0
    reads_with_chimera: int = 0

    def summary(self) -> str:
        return (
            f"ONT: {self.total_reads} reads, "
            f"Length: {self.mean_read_length:.0f}bp (max {self.max_read_length}), "
            f"Truncated: {self.truncated_reads} (branch: {self.branch_truncated}, pore: {self.pore_clogged}), "
            f"Chimera: {self.reads_with_chimera}"
        )


class ONTLibraryGenerator:
    """Nanopore文库生成器"""
    
    def __init__(
        self,
        config: SimConfig,
        rng: Optional[np.random.Generator] = None
    ):
        self.config = config
        self.rng = rng if rng is not None else np.random.default_rng()
        self._read_counter = 0

        # 初始化Nanopore错误模型（高错误率 + 同聚物偏差）
        self.error_model = NanoporeErrorModel(
            substitution_rate=0.03,        # 3% 替换
            insertion_rate=0.02,           # 2% 插入
            deletion_rate=0.03,            # 3% 删除
            homopolymer_error_rate=0.15,   # 同聚物区域额外15%
            rng=self.rng
        )
    
    def check_truncation(
        self,
        mol: LinearMolecule,
        read_start: int,
        read_length: int
    ) -> Tuple[bool, int, str]:
        """
        检查是否因分支锚点而截断（解耦版本）

        - 未去支化锚点：高截断概率 (p_trunc_per_anchor)
        - 已去支化锚点：低残余截断概率 (p_trunc_debranched)

        这解耦了"是否去支化"与"是否截断"的强关联。

        Args:
            mol: 线性分子
            read_start: 读取起点
            read_length: 原始读取长度

        Returns:
            (是否截断, 截断后长度, 截断原因)
        """
        p_trunc_active = self.config.constants.p_trunc_per_anchor
        p_trunc_debranched = self.config.constants.p_trunc_debranched

        # 合并所有潜在截断点：(位置, 概率, 类型)
        truncation_candidates = []

        for anchor in mol.active_branch_anchors:
            if read_start <= anchor < read_start + read_length:
                truncation_candidates.append((anchor, p_trunc_active, "active_anchor"))

        for anchor in mol.debranched_anchors:
            if read_start <= anchor < read_start + read_length:
                truncation_candidates.append((anchor, p_trunc_debranched, "debranched_anchor"))

        if not truncation_candidates:
            return False, read_length, ""

        # 按位置排序，检查每个潜在截断点
        truncation_candidates.sort(key=lambda x: x[0])

        for anchor, p_trunc, anchor_type in truncation_candidates:
            if self.rng.random() < p_trunc:
                truncated_len = anchor - read_start
                if truncated_len >= self.config.constants.ont_min_length:
                    return True, truncated_len, f"{anchor_type}@{anchor}"
                # 太短了，继续检查下一个锚点

        return False, read_length, ""

    def check_pore_clog(self, read_length: int) -> Tuple[bool, int, str]:
        """
        检查是否因孔道堵塞而提前终止（与分支截断独立）

        使用指数分布模型：P(survive to L) = exp(-lambda_clog * L)
        从指数分布采样截断位置。

        Args:
            read_length: 原始读取长度

        Returns:
            (是否截断, 截断后长度, 截断原因)
        """
        lambda_clog = self.config.constants.lambda_clog

        if lambda_clog <= 0:
            return False, read_length, ""

        # 从指数分布采样截断位置
        # 指数分布的参数是 1/lambda (scale)
        clog_position = int(self.rng.exponential(1.0 / lambda_clog))

        if clog_position >= read_length:
            # 没有发生堵塞
            return False, read_length, ""

        # 发生堵塞
        if clog_position < self.config.constants.ont_min_length:
            # 堵塞位置太近，产生的读长太短
            return True, 0, "pore_clog_too_short"

        return True, clog_position, f"pore_clog@{clog_position}"

    def generate_read(
        self,
        mol: LinearMolecule,
        mol_seq: str,
        read_start: int = 0,
        read_length: Optional[int] = None
    ) -> Optional[SequencedRead]:
        """
        从分子生成ONT read
        
        Args:
            mol: 线性分子
            mol_seq: 分子序列
            read_start: 读取起点（0-based）
            read_length: 读取长度（None表示读到末端）
        
        Returns:
            SequencedRead或None（如果太短）
        """
        mol_len = len(mol_seq)
        min_length = self.config.constants.ont_min_length

        if read_start < 0 or read_start >= mol_len:
            return None

        if read_length is None:
            read_length = mol_len - read_start
        else:
            read_length = min(read_length, mol_len - read_start)

        if read_length < min_length:
            return None
        
        self._read_counter += 1
        read_id = f"ont_read_{self._read_counter}"
        
        # ONT通常读取全长（可能从任意位置开始）
        # 默认从头开始读取，可显式指定read_start/read_length

        # 检查分支锚点截断
        branch_trunc, branch_len, branch_reason = self.check_truncation(
            mol, read_start, read_length
        )

        # 检查孔道堵塞（独立机制）
        clog_trunc, clog_len, clog_reason = self.check_pore_clog(read_length)

        # 应用先发生的截断（取较短的结果）
        truncated = False
        final_length = read_length
        trunc_reason = ""

        if branch_trunc and clog_trunc:
            # 两种机制都触发 - 取较短的
            if branch_len <= clog_len:
                final_length = branch_len
                truncated = True
                trunc_reason = branch_reason
            else:
                final_length = clog_len
                truncated = True
                trunc_reason = clog_reason
        elif branch_trunc:
            final_length = branch_len
            truncated = True
            trunc_reason = branch_reason
        elif clog_trunc:
            final_length = clog_len
            truncated = True
            trunc_reason = clog_reason

        if final_length < min_length:
            return None
        
        # 提取序列
        read_seq = mol_seq[read_start:read_start + final_length]

        # 应用Nanopore错误模型（高indel率 + 同聚物错误）
        read_seq, read_qual = self.error_model.apply(read_seq)
        
        # 提取segments
        segments = self._extract_segments(mol, read_start, final_length)
        
        # 收集来源eccDNA
        ecc_ids = list(set(s.ecc_id for s in segments))
        
        # 检查chimera
        has_chimera = False
        chimera_bps = []
        if mol.has_chimera:
            for pos in mol.chimera_positions:
                if read_start <= pos < read_start + final_length:
                    has_chimera = True
                    chimera_bps.append(pos - read_start)
        
        # 检查是否覆盖分支锚点
        covered_anchor = any(
            read_start <= a < read_start + final_length
            for a in mol.active_branch_anchors
        )
        
        return SequencedRead(
            read_id=read_id,
            platform="ONT",
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
            truncated=truncated,
            truncation_reason=trunc_reason,
            junction_covered_possible=mol.junction_covered_possible,
            # Background DNA fields
            is_background=mol.is_background,
            background_chrom=mol.background_chrom,
            background_start=mol.background_start,
            background_end=mol.background_end
        )
    
    def _extract_segments(
        self,
        mol: LinearMolecule,
        start: int,
        length: int
    ) -> List[Segment]:
        """提取segments（委托给共享函数）"""
        return extract_segments(mol, start, length)

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
    ) -> Tuple[List[SequencedRead], ONTLibraryStats]:
        """
        生成ONT reads
        
        Args:
            pool: 线性分子池
            target_reads: 目标read数
        
        Returns:
            (reads列表, 统计信息)
        """
        stats = ONTLibraryStats()
        reads = []
        read_lengths = []
        
        attempts = 0
        max_attempts = target_reads * 3  # 防止无限循环
        
        while stats.total_reads < target_reads and attempts < max_attempts:
            attempts += 1
            
            # 抽样分子
            mol = pool.sample(1, self.rng)[0]
            mol_seq = pool.get_sequence(mol)
            
            # 生成read
            read = self.generate_read(mol, mol_seq)
            
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
            
            if read.truncated:
                stats.truncated_reads += 1
                # 区分截断类型
                if "branch_anchor" in read.truncation_reason:
                    stats.branch_truncated += 1
                elif "pore_clog" in read.truncation_reason:
                    stats.pore_clogged += 1

            if read.has_inter_chimera:
                stats.reads_with_chimera += 1

        # 检查是否因达到最大尝试次数而退出
        if attempts >= max_attempts and stats.total_reads < target_reads:
            import warnings
            warnings.warn(
                f"ONT sampling reached max attempts ({max_attempts}). "
                f"Generated {stats.total_reads}/{target_reads} reads. "
                f"Consider checking molecule pool or increasing max_attempts."
            )

        # 计算长度统计
        if read_lengths:
            stats.mean_read_length = np.mean(read_lengths)
            stats.max_read_length = max(read_lengths)
        
        return reads, stats
