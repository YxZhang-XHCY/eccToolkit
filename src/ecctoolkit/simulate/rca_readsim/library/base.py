"""
文库生成器基础模块

包含各平台共享的工具函数和基类。
"""

from typing import List, Optional
from abc import ABC, abstractmethod
import numpy as np

from ..models import (
    LinearMolecule, SequencedRead, Segment,
    SegmentType, get_ecc_length
)
from ..config import SimConfig
from ..debranch import LinearMoleculePool


def extract_segments(
    mol: LinearMolecule,
    start: int,
    length: int
) -> List[Segment]:
    """
    从分子中提取指定区域的segments

    Args:
        mol: 线性分子
        start: 起始位置
        length: 提取长度

    Returns:
        提取的Segment列表
    """
    result = []
    current_pos = 0
    remaining = length
    extract_start = start

    for seg in mol.segments:
        seg_end = current_pos + seg.length

        if seg_end <= extract_start:
            current_pos = seg_end
            continue

        local_start = max(0, extract_start - current_pos)
        local_end = min(seg.length, local_start + remaining)
        local_length = local_end - local_start

        if local_length > 0:
            ecc_len = get_ecc_length(seg.ecc_id)

            new_seg = Segment(
                ecc_id=seg.ecc_id,
                ecc_offset=(seg.ecc_offset + local_start) % ecc_len,
                length=local_length,
                strand=seg.strand,
                segment_type=seg.segment_type,
                parent_offset=current_pos + local_start
            )
            result.append(new_seg)
            remaining -= local_length
            extract_start = seg_end

        current_pos = seg_end

        if remaining <= 0:
            break

    return result


def compute_repeat_count_by_source(segments: List[Segment]) -> dict:
    """
    计算read中各eccDNA来源的repeat比例

    Returns:
        {ecc_id: repeat_float}
    """
    if not segments:
        return {}

    def _collect_lengths(allowed_types: set) -> dict:
        lengths = {}
        for seg in segments:
            if seg.segment_type in allowed_types:
                lengths[seg.ecc_id] = lengths.get(seg.ecc_id, 0) + seg.length
        return lengths

    lengths = _collect_lengths({SegmentType.TRUNK, SegmentType.CHIMERA})
    if not lengths:
        lengths = _collect_lengths({SegmentType.BRANCH})

    repeat_by_ecc = {}
    for ecc_id, total_len in lengths.items():
        ecc_len = get_ecc_length(ecc_id)
        if ecc_len > 0:
            repeat_by_ecc[ecc_id] = total_len / ecc_len
    return repeat_by_ecc


def compute_repeat_count(segments: List[Segment]) -> int:
    """
    计算read覆盖的repeat次数（取各来源的最大整数repeat）
    """
    repeat_by_ecc = compute_repeat_count_by_source(segments)
    if not repeat_by_ecc:
        return 0

    max_repeat = 1
    for repeat in repeat_by_ecc.values():
        max_repeat = max(max_repeat, max(1, int(repeat)))

    return max_repeat


class BaseLibraryGenerator(ABC):
    """文库生成器基类"""

    def __init__(
        self,
        config: SimConfig,
        rng: Optional[np.random.Generator] = None
    ):
        self.config = config
        self.rng = rng if rng is not None else np.random.default_rng()
        self._read_counter = 0

    def _extract_segments(
        self,
        mol: LinearMolecule,
        start: int,
        length: int
    ) -> List[Segment]:
        """提取segments（委托给模块函数）"""
        return extract_segments(mol, start, length)

    def _compute_repeat_count(self, segments: List[Segment]) -> int:
        """计算repeat次数（委托给模块函数）"""
        return compute_repeat_count(segments)

    def _compute_repeat_count_by_source(self, segments: List[Segment]) -> dict:
        """计算各eccDNA来源repeat比例（委托给模块函数）"""
        return compute_repeat_count_by_source(segments)

    @abstractmethod
    def generate(
        self,
        pool: LinearMoleculePool,
        target_reads: int
    ):
        """生成reads（子类实现）"""
        pass
