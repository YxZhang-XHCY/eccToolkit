"""
核心数据结构定义

设计原则：
1. 全程维护坐标映射，支持高精度truth追踪
2. 不可变数据用dataclass(frozen=True)
3. 复杂结构用普通dataclass便于构建
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict, Any
from enum import Enum
import numpy as np


class Strand(Enum):
    """链方向"""
    FORWARD = "+"
    REVERSE = "-"


class SegmentType(Enum):
    """片段类型"""
    TRUNK = "trunk"
    BRANCH = "branch"
    CHIMERA = "chimera"
    BACKGROUND = "background"  # 背景线性DNA


# =============================================================================
# 输入数据结构
# =============================================================================

@dataclass
class EccDNA:
    """eccDNA模板"""
    id: str
    seq: str
    weight: float = 1.0
    
    @property
    def length(self) -> int:
        return len(self.seq)
    
    def get_circular_substr(self, offset: int, length: int) -> str:
        """从环状模板提取子序列（自动绕环）"""
        offset = offset % self.length
        if offset + length <= self.length:
            return self.seq[offset:offset + length]
        else:
            # 需要绕环
            result = self.seq[offset:]
            remaining = length - len(result)
            full_copies = remaining // self.length
            result += self.seq * full_copies
            result += self.seq[:remaining % self.length]
            return result


# =============================================================================
# RCA分子图结构
# =============================================================================

@dataclass
class Segment:
    """
    序列片段的精确坐标描述
    用于truth追踪的核心单元
    """
    ecc_id: str                    # 来源eccDNA
    ecc_offset: int                # 在eccDNA上的起始offset (0-based)
    length: int                    # 片段长度
    strand: Strand = Strand.FORWARD
    segment_type: SegmentType = SegmentType.TRUNK
    
    # 可选：在更大结构中的位置信息
    parent_offset: Optional[int] = None  # 在父结构（如主干）中的位置
    
    def to_dict(self) -> Dict[str, Any]:
        """转换为可序列化的dict"""
        return {
            "ecc_id": self.ecc_id,
            "ecc_offset": int(self.ecc_offset),
            "length": int(self.length),
            "strand": self.strand.value,
            "segment_type": self.segment_type.value
        }


@dataclass
class BranchNode:
    """分支节点"""
    branch_id: int                 # 分支编号
    anchor_pos: int                # 在主干上的锚点位置 (0-based)
    segment: Segment               # 分支序列信息
    debranched: bool = False       # 是否已去支化
    
    @property
    def length(self) -> int:
        return self.segment.length


@dataclass
class ChimeraJunction:
    """跨分子嵌合断点"""
    position: int                  # 在受体主干上的断点位置
    donor_ecc_id: str              # 供体eccDNA
    donor_segment: Segment         # 插入的供体片段
    

@dataclass
class RCAMoleculeGraph:
    """
    RCA分子图：主干 + 分支 + 嵌合
    
    这是compartment内单个eccDNA实例扩增后的完整结构表示
    """
    instance_id: str               # 唯一实例ID
    source_ecc: EccDNA             # 来源eccDNA模板
    compartment_id: int            # 所属compartment
    
    # 主干信息
    repeat_count: int              # 重复次数R
    trunk_length: int              # 主干总长度 = R * L
    
    # 分支列表
    branches: List[BranchNode] = field(default_factory=list)
    
    # 嵌合事件
    chimera_junctions: List[ChimeraJunction] = field(default_factory=list)
    
    def get_trunk_segment(self) -> Segment:
        """获取主干的segment描述"""
        return Segment(
            ecc_id=self.source_ecc.id,
            ecc_offset=0,
            length=self.trunk_length,
            strand=Strand.FORWARD,
            segment_type=SegmentType.TRUNK
        )
    
    def get_active_branches(self) -> List[BranchNode]:
        """获取未去支化的分支"""
        return [b for b in self.branches if not b.debranched]
    
    def get_debranched_molecules(self) -> List['LinearMolecule']:
        """获取去支化后产生的独立分子"""
        molecules = []
        for b in self.branches:
            if b.debranched:
                mol = LinearMolecule(
                    molecule_id=f"{self.instance_id}_branch_{b.branch_id}",
                    segments=[b.segment],
                    source_graph_id=self.instance_id,
                    is_from_branch=True
                )
                molecules.append(mol)
        return molecules


# =============================================================================
# 线性分子（可测序）
# =============================================================================

@dataclass
class ActiveBranchInfo:
    """未去支化分支的信息，用于模拟假chimera"""
    anchor_pos: int           # 在主干上的锚点位置
    branch_segment: Segment   # 分支的segment信息
    branch_id: int            # 分支ID


@dataclass
class LinearMolecule:
    """
    线性化后的可测序分子

    由一个或多个Segment组成，保留完整的来源追踪
    """
    molecule_id: str
    segments: List[Segment]        # 按顺序排列的片段
    source_graph_id: str           # 来源RCA分子图ID
    is_from_branch: bool = False   # 是否来自去支化的分支

    # RCA扩增信息（用于正确的采样权重）
    repeat_count: int = 1          # 来源RCA分子的重复次数（小环更高）
    source_ecc_length: int = 0     # 来源eccDNA的模板长度

    # 背景DNA标记
    is_background: bool = False    # 是否为背景线性DNA（非eccDNA来源）
    background_chrom: Optional[str] = None   # 背景DNA来源染色体
    background_start: Optional[int] = None   # 背景DNA在染色体上的起始位置
    background_end: Optional[int] = None     # 背景DNA在染色体上的结束位置

    # 结构标记
    has_chimera: bool = False
    chimera_positions: List[int] = field(default_factory=list)

    # 未去支化分支的锚点（可能影响测序）
    active_branch_anchors: List[int] = field(default_factory=list)

    # 已去支化分支的锚点（解耦：仍可能有残余结构影响）
    debranched_anchors: List[int] = field(default_factory=list)

    # 未去支化分支的详细信息（用于模拟假chimera）
    active_branches: List[ActiveBranchInfo] = field(default_factory=list)

    # 是否存在覆盖环状junction的可能
    junction_covered_possible: bool = False
    
    @property
    def total_length(self) -> int:
        return sum(s.length for s in self.segments)
    
    def get_sequence(self, ecc_db: Dict[str, EccDNA]) -> str:
        """根据segments重建完整序列"""
        result = []
        for seg in self.segments:
            ecc = ecc_db[seg.ecc_id]
            subseq = ecc.get_circular_substr(seg.ecc_offset, seg.length)
            if seg.strand == Strand.REVERSE:
                subseq = reverse_complement(subseq)
            result.append(subseq)
        return "".join(result)
    
    def extract_subregion(self, start: int, length: int) -> List[Segment]:
        """
        提取子区域的segments追踪信息

        这是truth追踪的核心方法

        Args:
            start: 起始位置
            length: 提取长度

        Returns:
            子区域对应的segment列表
        """
        sub_segments = []
        current_pos = 0
        remaining = length
        extract_start = start
        
        for seg in self.segments:
            seg_end = current_pos + seg.length
            
            # 跳过start之前的segment
            if seg_end <= extract_start:
                current_pos = seg_end
                continue
            
            # 计算在当前segment中的截取范围
            local_start = max(0, extract_start - current_pos)
            local_end = min(seg.length, local_start + remaining)
            local_length = local_end - local_start
            
            if local_length > 0:
                # 创建子segment
                new_seg = Segment(
                    ecc_id=seg.ecc_id,
                    ecc_offset=(seg.ecc_offset + local_start) % get_ecc_length(seg.ecc_id),
                    length=local_length,
                    strand=seg.strand,
                    segment_type=seg.segment_type,
                    parent_offset=current_pos + local_start
                )
                sub_segments.append(new_seg)
                remaining -= local_length
                extract_start = seg_end
            
            current_pos = seg_end
            
            if remaining <= 0:
                break
        
        return sub_segments


# =============================================================================
# Read结构
# =============================================================================

@dataclass
class SequencedRead:
    """
    最终的测序read，包含完整truth信息
    """
    read_id: str
    platform: str                  # NGS/HiFi/ONT
    sequence: str
    quality: Optional[str] = None  # FASTQ质量字符串

    # Truth信息
    source_molecule_id: str = ""
    source_ecc_ids: List[str] = field(default_factory=list)
    segments: List[Segment] = field(default_factory=list)

    # 背景DNA标记
    is_background: bool = False    # 是否来自背景线性DNA
    background_chrom: Optional[str] = None   # 背景DNA来源染色体
    background_start: Optional[int] = None   # 背景DNA起始位置
    background_end: Optional[int] = None     # 背景DNA结束位置

    # 结构事件标记
    repeat_count_truth: int = 0    # 覆盖的repeat数
    repeat_count_by_source: Dict[str, float] = field(default_factory=dict)
    has_inter_chimera: bool = False
    chimera_breakpoints: List[int] = field(default_factory=list)
    covered_branch_anchor: bool = False
    has_branch_chimera: bool = False  # 支链导致的假chimera
    truncated: bool = False
    truncation_reason: str = ""
    junction_covered_possible: bool = False

    # Paired-end专用
    is_paired: bool = False
    mate_id: Optional[str] = None
    read_number: int = 1           # 1 or 2
    insert_size: int = 0
    
    def to_fastq(self) -> str:
        """转换为FASTQ格式"""
        qual = self.quality if self.quality else "I" * len(self.sequence)
        return f"@{self.read_id}\n{self.sequence}\n+\n{qual}\n"
    
    def to_truth_dict(self) -> Dict[str, Any]:
        """转换为truth记录"""
        result = {
            "read_id": self.read_id,
            "platform": self.platform,
            "source_type": "background" if self.is_background else "eccDNA",
            "source_molecule_id": self.source_molecule_id,
            "source_ecc_id": self.source_ecc_ids[0] if len(self.source_ecc_ids) == 1 else "multi",
            "source_ecc_ids": self.source_ecc_ids,
            "repeat_count_truth": self.repeat_count_truth,
            "repeat_count_by_source": self.repeat_count_by_source,
            "segments": [s.to_dict() for s in self.segments],
            "has_inter_chimera": int(self.has_inter_chimera),
            "has_branch_chimera": int(self.has_branch_chimera),
            "chimera_breakpoints": self.chimera_breakpoints,
            "covered_branch_anchor": int(self.covered_branch_anchor),
            "truncated": int(self.truncated),
            "truncation_reason": self.truncation_reason,
            "junction_covered_possible": int(self.junction_covered_possible),
            "is_paired": int(self.is_paired),
            "mate_id": self.mate_id or "",
            "read_number": self.read_number,
            "insert_size": self.insert_size
        }
        # 添加背景DNA特有字段
        if self.is_background:
            result["background_chrom"] = self.background_chrom or ""
            result["background_start"] = self.background_start or 0
            result["background_end"] = self.background_end or 0
        return result


# =============================================================================
# Compartment结构
# =============================================================================

@dataclass
class Compartment:
    """空间compartment"""
    compartment_id: int
    ecc_instances: List[Tuple[str, int]]  # (ecc_id, instance_index)
    size: int = 0                          # N
    
    def __post_init__(self):
        self.size = len(self.ecc_instances)


# =============================================================================
# 辅助函数
# =============================================================================

_ecc_length_cache: Dict[str, int] = {}

def register_ecc_length(ecc_id: str, length: int):
    """注册eccDNA长度（用于segment offset计算）"""
    _ecc_length_cache[ecc_id] = length

def get_ecc_length(ecc_id: str, warn_on_missing: bool = True) -> int:
    """
    获取eccDNA长度

    Args:
        ecc_id: eccDNA标识符
        warn_on_missing: 如果为True且ecc_id未注册，发出警告

    Returns:
        eccDNA长度，未注册时返回默认值1000
    """
    if ecc_id in _ecc_length_cache:
        return _ecc_length_cache[ecc_id]
    else:
        if warn_on_missing:
            import warnings
            warnings.warn(
                f"ecc_id '{ecc_id}' not registered in length cache. "
                f"Using default length 1000. This may cause incorrect offset calculations. "
                f"Call register_ecc_length() to register eccDNA lengths.",
                stacklevel=2
            )
        return 1000  # 默认值


def compute_junction_covered_possible(segments: List[Segment]) -> bool:
    """
    判断segments中是否存在跨环状junction的可能

    规则：若某段长度超过该eccDNA从ecc_offset到末端的剩余长度，则该段包含wrap
    """
    for seg in segments:
        if seg.segment_type == SegmentType.BACKGROUND:
            continue
        ecc_len = get_ecc_length(seg.ecc_id, warn_on_missing=False)
        if ecc_len <= 0:
            continue
        offset = seg.ecc_offset % ecc_len
        if seg.length > ecc_len - offset:
            return True
    return False

def reverse_complement(seq: str) -> str:
    """反向互补"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return "".join(complement.get(base, 'N') for base in reversed(seq))
