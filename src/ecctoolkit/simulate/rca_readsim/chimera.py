"""
Inter-molecule Chimera模块

实现§7：跨eccDNA拼接
- compartment级inter-chimera概率
- chimera事件注入
"""

from typing import List, Dict, Tuple, Optional
import numpy as np
import math
from dataclasses import dataclass

from .models import (
    EccDNA, RCAMoleculeGraph, ChimeraJunction, Segment,
    Strand, SegmentType
)
from .config import SimConfig


@dataclass
class ChimeraStats:
    """Chimera统计"""
    total_chimera_events: int = 0
    molecules_with_chimera: int = 0
    chimera_by_compartment_size: Dict[int, int] = None
    
    def __post_init__(self):
        if self.chimera_by_compartment_size is None:
            self.chimera_by_compartment_size = {}
    
    def summary(self) -> str:
        return (
            f"Chimera events: {self.total_chimera_events}, "
            f"Molecules with chimera: {self.molecules_with_chimera}"
        )


class ChimeraInjector:
    """Chimera注入器"""
    
    def __init__(
        self,
        config: SimConfig,
        rng: Optional[np.random.Generator] = None,
        background_donors: Optional[List[RCAMoleculeGraph]] = None
    ):
        self.config = config
        self.rng = rng if rng is not None else np.random.default_rng()
        self.background_donors = background_donors or []
        self.background_ratio = max(0.0, min(1.0, self.config.constants.chimera_background_ratio))

    def _sample_breakpoint(self, trunk_length: int) -> int:
        """采样chimera断点位置"""
        if trunk_length <= 1:
            return 0

        mode = self.config.constants.chimera_breakpoint_mode
        if mode == "beta":
            alpha = self.config.constants.chimera_breakpoint_alpha
            beta = self.config.constants.chimera_breakpoint_beta
            value = self.rng.beta(alpha, beta)
            pos = int(value * trunk_length)
            return min(max(0, pos), trunk_length - 1)

        return int(self.rng.integers(0, trunk_length))

    def _use_background_donor(self) -> bool:
        if not self.background_donors:
            return False
        return self.background_ratio > 0 and self.rng.random() < self.background_ratio

    @staticmethod
    def _extract_donor_info(donor) -> Tuple[str, int, int]:
        """返回 (ecc_id, ecc_length, trunk_length)"""
        if hasattr(donor, "source_ecc"):
            return donor.source_ecc.id, donor.source_ecc.length, donor.trunk_length
        if isinstance(donor, dict):
            return donor["ecc_id"], donor["ecc_length"], donor.get("trunk_length", donor["ecc_length"])
        if isinstance(donor, tuple):
            if len(donor) == 2:
                ecc_id, ecc_length = donor
                return ecc_id, ecc_length, ecc_length
            ecc_id, ecc_length, trunk_length = donor[:3]
            return ecc_id, ecc_length, trunk_length
        raise ValueError("Unsupported donor type for chimera injection")
    
    def compute_inter_chimera_prob(self, compartment_size: int) -> float:
        """
        计算compartment内inter-chimera概率
        
        p_inter(N) = 1 - exp(-β * (N-1))
        
        Args:
            compartment_size: compartment大小N
        
        Returns:
            chimera发生概率
        """
        if compartment_size <= 1:
            return 0.0
        
        beta = self.config.beta
        return 1.0 - math.exp(-beta * (compartment_size - 1))
    
    def compute_chimera_length(self, donor_trunk_length: int) -> int:
        """
        计算chimera片段长度

        len_chim = min(0.1 * T_donor, max_len)

        修复：确保chimera长度不超过donor长度，避免凭空创造序列

        Args:
            donor_trunk_length: 供体主干长度

        Returns:
            chimera片段长度（保证 <= donor_trunk_length）
        """
        if donor_trunk_length <= 0:
            return 0

        ratio = self.config.constants.chimera_len_ratio
        max_len = self.config.constants.chimera_max_len

        # 基础长度计算
        length = min(int(ratio * donor_trunk_length), max_len)

        # 确保不超过donor长度（关键修复：避免凭空创造序列）
        length = min(length, donor_trunk_length)

        # 最小长度限制，但不能超过donor长度
        min_chimera = min(100, donor_trunk_length)
        return max(min_chimera, length)
    
    def inject_chimera(
        self,
        recipient: RCAMoleculeGraph,
        donor: RCAMoleculeGraph
    ) -> bool:
        """
        向受体分子注入chimera事件

        Args:
            recipient: 受体分子图
            donor: 供体分子图

        Returns:
            是否成功注入
        """
        T_rec = recipient.trunk_length
        donor_ecc_id, donor_ecc_len, T_don = self._extract_donor_info(donor)

        # 验证：受体和供体都需要有正长度
        if T_rec <= 0 or T_don <= 0:
            return False

        # 在受体主干上选断点
        breakpoint = self._sample_breakpoint(T_rec)

        # 计算chimera长度
        chim_len = self.compute_chimera_length(T_don)

        # 验证：chimera长度必须>0
        if chim_len <= 0:
            return False

        # 从供体选取片段起点
        # 修复：正确处理 T_don <= chim_len 的情况
        if T_don > chim_len:
            # 正常情况：donor有多余空间，随机选择起点
            donor_start = self.rng.integers(0, T_don - chim_len + 1)
        else:
            # donor刚好等于或小于chimera长度，从头开始
            donor_start = 0
            chim_len = T_don  # 使用整个donor

        donor_offset = donor_start % donor_ecc_len
        
        # 创建chimera segment
        chim_segment = Segment(
            ecc_id=donor_ecc_id,
            ecc_offset=donor_offset,
            length=chim_len,
            strand=Strand.FORWARD,
            segment_type=SegmentType.CHIMERA,
            parent_offset=breakpoint
        )
        
        # 创建chimera junction
        junction = ChimeraJunction(
            position=breakpoint,
            donor_ecc_id=donor_ecc_id,
            donor_segment=chim_segment
        )
        
        recipient.chimera_junctions.append(junction)
        return True
    
    def process_compartment_molecules(
        self,
        molecules: List[RCAMoleculeGraph],
        compartment_size: int
    ) -> int:
        """
        处理单个compartment中的chimera事件

        支持两种模式：
        - per_molecule: 每个分子以 P=1-exp(-β*(N-1)) 概率发生嵌合
        - per_pair: K ~ Poisson(alpha * #pairs)，#pairs = N*(N-1)/2
          然后随机选择 K 个分子对进行嵌合

        Args:
            molecules: 该compartment中的分子列表
            compartment_size: compartment大小

        Returns:
            发生的chimera事件数
        """
        if compartment_size <= 1 or len(molecules) <= 1:
            return 0

        chimera_mode = self.config.constants.chimera_mode
        chimera_count = 0

        if chimera_mode == "per_pair":
            # 推荐模式：基于分子对数量的 Poisson 模型
            n_pairs = compartment_size * (compartment_size - 1) // 2
            alpha = self.config.constants.chimera_alpha
            expected_events = alpha * n_pairs
            k_events = self.rng.poisson(expected_events)

            # 随机选择 k_events 个嵌合事件
            for _ in range(k_events):
                if self._use_background_donor():
                    recipient = self.rng.choice(molecules)
                    donor = self.rng.choice(self.background_donors)
                else:
                    # 随机选择受体和供体（不同分子）
                    indices = self.rng.choice(len(molecules), size=2, replace=False)
                    recipient = molecules[indices[0]]
                    donor = molecules[indices[1]]
                if self.inject_chimera(recipient, donor):
                    chimera_count += 1

        else:  # per_molecule (原始模式)
            p_inter = self.compute_inter_chimera_prob(compartment_size)
            for mol in molecules:
                # 以概率p_inter发生chimera
                if self.rng.random() < p_inter:
                    if self._use_background_donor():
                        donor = self.rng.choice(self.background_donors)
                        if self.inject_chimera(mol, donor):
                            chimera_count += 1
                        continue

                    # 从其他分子中随机选择donor
                    other_mols = [m for m in molecules if m.instance_id != mol.instance_id]
                    if other_mols:
                        donor = self.rng.choice(other_mols)
                        if self.inject_chimera(mol, donor):
                            chimera_count += 1

        return chimera_count
    
    def process_all(
        self,
        molecules: List[RCAMoleculeGraph]
    ) -> Tuple[List[RCAMoleculeGraph], ChimeraStats]:
        """
        处理所有分子的chimera注入
        
        Args:
            molecules: 所有分子图
        
        Returns:
            (处理后的分子列表, 统计信息)
        """
        stats = ChimeraStats()
        
        # 按compartment分组
        by_compartment: Dict[int, List[RCAMoleculeGraph]] = {}
        for mol in molecules:
            comp_id = mol.compartment_id
            if comp_id not in by_compartment:
                by_compartment[comp_id] = []
            by_compartment[comp_id].append(mol)
        
        # 处理每个compartment
        for comp_id, comp_mols in by_compartment.items():
            comp_size = len(comp_mols)
            events = self.process_compartment_molecules(comp_mols, comp_size)
            stats.total_chimera_events += events
            
            if events > 0:
                stats.chimera_by_compartment_size[comp_size] = \
                    stats.chimera_by_compartment_size.get(comp_size, 0) + events
        
        # 统计有chimera的分子数
        stats.molecules_with_chimera = sum(
            1 for mol in molecules if mol.chimera_junctions
        )
        
        return molecules, stats
