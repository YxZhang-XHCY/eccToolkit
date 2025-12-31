"""
RCA分子图生成引擎

实现§6：RCA分子图生成（支链/超支化）
- 主干（trunk）生成
- 分支数量与位置采样
- 分支内容生成
- 形成分子图表示
"""

from typing import List, Dict, Tuple, Optional
import numpy as np
import math
from dataclasses import dataclass

from .models import (
    EccDNA, Compartment, RCAMoleculeGraph, BranchNode, Segment,
    Strand, SegmentType, register_ecc_length
)
from .config import SimConfig
from .kinetics import simulate_repeat_count


@dataclass
class RCAStats:
    """RCA生成统计"""
    total_molecules: int = 0
    total_branches: int = 0
    repeat_count_distribution: Dict[int, int] = None
    branch_count_distribution: Dict[int, int] = None
    
    def __post_init__(self):
        if self.repeat_count_distribution is None:
            self.repeat_count_distribution = {}
        if self.branch_count_distribution is None:
            self.branch_count_distribution = {}
    
    def summary(self) -> str:
        avg_repeat = 0
        avg_branch = 0
        if self.repeat_count_distribution:
            total = sum(self.repeat_count_distribution.values())
            avg_repeat = sum(k*v for k,v in self.repeat_count_distribution.items()) / total
        if self.branch_count_distribution:
            total = sum(self.branch_count_distribution.values())
            avg_branch = sum(k*v for k,v in self.branch_count_distribution.items()) / total
        return (
            f"Molecules: {self.total_molecules}, "
            f"Total branches: {self.total_branches}, "
            f"Avg repeat: {avg_repeat:.1f}, "
            f"Avg branches/mol: {avg_branch:.1f}"
        )


def _stochastic_round(
    value: float,
    rng: np.random.Generator,
    min_value: int,
    max_value: int,
) -> int:
    if value <= 0:
        return min_value
    base = int(math.floor(value))
    remainder = value - base
    if remainder > 0 and rng.random() < remainder:
        base += 1
    return max(min_value, min(max_value, base))


class RCAEngine:
    """RCA扩增模拟引擎"""
    
    def __init__(self, config: SimConfig, rng: Optional[np.random.Generator] = None):
        self.config = config
        self.rng = rng if rng is not None else np.random.default_rng()
        self._instance_counter = 0
    
    def sample_base_repeat(self, length: int) -> int:
        """
        采样基础重复次数（未考虑资源竞争）

        ln R_base(L) ~ N(μ_R - k_len * ln(L/L_ref), σ_R^2)
        R_base = clamp(round(exp(ln_R_base)), R_min, R_max)

        改进：
        1. 对 adjusted_mu 做 clamp，防止极端 L 导致均值爆炸/塌缩
        2. 对最终 R 做上下限截断

        Args:
            length: eccDNA长度（必须 > 0）

        Returns:
            基础重复次数（整数）

        Raises:
            ValueError: 如果 length <= 0
        """
        # 输入验证：防止 log(0) 或 log(负数) 导致崩溃
        if length <= 0:
            raise ValueError(f"eccDNA length must be positive, got {length}")

        mu_R = self.config.mu_R
        sigma_R = self.config.sigma_R
        k_len = self.config.k_len
        L_ref = self.config.L_ref
        R_min = self.config.R_min
        R_max = self.config.R_max

        # 长度修正后的均值
        raw_adjusted_mu = mu_R - k_len * math.log(length / L_ref)

        # clamp adjusted_mu 到合理范围 [ln(R_min), ln(R_max)]
        # 防止极端 L 导致期望值爆炸（极小 L）或塌缩（极大 L）
        mu_min = math.log(R_min) if R_min > 0 else -10
        mu_max = math.log(R_max) + 2 * sigma_R  # 允许一些超出空间
        adjusted_mu = max(mu_min, min(mu_max, raw_adjusted_mu))

        # 采样ln(R)
        ln_R = self.rng.normal(adjusted_mu, sigma_R)

        # 转换为整数repeat，并做上下限截断
        R_raw = int(round(math.exp(ln_R)))
        R = max(R_min, min(R_max, R_raw))

        return R

    def sample_repeat_count(self, length: int) -> int:
        """采样重复次数（支持经验模型与动力学模型）"""
        if self.config.rca.mode == "kinetic":
            return simulate_repeat_count(length, self.config, self.rng)
        return self.sample_base_repeat(length)
    
    def sample_branch_count(self, trunk_length: int, ecc_length: int, repeat_count: int) -> int:
        """
        采样分支数量（大环产生更多分支）

        支持两种模式：
        - per_repeat: E[B] = B_rate * (L/L_ref)^k_branch * R / (1 + R/B_sat)
          语义：每个重复周期产生分支，但有饱和效应
          B_sat=inf 时退化为线性
          优点：避免 L 被重复计算，防止高 R 时分支数爆炸

        - per_kb: E[B] = B_rate * (L/L_ref)^k_branch * T/1000
          语义：每 kb 产生一定数量分支
          缺点：L 同时出现在 T 和 (L/L_ref)^k_branch 中

        Args:
            trunk_length: 主干长度 T = R * L
            ecc_length: eccDNA模板长度 L
            repeat_count: 重复次数 R

        Returns:
            分支数量
        """
        B_rate = self.config.B_rate
        k_branch = self.config.k_branch
        L_ref = self.config.L_ref
        B_sat = self.config.B_sat
        branch_mode = self.config.branch.branch_mode

        # 长度修正：大环分支率更高
        adjusted_B_rate = B_rate * math.pow(ecc_length / L_ref, k_branch)

        if branch_mode == "per_repeat":
            # 推荐模式：分支数与重复次数成正比，带饱和效应
            # E[B] = adjusted_B_rate * R / (1 + R/B_sat)
            # 当 B_sat=inf 时，分母=1，退化为线性
            # 当 R >> B_sat 时，E[B] → adjusted_B_rate * B_sat（渐近上限）
            if math.isinf(B_sat):
                expected_branches = adjusted_B_rate * repeat_count
            else:
                expected_branches = adjusted_B_rate * repeat_count / (1 + repeat_count / B_sat)
        else:  # per_kb (原始模式)
            expected_branches = adjusted_B_rate * trunk_length / 1000.0

        return self.rng.poisson(expected_branches)
    
    def sample_branch_length(self, trunk_length: int, ecc_length: int) -> int:
        """
        采样分支长度
        
        ln(r) ~ N(μ_br, σ_br^2)
        len_branch = min(base_len, max(1, floor(r * base_len)))
        
        Args:
            trunk_length: 主干长度T
            ecc_length: eccDNA模板长度L
        
        Returns:
            分支长度
        """
        mu_br = self.config.mu_br
        sigma_br = self.config.sigma_br
        length_mode = self.config.branch.branch_length_mode
        
        ln_r = self.rng.normal(mu_br, sigma_br)
        r = math.exp(ln_r)
        
        if length_mode == "ratio_ecc":
            base_len = min(trunk_length, ecc_length)
        else:
            base_len = trunk_length

        branch_len = min(base_len, max(1, int(r * base_len)))
        return branch_len
    
    def generate_molecule(
        self,
        ecc: EccDNA,
        compartment_id: int,
        effective_repeat: int
    ) -> RCAMoleculeGraph:
        """
        生成单个RCA分子图
        
        Args:
            ecc: eccDNA模板
            compartment_id: compartment ID
            effective_repeat: 有效重复次数（已考虑资源竞争）
        
        Returns:
            RCA分子图
        """
        self._instance_counter += 1
        instance_id = f"mol_{compartment_id}_{self._instance_counter}"
        
        L = ecc.length
        R = effective_repeat
        T = R * L  # 主干长度

        # 生成分支（传入eccDNA长度和重复次数）
        num_branches = self.sample_branch_count(T, L, R)
        branches = []
        
        for b_idx in range(num_branches):
            # 分支锚点：在主干上均匀采样
            anchor_pos = self.rng.integers(0, T)
            
            # 分支长度
            branch_len = self.sample_branch_length(T, L)
            
            # 分支起始offset（继承锚点在模板上的位置）
            branch_offset = anchor_pos % L
            
            # 分支方向（v1默认正向，可选翻转）
            strand = Strand.FORWARD
            if self.rng.random() < self.config.constants.p_flip:
                strand = Strand.REVERSE
            
            # 创建分支segment
            branch_segment = Segment(
                ecc_id=ecc.id,
                ecc_offset=branch_offset,
                length=branch_len,
                strand=strand,
                segment_type=SegmentType.BRANCH,
                parent_offset=anchor_pos
            )
            
            branch = BranchNode(
                branch_id=b_idx,
                anchor_pos=anchor_pos,
                segment=branch_segment,
                debranched=False
            )
            branches.append(branch)
        
        # 创建分子图
        molecule = RCAMoleculeGraph(
            instance_id=instance_id,
            source_ecc=ecc,
            compartment_id=compartment_id,
            repeat_count=R,
            trunk_length=T,
            branches=branches
        )
        
        return molecule
    
    def process_compartment(
        self,
        compartment: Compartment,
        ecc_db: Dict[str, EccDNA]
    ) -> List[RCAMoleculeGraph]:
        """
        处理单个compartment中的所有eccDNA实例

        支持两种资源竞争模式：
        - power_law: 每个分子独立计算 R_eff = R_base / N^gamma
        - total_pool: 先采样微区总资源，再按权重分配给各分子

        Args:
            compartment: compartment对象
            ecc_db: eccDNA字典 {id: EccDNA}

        Returns:
            RCA分子图列表
        """
        from .compartment import CompartmentGenerator

        comp_gen = CompartmentGenerator(self.config, self.rng)
        N = compartment.size
        competition_mode = self.config.constants.competition_mode

        # 收集所有 eccDNA 信息
        ecc_list = [ecc_db[ecc_id] for ecc_id, _ in compartment.ecc_instances]
        base_repeats = [self.sample_repeat_count(ecc.length) for ecc in ecc_list]
        lengths = [ecc.length for ecc in ecc_list]

        # 根据模式计算有效重复数
        if self.config.rca.mode == "kinetic":
            resource_mode = self.config.kinetics.resource_mode
            if resource_mode == "none":
                effective_repeats = base_repeats
            elif resource_mode == "noise":
                sigma = self.config.kinetics.resource_noise_sigma
                if sigma <= 0:
                    effective_repeats = base_repeats
                else:
                    noise = self.rng.lognormal(mean=0.0, sigma=sigma, size=len(base_repeats))
                    effective_repeats = [
                        _stochastic_round(r * n, self.rng, self.config.R_min, self.config.R_max)
                        for r, n in zip(base_repeats, noise)
                    ]
            elif resource_mode == "total_pool":
                effective_repeats = comp_gen.allocate_repeats_total_pool(base_repeats, lengths)
            elif resource_mode == "total_pool_uniform":
                effective_repeats = comp_gen.allocate_repeats_total_pool(
                    base_repeats, lengths, weights=[1.0] * len(base_repeats)
                )
            else:
                # inherit default competition mode
                competition_mode = self.config.constants.competition_mode
                if competition_mode == "total_pool":
                    effective_repeats = comp_gen.allocate_repeats_total_pool(base_repeats, lengths)
                else:
                    effective_repeats = [
                        comp_gen.compute_effective_repeat(r, N) for r in base_repeats
                    ]
        else:
            if competition_mode == "total_pool":
                effective_repeats = comp_gen.allocate_repeats_total_pool(base_repeats, lengths)
            else:  # power_law (default)
                effective_repeats = [
                    comp_gen.compute_effective_repeat(r, N) for r in base_repeats
                ]

        # 生成分子图
        molecules = []
        for ecc, eff_repeat in zip(ecc_list, effective_repeats):
            mol = self.generate_molecule(ecc, compartment.compartment_id, eff_repeat)
            molecules.append(mol)

        return molecules
    
    def generate_all(
        self,
        compartments: List[Compartment],
        ecc_db: Dict[str, EccDNA]
    ) -> Tuple[List[RCAMoleculeGraph], RCAStats]:
        """
        生成所有compartment的RCA分子图
        
        Args:
            compartments: compartment列表
            ecc_db: eccDNA字典
        
        Returns:
            (分子图列表, 统计信息)
        """
        # 注册eccDNA长度（用于segment计算）
        for ecc_id, ecc in ecc_db.items():
            register_ecc_length(ecc_id, ecc.length)
        
        all_molecules = []
        stats = RCAStats()
        
        for comp in compartments:
            molecules = self.process_compartment(comp, ecc_db)
            all_molecules.extend(molecules)
            
            for mol in molecules:
                # 更新统计
                R = mol.repeat_count
                B = len(mol.branches)
                stats.repeat_count_distribution[R] = stats.repeat_count_distribution.get(R, 0) + 1
                stats.branch_count_distribution[B] = stats.branch_count_distribution.get(B, 0) + 1
                stats.total_branches += B
        
        stats.total_molecules = len(all_molecules)
        return all_molecules, stats


def get_trunk_sequence(mol: RCAMoleculeGraph) -> str:
    """
    获取分子主干序列
    
    Args:
        mol: RCA分子图
    
    Returns:
        主干序列字符串
    """
    ecc = mol.source_ecc
    return ecc.get_circular_substr(0, mol.trunk_length)


def get_branch_sequence(mol: RCAMoleculeGraph, branch: BranchNode) -> str:
    """
    获取分支序列
    
    Args:
        mol: RCA分子图
        branch: 分支节点
    
    Returns:
        分支序列字符串
    """
    ecc = mol.source_ecc
    seq = ecc.get_circular_substr(branch.segment.ecc_offset, branch.segment.length)
    
    if branch.segment.strand == Strand.REVERSE:
        from .models import reverse_complement
        seq = reverse_complement(seq)
    
    return seq
