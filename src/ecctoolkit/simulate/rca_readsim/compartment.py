"""
Compartment划分模块

实现§5：空间/稀释模型
- 采样compartment大小（Poisson）
- 将eccDNA实例按权重分配到compartments
- 计算资源竞争对RCA产量的影响
"""

from typing import List, Dict, Tuple, Optional
import numpy as np
from dataclasses import dataclass, field

from .models import EccDNA, Compartment
from .config import SimConfig


@dataclass
class CompartmentStats:
    """Compartment统计信息"""
    total_compartments: int = 0
    total_instances: int = 0
    size_distribution: Dict[int, int] = field(default_factory=dict)
    ecc_instance_counts: Dict[str, int] = field(default_factory=dict)
    
    def summary(self) -> str:
        sizes = list(self.size_distribution.keys())
        if not sizes:
            return "No compartments generated"
        return (
            f"Compartments: {self.total_compartments}, "
            f"Total instances: {self.total_instances}, "
            f"Size range: {min(sizes)}-{max(sizes)}, "
            f"Mean size: {self.total_instances/self.total_compartments:.2f}"
        )


class CompartmentGenerator:
    """Compartment生成器"""
    
    def __init__(self, config: SimConfig, rng: Optional[np.random.Generator] = None):
        self.config = config
        self.rng = rng if rng is not None else np.random.default_rng()
        
    def generate(
        self, 
        eccdnas: List[EccDNA], 
        target_instances: int
    ) -> Tuple[List[Compartment], CompartmentStats]:
        """
        生成compartments并分配eccDNA实例
        
        Args:
            eccdnas: eccDNA模板列表
            target_instances: 目标总实例数（用于控制输出规模）
        
        Returns:
            (compartments列表, 统计信息)
        """
        lambda_comp = self.config.lambda_comp
        
        # 计算权重分布
        weights = np.array([e.weight for e in eccdnas])
        weight_sum = weights.sum()
        if weight_sum > 0:
            weights = weights / weight_sum
        else:
            weights = np.ones(len(eccdnas)) / len(eccdnas)
        ecc_ids = [e.id for e in eccdnas]
        
        compartments = []
        stats = CompartmentStats()
        total_instances = 0
        comp_id = 0
        
        # 持续生成compartments直到达到目标实例数
        while total_instances < target_instances:
            # 采样compartment大小
            N = self.rng.poisson(lambda_comp)
            
            # 跳过空compartment
            if N == 0:
                continue
            
            # 限制不要超过目标太多
            if total_instances + N > target_instances * 1.1:
                N = max(1, target_instances - total_instances)
            
            # 按权重分配eccDNA实例
            instances = []
            sampled_ids = self.rng.choice(
                ecc_ids, 
                size=N, 
                replace=True, 
                p=weights
            )
            
            for idx, eid in enumerate(sampled_ids):
                instances.append((eid, idx))
                stats.ecc_instance_counts[eid] = stats.ecc_instance_counts.get(eid, 0) + 1
            
            comp = Compartment(
                compartment_id=comp_id,
                ecc_instances=instances,
                size=N
            )
            compartments.append(comp)
            
            # 更新统计
            stats.size_distribution[N] = stats.size_distribution.get(N, 0) + 1
            total_instances += N
            comp_id += 1
        
        stats.total_compartments = len(compartments)
        stats.total_instances = total_instances
        
        return compartments, stats
    
    def compute_effective_repeat(
        self,
        base_repeat: int,
        compartment_size: int
    ) -> int:
        """
        计算资源竞争后的有效重复次数（power_law 模式）

        R = max(1, floor(R_base / N^gamma))

        总产量标度: N * R_eff ≈ N^(1-gamma)
        - gamma=1: 总产量守恒（推荐）
        - gamma<1: 总产量随N增长
        - gamma>1: 总产量随N下降

        Args:
            base_repeat: 基础重复次数（无竞争时）
            compartment_size: compartment大小N

        Returns:
            有效重复次数
        """
        gamma = self.config.gamma
        if compartment_size <= 1:
            return base_repeat

        effective = base_repeat / (compartment_size ** gamma)
        return max(1, int(effective))

    def allocate_repeats_total_pool(
        self,
        base_repeats: List[int],
        lengths: List[int],
        weights: Optional[List[float]] = None
    ) -> List[int]:
        """
        资源竞争的总资源池模式

        先采样微区总产量 Y_comp ~ Gamma，再按权重分配给各分子。
        这确保了微区总产量有上限，不会随 N 无限增长。

        权重 w_i ∝ R_base_i（基础重复数，体现模板效率差异）

        Args:
            base_repeats: 各分子的基础重复数列表
            lengths: 各分子的模板长度列表（用于加权，可选）

        Returns:
            各分子的有效重复数列表
        """
        N = len(base_repeats)
        if N == 0:
            return []
        if N == 1:
            return base_repeats.copy()

        # Gamma 分布参数：mean = alpha/beta, var = alpha/beta^2
        # cv = std/mean = 1/sqrt(alpha) => alpha = 1/cv^2
        Y_mean = self.config.constants.Y_comp_mean
        Y_cv = self.config.constants.Y_comp_cv
        Y_max = self.config.constants.Y_comp_max

        if Y_cv > 0:
            alpha = 1.0 / (Y_cv ** 2)
            beta = alpha / Y_mean
            Y_comp = self.rng.gamma(alpha, 1.0 / beta)
            # 应用上限：防止极端采样导致不合理的高重复数
            Y_comp = min(Y_comp, Y_max)
        else:
            Y_comp = Y_mean

        # 权重：默认按基础重复数，也可显式传入
        if weights is None:
            weights = base_repeats
        weights = np.array(weights, dtype=float)
        weight_sum = weights.sum()
        if weight_sum <= 0:
            weights = np.ones(N)
            weight_sum = N

        # 分配总产量
        allocations = Y_comp * (weights / weight_sum)

        # 转换为整数，确保至少为1
        effective_repeats = [max(1, int(round(a))) for a in allocations]

        return effective_repeats


def estimate_total_instances(
    eccdnas: List[EccDNA],
    config: SimConfig,
    target_reads: int,
    platform: str = "NGS"
) -> int:
    """
    估算达到目标reads数所需的eccDNA实例数
    
    这是一个粗略估计，用于初始化compartment生成
    
    Args:
        eccdnas: eccDNA列表
        config: 配置
        target_reads: 目标read数
        platform: 平台类型
    
    Returns:
        估计所需的eccDNA实例数
    """
    # 平均模板长度
    avg_length = np.mean([e.length for e in eccdnas])
    
    # 估计平均repeat次数（忽略竞争）
    import math
    ln_R_mean = config.mu_R - config.k_len * math.log(avg_length / config.L_ref)
    avg_repeat = math.exp(ln_R_mean)
    
    # 估计每个实例产生的平均分子长度
    avg_molecule_len = avg_length * avg_repeat
    
    # 根据平台估计每个分子产生的reads数
    if platform == "NGS":
        # 一个分子可以产生多个insert
        avg_reads_per_molecule = max(1, avg_molecule_len / config.I_mu)
    elif platform == "HiFi":
        # HiFi通常一个分子一条read
        avg_reads_per_molecule = 1.0
    else:  # ONT
        avg_reads_per_molecule = 1.0
    
    # 估算所需实例数（加一些余量）
    estimated_instances = int(target_reads / avg_reads_per_molecule * 1.5)
    
    return max(100, estimated_instances)  # 至少100个实例
