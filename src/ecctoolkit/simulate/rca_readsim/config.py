"""
配置模块：12个核心参数 + 固定常数

参数分组：
A. RCA主干重复次数（3）: mu_R, sigma_R, k_len
B. 支链强度与分支长度（3）: B_rate, mu_br, sigma_br
C. 空间/稀释（1）: lambda_comp
D. 去支化（1）: D_eff
E. NGS文库（1）: I_mu
F. 输出规模（3）: Cov_NGS, Cov_HiFi, Cov_ONT
"""

from dataclasses import dataclass, field
from typing import Optional, Union, Literal
import math
import yaml
import json


def _serialize_float(value: float) -> Union[float, str]:
    """序列化浮点数，将inf转换为字符串以兼容JSON"""
    if math.isinf(value):
        return "inf" if value > 0 else "-inf"
    if math.isnan(value):
        return "nan"
    return value


def _deserialize_float(value: Union[float, str]) -> float:
    """反序列化浮点数，将字符串inf转换回float"""
    if isinstance(value, str):
        if value == "inf":
            return float('inf')
        elif value == "-inf":
            return float('-inf')
        elif value == "nan":
            return float('nan')
        else:
            return float(value)
    return value


@dataclass
class RCAParams:
    """A. RCA主干重复次数参数"""
    mu_R: float = 3.0          # ln(R)均值 @ L_ref
    sigma_R: float = 0.5       # ln(R)标准差
    k_len: float = 0.5         # 长度修正指数（>0使小环效率更高）
    R_min: int = 1             # 重复数下限
    R_max: int = 10000         # 重复数上限（防止极端L时爆炸）
    mode: Literal["empirical", "kinetic"] = "kinetic"  # 重复数模型


@dataclass
class RCAKineticsParams:
    """RCA动力学参数（Gillespie/CTMC，单位默认: min）"""
    reaction_time_min: float = 960.0  # 16h
    v_nt_per_min: float = 2280.0
    k_on: float = 0.1
    k_start: float = 1.0
    k_pause: float = 0.02
    k_resume: float = 0.05
    k_off: float = 0.01
    k_off_pause: Optional[float] = None
    k_inact_free: float = 0.0005
    k_inact_bound: float = 0.0005
    len_exp_on: float = 0.5
    len_exp_pause: float = 0.5
    len_exp_off: float = 0.5
    len_exp_v: float = 0.0
    len_scale_mode: Literal["power", "exp", "saturating"] = "saturating"
    len_scale_tau: float = 1.0
    length_penalty_mode: Literal["full", "pause_off", "speed_only", "off_only", "none"] = "pause_off"
    rebinding_target: Literal["same_molecule", "random_in_compartment"] = "same_molecule"
    resource_mode: Literal["inherit", "none", "noise", "total_pool", "total_pool_uniform"] = "noise"
    resource_noise_sigma: float = 0.3
    dntp: Optional[float] = None
    k_dntp: float = 1.0
    ppi: Optional[float] = None
    k_ppi: float = 1.0
    off_pause_multiplier: float = 2.0
    max_events: int = 20000


@dataclass
class BranchParams:
    """B. 支链参数"""
    # B_rate 单位说明：
    # - per_repeat模式: 每个重复周期在L_ref处的期望分支数
    # - per_kb模式: 每kb在L_ref处的期望分支数
    B_rate: float = 0.5        # 分支强度（见上述单位说明）
    mu_br: float = -1.5        # 分支长度比例的log均值 (e^-1.5 ≈ 0.22)
    sigma_br: float = 0.5      # 分支长度比例的log标准差
    k_branch: float = 0.3      # 长度修正指数（>0使大环产生更多分支）
    # 分支数计算模式：
    # "per_repeat": E[B] = B_rate * (L/L_ref)^k_branch * R / (1 + R/B_sat)
    #     - B_rate 单位: 分支数/重复周期 @ L_ref
    #     - B_sat: 饱和参数，防止高R时分支数爆炸
    #     - 推荐：避免L被重复计算
    # "per_kb": E[B] = B_rate * (L/L_ref)^k_branch * T/1000
    #     - B_rate 单位: 分支数/kb @ L_ref
    #     - 遗留模式：L同时出现在T和幂律中
    branch_mode: str = "per_repeat"
    # 分支饱和参数：E[B] = base * R / (1 + R/B_sat)
    # B_sat = inf (默认) 表示无饱和（线性增长）
    # B_sat = 50 表示当 R=50 时，实际分支数约为线性预测的一半
    B_sat: float = float('inf')  # 分支饱和阈值（inf=无饱和）
    # 分支长度模式：
    # "ratio_trunk": len = r * T (默认，历史行为)
    # "ratio_ecc": len = r * L (更偏局部、避免过长分支)
    branch_length_mode: str = "ratio_trunk"


@dataclass
class SpatialParams:
    """C. 空间/稀释参数"""
    lambda_comp: float = 3.0   # compartment中eccDNA数的Poisson均值


@dataclass
class DebranchParams:
    """D. 去支化参数"""
    D_eff: float = 0.8         # 去支化基础效率 (0-1)
    # 长度依赖性：D_eff_actual = D_eff * exp(-k_debranch * (L_branch / L_ref))
    # k_debranch > 0 表示长分支更难去除
    # k_debranch = 0 表示长度无关（默认，保持向后兼容）
    k_debranch: float = 0.0    # 长度修正指数（>=0，默认0表示长度无关）
    L_ref_debranch: int = 500  # 去支化参考长度（bp）


@dataclass
class NGSLibraryParams:
    """E. NGS文库参数"""
    I_mu: float = 400.0        # insert size均值


@dataclass
class OutputScaleParams:
    """
    F. 输出规模参数

    支持两种模式：
    - 按read数：直接指定条数
    - 按覆盖度：指定目标覆盖倍数，自动换算
    """
    Cov_NGS: Union[int, float] = 10000      # NGS输出规模
    Cov_HiFi: Union[int, float] = 1000      # HiFi输出规模
    Cov_ONT: Union[int, float] = 1000       # ONT输出规模

    # 输出模式
    mode: Literal["read_count", "coverage"] = "read_count"


@dataclass
class BackgroundDNAParams:
    """
    G. 背景线性DNA参数

    用于模拟真实测序中的线性基因组DNA污染/背景。
    这些线性DNA不经过RCA扩增，直接进行测序。
    """
    enabled: bool = False                   # 是否启用背景DNA
    fasta_path: Optional[str] = None        # 参考基因组FASTA路径
    ratio: float = 0.1                      # 背景DNA占总reads的比例 (0-1)
    min_length: int = 200                   # 线性片段最小长度
    max_length: int = 10000                 # 线性片段最大长度
    # 长度分布模式
    # "uniform": 均匀分布
    # "lognormal": 对数正态分布（更接近真实片段化）
    length_distribution: str = "lognormal"
    length_mean: float = 2000.0             # lognormal模式下的均值
    length_std: float = 1500.0              # lognormal模式下的标准差


# =============================================================================
# 固定常数（v1不作为参数开放）
# =============================================================================

@dataclass
class FixedConstants:
    """固定常数，可通过profile修改但不算入12参数"""
    
    # 参考长度
    L_ref: int = 1000                    # 参考长度（bp）
    
    # Inter-chimera
    # =========================================================================
    # 嵌合体率校准指南：
    # - 文献报道实际嵌合体率通常在 0.1%-2% 范围内
    # - per_pair模式：E[K] = alpha * N*(N-1)/2
    #   - N=3, alpha=0.1 → E[K]=0.3 事件，约10%分子含嵌合体（合理）
    #   - N=5, alpha=0.1 → E[K]=1.0 事件，约33%分子含嵌合体（偏高）
    #   - 建议：若lambda_comp > 5，考虑降低alpha至0.02-0.05
    # - per_molecule模式：P = 1-exp(-beta*(N-1))
    #   - N=3, beta=0.2 → P=33%（可能过高）
    #   - 建议：beta=0.02-0.05 更接近实际
    # =========================================================================
    beta: float = 0.2                    # per_molecule模式：嵌合强度（遗留默认值）
    chimera_max_len: int = 2000          # chimera片段最大长度
    chimera_len_ratio: float = 0.1       # chimera长度占donor的比例
    # 嵌合体模式：
    # "per_molecule": 每个分子以 P=1-exp(-β*(N-1)) 概率发生嵌合（遗留模式）
    # "per_pair": K ~ Poisson(alpha * #pairs)，#pairs = N*(N-1)/2（推荐）
    chimera_mode: str = "per_pair"
    chimera_alpha: float = 0.1           # per_pair模式：每个pair的嵌合强度
    # Chimera断点分布
    # "uniform": 均匀分布
    # "beta": Beta(alpha, beta)，可产生端点偏好（alpha,beta<1）
    chimera_breakpoint_mode: str = "uniform"
    chimera_breakpoint_alpha: float = 0.5
    chimera_breakpoint_beta: float = 0.5
    # Chimera供体来自背景DNA的比例（需要提供reference）
    chimera_background_ratio: float = 0.0
    chimera_background_pool: int = 200   # 背景供体池大小（仅用于chimera）
    
    # 资源竞争
    gamma: float = 1.0                   # 资源竞争指数 (power_law模式使用)
    # 资源竞争模式：
    # "total_pool": 微区总资源 Y_comp ~ Gamma，按权重分配（推荐，更真实）
    # "power_law": R_eff = R_base / N^gamma（遗留模式，对N敏感）
    competition_mode: str = "total_pool"  # 默认使用更稳健的total_pool
    # total_pool 模式参数
    Y_comp_mean: float = 50.0            # 微区平均总重复产出
    Y_comp_cv: float = 0.3               # 变异系数 (std/mean)
    Y_comp_max: float = 500.0            # 微区总产出上限（防止极端采样）
    
    # NGS
    I_sigma_ratio: float = 0.1           # I_sigma = I_mu * ratio
    I_sigma_fixed: Optional[float] = None  # 或者固定值（优先）
    read_length: int = 150               # Illumina read长度
    
    # 分支方向翻转
    p_flip: float = 0.0                  # 分支方向翻转概率（v1默认0）
    
    # ONT截断（解耦去支化与截断的强关联）
    p_trunc_per_anchor: float = 0.3      # 每个未去支化锚点的截断概率
    p_trunc_debranched: float = 0.05     # 已去支化锚点的残余截断概率（结构损伤）
    
    # HiFi
    hifi_min_length: int = 15000         # HiFi最小分子长度
    hifi_max_length: int = 25000         # HiFi最大分子长度
    hifi_target_size_mean: int = 20000   # HiFi打断目标片段均值
    hifi_target_size_std: int = 1200     # HiFi打断目标片段标准差

    # ONT
    ont_min_length: int = 200            # ONT最小read长度
    lambda_clog: float = 0.00005         # ONT孔道堵塞率/碱基（20kb读长约63%堵塞概率）

    # 去支化导致的nick断裂（T7处理后的效应）
    p_nick_per_debranch: float = 0.7     # 每个去支化锚点产生nick的概率
    p_nick_to_break: float = 0.2         # nick转为双链断裂的概率

    # 未去支化分支导致的假chimera（建库打断时）
    # 当打断点落在未去支化锚点附近时，产生主干+分支的假chimera片段
    p_branch_chimera: float = 0.3        # 覆盖锚点的insert产生假chimera的概率
    branch_chimera_window: int = 100     # 锚点附近多少bp内认为是"覆盖"

    # 输出规模足够大时，强制每个eccDNA至少出现一次
    force_ecc_coverage: bool = True
    force_ecc_min_reads_factor: float = 2.0

    @property
    def I_sigma(self) -> float:
        """获取insert size标准差"""
        if self.I_sigma_fixed is not None:
            return self.I_sigma_fixed
        return 50.0  # 默认固定50bp，在SimConfig中会用I_mu计算


# =============================================================================
# 完整配置
# =============================================================================

@dataclass
class SimConfig:
    """完整模拟配置"""

    # 12个核心参数（分组）
    rca: RCAParams = field(default_factory=RCAParams)
    kinetics: RCAKineticsParams = field(default_factory=RCAKineticsParams)
    branch: BranchParams = field(default_factory=BranchParams)
    spatial: SpatialParams = field(default_factory=SpatialParams)
    debranch: DebranchParams = field(default_factory=DebranchParams)
    ngs_library: NGSLibraryParams = field(default_factory=NGSLibraryParams)
    output_scale: OutputScaleParams = field(default_factory=OutputScaleParams)

    # 背景线性DNA参数
    background: BackgroundDNAParams = field(default_factory=BackgroundDNAParams)

    # 固定常数
    constants: FixedConstants = field(default_factory=FixedConstants)

    # 随机种子（不计入12参数）
    seed: Optional[int] = None
    
    # =========================================================================
    # 便捷属性访问（扁平化）
    # =========================================================================
    
    @property
    def mu_R(self) -> float:
        return self.rca.mu_R
    
    @property
    def sigma_R(self) -> float:
        return self.rca.sigma_R
    
    @property
    def k_len(self) -> float:
        return self.rca.k_len

    @property
    def R_min(self) -> int:
        return self.rca.R_min

    @property
    def R_max(self) -> int:
        return self.rca.R_max

    @property
    def B_rate(self) -> float:
        return self.branch.B_rate
    
    @property
    def mu_br(self) -> float:
        return self.branch.mu_br
    
    @property
    def sigma_br(self) -> float:
        return self.branch.sigma_br

    @property
    def k_branch(self) -> float:
        return self.branch.k_branch

    @property
    def B_sat(self) -> float:
        return self.branch.B_sat

    @property
    def lambda_comp(self) -> float:
        return self.spatial.lambda_comp
    
    @property
    def D_eff(self) -> float:
        return self.debranch.D_eff

    @property
    def k_debranch(self) -> float:
        return self.debranch.k_debranch

    @property
    def L_ref_debranch(self) -> int:
        return self.debranch.L_ref_debranch

    @property
    def I_mu(self) -> float:
        return self.ngs_library.I_mu
    
    @property
    def I_sigma(self) -> float:
        if self.constants.I_sigma_fixed is not None:
            return self.constants.I_sigma_fixed
        return self.I_mu * self.constants.I_sigma_ratio
    
    @property
    def L_ref(self) -> int:
        return self.constants.L_ref
    
    @property
    def beta(self) -> float:
        return self.constants.beta
    
    @property
    def gamma(self) -> float:
        return self.constants.gamma
    
    # =========================================================================
    # 序列化/反序列化
    # =========================================================================
    
    def to_dict(self) -> dict:
        """转换为字典（完整版，支持 round-trip 序列化）"""
        return {
            "rca": {
                "mu_R": self.rca.mu_R,
                "sigma_R": self.rca.sigma_R,
                "k_len": self.rca.k_len,
                "R_min": self.rca.R_min,
                "R_max": self.rca.R_max,
                "mode": self.rca.mode
            },
            "kinetics": {
                "reaction_time_min": self.kinetics.reaction_time_min,
                "v_nt_per_min": self.kinetics.v_nt_per_min,
                "k_on": self.kinetics.k_on,
                "k_start": self.kinetics.k_start,
                "k_pause": self.kinetics.k_pause,
                "k_resume": self.kinetics.k_resume,
                "k_off": self.kinetics.k_off,
                "k_off_pause": self.kinetics.k_off_pause,
                "k_inact_free": self.kinetics.k_inact_free,
                "k_inact_bound": self.kinetics.k_inact_bound,
                "len_exp_on": self.kinetics.len_exp_on,
                "len_exp_pause": self.kinetics.len_exp_pause,
                "len_exp_off": self.kinetics.len_exp_off,
                "len_exp_v": self.kinetics.len_exp_v,
                "len_scale_mode": self.kinetics.len_scale_mode,
                "len_scale_tau": self.kinetics.len_scale_tau,
                "length_penalty_mode": self.kinetics.length_penalty_mode,
                "rebinding_target": self.kinetics.rebinding_target,
                "resource_mode": self.kinetics.resource_mode,
                "resource_noise_sigma": self.kinetics.resource_noise_sigma,
                "dntp": self.kinetics.dntp,
                "k_dntp": self.kinetics.k_dntp,
                "ppi": self.kinetics.ppi,
                "k_ppi": self.kinetics.k_ppi,
                "off_pause_multiplier": self.kinetics.off_pause_multiplier,
                "max_events": self.kinetics.max_events
            },
            "branch": {
                "B_rate": self.branch.B_rate,
                "mu_br": self.branch.mu_br,
                "sigma_br": self.branch.sigma_br,
                "k_branch": self.branch.k_branch,
                "branch_mode": self.branch.branch_mode,
                "B_sat": _serialize_float(self.branch.B_sat),
                "branch_length_mode": self.branch.branch_length_mode
            },
            "spatial": {
                "lambda_comp": self.spatial.lambda_comp
            },
            "debranch": {
                "D_eff": self.debranch.D_eff,
                "k_debranch": self.debranch.k_debranch,
                "L_ref_debranch": self.debranch.L_ref_debranch
            },
            "ngs_library": {
                "I_mu": self.ngs_library.I_mu
            },
            "output_scale": {
                "Cov_NGS": self.output_scale.Cov_NGS,
                "Cov_HiFi": self.output_scale.Cov_HiFi,
                "Cov_ONT": self.output_scale.Cov_ONT,
                "mode": self.output_scale.mode
            },
            "background": {
                "enabled": self.background.enabled,
                "fasta_path": self.background.fasta_path,
                "ratio": self.background.ratio,
                "min_length": self.background.min_length,
                "max_length": self.background.max_length,
                "length_distribution": self.background.length_distribution,
                "length_mean": self.background.length_mean,
                "length_std": self.background.length_std
            },
            "constants": {
                "L_ref": self.constants.L_ref,
                "beta": self.constants.beta,
                "chimera_max_len": self.constants.chimera_max_len,
                "chimera_len_ratio": self.constants.chimera_len_ratio,
                "chimera_mode": self.constants.chimera_mode,
                "chimera_alpha": self.constants.chimera_alpha,
                "chimera_breakpoint_mode": self.constants.chimera_breakpoint_mode,
                "chimera_breakpoint_alpha": self.constants.chimera_breakpoint_alpha,
                "chimera_breakpoint_beta": self.constants.chimera_breakpoint_beta,
                "chimera_background_ratio": self.constants.chimera_background_ratio,
                "chimera_background_pool": self.constants.chimera_background_pool,
                "gamma": self.constants.gamma,
                "competition_mode": self.constants.competition_mode,
                "Y_comp_mean": self.constants.Y_comp_mean,
                "Y_comp_cv": self.constants.Y_comp_cv,
                "Y_comp_max": self.constants.Y_comp_max,
                "I_sigma_ratio": self.constants.I_sigma_ratio,
                "I_sigma_fixed": self.constants.I_sigma_fixed,
                "read_length": self.constants.read_length,
                "p_flip": self.constants.p_flip,
                "p_trunc_per_anchor": self.constants.p_trunc_per_anchor,
                "p_trunc_debranched": self.constants.p_trunc_debranched,
                "hifi_min_length": self.constants.hifi_min_length,
                "hifi_max_length": self.constants.hifi_max_length,
                "hifi_target_size_mean": self.constants.hifi_target_size_mean,
                "hifi_target_size_std": self.constants.hifi_target_size_std,
                "ont_min_length": self.constants.ont_min_length,
                "lambda_clog": self.constants.lambda_clog,
                "p_nick_per_debranch": self.constants.p_nick_per_debranch,
                "p_nick_to_break": self.constants.p_nick_to_break,
                "p_branch_chimera": self.constants.p_branch_chimera,
                "branch_chimera_window": self.constants.branch_chimera_window,
                "force_ecc_coverage": self.constants.force_ecc_coverage,
                "force_ecc_min_reads_factor": self.constants.force_ecc_min_reads_factor
            },
            "seed": self.seed
        }
    
    @classmethod
    def from_dict(cls, d: dict) -> 'SimConfig':
        """从字典创建"""
        config = cls()

        if "rca" in d:
            config.rca = RCAParams(**d["rca"])
        if "kinetics" in d:
            config.kinetics = RCAKineticsParams(**d["kinetics"])
        if "branch" in d:
            branch_dict = d["branch"].copy()
            # 处理B_sat的特殊序列化值
            if "B_sat" in branch_dict:
                branch_dict["B_sat"] = _deserialize_float(branch_dict["B_sat"])
            config.branch = BranchParams(**branch_dict)
        if "spatial" in d:
            config.spatial = SpatialParams(**d["spatial"])
        if "debranch" in d:
            config.debranch = DebranchParams(**d["debranch"])
        if "ngs_library" in d:
            config.ngs_library = NGSLibraryParams(**d["ngs_library"])
        if "output_scale" in d:
            config.output_scale = OutputScaleParams(**d["output_scale"])
        if "background" in d:
            config.background = BackgroundDNAParams(**d["background"])
        if "constants" in d:
            config.constants = FixedConstants(**d["constants"])
        if "seed" in d:
            config.seed = d["seed"]

        return config
    
    @classmethod
    def from_yaml(cls, path: str) -> 'SimConfig':
        """从YAML文件加载"""
        with open(path, 'r') as f:
            d = yaml.safe_load(f)
        return cls.from_dict(d)
    
    def to_yaml(self, path: str):
        """保存为YAML文件"""
        with open(path, 'w') as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False)
    
    @classmethod
    def from_json(cls, path: str) -> 'SimConfig':
        """从JSON文件加载"""
        with open(path, 'r') as f:
            d = json.load(f)
        return cls.from_dict(d)
    
    def to_json(self, path: str):
        """保存为JSON文件"""
        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)
    
    # =========================================================================
    # 验证
    # =========================================================================
    
    def validate(self) -> list:
        """验证配置有效性，返回警告列表"""
        warnings = []
        
        # RCA参数
        if self.rca.mode == "empirical":
            if self.k_len <= 0:
                warnings.append("k_len应>0以保证小环效率更高")
            if self.sigma_R <= 0:
                warnings.append("sigma_R必须>0")
        elif self.rca.mode == "kinetic":
            if self.kinetics.reaction_time_min <= 0:
                warnings.append("kinetics.reaction_time_min必须>0")
            if self.kinetics.v_nt_per_min <= 0:
                warnings.append("kinetics.v_nt_per_min必须>0")
            if self.kinetics.max_events <= 0:
                warnings.append("kinetics.max_events必须>0")
            if self.kinetics.len_scale_tau <= 0:
                warnings.append("kinetics.len_scale_tau必须>0")
            if self.kinetics.resource_noise_sigma < 0:
                warnings.append("kinetics.resource_noise_sigma必须>=0")
            if min(
                self.kinetics.k_on,
                self.kinetics.k_start,
                self.kinetics.k_pause,
                self.kinetics.k_resume,
                self.kinetics.k_off,
                self.kinetics.k_inact_free,
                self.kinetics.k_inact_bound,
            ) < 0:
                warnings.append("kinetics速率参数不能为负")
            if self.kinetics.off_pause_multiplier < 1:
                warnings.append("kinetics.off_pause_multiplier应>=1")
        else:
            warnings.append(f"未知rca.mode: {self.rca.mode}")
        
        # 分支参数
        if self.B_rate < 0:
            warnings.append("B_rate不能为负")
        if self.sigma_br <= 0:
            warnings.append("sigma_br必须>0")
        if self.branch.branch_length_mode not in ["ratio_trunk", "ratio_ecc"]:
            warnings.append(f"未知branch_length_mode: {self.branch.branch_length_mode}")
        
        # 空间参数
        if self.lambda_comp <= 0:
            warnings.append("lambda_comp必须>0")
        
        # 去支化
        if not 0 <= self.D_eff <= 1:
            warnings.append("D_eff必须在[0,1]范围内")
        
        # NGS
        if self.I_mu < 2 * self.constants.read_length:
            warnings.append(f"I_mu({self.I_mu})应>=2*read_length({self.constants.read_length})")

        # 背景DNA
        if self.background.enabled:
            if not 0 < self.background.ratio < 1:
                warnings.append("background.ratio必须在(0,1)范围内")
            if self.background.fasta_path is None:
                warnings.append("启用背景DNA时必须提供fasta_path")
            if self.background.min_length >= self.background.max_length:
                warnings.append("background.min_length必须小于max_length")

        # Chimera扩展参数
        if self.constants.chimera_breakpoint_mode not in ["uniform", "beta"]:
            warnings.append(f"未知chimera_breakpoint_mode: {self.constants.chimera_breakpoint_mode}")
        if not 0 <= self.constants.chimera_background_ratio <= 1:
            warnings.append("chimera_background_ratio必须在[0,1]范围内")
        if self.constants.chimera_background_pool < 0:
            warnings.append("chimera_background_pool必须>=0")
        if self.constants.chimera_background_ratio > 0 and not self.background.fasta_path:
            warnings.append("chimera_background_ratio>0时必须提供reference/背景fasta_path")

        # 强制覆盖
        if self.constants.force_ecc_coverage and self.constants.force_ecc_min_reads_factor <= 0:
            warnings.append("force_ecc_min_reads_factor必须>0")

        return warnings


# =============================================================================
# 预设配置
# =============================================================================

def get_default_config() -> SimConfig:
    """获取默认配置"""
    return SimConfig()


def get_high_branch_config() -> SimConfig:
    """高分支率配置"""
    config = SimConfig()
    config.branch.B_rate = 2.0
    config.branch.mu_br = -1.0
    return config


def get_dilute_config() -> SimConfig:
    """高稀释配置（减少chimera）"""
    config = SimConfig()
    config.spatial.lambda_comp = 1.0
    return config


def get_dense_config() -> SimConfig:
    """高密度配置（增加chimera）"""
    config = SimConfig()
    config.spatial.lambda_comp = 10.0
    return config
