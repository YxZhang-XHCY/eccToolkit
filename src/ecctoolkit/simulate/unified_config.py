"""
统一模拟配置模块

整合 sim-region（区域生成）和 readsim（读段模拟）的所有配置参数。
"""

from dataclasses import dataclass, field
from typing import Optional, Literal, Union
from pathlib import Path
import yaml


@dataclass
class LengthDistributionConfig:
    """eccDNA 长度分布配置"""
    mode: float = 400.0           # Lognormal 峰值 (bp)
    sigma: float = 0.8            # Lognormal 标准差
    tail_weight: float = 0.05     # 尾部分布权重 (0-1)
    tail_min: int = 5000          # 尾部最小长度 (bp)
    tail_max: int = 500000        # 尾部最大长度 (bp)
    min_length: int = 100         # 全局最小长度 (bp)
    max_length: int = 500000      # 全局最大长度 (bp)


@dataclass
class ClassificationConfig:
    """eccDNA 分类阈值配置"""
    identity: float = 99.0             # 最小 identity 阈值 (%)
    min_coverage: float = 90.0         # self-hit 覆盖率阈值 (%)
    length_consistency: float = 99.0   # 长度一致性阈值 (%)
    multi_coverage: float = 95.0       # multi-hit 覆盖率阈值 (%)
    max_secondary: int = 50            # 最多副比对数
    no_hit_policy: Literal["skip", "unique"] = "skip"


@dataclass
class MinimapConfig:
    """minimap2 预设配置（与 CircleSeeker HiFi 检测保持一致）"""
    split_by_length: bool = False      # 是否按长度分流
    split_length: int = 10000          # 分流阈值 (bp)
    preset_short: str = "map-hifi"     # 短序列预设（HiFi 一致性）
    preset_long: str = "map-hifi"      # 长序列预设（HiFi 一致性）


@dataclass
class RegionConfig:
    """区域模拟配置 (sim-region)"""
    # 目标数量
    num_unique: int = 1000             # UeccDNA 数量
    num_multi: int = 0                 # MeccDNA 数量
    num_chimeric: int = 0              # CeccDNA 数量
    max_length_multi: int = 10000      # MeccDNA 最大长度 (bp)

    # 长度分布
    length_distribution: LengthDistributionConfig = field(
        default_factory=LengthDistributionConfig
    )

    # 分类阈值
    classification: ClassificationConfig = field(
        default_factory=ClassificationConfig
    )

    # minimap2 设置
    minimap: MinimapConfig = field(default_factory=MinimapConfig)

    # 候选倍数
    candidate_multiplier_u: float = 2.0
    candidate_multiplier_m: float = 10.0

    # 输出选项
    keep_tmp: bool = False


@dataclass
class NGSConfig:
    """NGS 测序参数"""
    platform: str = "HS25"
    read_length: int = 150
    insert_mean: float = 400.0
    insert_std: float = 125.0


@dataclass
class ONTConfig:
    """ONT 测序参数

    默认长度参数基于真实 ONT R9.4.1 PromethION 数据:
    - mean_length: 3500 bp
    - std_length: 3000 bp
    """
    model: str = "R94"
    mean_length: float = 3500.0
    std_length: float = 3000.0


@dataclass
class HiFiConfig:
    """HiFi 测序参数"""
    mode: Literal["auto", "sampling", "simple"] = "auto"
    sample_fastq: Optional[str] = None
    profile_id: Optional[str] = None
    profile_root: Optional[str] = None
    len_min: int = 5000
    len_peak_min: int = 10000
    len_peak_max: int = 25000
    len_max: int = 60000
    qmin: int = 20
    qmean: int = 30
    qsd: float = 0.0
    total_reads: Optional[int] = None


@dataclass
class ReadsimParams:
    """读段模拟参数"""
    # 分子池生成
    circular_number: int = 5000
    linear_number: int = 5000
    simple_ratio: Optional[float] = None
    simple_template: Optional[str] = None
    chimeric_template: Optional[str] = None

    # RCA 扩增
    amp: int = 50000           # 最小扩增长度 (bp)
    min_repeats: int = 5       # 每个 eccDNA 最少重复次数
    meancov: float = 25.0

    # 平台参数
    ngs: NGSConfig = field(default_factory=NGSConfig)
    ont: ONTConfig = field(default_factory=ONTConfig)
    hifi: HiFiConfig = field(default_factory=HiFiConfig)


@dataclass
class PlatformConfig:
    """平台选择配置"""
    ngs: bool = True
    hifi: bool = True
    ont: bool = True


@dataclass
class ReadsimConfig:
    """读段模拟配置"""
    # mode: "on" 启用读段模拟, "none" 禁用
    mode: Literal["on", "none"] = "on"
    platforms: PlatformConfig = field(default_factory=PlatformConfig)
    params: ReadsimParams = field(default_factory=ReadsimParams)

    @property
    def enabled(self) -> bool:
        """读段模拟是否启用"""
        return self.mode != "none"


@dataclass
class UnifiedSimulateConfig:
    """统一模拟配置"""

    # 全局设置
    sample: str = "sim_ecc"
    seed: Optional[int] = 42
    threads: int = 8

    # 区域模拟
    region: RegionConfig = field(default_factory=RegionConfig)

    # 读段模拟
    readsim: ReadsimConfig = field(default_factory=ReadsimConfig)

    # =========================================================================
    # 序列化/反序列化
    # =========================================================================

    def to_dict(self) -> dict:
        """转换为字典"""
        return {
            "sample": self.sample,
            "seed": self.seed,
            "threads": self.threads,
            "region": {
                "num_unique": self.region.num_unique,
                "num_multi": self.region.num_multi,
                "num_chimeric": self.region.num_chimeric,
                "max_length_multi": self.region.max_length_multi,
                "length_distribution": {
                    "mode": self.region.length_distribution.mode,
                    "sigma": self.region.length_distribution.sigma,
                    "tail_weight": self.region.length_distribution.tail_weight,
                    "tail_min": self.region.length_distribution.tail_min,
                    "tail_max": self.region.length_distribution.tail_max,
                    "min_length": self.region.length_distribution.min_length,
                    "max_length": self.region.length_distribution.max_length,
                },
                "classification": {
                    "identity": self.region.classification.identity,
                    "min_coverage": self.region.classification.min_coverage,
                    "length_consistency": self.region.classification.length_consistency,
                    "multi_coverage": self.region.classification.multi_coverage,
                    "max_secondary": self.region.classification.max_secondary,
                    "no_hit_policy": self.region.classification.no_hit_policy,
                },
                "minimap": {
                    "split_by_length": self.region.minimap.split_by_length,
                    "split_length": self.region.minimap.split_length,
                    "preset_short": self.region.minimap.preset_short,
                    "preset_long": self.region.minimap.preset_long,
                },
                "candidate_multiplier_u": self.region.candidate_multiplier_u,
                "candidate_multiplier_m": self.region.candidate_multiplier_m,
                "keep_tmp": self.region.keep_tmp,
            },
            "readsim": {
                "mode": self.readsim.mode,
                "platforms": {
                    "ngs": self.readsim.platforms.ngs,
                    "hifi": self.readsim.platforms.hifi,
                    "ont": self.readsim.platforms.ont,
                },
                "params": {
                    "circular_number": self.readsim.params.circular_number,
                    "linear_number": self.readsim.params.linear_number,
                    "simple_ratio": self.readsim.params.simple_ratio,
                    "simple_template": self.readsim.params.simple_template,
                    "chimeric_template": self.readsim.params.chimeric_template,
                    "amp": self.readsim.params.amp,
                    "min_repeats": self.readsim.params.min_repeats,
                    "meancov": self.readsim.params.meancov,
                    "ngs": {
                        "platform": self.readsim.params.ngs.platform,
                        "read_length": self.readsim.params.ngs.read_length,
                        "insert_mean": self.readsim.params.ngs.insert_mean,
                        "insert_std": self.readsim.params.ngs.insert_std,
                    },
                    "ont": {
                        "model": self.readsim.params.ont.model,
                        "mean_length": self.readsim.params.ont.mean_length,
                        "std_length": self.readsim.params.ont.std_length,
                    },
                    "hifi": {
                        "mode": self.readsim.params.hifi.mode,
                        "sample_fastq": self.readsim.params.hifi.sample_fastq,
                        "profile_id": self.readsim.params.hifi.profile_id,
                        "profile_root": self.readsim.params.hifi.profile_root,
                        "len_min": self.readsim.params.hifi.len_min,
                        "len_peak_min": self.readsim.params.hifi.len_peak_min,
                        "len_peak_max": self.readsim.params.hifi.len_peak_max,
                        "len_max": self.readsim.params.hifi.len_max,
                        "qmin": self.readsim.params.hifi.qmin,
                        "qmean": self.readsim.params.hifi.qmean,
                        "qsd": self.readsim.params.hifi.qsd,
                        "total_reads": self.readsim.params.hifi.total_reads,
                    },
                },
            },
        }

    @classmethod
    def from_dict(cls, d: dict) -> "UnifiedSimulateConfig":
        """从字典创建配置"""
        config = cls()

        # 全局设置
        if "sample" in d:
            config.sample = d["sample"]
        if "seed" in d:
            config.seed = d["seed"]
        if "threads" in d:
            config.threads = d["threads"]

        # 区域模拟配置
        if "region" in d:
            r = d["region"]
            config.region.num_unique = r.get("num_unique", config.region.num_unique)
            config.region.num_multi = r.get("num_multi", config.region.num_multi)
            config.region.num_chimeric = r.get("num_chimeric", config.region.num_chimeric)
            config.region.max_length_multi = r.get(
                "max_length_multi", config.region.max_length_multi
            )

            if "length_distribution" in r:
                ld = r["length_distribution"]
                config.region.length_distribution.mode = ld.get(
                    "mode", config.region.length_distribution.mode
                )
                config.region.length_distribution.sigma = ld.get(
                    "sigma", config.region.length_distribution.sigma
                )
                config.region.length_distribution.tail_weight = ld.get(
                    "tail_weight", config.region.length_distribution.tail_weight
                )
                config.region.length_distribution.tail_min = ld.get(
                    "tail_min", config.region.length_distribution.tail_min
                )
                config.region.length_distribution.tail_max = ld.get(
                    "tail_max", config.region.length_distribution.tail_max
                )
                config.region.length_distribution.min_length = ld.get(
                    "min_length", config.region.length_distribution.min_length
                )
                config.region.length_distribution.max_length = ld.get(
                    "max_length", config.region.length_distribution.max_length
                )

            if "classification" in r:
                c = r["classification"]
                config.region.classification.identity = c.get(
                    "identity", config.region.classification.identity
                )
                config.region.classification.min_coverage = c.get(
                    "min_coverage", config.region.classification.min_coverage
                )
                config.region.classification.length_consistency = c.get(
                    "length_consistency", config.region.classification.length_consistency
                )
                config.region.classification.multi_coverage = c.get(
                    "multi_coverage", config.region.classification.multi_coverage
                )
                config.region.classification.max_secondary = c.get(
                    "max_secondary", config.region.classification.max_secondary
                )
                config.region.classification.no_hit_policy = c.get(
                    "no_hit_policy", config.region.classification.no_hit_policy
                )

            if "minimap" in r:
                m = r["minimap"]
                config.region.minimap.split_by_length = m.get(
                    "split_by_length", config.region.minimap.split_by_length
                )
                config.region.minimap.split_length = m.get(
                    "split_length", config.region.minimap.split_length
                )
                config.region.minimap.preset_short = m.get(
                    "preset_short", config.region.minimap.preset_short
                )
                config.region.minimap.preset_long = m.get(
                    "preset_long", config.region.minimap.preset_long
                )

            config.region.candidate_multiplier_u = r.get(
                "candidate_multiplier_u", config.region.candidate_multiplier_u
            )
            config.region.candidate_multiplier_m = r.get(
                "candidate_multiplier_m", config.region.candidate_multiplier_m
            )
            config.region.keep_tmp = r.get("keep_tmp", config.region.keep_tmp)

        # 读段模拟配置
        if "readsim" in d:
            rs = d["readsim"]
            config.readsim.mode = rs.get("mode", config.readsim.mode)

            if "platforms" in rs:
                p = rs["platforms"]
                config.readsim.platforms.ngs = p.get("ngs", config.readsim.platforms.ngs)
                config.readsim.platforms.hifi = p.get("hifi", config.readsim.platforms.hifi)
                config.readsim.platforms.ont = p.get("ont", config.readsim.platforms.ont)

            # 读段模拟参数
            if "params" in rs:
                p = rs["params"]
                config.readsim.params.circular_number = p.get(
                    "circular_number", config.readsim.params.circular_number
                )
                config.readsim.params.linear_number = p.get(
                    "linear_number", config.readsim.params.linear_number
                )
                config.readsim.params.simple_ratio = p.get(
                    "simple_ratio", config.readsim.params.simple_ratio
                )
                config.readsim.params.simple_template = p.get(
                    "simple_template", config.readsim.params.simple_template
                )
                config.readsim.params.chimeric_template = p.get(
                    "chimeric_template", config.readsim.params.chimeric_template
                )
                config.readsim.params.amp = p.get("amp", config.readsim.params.amp)
                config.readsim.params.min_repeats = p.get(
                    "min_repeats", config.readsim.params.min_repeats
                )
                config.readsim.params.meancov = p.get(
                    "meancov", config.readsim.params.meancov
                )

                if "ngs" in p:
                    ngs = p["ngs"]
                    config.readsim.params.ngs.platform = ngs.get(
                        "platform", config.readsim.params.ngs.platform
                    )
                    config.readsim.params.ngs.read_length = ngs.get(
                        "read_length", config.readsim.params.ngs.read_length
                    )
                    config.readsim.params.ngs.insert_mean = ngs.get(
                        "insert_mean", config.readsim.params.ngs.insert_mean
                    )
                    config.readsim.params.ngs.insert_std = ngs.get(
                        "insert_std", config.readsim.params.ngs.insert_std
                    )

                if "ont" in p:
                    ont = p["ont"]
                    config.readsim.params.ont.model = ont.get(
                        "model", config.readsim.params.ont.model
                    )
                    config.readsim.params.ont.mean_length = ont.get(
                        "mean_length", config.readsim.params.ont.mean_length
                    )
                    config.readsim.params.ont.std_length = ont.get(
                        "std_length", config.readsim.params.ont.std_length
                    )

                if "hifi" in p:
                    hifi = p["hifi"]
                    config.readsim.params.hifi.mode = hifi.get(
                        "mode", config.readsim.params.hifi.mode
                    )
                    config.readsim.params.hifi.sample_fastq = hifi.get(
                        "sample_fastq", config.readsim.params.hifi.sample_fastq
                    )
                    config.readsim.params.hifi.profile_id = hifi.get(
                        "profile_id", config.readsim.params.hifi.profile_id
                    )
                    config.readsim.params.hifi.profile_root = hifi.get(
                        "profile_root", config.readsim.params.hifi.profile_root
                    )
                    config.readsim.params.hifi.len_min = hifi.get(
                        "len_min", config.readsim.params.hifi.len_min
                    )
                    config.readsim.params.hifi.len_peak_min = hifi.get(
                        "len_peak_min", config.readsim.params.hifi.len_peak_min
                    )
                    config.readsim.params.hifi.len_peak_max = hifi.get(
                        "len_peak_max", config.readsim.params.hifi.len_peak_max
                    )
                    config.readsim.params.hifi.len_max = hifi.get(
                        "len_max", config.readsim.params.hifi.len_max
                    )
                    config.readsim.params.hifi.qmin = hifi.get(
                        "qmin", config.readsim.params.hifi.qmin
                    )
                    config.readsim.params.hifi.qmean = hifi.get(
                        "qmean", config.readsim.params.hifi.qmean
                    )
                    config.readsim.params.hifi.qsd = hifi.get(
                        "qsd", config.readsim.params.hifi.qsd
                    )
                    config.readsim.params.hifi.total_reads = hifi.get(
                        "total_reads", config.readsim.params.hifi.total_reads
                    )

        return config

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "UnifiedSimulateConfig":
        """从 YAML 文件加载配置"""
        with open(path, "r") as f:
            d = yaml.safe_load(f)
        return cls.from_dict(d if d else {})

    def to_yaml(self, path: Union[str, Path]) -> None:
        """保存为 YAML 文件"""
        with open(path, "w") as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False, sort_keys=False)

    def validate(self) -> list:
        """验证配置有效性，返回警告列表"""
        warnings = []

        # 区域模拟验证
        if self.region.num_unique < 0:
            warnings.append("region.num_unique 不能为负")
        if self.region.num_multi < 0:
            warnings.append("region.num_multi 不能为负")
        if self.region.num_chimeric < 0:
            warnings.append("region.num_chimeric 不能为负")
        if (
            self.region.num_unique == 0
            and self.region.num_multi == 0
            and self.region.num_chimeric == 0
        ):
            warnings.append("至少需要生成一种 eccDNA 类型")

        # 长度分布验证
        ld = self.region.length_distribution
        if ld.min_length >= ld.max_length:
            warnings.append("min_length 必须小于 max_length")
        if not 0 <= ld.tail_weight <= 1:
            warnings.append("tail_weight 必须在 [0, 1] 范围内")
        if self.region.max_length_multi <= 0:
            warnings.append("region.max_length_multi 必须大于 0")
        if self.region.max_length_multi < ld.min_length:
            warnings.append("region.max_length_multi 必须大于等于 min_length")

        # 分类阈值验证
        c = self.region.classification
        if not 0 < c.identity <= 100:
            warnings.append("identity 必须在 (0, 100] 范围内")
        if not 0 < c.min_coverage <= 100:
            warnings.append("min_coverage 必须在 (0, 100] 范围内")

        # readsim 验证
        if self.readsim.mode not in ["on", "none"]:
            warnings.append(f"未知 readsim.mode: {self.readsim.mode}，仅支持 'on' 或 'none'")

        # 线程数验证
        if self.threads <= 0:
            warnings.append("threads 必须大于 0")

        return warnings
