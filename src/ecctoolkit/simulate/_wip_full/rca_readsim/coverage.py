"""
覆盖度分析与补充模块

实现：
1. 计算每个 eccDNA 的位点覆盖度分布
2. 识别覆盖度不足的区域
3. 生成补充 reads 以达到目标覆盖度
"""

import logging
from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np

from .models import EccDNA, SequencedRead, Segment, SegmentType, Strand

logger = logging.getLogger(__name__)


@dataclass
class CoverageStats:
    """单个 eccDNA 的覆盖度统计"""
    ecc_id: str
    length: int
    min_coverage: int = 0
    max_coverage: int = 0
    mean_coverage: float = 0.0
    median_coverage: float = 0.0
    zero_coverage_bases: int = 0  # 没有被覆盖的碱基数
    low_coverage_bases: int = 0   # 覆盖度 < 目标值的碱基数
    coverage_array: Optional[np.ndarray] = None

    @property
    def zero_coverage_ratio(self) -> float:
        return self.zero_coverage_bases / self.length if self.length > 0 else 0.0

    @property
    def is_fully_covered(self) -> bool:
        return self.zero_coverage_bases == 0

    def summary(self) -> str:
        return (
            f"{self.ecc_id}: len={self.length}, "
            f"cov={self.min_coverage}-{self.max_coverage} (mean={self.mean_coverage:.1f}), "
            f"uncovered={self.zero_coverage_bases} ({self.zero_coverage_ratio:.1%})"
        )


@dataclass
class CoverageReport:
    """整体覆盖度报告"""
    total_eccdnas: int = 0
    fully_covered: int = 0
    partially_covered: int = 0
    target_coverage: int = 1
    stats_by_ecc: Dict[str, CoverageStats] = field(default_factory=dict)

    @property
    def fully_covered_ratio(self) -> float:
        return self.fully_covered / self.total_eccdnas if self.total_eccdnas > 0 else 0.0

    def get_undercovered_eccs(self, min_coverage: int = 1) -> List[str]:
        """获取覆盖度不足的 eccDNA ID 列表"""
        result = []
        for ecc_id, stats in self.stats_by_ecc.items():
            if stats.min_coverage < min_coverage:
                result.append(ecc_id)
        return result

    def summary(self) -> str:
        return (
            f"Coverage Report: {self.fully_covered}/{self.total_eccdnas} "
            f"({self.fully_covered_ratio:.1%}) fully covered at {self.target_coverage}X"
        )


class CoverageAnalyzer:
    """覆盖度分析器"""

    def __init__(self, eccdnas: List[EccDNA]):
        """
        Args:
            eccdnas: eccDNA 列表
        """
        self.eccdnas = {e.id: e for e in eccdnas}
        self.ecc_lengths = {e.id: e.length for e in eccdnas}
        # 每个 eccDNA 的覆盖度数组
        self.coverage_arrays: Dict[str, np.ndarray] = {}
        for ecc in eccdnas:
            self.coverage_arrays[ecc.id] = np.zeros(ecc.length, dtype=np.int32)

    def add_read(self, read: SequencedRead) -> None:
        """添加一条 read 的覆盖度信息"""
        if read.is_background:
            return

        for seg in read.segments:
            self._add_segment(seg)

    def add_reads(self, reads: List[SequencedRead]) -> None:
        """批量添加 reads"""
        for read in reads:
            self.add_read(read)

    def _add_segment(self, seg: Segment) -> None:
        """添加一个 segment 的覆盖度"""
        ecc_id = seg.ecc_id
        if ecc_id not in self.coverage_arrays:
            return

        arr = self.coverage_arrays[ecc_id]
        ecc_len = len(arr)

        start = int(seg.ecc_offset)
        length = int(seg.length)

        # 处理环状结构：segment 可能跨越 eccDNA 的起点
        if start + length <= ecc_len:
            # 不跨越
            arr[start:start + length] += 1
        else:
            # 跨越起点
            arr[start:] += 1
            remaining = (start + length) - ecc_len
            # 可能绕多圈
            full_circles = remaining // ecc_len
            if full_circles > 0:
                arr[:] += full_circles
                remaining = remaining % ecc_len
            if remaining > 0:
                arr[:remaining] += 1

    def add_from_truth_records(self, records: List[dict]) -> None:
        """从 truth 记录添加覆盖度信息"""
        for rec in records:
            if rec.get("source_type") == "background":
                continue

            segments = rec.get("segments", [])
            if isinstance(segments, str):
                import json
                segments = json.loads(segments)

            for seg_dict in segments:
                ecc_id = seg_dict.get("ecc_id")
                if ecc_id not in self.coverage_arrays:
                    continue

                arr = self.coverage_arrays[ecc_id]
                ecc_len = len(arr)

                start = int(seg_dict.get("ecc_offset", 0))
                length = int(seg_dict.get("length", 0))

                if start + length <= ecc_len:
                    arr[start:start + length] += 1
                else:
                    arr[start:] += 1
                    remaining = (start + length) - ecc_len
                    full_circles = remaining // ecc_len
                    if full_circles > 0:
                        arr[:] += full_circles
                        remaining = remaining % ecc_len
                    if remaining > 0:
                        arr[:remaining] += 1

    def compute_stats(self, target_coverage: int = 1) -> CoverageReport:
        """计算覆盖度统计"""
        report = CoverageReport(
            total_eccdnas=len(self.eccdnas),
            target_coverage=target_coverage
        )

        for ecc_id, arr in self.coverage_arrays.items():
            stats = CoverageStats(
                ecc_id=ecc_id,
                length=len(arr),
                min_coverage=int(arr.min()),
                max_coverage=int(arr.max()),
                mean_coverage=float(arr.mean()),
                median_coverage=float(np.median(arr)),
                zero_coverage_bases=int(np.sum(arr == 0)),
                low_coverage_bases=int(np.sum(arr < target_coverage)),
                coverage_array=arr.copy()
            )
            report.stats_by_ecc[ecc_id] = stats

            if stats.min_coverage >= target_coverage:
                report.fully_covered += 1
            else:
                report.partially_covered += 1

        return report

    def get_uncovered_regions(
        self,
        ecc_id: str,
        min_coverage: int = 1
    ) -> List[Tuple[int, int]]:
        """
        获取覆盖度不足的区域

        Returns:
            [(start, end), ...] 列表，表示覆盖度 < min_coverage 的区域
        """
        if ecc_id not in self.coverage_arrays:
            return []

        arr = self.coverage_arrays[ecc_id]
        low_cov = arr < min_coverage

        regions = []
        in_region = False
        start = 0

        for i, is_low in enumerate(low_cov):
            if is_low and not in_region:
                start = i
                in_region = True
            elif not is_low and in_region:
                regions.append((start, i))
                in_region = False

        if in_region:
            regions.append((start, len(arr)))

        return regions


class SupplementaryReadGenerator:
    """补充 reads 生成器"""

    def __init__(
        self,
        eccdnas: List[EccDNA],
        ecc_db: Dict[str, EccDNA],
        rng: Optional[np.random.Generator] = None
    ):
        """
        Args:
            eccdnas: eccDNA 列表
            ecc_db: eccDNA 数据库 {id: EccDNA}
            rng: 随机数生成器
        """
        self.eccdnas = {e.id: e for e in eccdnas}
        self.ecc_db = ecc_db
        self.rng = rng if rng is not None else np.random.default_rng()

    def _get_sequence(self, ecc_id: str) -> Optional[str]:
        """获取 eccDNA 序列"""
        ecc = self.ecc_db.get(ecc_id)
        if ecc is None:
            return None
        # ecc_db 可能存储 EccDNA 对象或字符串
        if hasattr(ecc, 'seq'):
            return ecc.seq
        elif hasattr(ecc, 'sequence'):
            return ecc.sequence
        elif isinstance(ecc, str):
            return ecc
        return None

    def generate_targeted_reads(
        self,
        analyzer: CoverageAnalyzer,
        platform: str,
        target_coverage: int = 1,
        read_length_mean: int = 15000,
        read_length_std: int = 2000
    ) -> List[SequencedRead]:
        """
        生成针对性的补充 reads

        针对覆盖度不足的区域生成额外的 reads

        Args:
            analyzer: 覆盖度分析器
            platform: 测序平台 (HiFi, ONT, NGS)
            target_coverage: 目标覆盖度
            read_length_mean: read 长度均值
            read_length_std: read 长度标准差

        Returns:
            补充 reads 列表
        """
        supplementary_reads = []
        read_counter = 0

        for ecc_id, arr in analyzer.coverage_arrays.items():
            ecc_seq = self._get_sequence(ecc_id)
            if ecc_seq is None:
                continue

            ecc_len = len(ecc_seq)

            # 找出需要补充的位点
            need_supplement = arr < target_coverage
            supplement_count = target_coverage - arr
            supplement_count[supplement_count < 0] = 0

            if supplement_count.sum() == 0:
                continue

            # 计算需要多少条 reads 来覆盖
            # 简单策略：在需要补充的区域随机放置 reads
            uncovered_regions = analyzer.get_uncovered_regions(ecc_id, target_coverage)

            if not uncovered_regions:
                continue

            logger.info(
                f"{ecc_id}: {len(uncovered_regions)} uncovered regions, "
                f"total {int(supplement_count.sum())} bases need supplement"
            )

            # 对于每个未覆盖区域，生成足够的 reads
            for region_start, region_end in uncovered_regions:
                region_len = region_end - region_start

                # 计算需要多少 reads
                # 考虑到 read 长度通常远大于未覆盖区域，一条 read 可能覆盖多个区域
                reads_needed = max(1, int(np.ceil(
                    region_len * target_coverage / read_length_mean
                )))

                for _ in range(reads_needed):
                    # 生成 read 长度
                    read_len = max(
                        1000,
                        int(self.rng.normal(read_length_mean, read_length_std))
                    )
                    read_len = min(read_len, ecc_len * 3)  # 不超过 3 圈

                    # 选择起始位置：让 read 中心落在未覆盖区域
                    region_center = (region_start + region_end) // 2
                    read_start = region_center - read_len // 2
                    read_start = read_start % ecc_len  # 处理负数

                    # 生成序列
                    if read_len <= ecc_len:
                        if read_start + read_len <= ecc_len:
                            seq = ecc_seq[read_start:read_start + read_len]
                        else:
                            seq = ecc_seq[read_start:] + ecc_seq[:read_start + read_len - ecc_len]
                    else:
                        # 多圈
                        seq = ""
                        pos = read_start
                        remaining = read_len
                        while remaining > 0:
                            chunk_len = min(remaining, ecc_len - pos)
                            seq += ecc_seq[pos:pos + chunk_len]
                            remaining -= chunk_len
                            pos = 0

                    # 创建 segment
                    segment = Segment(
                        ecc_id=ecc_id,
                        ecc_offset=read_start,
                        length=read_len,
                        strand=Strand.FORWARD,
                        segment_type=SegmentType.TRUNK
                    )

                    # 创建 read
                    read_counter += 1
                    read = SequencedRead(
                        read_id=f"supplement_{platform.lower()}_{read_counter}",
                        sequence=seq,
                        quality="I" * len(seq),  # 高质量
                        platform=platform,
                        source_molecule_id=f"supplement_{ecc_id}",
                        source_ecc_ids=[ecc_id],
                        segments=[segment],
                        repeat_count_truth=read_len / ecc_len,
                        repeat_count_by_source={ecc_id: read_len / ecc_len}
                    )

                    supplementary_reads.append(read)

                    # 更新覆盖度数组
                    analyzer._add_segment(segment)

        return supplementary_reads

    def ensure_minimum_coverage(
        self,
        analyzer: CoverageAnalyzer,
        platform: str,
        target_coverage: int = 1,
        read_length_mean: int = 15000,
        read_length_std: int = 2000,
        max_iterations: int = 10
    ) -> Tuple[List[SequencedRead], CoverageReport]:
        """
        确保所有 eccDNA 达到最低覆盖度

        迭代生成补充 reads 直到所有位点达到目标覆盖度

        Args:
            analyzer: 覆盖度分析器
            platform: 测序平台
            target_coverage: 目标覆盖度
            read_length_mean: read 长度均值
            read_length_std: read 长度标准差
            max_iterations: 最大迭代次数

        Returns:
            (补充 reads 列表, 最终覆盖度报告)
        """
        all_supplementary = []

        for iteration in range(max_iterations):
            report = analyzer.compute_stats(target_coverage)

            if report.fully_covered == report.total_eccdnas:
                logger.info(
                    f"All {report.total_eccdnas} eccDNAs fully covered "
                    f"at {target_coverage}X after {iteration} iterations"
                )
                break

            logger.info(
                f"Iteration {iteration + 1}: "
                f"{report.fully_covered}/{report.total_eccdnas} fully covered, "
                f"generating supplementary reads..."
            )

            new_reads = self.generate_targeted_reads(
                analyzer, platform, target_coverage,
                read_length_mean, read_length_std
            )

            if not new_reads:
                logger.warning("No supplementary reads generated, stopping")
                break

            all_supplementary.extend(new_reads)
            logger.info(f"Generated {len(new_reads)} supplementary reads")

        final_report = analyzer.compute_stats(target_coverage)

        if final_report.fully_covered < final_report.total_eccdnas:
            undercovered = final_report.get_undercovered_eccs(target_coverage)
            logger.warning(
                f"After {max_iterations} iterations, "
                f"{len(undercovered)} eccDNAs still undercovered: "
                f"{undercovered[:5]}{'...' if len(undercovered) > 5 else ''}"
            )

        return all_supplementary, final_report


def verify_and_supplement_coverage(
    reads: List[SequencedRead],
    eccdnas: List[EccDNA],
    ecc_db: Dict[str, EccDNA],
    platform: str,
    target_coverage: int = 1,
    read_length_mean: int = 15000,
    read_length_std: int = 2000,
    rng: Optional[np.random.Generator] = None
) -> Tuple[List[SequencedRead], CoverageReport]:
    """
    验证覆盖度并生成补充 reads

    这是主要的入口函数

    Args:
        reads: 已生成的 reads
        eccdnas: eccDNA 列表
        ecc_db: eccDNA 数据库 {id: EccDNA}
        platform: 测序平台
        target_coverage: 目标覆盖度
        read_length_mean: read 长度均值
        read_length_std: read 长度标准差
        rng: 随机数生成器

    Returns:
        (所有 reads（原始 + 补充）, 覆盖度报告)
    """
    # 分析当前覆盖度
    analyzer = CoverageAnalyzer(eccdnas)
    analyzer.add_reads(reads)

    initial_report = analyzer.compute_stats(target_coverage)
    logger.info(f"Initial coverage: {initial_report.summary()}")

    if initial_report.fully_covered == initial_report.total_eccdnas:
        logger.info("All eccDNAs fully covered, no supplement needed")
        return reads, initial_report

    # 生成补充 reads
    generator = SupplementaryReadGenerator(eccdnas, ecc_db, rng)
    supplementary, final_report = generator.ensure_minimum_coverage(
        analyzer, platform, target_coverage,
        read_length_mean, read_length_std
    )

    logger.info(
        f"Coverage supplementation: added {len(supplementary)} reads, "
        f"final: {final_report.summary()}"
    )

    # 合并所有 reads
    all_reads = list(reads) + supplementary

    return all_reads, final_report


def downsample_reads(
    reads: List[SequencedRead],
    target_count: int,
    protected_ids: Optional[Set[str]] = None,
    rng: Optional[np.random.Generator] = None
) -> List[SequencedRead]:
    """
    智能下采样 reads

    优先保留被保护的 reads（如补充 reads），按长度加权随机移除其他 reads

    Args:
        reads: 原始 reads 列表
        target_count: 目标 read 数量
        protected_ids: 被保护的 read IDs（不会被移除）
        rng: 随机数生成器

    Returns:
        下采样后的 reads 列表
    """
    if len(reads) <= target_count:
        return reads

    if rng is None:
        rng = np.random.default_rng()

    if protected_ids is None:
        protected_ids = set()

    # 分离被保护的和可移除的 reads
    protected = [r for r in reads if r.read_id in protected_ids]
    removable = [r for r in reads if r.read_id not in protected_ids]

    # 如果被保护的 reads 已经超过目标，直接返回保护的
    if len(protected) >= target_count:
        logger.warning(
            f"Downsample: protected reads ({len(protected)}) >= target ({target_count}), "
            "returning protected reads only"
        )
        return protected

    # 计算需要保留的可移除 reads 数量
    keep_count = target_count - len(protected)

    if keep_count >= len(removable):
        return reads

    # 按长度加权采样保留
    lengths = np.array([len(r.sequence) for r in removable], dtype=float)
    weights = lengths / lengths.sum()

    keep_indices = rng.choice(
        len(removable),
        size=keep_count,
        replace=False,
        p=weights
    )

    kept_removable = [removable[i] for i in sorted(keep_indices)]

    result = protected + kept_removable

    logger.info(
        f"Downsample: {len(reads)} -> {len(result)} reads "
        f"(protected: {len(protected)}, kept: {len(kept_removable)})"
    )

    return result


def compute_total_coverage(
    reads: List[SequencedRead],
    eccdnas: List[EccDNA]
) -> float:
    """
    计算总体平均覆盖度

    Args:
        reads: reads 列表
        eccdnas: eccDNA 列表

    Returns:
        平均覆盖度 (总 read bases / 总 eccDNA length)
    """
    total_ecc_length = sum(e.length for e in eccdnas)
    if total_ecc_length == 0:
        return 0.0

    # 只计算 eccDNA reads，不包括背景
    total_read_bases = sum(
        len(r.sequence) for r in reads
        if not r.is_background
    )

    return total_read_bases / total_ecc_length


class AdaptiveCoverageSampler:
    """
    两阶段自适应覆盖度采样器

    实现流程：
    1. 初始采样 (initial_ratio * target_coverage)
    2. 分析覆盖度，补充到 min_coverage
    3. 计算差值，继续采样或下采样到目标
    """

    def __init__(
        self,
        eccdnas: List[EccDNA],
        ecc_db: Dict[str, EccDNA],
        target_coverage: float = 25.0,
        min_coverage: float = 3.0,
        rng: Optional[np.random.Generator] = None
    ):
        """
        Args:
            eccdnas: eccDNA 列表
            ecc_db: eccDNA 数据库 {id: EccDNA}
            target_coverage: 目标平均覆盖度
            min_coverage: 每个 eccDNA 的最小覆盖度
            rng: 随机数生成器
        """
        self.eccdnas = eccdnas
        self.ecc_db = ecc_db
        self.target_coverage = target_coverage
        self.min_coverage = min_coverage
        self.rng = rng if rng is not None else np.random.default_rng()

        # 计算总 eccDNA 长度
        self.total_ecc_length = sum(e.length for e in eccdnas)

    def compute_target_reads(self, coverage: float, read_length_mean: int) -> int:
        """根据目标覆盖度计算需要的 read 数量"""
        if self.total_ecc_length == 0 or read_length_mean == 0:
            return 0
        return max(1, int(coverage * self.total_ecc_length / read_length_mean))

    def _smart_downsample(
        self,
        reads: List[SequencedRead],
        protected_ids: Set[str],
        target_coverage: float,
        read_length_mean: int
    ) -> Tuple[List[SequencedRead], int]:
        """
        智能下采样：优先从高覆盖度 eccDNA 移除 reads

        Args:
            reads: 当前 reads 列表
            protected_ids: 被保护的 read IDs（补充 reads，不会被移除）
            target_coverage: 目标覆盖度
            read_length_mean: 平均 read 长度

        Returns:
            (下采样后的 reads, 移除的 reads 数量)
        """
        # 计算每个 eccDNA 的当前覆盖度
        ecc_coverage = {ecc.id: 0.0 for ecc in self.eccdnas}
        ecc_reads = {ecc.id: [] for ecc in self.eccdnas}

        for read in reads:
            if read.read_id in protected_ids:
                continue  # 保护的 reads 不参与下采样
            for ecc_id in read.source_ecc_ids:
                if ecc_id in ecc_coverage:
                    # 按比例分配覆盖度（read 可能跨多个 eccDNA）
                    contribution = len(read.sequence) / len(read.source_ecc_ids)
                    ecc = self.ecc_db.get(ecc_id)
                    if ecc and ecc.length > 0:
                        ecc_coverage[ecc_id] += contribution / ecc.length
                        ecc_reads[ecc_id].append(read)

        # 按覆盖度降序排列 eccDNA
        sorted_eccs = sorted(
            ecc_coverage.keys(),
            key=lambda x: ecc_coverage[x],
            reverse=True
        )

        current_coverage = compute_total_coverage(reads, self.eccdnas)
        reads_to_remove = set()
        tolerance = 0.5

        # 从高覆盖度 eccDNA 移除 reads，直到接近目标
        for ecc_id in sorted_eccs:
            if current_coverage <= target_coverage + tolerance:
                break

            # 从这个 eccDNA 的 reads 中选择要移除的
            available_reads = [r for r in ecc_reads[ecc_id] if r.read_id not in reads_to_remove]
            if not available_reads:
                continue

            # 优先移除较长的 reads（贡献更多覆盖度）
            available_reads.sort(key=lambda r: len(r.sequence), reverse=True)

            for read in available_reads:
                if current_coverage <= target_coverage + tolerance:
                    break
                reads_to_remove.add(read.read_id)
                # 估算移除后的覆盖度
                current_coverage -= len(read.sequence) / self.total_ecc_length

        # 构建最终的 reads 列表
        final_reads = [r for r in reads if r.read_id not in reads_to_remove]

        return final_reads, len(reads_to_remove)

    def run_adaptive_sampling(
        self,
        sample_func,
        platform: str,
        read_length_mean: int,
        read_length_std: int = 2000
    ) -> Tuple[List[SequencedRead], CoverageReport, Dict]:
        """
        执行两阶段自适应采样

        流程：
        1. 直接采样目标覆盖度的 reads（确定性计算）
        2. 分析覆盖度，补充未达到 min_coverage 的 eccDNA

        Args:
            sample_func: 采样函数，签名 (count: int) -> List[SequencedRead]
            platform: 测序平台 (NGS, HiFi, ONT)
            read_length_mean: read 平均长度
            read_length_std: read 长度标准差

        Returns:
            (final_reads, coverage_report, stats_dict)
        """
        stats = {
            "target_coverage": self.target_coverage,
            "min_coverage": self.min_coverage,
            "phase1_reads": 0,
            "phase2_supplement_reads": 0,
            "final_reads": 0,
            "final_coverage": 0.0,
        }

        # 阶段1：直接采样目标覆盖度
        target_reads = self.compute_target_reads(self.target_coverage, read_length_mean)

        logger.info(
            f"Adaptive sampling phase 1: sampling for {self.target_coverage}X "
            f"(~{target_reads} reads)"
        )

        phase1_reads = sample_func(target_reads)
        stats["phase1_reads"] = len(phase1_reads)

        current_coverage = compute_total_coverage(phase1_reads, self.eccdnas)
        logger.info(
            f"After phase 1: {len(phase1_reads)} reads, "
            f"coverage {current_coverage:.1f}X"
        )

        # 阶段2：分析覆盖度并补充到 min_coverage
        logger.info(
            f"Adaptive sampling phase 2: ensure each eccDNA >= {self.min_coverage}X"
        )

        analyzer = CoverageAnalyzer(self.eccdnas)
        analyzer.add_reads(phase1_reads)

        initial_report = analyzer.compute_stats(int(self.min_coverage))
        logger.info(f"Coverage check: {initial_report.summary()}")

        # 生成补充 reads（仅针对覆盖不足的 eccDNA）
        supplement_reads = []

        if initial_report.fully_covered < initial_report.total_eccdnas:
            generator = SupplementaryReadGenerator(self.eccdnas, self.ecc_db, self.rng)
            supplement_reads, _ = generator.ensure_minimum_coverage(
                analyzer,
                platform,
                target_coverage=int(self.min_coverage),
                read_length_mean=read_length_mean,
                read_length_std=read_length_std
            )
            logger.info(f"Added {len(supplement_reads)} supplement reads for undercovered eccDNAs")

        stats["phase2_supplement_reads"] = len(supplement_reads)

        # 合并 Phase 1 和 Phase 2 的 reads
        current_reads = phase1_reads + supplement_reads
        supplement_ids = {r.read_id for r in supplement_reads}
        current_coverage = compute_total_coverage(current_reads, self.eccdnas)

        logger.info(
            f"After phase 2: {len(current_reads)} reads, coverage {current_coverage:.1f}X"
        )

        # 阶段3：迭代调整到目标覆盖度
        tolerance = 1.0  # 允许 ±1X 误差
        stats["phase3_sampled"] = 0
        stats["phase3_downsampled"] = 0

        logger.info(f"Adaptive sampling phase 3: adjusting to target {self.target_coverage}X (±{tolerance}X)")

        max_iterations = 5
        for iteration in range(max_iterations):
            diff = current_coverage - self.target_coverage

            if abs(diff) <= tolerance:
                logger.info(f"  Target reached: {current_coverage:.1f}X")
                break

            if diff > tolerance:
                # 高于目标：智能下采样（优先从高覆盖度 eccDNA 移除）
                logger.info(f"  Iteration {iteration+1}: {current_coverage:.1f}X > target, downsampling...")
                current_reads, removed = self._smart_downsample(
                    current_reads,
                    supplement_ids,
                    self.target_coverage,
                    read_length_mean
                )
                stats["phase3_downsampled"] += removed
                current_coverage = compute_total_coverage(current_reads, self.eccdnas)
                logger.info(f"    Removed {removed} reads -> {current_coverage:.1f}X")

            elif diff < -tolerance:
                # 低于目标：继续随机采样
                remaining = self.target_coverage - current_coverage
                # 保守估计，避免过度采样
                additional_needed = max(1, self.compute_target_reads(remaining * 0.8, read_length_mean))
                logger.info(f"  Iteration {iteration+1}: {current_coverage:.1f}X < target, sampling {additional_needed} more...")
                additional_reads = sample_func(additional_needed)
                if not additional_reads:
                    logger.warning("    No more reads available")
                    break
                stats["phase3_sampled"] += len(additional_reads)
                current_reads = current_reads + additional_reads
                current_coverage = compute_total_coverage(current_reads, self.eccdnas)
                logger.info(f"    Added {len(additional_reads)} reads -> {current_coverage:.1f}X")

        if abs(current_coverage - self.target_coverage) > tolerance:
            logger.info(f"  Final coverage: {current_coverage:.1f}X (target: {self.target_coverage}X)")

        # 最终统计
        final_coverage = compute_total_coverage(current_reads, self.eccdnas)
        final_analyzer = CoverageAnalyzer(self.eccdnas)
        final_analyzer.add_reads(current_reads)
        final_report = final_analyzer.compute_stats(int(self.min_coverage))

        stats["final_reads"] = len(current_reads)
        stats["final_coverage"] = final_coverage

        logger.info(
            f"Adaptive sampling complete: {len(current_reads)} reads, "
            f"coverage {final_coverage:.2f}X, {final_report.summary()}"
        )

        return current_reads, final_report, stats
