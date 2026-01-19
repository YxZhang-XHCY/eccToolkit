#!/usr/bin/env python3
"""
Region Library Builder - 构建高质量的 eccDNA 区域库

设计思路：
1. 生成 5X 目标数量的候选区域
2. 用 TideHunter 检测模拟的串联重复序列
3. 如果 consensus 长度 ≈ 区域长度 → 干净区域
4. 如果 consensus 长度 << 区域长度 → 有内部串联重复，排除
5. 从干净区域中抽取需要的数量
"""

import logging
import subprocess
import tempfile
import os
import shutil
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Set
from pathlib import Path
import numpy as np

from .eccDNA_simulator import (
    SimulationConfig,
    EccRegion,
    RegionGenerator,
    Minimap2Manager,
    LengthSampler,
)


@dataclass
class LibraryConfig:
    """区域库配置"""
    # 候选倍数（生成 N 倍目标数量）
    candidate_multiplier: int = 5

    # TideHunter 检测参数
    tidehunter_copy_count: int = 10  # 模拟的串联重复次数
    tidehunter_min_copy: int = 2     # TideHunter -c 参数

    # 内部串联重复阈值
    # 如果 TideHunter consensus 长度 / 区域长度 < threshold → 有内部重复
    consensus_ratio_threshold: float = 0.8

    # 并行线程数
    threads: int = 12

    # 输出目录
    output_dir: str = "region_library"


@dataclass
class RegionLibrary:
    """验证过的区域库"""
    unique_regions: List[EccRegion] = field(default_factory=list)
    multi_regions: List[EccRegion] = field(default_factory=list)
    chimeric_regions: List[EccRegion] = field(default_factory=list)

    # 统计信息
    stats: Dict = field(default_factory=dict)


class TideHunterChecker:
    """使用 TideHunter 检测内部串联重复"""

    def __init__(self, config: LibraryConfig):
        self.config = config
        self.work_dir: Optional[str] = None

    def __enter__(self):
        self.work_dir = tempfile.mkdtemp(prefix="th_check_")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.work_dir and os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)

    def check_regions(
        self,
        regions: List[EccRegion],
        sequences: Dict[str, str],
    ) -> Tuple[List[EccRegion], List[EccRegion], Dict]:
        """
        检查区域是否有内部串联重复

        方法：模拟串联重复序列，用 TideHunter 检测
        如果 consensus 长度 ≈ 区域长度 → 干净
        如果 consensus 长度 << 区域长度 → 有内部重复

        Returns:
            (clean_regions, problematic_regions, stats)
        """
        if not regions:
            return [], [], {}

        assert self.work_dir is not None

        # 生成模拟的串联重复序列
        fasta_path = os.path.join(self.work_dir, "tandem_reads.fa")
        region_lengths = {}

        with open(fasta_path, 'w') as f:
            for r in regions:
                seq = r.get_sequence(sequences)
                region_lengths[r.region_id] = len(seq)
                # 生成串联重复序列（模拟 eccDNA read）
                tandem_seq = seq * self.config.tidehunter_copy_count
                f.write(f">{r.region_id}\n{tandem_seq}\n")

        # 运行 TideHunter
        th_output = self._run_tidehunter(fasta_path)

        # 解析结果，找出有问题的区域
        problematic_ids = self._parse_tidehunter_output(th_output, region_lengths)

        # 分类
        clean = []
        problematic = []
        for r in regions:
            if r.region_id in problematic_ids:
                problematic.append(r)
            else:
                clean.append(r)

        stats = {
            'total': len(regions),
            'clean': len(clean),
            'problematic': len(problematic),
            'problematic_ratio': len(problematic) / len(regions) if regions else 0,
        }

        return clean, problematic, stats

    def _run_tidehunter(self, fasta_path: str) -> str:
        """运行 TideHunter"""
        output_path = os.path.join(self.work_dir, "th_output.tsv")

        cmd = [
            "TideHunter",
            "-f", "2",  # 输出格式：tabular
            "-c", str(self.config.tidehunter_min_copy),
            "-o", output_path,
            fasta_path,
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,
            )
            if os.path.exists(output_path):
                with open(output_path) as f:
                    return f.read()
            return ""
        except subprocess.TimeoutExpired:
            logging.warning("TideHunter timeout")
            return ""
        except FileNotFoundError:
            logging.error("TideHunter not found!")
            raise

    def _parse_tidehunter_output(
        self,
        th_output: str,
        region_lengths: Dict[str, int],
    ) -> Set[str]:
        """
        解析 TideHunter 输出，找出有问题的区域

        TideHunter -f 2 输出格式（tab分隔）：
        read_name rep_id copy_num read_len start end cons_len ... consensus_seq
        """
        problematic_ids = set()
        best_consensus = {}  # region_id -> max consensus length

        for line in th_output.strip().split('\n'):
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 8:
                continue

            try:
                region_id = parts[0]
                consensus_len = len(parts[-1]) if parts[-1] else int(parts[6])

                # 记录最长的 consensus
                if region_id not in best_consensus or consensus_len > best_consensus[region_id]:
                    best_consensus[region_id] = consensus_len
            except (ValueError, IndexError):
                continue

        # 检查每个区域
        for region_id, consensus_len in best_consensus.items():
            region_len = region_lengths.get(region_id, 0)
            if region_len > 0:
                ratio = consensus_len / region_len
                if ratio < self.config.consensus_ratio_threshold:
                    # consensus 远小于区域长度 → 有内部串联重复
                    problematic_ids.add(region_id)
                    logging.debug(
                        f"Region {region_id} has internal repeat: "
                        f"consensus={consensus_len}bp, region={region_len}bp, "
                        f"ratio={ratio:.2f}"
                    )

        # 检查没有 TideHunter 输出的区域（可能太短或其他原因）
        for region_id in region_lengths:
            if region_id not in best_consensus:
                # TideHunter 没有检测到任何重复 - 可能有问题
                logging.debug(f"Region {region_id} not detected by TideHunter")
                # 不标记为问题，因为可能只是序列太短

        return problematic_ids


class RegionLibraryBuilder:
    """区域库构建器"""

    def __init__(
        self,
        sim_config: SimulationConfig,
        lib_config: LibraryConfig,
    ):
        self.sim_config = sim_config
        self.lib_config = lib_config
        self.rng = np.random.default_rng(sim_config.seed)

        # 组件
        self.generator: Optional[RegionGenerator] = None
        self.aligner: Optional[Minimap2Manager] = None
        self.sampler: Optional[LengthSampler] = None
        self.sequences: Optional[Dict[str, str]] = None

    def build(
        self,
        target_unique: int,
        target_multi: int,
        target_chimeric: int,
    ) -> RegionLibrary:
        """构建区域库"""
        logging.info("=" * 60)
        logging.info("Building Region Library")
        logging.info("=" * 60)
        logging.info(f"Target: {target_unique} U, {target_multi} M, {target_chimeric} C")
        logging.info(f"Candidate multiplier: {self.lib_config.candidate_multiplier}x")

        # 初始化组件
        self._init_components()

        library = RegionLibrary()

        # 1. 生成并验证 Unique 区域
        candidate_u = target_unique * self.lib_config.candidate_multiplier
        logging.info(f"\n[1/3] Generating {candidate_u} Unique candidates (target: {target_unique})...")
        library.unique_regions, u_stats = self._build_unique_pool(candidate_u, target_unique)
        library.stats['unique'] = u_stats

        # 2. 生成并验证 Multi 区域
        candidate_m = target_multi * self.lib_config.candidate_multiplier
        logging.info(f"\n[2/3] Generating {candidate_m} Multi candidates (target: {target_multi})...")
        library.multi_regions, m_stats = self._build_multi_pool(candidate_m, target_multi)
        library.stats['multi'] = m_stats

        # 3. 生成并验证 Chimeric 区域
        candidate_c = target_chimeric * self.lib_config.candidate_multiplier
        logging.info(f"\n[3/3] Generating {candidate_c} Chimeric candidates (target: {target_chimeric})...")
        library.chimeric_regions, c_stats = self._build_chimeric_pool(
            candidate_c, target_chimeric, library.unique_regions
        )
        library.stats['chimeric'] = c_stats

        # 输出统计
        self._print_summary(library)

        return library

    def _init_components(self):
        """初始化组件"""
        from pyfaidx import Fasta

        logging.info(f"Loading reference: {self.sim_config.reference}")
        fasta = Fasta(self.sim_config.reference)
        self.sequences = {name: str(fasta[name][:]) for name in fasta.keys()}

        self.sampler = LengthSampler(self.sim_config, self.rng)
        self.generator = RegionGenerator(self.sim_config, self.sampler, self.sequences, self.rng)
        self.aligner = Minimap2Manager(self.sim_config)

    def _build_unique_pool(self, num_candidates: int, target: int) -> Tuple[List[EccRegion], Dict]:
        """构建 Unique 区域池"""
        assert self.generator is not None
        assert self.aligner is not None
        assert self.sequences is not None

        stats = {'generated': 0, 'classified_unique': 0, 'th_clean': 0, 'th_problematic': 0, 'final': 0}

        # 批量生成候选
        batch_size = min(20000, num_candidates)
        all_unique: List[EccRegion] = []
        total_generated = 0

        while len(all_unique) < num_candidates and total_generated < num_candidates * 2:
            remaining = num_candidates - len(all_unique)
            current_batch = min(batch_size, remaining + 1000)

            candidates = self.generator.generate_candidates(current_batch, id_prefix="U_cand")
            total_generated += len(candidates)

            if not candidates:
                break

            # 分类
            unique, _ = self.aligner.classify_regions(candidates, self.sequences)
            stats['classified_unique'] += len(unique)
            all_unique.extend(unique)

            logging.info(f"  Batch: generated {len(candidates)}, got {len(unique)} unique, total {len(all_unique)}")

        stats['generated'] = total_generated

        # TideHunter 检查
        logging.info(f"  Running TideHunter check on {len(all_unique)} unique regions...")
        with TideHunterChecker(self.lib_config) as checker:
            clean, problematic, th_stats = checker.check_regions(all_unique, self.sequences)

        stats['th_clean'] = len(clean)
        stats['th_problematic'] = len(problematic)

        if problematic:
            logging.info(f"  Removed {len(problematic)} regions with internal tandem repeats")

        # 取目标数量
        final = clean[:target] if len(clean) >= target else clean

        # 重新编号
        for i, r in enumerate(final):
            r.region_id = f"UeccDNA_{i+1:06d}"
            r.ecc_type = 'U'

        stats['final'] = len(final)
        logging.info(f"  Final: {len(final)} / {target} unique regions")

        return final, stats

    def _build_multi_pool(self, num_candidates: int, target: int) -> Tuple[List[EccRegion], Dict]:
        """构建 Multi 区域池"""
        assert self.generator is not None
        assert self.aligner is not None
        assert self.sequences is not None

        stats = {'generated': 0, 'classified_multi': 0, 'th_clean': 0, 'th_problematic': 0, 'final': 0}

        # Multi 难找，需要更大的候选池
        batch_size = 20000
        all_multi: List[EccRegion] = []
        total_generated = 0
        max_attempts = num_candidates * 20

        while len(all_multi) < num_candidates and total_generated < max_attempts:
            remaining = num_candidates - len(all_multi)
            current_batch = min(batch_size, remaining * 5)

            candidates = self.generator.generate_candidates(current_batch, id_prefix="M_cand")
            total_generated += len(candidates)

            if not candidates:
                break

            # 分类
            _, multi = self.aligner.classify_regions(candidates, self.sequences)

            # 过滤超长的 Multi
            multi = [r for r in multi if r.length <= self.sim_config.max_length_multi]
            stats['classified_multi'] += len(multi)
            all_multi.extend(multi)

            logging.info(f"  Batch: generated {len(candidates)}, got {len(multi)} multi, total {len(all_multi)}")

        stats['generated'] = total_generated

        # TideHunter 检查 - 这是关键步骤！
        logging.info(f"  Running TideHunter check on {len(all_multi)} multi regions...")
        with TideHunterChecker(self.lib_config) as checker:
            clean, problematic, th_stats = checker.check_regions(all_multi, self.sequences)

        stats['th_clean'] = len(clean)
        stats['th_problematic'] = len(problematic)

        if problematic:
            logging.info(f"  *** Removed {len(problematic)} multi regions with internal tandem repeats ***")

        # 取目标数量
        final = clean[:target] if len(clean) >= target else clean

        # 重新编号
        for i, r in enumerate(final):
            r.region_id = f"MeccDNA_{i+1:06d}"
            r.ecc_type = 'M'

        stats['final'] = len(final)
        logging.info(f"  Final: {len(final)} / {target} multi regions")

        return final, stats

    def _build_chimeric_pool(
        self,
        num_candidates: int,
        target: int,
        unique_pool: List[EccRegion],
    ) -> Tuple[List[EccRegion], Dict]:
        """构建 Chimeric 区域池"""
        assert self.sequences is not None

        stats = {'generated': 0, 'th_clean': 0, 'th_problematic': 0, 'final': 0}

        if len(unique_pool) < 10:
            logging.warning("Not enough unique regions for chimeric generation")
            return [], stats

        # 按染色体分组
        by_chrom: Dict[str, List[EccRegion]] = {}
        for r in unique_pool:
            by_chrom.setdefault(r.chrom, []).append(r)

        chroms = list(by_chrom.keys())

        chimeric_regions: List[EccRegion] = []

        for i in range(num_candidates):
            # 随机选择 2-4 个片段
            num_frags = self.rng.choice([2, 2, 2, 3, 3, 4])

            # 选择片段
            if self.rng.random() < 0.3 and len(chroms) > 0:
                # 同一染色体
                chrom = self.rng.choice(chroms)
                if len(by_chrom[chrom]) < num_frags:
                    continue
                frags = list(self.rng.choice(by_chrom[chrom], size=num_frags, replace=False))
            else:
                # 不同染色体
                frags = []
                available = list(unique_pool)
                self.rng.shuffle(available)
                for r in available:
                    if len(frags) >= num_frags:
                        break
                    frags.append(r)

            if len(frags) < 2:
                continue

            # 创建 chimeric region
            fragments = [(f.chrom, f.start, f.end) for f in frags]
            total_length = sum(f.length for f in frags)

            cecc = EccRegion(
                chrom=frags[0].chrom,
                start=frags[0].start,
                end=frags[0].end,
                strand='+',
                length=total_length,
                region_id=f"C_cand_{i+1:06d}",
                ecc_type='C',
                fragments=fragments,
            )
            chimeric_regions.append(cecc)

        stats['generated'] = len(chimeric_regions)

        # TideHunter 检查
        logging.info(f"  Running TideHunter check on {len(chimeric_regions)} chimeric regions...")
        with TideHunterChecker(self.lib_config) as checker:
            clean, problematic, th_stats = checker.check_regions(chimeric_regions, self.sequences)

        stats['th_clean'] = len(clean)
        stats['th_problematic'] = len(problematic)

        # 取目标数量
        final = clean[:target] if len(clean) >= target else clean

        # 重新编号
        for i, r in enumerate(final):
            r.region_id = f"CeccDNA_{i+1:06d}"

        stats['final'] = len(final)
        logging.info(f"  Final: {len(final)} / {target} chimeric regions")

        return final, stats

    def _print_summary(self, library: RegionLibrary):
        """打印摘要"""
        logging.info("\n" + "=" * 60)
        logging.info("Region Library Summary")
        logging.info("=" * 60)

        for ecc_type, regions in [
            ('Unique', library.unique_regions),
            ('Multi', library.multi_regions),
            ('Chimeric', library.chimeric_regions),
        ]:
            if regions:
                lengths = [r.length for r in regions]
                logging.info(f"\n{ecc_type}: {len(regions)} regions")
                logging.info(f"  Length: min={min(lengths)}, max={max(lengths)}, mean={np.mean(lengths):.0f}")

        # 统计被过滤的比例
        for ecc_type in ['unique', 'multi', 'chimeric']:
            if ecc_type in library.stats:
                s = library.stats[ecc_type]
                if s.get('th_problematic', 0) > 0:
                    total = s.get('th_clean', 0) + s.get('th_problematic', 0)
                    pct = s['th_problematic'] / total * 100 if total > 0 else 0
                    logging.info(f"\n{ecc_type.capitalize()}: {s['th_problematic']} / {total} ({pct:.1f}%) had internal repeats")

        logging.info("\n" + "=" * 60)

    def save_library(self, library: RegionLibrary, output_dir: str):
        """保存区域库到文件"""
        os.makedirs(output_dir, exist_ok=True)

        # BED 文件
        for ecc_type, regions, filename in [
            ('U', library.unique_regions, 'unique.bed'),
            ('M', library.multi_regions, 'multi.bed'),
            ('C', library.chimeric_regions, 'chimeric.bed'),
        ]:
            if regions:
                path = os.path.join(output_dir, filename)
                with open(path, 'w') as f:
                    f.write("#chrom\tstart\tend\tname\tlength\tstrand\ttype\tfragments\n")
                    for r in regions:
                        frag_str = ';'.join([f"{c}:{s}-{e}" for c, s, e in r.fragments]) if r.fragments else '.'
                        f.write(f"{r.chrom}\t{r.start}\t{r.end}\t{r.region_id}\t{r.length}\t{r.strand}\t{r.ecc_type}\t{frag_str}\n")
                logging.info(f"Saved {len(regions)} regions to {path}")

        # FASTA 文件
        if self.sequences:
            for ecc_type, regions, filename in [
                ('U', library.unique_regions, 'unique.fa'),
                ('M', library.multi_regions, 'multi.fa'),
                ('C', library.chimeric_regions, 'chimeric.fa'),
            ]:
                if regions:
                    path = os.path.join(output_dir, filename)
                    with open(path, 'w') as f:
                        for r in regions:
                            seq = r.get_sequence(self.sequences)
                            if r.fragments:
                                frag_str = ';'.join([f"{c}:{s}-{e}" for c, s, e in r.fragments])
                                header = f">{r.region_id} {frag_str} length={r.length}"
                            else:
                                header = f">{r.region_id} {r.chrom}:{r.start}-{r.end} length={r.length}"
                            f.write(header + "\n")
                            for i in range(0, len(seq), 80):
                                f.write(seq[i:i+80] + "\n")
                    logging.info(f"Saved {len(regions)} sequences to {path}")


def main():
    """命令行入口"""
    import argparse

    parser = argparse.ArgumentParser(description="Build eccDNA Region Library")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome FASTA")
    parser.add_argument("-o", "--output", default="region_library", help="Output directory")
    parser.add_argument("--num-unique", type=int, default=1000, help="Target number of unique regions")
    parser.add_argument("--num-multi", type=int, default=100, help="Target number of multi regions")
    parser.add_argument("--num-chimeric", type=int, default=100, help="Target number of chimeric regions")
    parser.add_argument("--multiplier", type=int, default=5, help="Candidate multiplier (default: 5x)")
    parser.add_argument("-t", "--threads", type=int, default=12, help="Number of threads")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
    )

    # 配置
    sim_config = SimulationConfig(
        reference=args.reference,
        output_prefix=os.path.join(args.output, "library"),
        seed=args.seed,
        threads=args.threads,
        num_unique=0,
        num_multi=0,
        num_chimeric=0,
    )

    lib_config = LibraryConfig(
        candidate_multiplier=args.multiplier,
        threads=args.threads,
        output_dir=args.output,
    )

    # 构建
    builder = RegionLibraryBuilder(sim_config, lib_config)
    library = builder.build(
        target_unique=args.num_unique,
        target_multi=args.num_multi,
        target_chimeric=args.num_chimeric,
    )
    builder.save_library(library, args.output)


if __name__ == "__main__":
    main()
