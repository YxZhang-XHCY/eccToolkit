#!/usr/bin/env python3
"""
eccDNA Region Simulator v4 (minimap2 edition) - Final
======================================================

从参考基因组中模拟 eccDNA 区域，支持三种类型：
  - UeccDNA: 唯一比对 (Unique) - 在基因组上只有一处全长(长度一致>=99%) >99% 匹配
  - MeccDNA: 多处比对 (Multi-mapped) - 在基因组上存在多处全长(长度一致>=99%) >99% 匹配
  - CeccDNA: 嵌合体 (Chimeric) - 由多个片段拼接而成

特性：
  - 使用 minimap2 进行快速比对（比 BLAST 快 10-100x）
  - 按长度分流：短序列用 -x sr；长序列用 -x asm5
  - 索引复用：默认在 reference 所在目录寻找/复用 .mmi
  - 严格 identity 阈值：minimap2 加 -c 输出 CIGAR
  - self-hit 用 overlap fraction 过滤（比坐标误差更稳）
  - no-hit 默认 skip（更保守；可 --no-hit-policy unique）

用法示例：
  ecc sim-region \\
      -r ColCEN.fasta \\
      -u 700 -m 20 -c 10 \\
      -o sim_ecc \\
      -t 12 \\
      --seed 42

Author: Claude
Version: 4.0 Final
"""

import argparse
import logging
import math
import os
import gzip
import subprocess
import shutil
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
from pathlib import Path
import numpy as np


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class SimulationConfig:
    """模拟配置参数"""
    # 输入输出
    reference: str
    output_prefix: str
    output_dir: Optional[str] = None  # 自定义输出目录，None 则用 prefix 创建

    # 目标数量
    num_unique: int = 1000
    num_multi: int = 0
    num_chimeric: int = 0

    # 随机种子
    seed: int = 42

    # Lognormal 主分布参数
    mode: float = 400.0
    sigma: float = 0.8

    # 尾部分布参数
    tail_weight: float = 0.05
    tail_min: int = 5000
    tail_max: int = 500000

    # 全局约束
    min_length: int = 100
    max_length: int = 500000
    max_length_multi: int = 10000  # MeccDNA hard length cap (bp)

    # 重试设置
    max_retry: int = 100

    # minimap2 设置
    threads: int = 4
    identity_threshold: float = 99.0  # percent
    min_coverage: float = 90.0        # percent
    length_consistency_threshold: float = 99.0  # percent
    # Multi-mapping classification threshold (CircleSeeker theta_m ≈ 95%).
    # A candidate is only classified as Multi if it has at least one non-self hit
    # with coverage >= this threshold.
    multi_coverage_threshold: float = 95.0  # percent
    max_secondary: int = 50

    # preset 策略（使用 map-hifi 与 CircleSeeker HiFi 检测保持一致）
    minimap_preset_short: str = "map-hifi"
    minimap_preset_long: str = "map-hifi"
    split_by_length: bool = False
    split_length: int = 10000  # >= split_length 走 long preset

    # 无命中处理：skip 更保守
    no_hit_policy: str = "skip"  # skip|unique

    # 索引复用设置
    index_dir: Optional[str] = None  # None => 使用 reference 所在目录
    reuse_index: bool = True
    force_rebuild_index: bool = False

    # 嵌合体设置
    chimeric_fragments_weights: Dict[int, float] = field(
        default_factory=lambda: {2: 0.80, 3: 0.15, 4: 0.04, 5: 0.01}
    )
    enforce_unique_chimeric_fragments: bool = True
    chimeric_unique_chrom_fraction: float = 0.7
    chimeric_same_chrom_probability: float = 0.2
    chimeric_fragment_repeat_probability: float = 0.05
    chimeric_fragment_max_repeats: int = 1
    chimeric_fragment_pool_multiplier: float = 5.0
    chimeric_fragment_min_length: int = 100
    chimeric_fragment_max_length: int = 5000

    # 候选倍数
    candidate_multiplier_u: float = 2.0
    candidate_multiplier_m: float = 10.0

    # 是否保留临时文件
    keep_tmp: bool = False

    # 排除的染色体
    exclude_chroms: List[str] = field(
        default_factory=lambda: ['chrM', 'chrUn', '_random', '_alt', '_hap']
    )

    # 自身命中过滤阈值
    self_overlap_fraction: float = 0.90

    # ML 训练数据生成
    generate_ml_data: bool = False  # 是否生成 ML 训练数据
    ml_identity_threshold: float = 99.0  # ML 数据 identity 阈值

    # 内部串联重复检测（使用 TideHunter）
    # 启用后会过滤掉内部含有短串联重复的区域，避免 TideHunter 检测出错误的 consensus
    filter_internal_repeats: bool = True
    internal_repeat_copy_count: int = 10  # 模拟的串联重复次数
    internal_repeat_ratio_threshold: float = 0.8  # consensus/region_length < 此值则排除

    @property
    def mu(self) -> float:
        """计算 lognormal 的 mu 参数"""
        return np.log(self.mode) + self.sigma ** 2


@dataclass
class ChromInfo:
    """染色体信息"""
    name: str
    length: int


@dataclass
class EccRegion:
    """单个 eccDNA 区域"""
    chrom: str
    start: int
    end: int
    region_id: str
    length: int
    strand: str = '+'
    component: str = 'main'  # 'main', 'tail', 'chimeric'
    ecc_type: str = 'unknown'  # 'U', 'M', 'C'
    fragments: List[Tuple[str, int, int]] = field(default_factory=list)
    # MeccDNA 的多处匹配位置: [(chrom, start, end, identity, coverage), ...]
    multi_hits: List[Tuple[str, int, int, float, float]] = field(default_factory=list)

    def to_bed(self) -> str:
        """转换为 BED 格式"""
        if self.ecc_type == 'C' and self.fragments:
            frag_str = ';'.join([f"{c}:{s}-{e}" for c, s, e in self.fragments])
            return f"{self.chrom}\t{self.start}\t{self.end}\t{self.region_id}\t{self.length}\t{self.strand}\t{self.ecc_type}\t{frag_str}"
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.region_id}\t{self.length}\t{self.strand}\t{self.ecc_type}\t."

    def get_sequence(self, sequences: Dict[str, str]) -> str:
        """获取序列"""
        if self.ecc_type == 'C' and self.fragments:
            return ''.join([sequences[c][s:e] for c, s, e in self.fragments])
        return sequences[self.chrom][self.start:self.end]


# =============================================================================
# FASTA Reader
# =============================================================================

def read_fasta(filepath: str, exclude_patterns: Optional[List[str]] = None) -> Tuple[Dict[str, str], List[ChromInfo]]:
    """读取 FASTA 文件"""
    exclude_patterns = exclude_patterns or []
    sequences: Dict[str, str] = {}
    chrom_info: List[ChromInfo] = []

    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'

    current_chrom = None
    current_seq: List[str] = []

    logging.info(f"Reading reference FASTA: {filepath}")

    with open_func(filepath, mode) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_chrom is not None:
                    seq = ''.join(current_seq)
                    sequences[current_chrom] = seq
                    chrom_info.append(ChromInfo(current_chrom, len(seq)))

                chrom_name = line[1:].split()[0]
                if any(pat in chrom_name for pat in exclude_patterns):
                    current_chrom = None
                    current_seq = []
                    logging.debug(f"Excluding chromosome: {chrom_name}")
                else:
                    current_chrom = chrom_name
                    current_seq = []
            else:
                if current_chrom is not None:
                    current_seq.append(line)

    if current_chrom is not None:
        seq = ''.join(current_seq)
        sequences[current_chrom] = seq
        chrom_info.append(ChromInfo(current_chrom, len(seq)))

    logging.info(f"Loaded {len(sequences)} chromosomes, total length: {sum(c.length for c in chrom_info):,} bp")
    return sequences, chrom_info


# =============================================================================
# Minimap2 Alignment Manager
# =============================================================================

class Minimap2Manager:
    """Minimap2 比对和分类管理"""

    def __init__(self, config: SimulationConfig, reference_path: str):
        self.config = config
        self.reference_path = reference_path

        self.work_dir: Optional[str] = None
        self.index_dir: Optional[str] = None
        self.reference_for_align: Optional[str] = None
        self.index_short: Optional[str] = None
        self.index_long: Optional[str] = None

    def setup(self):
        """设置 minimap2 环境"""
        if not self._check_minimap2_installed():
            raise RuntimeError("minimap2 not found. Please install: conda install -c bioconda minimap2")

        output_dir = Path(self.config.output_prefix).resolve().parent
        output_base = Path(self.config.output_prefix).name

        # work_dir: query/paf 临时文件（放输出目录）
        self.work_dir = str(output_dir / f"{output_base}_minimap2_tmp")
        if os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)
        os.makedirs(self.work_dir)
        logging.info(f"Created working directory: {self.work_dir}")

        # index_dir 默认：reference 所在目录
        ref_dir = Path(self.reference_path).resolve().parent
        chosen_index_dir = Path(self.config.index_dir).resolve() if self.config.index_dir else ref_dir

        # 检查目录是否可写
        try:
            os.makedirs(chosen_index_dir, exist_ok=True)
            testfile = chosen_index_dir / ".mmi_write_test.tmp"
            with open(testfile, "w") as _:
                pass
            os.remove(testfile)
        except Exception as e:
            logging.warning(
                f"Index dir '{chosen_index_dir}' is not writable ({e}); fallback to output dir '{output_dir}'."
            )
            chosen_index_dir = output_dir
            os.makedirs(chosen_index_dir, exist_ok=True)

        self.index_dir = str(chosen_index_dir)
        logging.info(f"Index directory: {self.index_dir} (reuse_index={self.config.reuse_index})")

        # 若 reference 是 .gz，解压到 index_dir
        self.reference_for_align = self._materialize_reference(self.reference_path, self.index_dir)

        ref_tag = self._safe_ref_tag(self.reference_path)

        # 构建/复用索引
        self.index_short = os.path.join(self.index_dir, f"{ref_tag}.{self.config.minimap_preset_short}.mmi")
        self._ensure_index(self.index_short, self.config.minimap_preset_short)

        if self.config.split_by_length:
            self.index_long = os.path.join(self.index_dir, f"{ref_tag}.{self.config.minimap_preset_long}.mmi")
            self._ensure_index(self.index_long, self.config.minimap_preset_long)

    @staticmethod
    def _safe_ref_tag(ref_path: str) -> str:
        """从路径提取安全的标签名"""
        p = Path(ref_path)
        name = p.name
        if name.endswith('.gz'):
            name = name[:-3]
        for suf in ('.fa', '.fasta', '.fna'):
            if name.endswith(suf):
                name = name[: -len(suf)]
                break
        return name if name else "reference"

    def _check_minimap2_installed(self) -> bool:
        """检查 minimap2 是否安装"""
        try:
            subprocess.run(['minimap2', '--version'], capture_output=True, check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

    def _materialize_reference(self, ref_path: str, dest_dir: str) -> str:
        """如果是 gzip 压缩的，解压到目标目录"""
        if not ref_path.endswith('.gz'):
            return ref_path

        out_path = os.path.join(dest_dir, f"{self._safe_ref_tag(ref_path)}.decompressed.fa")

        try:
            gz_mtime = os.path.getmtime(ref_path)
            if os.path.exists(out_path):
                fa_mtime = os.path.getmtime(out_path)
                if fa_mtime >= gz_mtime and os.path.getsize(out_path) > 0:
                    logging.info(f"Reusing decompressed reference: {out_path}")
                    return out_path
        except OSError:
            pass

        logging.info(f"Reference is gzipped; decompressing to: {out_path}")
        with gzip.open(ref_path, 'rt') as fin, open(out_path, 'w') as fout:
            shutil.copyfileobj(fin, fout)
        return out_path

    def _ensure_index(self, index_path: str, preset: str):
        """确保索引存在，必要时构建"""
        assert self.reference_for_align is not None

        if self.config.reuse_index and (not self.config.force_rebuild_index) and os.path.exists(index_path):
            try:
                idx_mtime = os.path.getmtime(index_path)
                ref_mtime = os.path.getmtime(self.reference_for_align)
                if idx_mtime >= ref_mtime and os.path.getsize(index_path) > 0:
                    logging.info(f"Reusing existing index: {index_path}")
                    return
                else:
                    logging.info(f"Index exists but older than reference; rebuilding: {index_path}")
            except OSError:
                logging.info(f"Index exists but stat failed; rebuilding: {index_path}")

        self._build_index(index_path, preset)

    def _build_index(self, index_path: str, preset: str):
        """构建 minimap2 索引"""
        assert self.reference_for_align is not None
        logging.info(f"Building minimap2 index ({preset}) -> {index_path}")

        cmd_try = ['minimap2', '-x', preset, '-d', index_path, self.reference_for_align]
        result = subprocess.run(cmd_try, capture_output=True, text=True)
        if result.returncode == 0:
            logging.info(f"Index built successfully: {index_path}")
            return

        logging.warning(
            f"Index build with '-x {preset}' failed; fallback to plain '-d'. stderr: {result.stderr.strip()}"
        )
        cmd_fallback = ['minimap2', '-d', index_path, self.reference_for_align]
        result2 = subprocess.run(cmd_fallback, capture_output=True, text=True)
        if result2.returncode != 0:
            raise RuntimeError(f"Failed to build minimap2 index: {result2.stderr}")
        logging.info(f"Index built successfully (fallback): {index_path}")

    @staticmethod
    def _overlap_fraction(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
        """计算两个区间的重叠比例"""
        a_len = max(1, a_end - a_start)
        ov = max(0, min(a_end, b_end) - max(a_start, b_start))
        return ov / a_len

    def _run_minimap2(self, index_path: str, preset: str, query_fa: str, paf_out: str):
        """运行 minimap2"""
        cmd = [
            'minimap2',
            '-x', preset,
            '--secondary=yes',
            '-N', str(self.config.max_secondary),
            '-c',  # 输出 CIGAR，使 identity 计算更准确
            '--cs',  # 输出 cs tag，用于更准确的 identity 计算
            '-t', str(self.config.threads),
            '-o', paf_out,
            index_path,
            query_fa
        ]
        logging.info(f"Running minimap2 preset={preset} threads={self.config.threads} ...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"minimap2 failed: {result.stderr}")

    def classify_regions(self, regions: List[EccRegion], sequences: Dict[str, str]) -> Tuple[List[EccRegion], List[EccRegion]]:
        """通过 minimap2 分类 eccDNA 为 Unique 和 Multi-mapped

        修改说明 (v4.1):
        - 将序列双倍化后比对，模拟 rolling-circle amplification 产物
        - 使用 de:f 标签（per-base divergence）计算 identity，与 CircleSeeker 一致
        - 基于原始长度计算覆盖率，正确处理环状 DNA 特性
        """
        if not regions:
            return [], []

        assert self.work_dir is not None
        assert self.index_short is not None

        logging.info(f"Classifying {len(regions)} regions via minimap2 (doubled sequences)...")

        # 按长度分流写入不同文件，序列双倍化
        query_short = os.path.join(self.work_dir, "query.short.fa")
        query_long = os.path.join(self.work_dir, "query.long.fa")

        # 记录原始长度（用于后续覆盖率计算）
        original_lengths: Dict[str, int] = {}

        with open(query_short, 'w') as fa_s, open(query_long, 'w') as fa_l:
            for r in regions:
                seq = r.get_sequence(sequences)
                original_len = len(seq)
                original_lengths[r.region_id] = original_len

                # 双倍化序列，模拟环状 DNA 的 rolling-circle amplification
                doubled_seq = seq + seq

                header = f">{r.region_id}\n"
                body = f"{doubled_seq}\n"

                # 按原始长度分流（不是双倍化后的长度）
                if self.config.split_by_length and original_len >= self.config.split_length:
                    fa_l.write(header)
                    fa_l.write(body)
                else:
                    fa_s.write(header)
                    fa_s.write(body)

        paf_paths: List[str] = []

        # 短序列比对
        if os.path.exists(query_short) and os.path.getsize(query_short) > 0:
            paf_short = os.path.join(self.work_dir, "align.short.paf")
            self._run_minimap2(self.index_short, self.config.minimap_preset_short, query_short, paf_short)
            paf_paths.append(paf_short)

        # 长序列比对
        if self.config.split_by_length:
            if self.index_long is None:
                raise RuntimeError("split_by_length=True but long index is not built")
            if os.path.exists(query_long) and os.path.getsize(query_long) > 0:
                paf_long = os.path.join(self.work_dir, "align.long.paf")
                self._run_minimap2(self.index_long, self.config.minimap_preset_long, query_long, paf_long)
                paf_paths.append(paf_long)

        # 解析 PAF 结果
        # 存储结构: {qname: [(tname, tstart, tend, identity, coverage, qstart, qend), ...]}
        hits_by_q = defaultdict(list)

        for paf in paf_paths:
            with open(paf, 'r') as f:
                for line in f:
                    if not line.strip():
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 12:
                        continue

                    qname = parts[0]
                    qlen = int(parts[1])  # 这是双倍化后的长度
                    qstart = int(parts[2])
                    qend = int(parts[3])
                    tname = parts[5]
                    tstart = int(parts[7])
                    tend = int(parts[8])
                    nmatch = int(parts[9])
                    alen = int(parts[10])

                    if alen <= 0 or qlen <= 0:
                        continue

                    # 获取原始长度
                    original_len = original_lengths.get(qname, qlen // 2)

                    qspan = qend - qstart
                    tspan = abs(tend - tstart)
                    if qspan <= 0 or tspan <= 0:
                        continue

                    # 使用 de:f 标签计算 identity（与 CircleSeeker 一致）
                    identity = None
                    for tag in parts[12:]:
                        if tag.startswith("de:f:"):
                            try:
                                divergence = float(tag[5:])
                                identity = (1.0 - divergence) * 100.0
                            except ValueError:
                                pass
                            break

                    # 回退到 nmatch/alen 计算
                    if identity is None:
                        identity = (nmatch / alen) * 100.0

                    # 基于原始长度计算覆盖率和长度一致性
                    # 对于双倍化序列，target span 应该接近原始长度
                    tspan_delta_pct = (abs(tspan - original_len) / original_len) * 100.0
                    length_tolerance = max(0.0, 100.0 - self.config.length_consistency_threshold)

                    # 覆盖率：target span / original_len
                    # 如果 tspan >= original_len，说明完整覆盖了环状 DNA
                    coverage = min(100.0, (tspan / original_len) * 100.0)

                    # 长度一致性检查：tspan 应该接近 original_len
                    length_consistent = tspan_delta_pct <= length_tolerance

                    if length_consistent and identity >= self.config.identity_threshold and coverage >= self.config.min_coverage:
                        hits_by_q[qname].append((tname, tstart, tend, identity, coverage, qstart, qend))

        # 分类
        region_dict = {r.region_id: r for r in regions}
        unique_regions: List[EccRegion] = []
        multi_regions: List[EccRegion] = []
        skipped_no_hit = 0

        def dedupe_hits(hits: List[Tuple], bin_size: int = 100) -> List[Tuple]:
            """去重相似的 hits（基于基因组位置）"""
            seen = set()
            out = []
            for hit in hits:
                chrom, s, e = hit[0], hit[1], hit[2]
                key = (chrom, s // bin_size, e // bin_size)
                if key in seen:
                    continue
                seen.add(key)
                out.append(hit)
            return out

        for rid, r in region_dict.items():
            hits = hits_by_q.get(rid, [])
            if not hits:
                if self.config.no_hit_policy == "unique":
                    r.ecc_type = 'U'
                    unique_regions.append(r)
                else:
                    skipped_no_hit += 1
                continue

            hits = dedupe_hits(hits)

            # 过滤自身比对
            non_self_hits = []
            for hit in hits:
                chrom, s, e, ident, cov = hit[0], hit[1], hit[2], hit[3], hit[4]
                if chrom == r.chrom:
                    ovf = self._overlap_fraction(r.start, r.end, s, e)
                    if ovf >= self.config.self_overlap_fraction:
                        continue
                non_self_hits.append((chrom, s, e, ident, cov))

            full_coverage_hits = [
                hit for hit in non_self_hits if hit[4] >= self.config.multi_coverage_threshold
            ]
            if full_coverage_hits:
                r.ecc_type = 'M'
                # 记录满足 Multi 判定的 full-coverage hits（不包括自身）
                r.multi_hits = full_coverage_hits
                multi_regions.append(r)
            else:
                r.ecc_type = 'U'
                unique_regions.append(r)

        if skipped_no_hit > 0:
            logging.warning(
                f"{skipped_no_hit} regions had no hits meeting thresholds and were skipped (no_hit_policy=skip). "
                f"Use --no-hit-policy unique if you want to count them as Unique."
            )

        logging.info(f"Classification complete: {len(unique_regions)} Unique, {len(multi_regions)} Multi-mapped")
        return unique_regions, multi_regions

    def cleanup(self, keep: bool = False):
        """清理临时文件"""
        if self.work_dir and os.path.exists(self.work_dir):
            if keep:
                logging.info(f"Keeping working directory: {self.work_dir}")
            else:
                shutil.rmtree(self.work_dir)
                logging.info(f"Cleaned up working directory: {self.work_dir}")


# =============================================================================
# Length Distribution Sampler
# =============================================================================

class LengthSampler:
    """混合分布长度采样器"""

    def __init__(self, config: SimulationConfig, rng: np.random.Generator):
        self.config = config
        self.rng = rng
        self.mu = config.mu
        self.sigma = config.sigma

        # 统计计数
        self.main_count = 0
        self.tail_count = 0
        self.rejection_count = 0

    def sample_lognormal(self) -> int:
        """从 lognormal 分布采样"""
        return int(round(self.rng.lognormal(mean=self.mu, sigma=self.sigma)))

    def sample_log_uniform(self) -> int:
        """从 log-uniform 分布采样"""
        log_min = np.log10(self.config.tail_min)
        log_max = np.log10(self.config.tail_max)
        return int(round(10 ** self.rng.uniform(log_min, log_max)))

    def sample_one(self) -> Tuple[int, str]:
        """从混合分布采样一个长度"""
        if self.rng.random() < self.config.tail_weight:
            return self.sample_log_uniform(), 'tail'
        return self.sample_lognormal(), 'main'

    def sample_with_constraints(self, max_length: int) -> Tuple[Optional[int], str]:
        """采样一个满足约束的长度"""
        for _ in range(self.config.max_retry):
            length, component = self.sample_one()
            if self.config.min_length <= length <= min(max_length, self.config.max_length):
                if component == 'main':
                    self.main_count += 1
                else:
                    self.tail_count += 1
                return length, component
            self.rejection_count += 1
        return None, 'failed'

    def sample_uniform(self, min_len: int, max_len: int) -> int:
        """均匀采样（用于嵌合体片段）"""
        return int(self.rng.integers(min_len, max_len + 1))

    def get_stats(self) -> dict:
        """获取采样统计信息"""
        total = self.main_count + self.tail_count
        return {
            'total_sampled': total,
            'main_count': self.main_count,
            'tail_count': self.tail_count,
            'main_fraction': self.main_count / total if total > 0 else 0,
            'tail_fraction': self.tail_count / total if total > 0 else 0,
            'rejection_count': self.rejection_count,
        }


# =============================================================================
# Internal Repeat Checker (TideHunter)
# =============================================================================

class InternalRepeatChecker:
    """
    使用 TideHunter 检测区域内部是否含有短串联重复

    方法：将区域序列重复 N 次模拟 eccDNA read，用 TideHunter 检测
    如果 TideHunter 检测到的 consensus 长度 << 区域长度，说明区域内部有短串联重复
    这种区域在真实检测时会导致 TideHunter 检测出错误的 consensus，应该排除
    """

    def __init__(self, config: SimulationConfig):
        self.config = config
        self.work_dir: Optional[str] = None

    def __enter__(self):
        self.work_dir = tempfile.mkdtemp(prefix="th_check_")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.work_dir and os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)

    def filter_regions(
        self,
        regions: List['EccRegion'],
        sequences: Dict[str, str],
    ) -> Tuple[List['EccRegion'], int]:
        """
        过滤有内部串联重复的区域

        Returns:
            (clean_regions, num_filtered)
        """
        if not regions or not self.config.filter_internal_repeats:
            return regions, 0

        assert self.work_dir is not None

        # 生成模拟的串联重复序列
        fasta_path = os.path.join(self.work_dir, "tandem_reads.fa")
        region_lengths = {}

        with open(fasta_path, 'w') as f:
            for r in regions:
                seq = r.get_sequence(sequences)
                region_lengths[r.region_id] = len(seq)
                # 生成串联重复序列（模拟 eccDNA read）
                tandem_seq = seq * self.config.internal_repeat_copy_count
                f.write(f">{r.region_id}\n{tandem_seq}\n")

        # 运行 TideHunter
        th_output = self._run_tidehunter(fasta_path)

        # 解析结果，找出有问题的区域
        problematic_ids = self._parse_output(th_output, region_lengths)

        # 过滤
        clean = [r for r in regions if r.region_id not in problematic_ids]
        num_filtered = len(regions) - len(clean)

        return clean, num_filtered

    def _run_tidehunter(self, fasta_path: str) -> str:
        """运行 TideHunter"""
        output_path = os.path.join(self.work_dir, "th_output.tsv")

        cmd = [
            "TideHunter",
            "-f", "2",  # 输出格式：tabular
            "-c", "2",  # 最小 copy number
            "-o", output_path,
            fasta_path,
        ]

        try:
            subprocess.run(
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
            logging.warning("TideHunter timeout during internal repeat check")
            return ""
        except FileNotFoundError:
            logging.warning("TideHunter not found, skipping internal repeat filter")
            return ""

    def _parse_output(
        self,
        th_output: str,
        region_lengths: Dict[str, int],
    ) -> Set[str]:
        """解析 TideHunter 输出，找出有内部重复的区域"""
        problematic_ids: Set[str] = set()
        best_consensus: Dict[str, int] = {}

        for line in th_output.strip().split('\n'):
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 8:
                continue

            try:
                region_id = parts[0]
                # TideHunter -f 2 格式：最后一列是 consensus 序列
                consensus_len = len(parts[-1]) if parts[-1] else 0

                # 记录最长的 consensus
                if region_id not in best_consensus or consensus_len > best_consensus[region_id]:
                    best_consensus[region_id] = consensus_len
            except (ValueError, IndexError):
                continue

        # 检查每个区域
        threshold = self.config.internal_repeat_ratio_threshold
        for region_id, consensus_len in best_consensus.items():
            region_len = region_lengths.get(region_id, 0)
            if region_len > 0:
                ratio = consensus_len / region_len
                if ratio < threshold:
                    # consensus 远小于区域长度 → 有内部串联重复
                    problematic_ids.add(region_id)
                    logging.debug(
                        f"Region {region_id} has internal repeat: "
                        f"consensus={consensus_len}bp, region={region_len}bp, ratio={ratio:.2f}"
                    )

        return problematic_ids


# =============================================================================
# Region Generator
# =============================================================================

class RegionGenerator:
    """eccDNA 区域生成器"""

    def __init__(self, sequences: Dict[str, str], chrom_info: List[ChromInfo], config: SimulationConfig, rng: np.random.Generator):
        self.sequences = sequences
        self.chrom_info = chrom_info
        self.config = config
        self.rng = rng

        total_length = sum(c.length for c in chrom_info)
        self.chrom_weights = [c.length / total_length for c in chrom_info]
        self.length_sampler = LengthSampler(config, rng)

    def choose_chromosome(self) -> ChromInfo:
        """按长度加权随机选择染色体"""
        idx = int(self.rng.choice(len(self.chrom_info), p=self.chrom_weights))
        return self.chrom_info[idx]

    def generate_one(self, region_id: str) -> Optional[EccRegion]:
        """生成单个候选 eccDNA"""
        for _ in range(self.config.max_retry):
            chrom = self.choose_chromosome()
            length, component = self.length_sampler.sample_with_constraints(chrom.length)
            if length is None:
                continue
            max_start = chrom.length - length
            if max_start < 0:
                continue
            start = int(self.rng.integers(0, max_start + 1))
            end = start + length
            return EccRegion(
                chrom=chrom.name,
                start=start,
                end=end,
                region_id=region_id,
                length=length,
                component=component
            )
        return None

    def generate_candidates(self, num_candidates: int, id_prefix: str = "eccDNA") -> List[EccRegion]:
        """生成候选 eccDNA 区域"""
        regions: List[EccRegion] = []
        failed = 0
        logging.info(f"Generating {num_candidates} candidate eccDNA regions...")
        for i in range(num_candidates):
            rid = f"{id_prefix}_{i+1:06d}"
            r = self.generate_one(rid)
            if r is not None:
                regions.append(r)
            else:
                failed += 1
            if (i + 1) % 10000 == 0:
                logging.info(f"  Progress: {i+1}/{num_candidates} ({100*(i+1)/num_candidates:.1f}%)")
        logging.info(f"Generated {len(regions)} candidates, {failed} failed")
        return regions

    def generate_fragment_candidates(
        self,
        num_candidates: int,
        min_length: int,
        max_length: int,
        id_prefix: str = "frag",
    ) -> List[EccRegion]:
        """生成嵌合体候选片段"""
        if min_length <= 0 or max_length <= 0 or min_length > max_length:
            raise ValueError("Invalid fragment length bounds.")

        regions: List[EccRegion] = []
        failed = 0
        for i in range(num_candidates):
            for _ in range(self.config.max_retry):
                chrom = self.choose_chromosome()
                frag_len = self.length_sampler.sample_uniform(min_length, max_length)
                frag_len = min(frag_len, chrom.length - 100)
                if frag_len < min_length:
                    continue
                max_start = chrom.length - frag_len
                if max_start < 0:
                    continue
                start = int(self.rng.integers(0, max_start + 1))
                end = start + frag_len
                rid = f"{id_prefix}_{i+1:06d}"
                regions.append(
                    EccRegion(
                        chrom=chrom.name,
                        start=start,
                        end=end,
                        region_id=rid,
                        length=frag_len,
                        component="chimeric",
                    )
                )
                break
            else:
                failed += 1
        if failed:
            logging.info(f"Fragment candidate generation: {failed} failed out of {num_candidates}")
        return regions

    # Minimum distance between IntraChr fragments to ensure they are distinct loci
    MIN_INTRACHR_FRAGMENT_DISTANCE = 1000

    def _check_fragment_distance(
        self,
        new_start: int,
        new_end: int,
        existing_fragments: List[Tuple[str, int, int]],
        same_chrom: str,
        min_distance: int = 1000,
    ) -> bool:
        """Check if new fragment has sufficient distance from existing fragments on same chromosome.

        For IntraChr CeccDNA, fragments must be far enough apart to be recognized
        as distinct genomic loci (not merged during detection).

        Also rejects fragments that significantly overlap with existing ones (>50% overlap).

        Returns:
            True if distance is sufficient and no significant overlap, False otherwise
        """
        for frag_chrom, frag_start, frag_end in existing_fragments:
            if frag_chrom != same_chrom:
                continue

            # Check for significant overlap (reject if >50% overlap)
            overlap_start = max(new_start, frag_start)
            overlap_end = min(new_end, frag_end)
            if overlap_start < overlap_end:
                overlap_len = overlap_end - overlap_start
                new_len = new_end - new_start
                frag_len = frag_end - frag_start
                overlap_ratio = overlap_len / min(new_len, frag_len)
                if overlap_ratio > 0.5:
                    return False

            # Calculate distance between regions
            if new_start < frag_start:
                distance = frag_start - new_end
            else:
                distance = new_start - frag_end
            if distance < min_distance:
                return False
        return True

    def generate_chimeric(self, num_chimeric: int, id_prefix: str = "CeccDNA") -> List[EccRegion]:
        """生成嵌合体 eccDNA

        For IntraChr chimeric eccDNA (all fragments on same chromosome),
        ensures fragments are separated by at least MIN_INTRACHR_FRAGMENT_DISTANCE
        to be recognized as distinct genomic loci during detection.
        """
        if num_chimeric == 0:
            return []

        logging.info(f"Generating {num_chimeric} chimeric eccDNA regions...")
        chimeric_regions: List[EccRegion] = []
        seen_fragments = set()
        duplicate_chimeric = 0

        weights = self.config.chimeric_fragments_weights
        options = list(weights.keys())
        probs = np.array([weights[k] for k in options], dtype=float)
        probs = probs / probs.sum()

        max_attempts_per_fragment = 50

        for i in range(num_chimeric):
            num_frags = int(self.rng.choice(options, p=probs))
            fragments: List[Tuple[str, int, int]] = []
            total_len = 0

            for frag_idx in range(num_frags):
                chrom = self.choose_chromosome()
                frag_len = self.length_sampler.sample_uniform(100, 5000)
                frag_len = min(frag_len, chrom.length - 100)
                if frag_len < 100:
                    continue

                # Try to find a valid position with sufficient distance from existing fragments
                found_valid = False
                for attempt in range(max_attempts_per_fragment):
                    max_start = chrom.length - frag_len
                    start = int(self.rng.integers(0, max_start + 1))
                    end = start + frag_len

                    # Check distance for IntraChr case
                    if self._check_fragment_distance(
                        start, end, fragments, chrom.name,
                        min_distance=self.MIN_INTRACHR_FRAGMENT_DISTANCE
                    ):
                        fragments.append((chrom.name, start, end))
                        total_len += frag_len
                        found_valid = True
                        break

                # Fallback: try a different chromosome to avoid overlap/duplicate issues
                if not found_valid:
                    # Try other chromosomes to find a valid fragment
                    fallback_found = False
                    for _ in range(10):  # Try up to 10 different chromosomes
                        alt_chrom = self.choose_chromosome()
                        alt_frag_len = self.length_sampler.sample_uniform(100, 5000)
                        alt_frag_len = min(alt_frag_len, alt_chrom.length - 100)
                        if alt_frag_len < 100:
                            continue
                        max_start = alt_chrom.length - alt_frag_len
                        start = int(self.rng.integers(0, max_start + 1))
                        end = start + alt_frag_len
                        # Check this fragment is valid (different chrom = always valid for distance)
                        if self._check_fragment_distance(start, end, fragments, alt_chrom.name, min_distance=0):
                            fragments.append((alt_chrom.name, start, end))
                            total_len += alt_frag_len
                            fallback_found = True
                            break
                    # If still not found, skip this fragment entirely (better than duplicate)
                    if not fallback_found:
                        logging.debug(f"Could not find valid fragment for chimeric {i+1}, skipping fragment {frag_idx+1}")

            if len(fragments) < 2:
                continue

            frag_key = tuple(fragments)
            if frag_key in seen_fragments:
                duplicate_chimeric += 1
                continue
            seen_fragments.add(frag_key)

            region_id = f"{id_prefix}_{i+1:06d}"
            first = fragments[0]
            chimeric_regions.append(
                EccRegion(
                    chrom=first[0],
                    start=first[1],
                    end=first[2],
                    region_id=region_id,
                    length=total_len,
                    component='chimeric',
                    ecc_type='C',
                    fragments=fragments
                )
            )

        frag_counts = defaultdict(int)
        for r in chimeric_regions:
            frag_counts[len(r.fragments)] += 1
        msg = f"Generated {len(chimeric_regions)} chimeric eccDNA; fragment distribution: {dict(frag_counts)}"
        if duplicate_chimeric:
            msg += f" (deduplicated {duplicate_chimeric})"
        logging.info(msg)
        return chimeric_regions


# =============================================================================
# Main Pipeline
# =============================================================================

class EccDNASimulator:
    """eccDNA 模拟主流程"""

    def __init__(self, config: SimulationConfig):
        self.config = config
        self.rng = np.random.default_rng(config.seed)
        self.sequences: Optional[Dict[str, str]] = None
        self.chrom_info: Optional[List[ChromInfo]] = None
        self.aligner: Optional[Minimap2Manager] = None
        self.generator: Optional[RegionGenerator] = None

    def run(self) -> List[EccRegion]:
        """运行完整模拟流程"""
        # 设置输出目录
        self._setup_output_dir()
        try:
            logging.info("=" * 60)
            logging.info("eccDNA Region Simulator v4 (minimap2) - Final")
            logging.info("=" * 60)

            self._log_config()

            # 1. 读取参考基因组
            self.sequences, self.chrom_info = read_fasta(self.config.reference, self.config.exclude_chroms)
            if not self.sequences:
                raise RuntimeError("No chromosomes loaded from reference!")

            # 2. 初始化生成器
            self.generator = RegionGenerator(self.sequences, self.chrom_info, self.config, self.rng)

            need_um = (self.config.num_unique > 0 or self.config.num_multi > 0)
            need_chimeric_alignment = (
                self.config.num_chimeric > 0 and self.config.enforce_unique_chimeric_fragments
            )
            need_alignment = need_um or need_chimeric_alignment

            unique_regions: List[EccRegion] = []
            multi_regions: List[EccRegion] = []
            chimeric_regions: List[EccRegion] = []

            # 3. 生成并分类 U/M + 嵌合体片段筛选
            if need_alignment:
                self.aligner = Minimap2Manager(self.config, self.config.reference)
                self.aligner.setup()
                try:
                    if need_um:
                        unique_regions, multi_regions = self._generate_and_classify()
                    if self.config.num_chimeric > 0:
                        if self.config.enforce_unique_chimeric_fragments:
                            chimeric_regions = self._generate_chimeric_unique(self.config.num_chimeric)
                        else:
                            chimeric_regions = self.generator.generate_chimeric(self.config.num_chimeric)
                finally:
                    self.aligner.cleanup(keep=self.config.keep_tmp)
            else:
                if self.config.num_chimeric > 0:
                    chimeric_regions = self.generator.generate_chimeric(self.config.num_chimeric)

            # 5. 输出前统一去重并整理输出
            unique_regions, multi_regions, chimeric_regions = self._dedupe_final_regions(
                unique_regions,
                multi_regions,
                chimeric_regions,
            )
            all_regions = self._finalize_regions(unique_regions, multi_regions, chimeric_regions)
            self._write_outputs(all_regions)

            # 6. 生成 ML 训练数据（可选）
            if self.config.generate_ml_data:
                self._generate_ml_training_data(all_regions)

            self._print_summary(all_regions)
            return all_regions
        except Exception:
            self._write_error_log()
            raise

    def _setup_output_dir(self):
        """设置输出目录"""
        prefix = self.config.output_prefix
        prefix_path = Path(prefix)
        
        if self.config.output_dir:
            # 用户指定了输出目录
            out_dir = Path(self.config.output_dir)
            # prefix 只取文件名部分
            prefix_name = prefix_path.name
        else:
            # 未指定，用 prefix 创建目录
            if prefix_path.is_absolute():
                # 绝对路径：/path/to/sim_ecc → /path/to/sim_ecc/sim_ecc.*
                out_dir = prefix_path
                prefix_name = prefix_path.name
            else:
                # 相对路径：out/sim_ecc → ./out/sim_ecc/sim_ecc.*
                #          sim_ecc → ./sim_ecc/sim_ecc.*
                out_dir = Path.cwd() / prefix
                prefix_name = prefix_path.name  # 只取最后的文件名部分，不是整个路径
        
        # 创建目录
        out_dir.mkdir(parents=True, exist_ok=True)
        
        # 更新 output_prefix 为完整路径
        self.config.output_prefix = str(out_dir / prefix_name)
        self.config.output_dir = str(out_dir)
        
        logging.info(f"Output directory: {out_dir}")
        logging.info(f"Output prefix: {self.config.output_prefix}")

    def _log_config(self):
        """记录配置信息"""
        logging.info(f"Random seed: {self.config.seed}")
        logging.info(f"Target: {self.config.num_unique} Unique, {self.config.num_multi} Multi, {self.config.num_chimeric} Chimeric")
        logging.info(f"Lognormal mode: {self.config.mode} bp, sigma: {self.config.sigma}, mu: {self.config.mu:.4f}")
        logging.info(f"Tail weight: {self.config.tail_weight:.1%}, range: {self.config.tail_min:,}-{self.config.tail_max:,} bp")
        logging.info(f"Length constraints: {self.config.min_length}-{self.config.max_length:,} bp")
        logging.info(f"MeccDNA max length: {self.config.max_length_multi:,} bp")
        length_tol = max(0.0, 100.0 - self.config.length_consistency_threshold)
        logging.info(
            f"Threads: {self.config.threads} | Identity >= {self.config.identity_threshold}% | "
            f"Coverage >= {self.config.min_coverage}% | Length consistency >= {self.config.length_consistency_threshold}% "
            f"(target span within +/-{length_tol:.1f}%)"
        )
        logging.info(
            f"Multi classification: non-self hits must have coverage >= {self.config.multi_coverage_threshold}%"
        )
        logging.info(f"Split by length: {self.config.split_by_length} (>= {self.config.split_length} bp -> {self.config.minimap_preset_long})")
        logging.info(f"Index dir: {self.config.index_dir or str(Path(self.config.reference).resolve().parent)}")
        logging.info(f"Reuse index: {self.config.reuse_index} | Force rebuild: {self.config.force_rebuild_index}")
        logging.info(f"Candidate multiplier: U={self.config.candidate_multiplier_u}x, M={self.config.candidate_multiplier_m}x")
        if self.config.num_chimeric > 0:
            logging.info(
                "Chimeric fragments: unique=%s, unique_chrom_fraction>=%.2f, same_chrom_prob=%.2f, repeat_prob=%.2f",
                self.config.enforce_unique_chimeric_fragments,
                self.config.chimeric_unique_chrom_fraction,
                self.config.chimeric_same_chrom_probability,
                self.config.chimeric_fragment_repeat_probability,
            )

    def _generate_and_classify(self) -> Tuple[List[EccRegion], List[EccRegion]]:
        """生成候选并通过 minimap2 分类，最多 3 轮以确保达到目标数量"""
        assert self.generator is not None
        assert self.aligner is not None
        assert self.sequences is not None

        target_u = self.config.num_unique
        target_m = self.config.num_multi
        max_rounds = 10  # 最多 10 轮

        all_unique: List[EccRegion] = []
        all_multi: List[EccRegion] = []
        seen_coords: set = set()
        deduped_total = 0

        def dedupe_regions(regions: List[EccRegion]) -> Tuple[List[EccRegion], int]:
            kept: List[EccRegion] = []
            dupes = 0
            for r in regions:
                key = (r.chrom, r.start, r.end)
                if key in seen_coords:
                    dupes += 1
                    continue
                seen_coords.add(key)
                kept.append(r)
            return kept, dupes

        for round_num in range(1, max_rounds + 1):
            remaining_u = max(0, target_u - len(all_unique))
            remaining_m = max(0, target_m - len(all_multi))

            if remaining_u == 0 and remaining_m == 0:
                logging.info(f"Round {round_num}: Target already reached (U={len(all_unique)}, M={len(all_multi)})")
                break

            logging.info(f"{'='*60}")
            logging.info(f"Round {round_num}/{max_rounds}: Need {remaining_u} more U, {remaining_m} more M")
            logging.info(f"{'='*60}")

            # 计算本轮需要生成的候选数量
            candidates_for_u = int(remaining_u * self.config.candidate_multiplier_u) if remaining_u > 0 else 0
            candidates_for_m = int(remaining_m * self.config.candidate_multiplier_m) if remaining_m > 0 else 0
            total_candidates = candidates_for_u + candidates_for_m

            if total_candidates == 0:
                break

            logging.info(
                f"Candidate strategy: {candidates_for_u} for U (x{self.config.candidate_multiplier_u}), "
                f"{candidates_for_m} for M (x{self.config.candidate_multiplier_m})"
            )

            batch_size = 20000
            round_unique: List[EccRegion] = []
            round_multi: List[EccRegion] = []
            candidates_generated = 0
            batch_num = 0
            # 使用固定 multiplier: U=2x, M=10x
            max_total_candidates = max(10000, total_candidates * 3)

            # 批次生成循环
            while candidates_generated < max_total_candidates:
                if len(round_unique) >= remaining_u and len(round_multi) >= remaining_m:
                    logging.info("Round target reached for both U and M")
                    break

                batch_num += 1
                still_need_u = max(0, remaining_u - len(round_unique))
                still_need_m = max(0, remaining_m - len(round_multi))

                # 使用固定 multiplier
                needed_for_u = int(still_need_u * self.config.candidate_multiplier_u) if still_need_u > 0 else 0
                needed_for_m = int(still_need_m * self.config.candidate_multiplier_m) if still_need_m > 0 else 0
                current_batch = min(batch_size, needed_for_u + needed_for_m)
                current_batch = max(current_batch, 1000)

                logging.info(f"  Batch {batch_num}: generating {current_batch} candidates (need {still_need_u} U, {still_need_m} M)")

                candidates = self.generator.generate_candidates(current_batch, id_prefix=f"r{round_num}_b{batch_num}")
                candidates_generated += len(candidates)
                if not candidates:
                    logging.warning("  No candidates generated in this batch")
                    continue

                unique, multi = self.aligner.classify_regions(candidates, self.sequences)
                cap = int(self.config.max_length_multi)
                before = len(multi)
                if before:
                    multi = [r for r in multi if r.length <= cap]
                filtered = before - len(multi)
                if filtered:
                    logging.info(f"  Batch {batch_num}: filtered {filtered} multi candidates > {cap} bp")
                unique, dup_u = dedupe_regions(unique)
                multi, dup_m = dedupe_regions(multi)
                if dup_u or dup_m:
                    logging.info(f"  Batch {batch_num}: deduplicated {dup_u + dup_m} regions")
                    deduped_total += dup_u + dup_m

                # 只保留需要的数量，避免浪费
                if still_need_u > 0:
                    unique_to_add = unique[:still_need_u]
                    round_unique.extend(unique_to_add)
                if still_need_m > 0:
                    multi_to_add = multi[:still_need_m]
                    round_multi.extend(multi_to_add)

                u_rate = (len(unique) / len(candidates) * 100) if candidates else 0
                m_rate = (len(multi) / len(candidates) * 100) if candidates else 0
                logging.info(f"  Batch result: {len(unique)} U ({u_rate:.1f}%), {len(multi)} M ({m_rate:.1f}%)")
                logging.info(f"  Round cumulative: {len(round_unique)} U, {len(round_multi)} M")

            # 本轮 TideHunter 过滤
            if self.config.filter_internal_repeats and (round_unique or round_multi):
                logging.info(f"Round {round_num}: Filtering internal tandem repeats (TideHunter)...")
                with InternalRepeatChecker(self.config) as checker:
                    if round_unique:
                        round_unique, filtered_u = checker.filter_regions(round_unique, self.sequences)
                        if filtered_u:
                            logging.info(f"  Filtered {filtered_u} Unique regions with internal repeats")
                    if round_multi:
                        round_multi, filtered_m = checker.filter_regions(round_multi, self.sequences)
                        if filtered_m:
                            logging.info(f"  Filtered {filtered_m} Multi regions with internal repeats")

            # 添加本轮结果
            all_unique.extend(round_unique)
            all_multi.extend(round_multi)

            logging.info(f"Round {round_num} complete: Total {len(all_unique)} U, {len(all_multi)} M")

            # 检查是否达到目标
            if len(all_unique) >= target_u and len(all_multi) >= target_m:
                logging.info(f"Target reached after round {round_num}")
                break

        # 截取最终结果
        unique_final = all_unique[:target_u]
        multi_final = all_multi[:target_m]

        if len(unique_final) < target_u:
            logging.warning(
                f"Only got {len(unique_final)} Unique regions after {max_rounds} rounds (target: {target_u}). "
                f"Try higher --multiplier-u or check reference genome."
            )
        if target_m > 0 and len(multi_final) < target_m:
            logging.warning(
                f"Only got {len(multi_final)} Multi regions after {max_rounds} rounds (target: {target_m}). "
                f"Try higher --multiplier-m / increase --max-secondary / slightly lower --identity."
            )
        if deduped_total:
            logging.info(f"Deduplicated {deduped_total} regions during candidate selection")
        return unique_final, multi_final

    def _estimate_chimeric_pool_target(self, num_chimeric: int) -> int:
        weights = self.config.chimeric_fragments_weights
        if not weights:
            expected = 2.0
        else:
            total_weight = float(sum(weights.values()))
            if total_weight <= 0:
                expected = 2.0
            else:
                expected = sum(k * w for k, w in weights.items()) / total_weight
        target = int(math.ceil(num_chimeric * expected * self.config.chimeric_fragment_pool_multiplier))
        return max(target, num_chimeric * 2)

    def _build_unique_fragment_pool(self, target_size: int) -> List[EccRegion]:
        assert self.generator is not None
        assert self.aligner is not None
        assert self.sequences is not None

        if target_size <= 0:
            return []

        min_len = int(self.config.chimeric_fragment_min_length)
        max_len = int(self.config.chimeric_fragment_max_length)
        batch_size = 20000
        pool: List[EccRegion] = []
        seen_coords = set()
        generated = 0
        batch_num = 0
        max_total_candidates = max(1000, target_size * 5)

        logging.info("Building unique fragment pool (target=%d)...", target_size)
        while len(pool) < target_size and generated < max_total_candidates:
            batch_num += 1
            remaining = target_size - len(pool)
            current_batch = min(batch_size, max(remaining, 1000))
            candidates = self.generator.generate_fragment_candidates(
                current_batch,
                min_length=min_len,
                max_length=max_len,
                id_prefix=f"cfrag_b{batch_num}",
            )
            generated += len(candidates)
            if not candidates:
                break
            unique, _ = self.aligner.classify_regions(candidates, self.sequences)
            for r in unique:
                key = (r.chrom, r.start, r.end)
                if key in seen_coords:
                    continue
                seen_coords.add(key)
                pool.append(r)
                if len(pool) >= target_size:
                    break
            logging.info("  pool size: %d (batch %d, uniques %d)", len(pool), batch_num, len(unique))

        if len(pool) < target_size:
            logging.warning(
                "Unique fragment pool undersized: %d / %d. Chimeric output may be reduced.",
                len(pool),
                target_size,
            )
        return pool

    def _pick_additional_fragment_indices(
        self,
        pool_by_chrom: dict,
        pool_size: int,
        used_indices: set,
        selected_chroms: set,
        remaining: int,
    ) -> Optional[List[int]]:
        extras: List[int] = []
        for _ in range(remaining):
            use_same = bool(selected_chroms) and (
                self.rng.random() < self.config.chimeric_same_chrom_probability
            )
            if use_same:
                chrom_candidates = list(selected_chroms)
            else:
                chrom_candidates = [c for c in pool_by_chrom if c not in selected_chroms]
                if not chrom_candidates:
                    chrom_candidates = list(selected_chroms)
            if not chrom_candidates:
                return None
            chrom = chrom_candidates[int(self.rng.integers(0, len(chrom_candidates)))]
            indices = [i for i in pool_by_chrom[chrom] if i not in used_indices]
            if not indices:
                indices = [i for i in range(pool_size) if i not in used_indices]
                if not indices:
                    return None
            idx = indices[int(self.rng.integers(0, len(indices)))]
            used_indices.add(idx)
            extras.append(idx)
        return extras

    def _select_fragments_from_pool(
        self,
        pool: List[EccRegion],
        pool_by_chrom: dict,
        distinct_needed: int,
    ) -> Optional[List[EccRegion]]:
        if distinct_needed <= 0 or len(pool) < distinct_needed:
            return None

        chroms = list(pool_by_chrom.keys())
        if not chroms:
            return None

        min_unique = int(math.ceil(distinct_needed * self.config.chimeric_unique_chrom_fraction))
        min_unique = max(1, min(min_unique, distinct_needed, len(chroms)))
        if min_unique > len(chroms):
            min_unique = len(chroms)

        max_attempts = max(50, self.config.max_retry)
        pool_size = len(pool)
        for _ in range(max_attempts):
            used_indices = set()
            selected: List[EccRegion] = []
            selected_chroms = set()

            chosen_chroms = list(self.rng.choice(chroms, size=min_unique, replace=False))
            for chrom in chosen_chroms:
                indices = [i for i in pool_by_chrom[chrom] if i not in used_indices]
                if not indices:
                    break
                idx = indices[int(self.rng.integers(0, len(indices)))]
                used_indices.add(idx)
                selected.append(pool[idx])
                selected_chroms.add(chrom)

            if len(selected_chroms) < min_unique:
                continue

            remaining = distinct_needed - len(selected)
            if remaining > 0:
                extras = self._pick_additional_fragment_indices(
                    pool_by_chrom,
                    pool_size,
                    used_indices,
                    selected_chroms,
                    remaining,
                )
                if extras is None:
                    continue
                selected.extend(pool[i] for i in extras)

            return selected
        return None

    def _generate_chimeric_unique(self, num_chimeric: int) -> List[EccRegion]:
        """生成 Chimeric eccDNA，最多 10 轮以确保达到目标数量"""
        assert self.generator is not None

        if num_chimeric == 0:
            return []

        max_rounds = 10  # 最多 10 轮
        chimeric_regions: List[EccRegion] = []
        seen_fragments: set = set()
        total_duplicate = 0
        total_failed = 0

        weights = self.config.chimeric_fragments_weights
        options = list(weights.keys())
        probs = np.array([weights[k] for k in options], dtype=float)
        probs = probs / probs.sum()

        for round_num in range(1, max_rounds + 1):
            remaining = num_chimeric - len(chimeric_regions)
            if remaining <= 0:
                break

            logging.info(f"{'='*60}")
            logging.info(f"Chimeric Round {round_num}/{max_rounds}: Need {remaining} more C")
            logging.info(f"{'='*60}")

            # 估算并构建 fragment pool
            target_pool = self._estimate_chimeric_pool_target(remaining)
            pool = self._build_unique_fragment_pool(target_pool)
            if not pool:
                logging.warning(f"Round {round_num}: No unique fragments available for chimeric generation.")
                continue

            pool_by_chrom = defaultdict(list)
            for idx, frag in enumerate(pool):
                pool_by_chrom[frag.chrom].append(idx)

            round_chimeric: List[EccRegion] = []
            duplicate_chimeric = 0
            failed = 0

            for i in range(remaining):
                num_frags = int(self.rng.choice(options, p=probs))
                repeat_count = 0
                if (
                    self.config.chimeric_fragment_max_repeats > 0
                    and num_frags > 1
                    and self.rng.random() < self.config.chimeric_fragment_repeat_probability
                ):
                    repeat_count = min(self.config.chimeric_fragment_max_repeats, num_frags - 1)

                distinct_needed = num_frags - repeat_count
                selected = self._select_fragments_from_pool(pool, pool_by_chrom, distinct_needed)
                if not selected:
                    failed += 1
                    continue

                if repeat_count > 0:
                    for _ in range(repeat_count):
                        rep_idx = int(self.rng.integers(0, len(selected)))
                        selected.append(selected[rep_idx])

                fragments = [(f.chrom, f.start, f.end) for f in selected]
                self.rng.shuffle(fragments)

                frag_key = tuple(fragments)
                if frag_key in seen_fragments:
                    duplicate_chimeric += 1
                    continue
                seen_fragments.add(frag_key)

                total_len = sum(end - start for _, start, end in fragments)
                region_id = f"CeccDNA_{len(chimeric_regions) + len(round_chimeric) + 1:06d}"
                first = fragments[0]
                round_chimeric.append(
                    EccRegion(
                        chrom=first[0],
                        start=first[1],
                        end=first[2],
                        region_id=region_id,
                        length=total_len,
                        component="chimeric",
                        ecc_type="C",
                        fragments=fragments,
                    )
                )

            total_duplicate += duplicate_chimeric
            total_failed += failed
            chimeric_regions.extend(round_chimeric)

            frag_counts = defaultdict(int)
            for r in round_chimeric:
                frag_counts[len(r.fragments)] += 1
            msg = f"Round {round_num}: Generated {len(round_chimeric)} chimeric; distribution: {dict(frag_counts)}"
            if duplicate_chimeric:
                msg += f" (deduplicated {duplicate_chimeric})"
            if failed:
                msg += f" (failed {failed})"
            logging.info(msg)
            logging.info(f"Total chimeric: {len(chimeric_regions)}/{num_chimeric}")

            if len(chimeric_regions) >= num_chimeric:
                logging.info(f"Chimeric target reached after round {round_num}")
                break

        # 最终统计
        if len(chimeric_regions) < num_chimeric:
            logging.warning(
                f"Only got {len(chimeric_regions)} Chimeric regions after {max_rounds} rounds (target: {num_chimeric}). "
                f"Try increasing fragment pool size or relaxing chimeric constraints."
            )

        frag_counts = defaultdict(int)
        for r in chimeric_regions:
            frag_counts[len(r.fragments)] += 1
        msg = f"Final: {len(chimeric_regions)} chimeric eccDNA; fragment distribution: {dict(frag_counts)}"
        if total_duplicate:
            msg += f" (total deduplicated {total_duplicate})"
        if total_failed:
            msg += f" (total failed {total_failed})"
        logging.info(msg)
        return chimeric_regions[:num_chimeric]

    def _dedupe_final_regions(
        self,
        unique: List[EccRegion],
        multi: List[EccRegion],
        chimeric: List[EccRegion],
    ) -> Tuple[List[EccRegion], List[EccRegion], List[EccRegion]]:
        """输出前统一去重"""
        seen = set()
        removed_u = removed_m = removed_c = 0

        def region_key(r: EccRegion):
            if r.ecc_type == 'C' and r.fragments:
                return ("C", tuple(r.fragments))
            return ("L", r.chrom, r.start, r.end)

        def dedupe_list(regions: List[EccRegion]):
            kept = []
            removed = 0
            for r in regions:
                key = region_key(r)
                if key in seen:
                    removed += 1
                    continue
                seen.add(key)
                kept.append(r)
            return kept, removed

        unique, removed_u = dedupe_list(unique)
        multi, removed_m = dedupe_list(multi)
        chimeric, removed_c = dedupe_list(chimeric)

        total_removed = removed_u + removed_m + removed_c
        if total_removed:
            logging.info(
                "Final deduplication removed "
                f"{total_removed} regions (U:{removed_u}, M:{removed_m}, C:{removed_c})"
            )

        return unique, multi, chimeric

    def _finalize_regions(self, unique: List[EccRegion], multi: List[EccRegion], chimeric: List[EccRegion]) -> List[EccRegion]:
        """整理并重新编号所有区域"""
        all_regions: List[EccRegion] = []
        for i, r in enumerate(unique):
            r.region_id = f"UeccDNA_{i+1:06d}"
            all_regions.append(r)
        for i, r in enumerate(multi):
            r.region_id = f"MeccDNA_{i+1:06d}"
            all_regions.append(r)
        for i, r in enumerate(chimeric):
            r.region_id = f"CeccDNA_{i+1:06d}"
            all_regions.append(r)
        return all_regions

    def _write_outputs(self, regions: List[EccRegion]):
        """写入输出文件"""
        u_regions = [r for r in regions if r.ecc_type == 'U']
        m_regions = [r for r in regions if r.ecc_type == 'M']
        c_regions = [r for r in regions if r.ecc_type == 'C']

        prefix = self.config.output_prefix

        # BED 文件
        self._write_bed(regions, f"{prefix}.all.bed")
        if u_regions:
            self._write_bed(u_regions, f"{prefix}.unique.bed")
        if m_regions:
            self._write_bed(m_regions, f"{prefix}.multi.bed")
        if c_regions:
            self._write_bed(c_regions, f"{prefix}.chimeric.bed")

        # FASTA 文件
        self._write_fasta(regions, f"{prefix}.all.fa")
        if u_regions:
            self._write_fasta(u_regions, f"{prefix}.unique.fa")
        if m_regions:
            self._write_fasta(m_regions, f"{prefix}.multi.fa")
        if c_regions:
            self._write_fasta(c_regions, f"{prefix}.chimeric.fa")

        # QC 日志和长度数据
        self._write_qc_log(regions)
        self._write_length_data(regions)
        
        # MeccDNA 多处匹配位置
        if m_regions:
            self._write_multi_hits(m_regions)

    def _write_bed(self, regions: List[EccRegion], output_path: str):
        """写入 BED 文件"""
        logging.info(f"Writing BED file: {output_path}")
        with open(output_path, 'w') as f:
            f.write("#chrom\tstart\tend\tname\tlength\tstrand\ttype\tfragments\n")
            for r in regions:
                f.write(r.to_bed() + "\n")

    def _write_fasta(self, regions: List[EccRegion], output_path: str):
        """写入 FASTA 文件"""
        assert self.sequences is not None
        logging.info(f"Writing FASTA file: {output_path}")
        with open(output_path, 'w') as f:
            for r in regions:
                seq = r.get_sequence(self.sequences)
                if r.ecc_type == 'C' and r.fragments:
                    frag_str = ';'.join([f"{c}:{s}-{e}" for c, s, e in r.fragments])
                    header = f">{r.region_id} {frag_str} length={r.length} type={r.ecc_type}"
                else:
                    header = f">{r.region_id} {r.chrom}:{r.start}-{r.end} length={r.length} type={r.ecc_type}"
                f.write(header + "\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")

    def _write_qc_log(self, regions: List[EccRegion]):
        """写入详细的 QC 日志"""
        output_path = f"{self.config.output_prefix}.qc.log"
        logging.info(f"Writing QC log: {output_path}")

        u_regions = [r for r in regions if r.ecc_type == 'U']
        m_regions = [r for r in regions if r.ecc_type == 'M']
        c_regions = [r for r in regions if r.ecc_type == 'C']

        with open(output_path, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("eccDNA Simulation QC Report v4 (minimap2) - Final\n")
            f.write("=" * 70 + "\n\n")

            # ===== 参数配置 =====
            f.write("## Simulation Parameters\n")
            f.write("-" * 40 + "\n")
            f.write(f"Reference: {self.config.reference}\n")
            f.write(f"Random seed: {self.config.seed}\n")
            f.write(f"Target: {self.config.num_unique} Unique, {self.config.num_multi} Multi, {self.config.num_chimeric} Chimeric\n")
            f.write(f"Lognormal mode: {self.config.mode} bp, sigma: {self.config.sigma}\n")
            f.write(f"Lognormal mu (calculated): {self.config.mu:.4f}\n")
            f.write(f"Tail weight: {self.config.tail_weight:.1%}, range: {self.config.tail_min:,}-{self.config.tail_max:,} bp\n")
            f.write(f"Length constraints: {self.config.min_length}-{self.config.max_length:,} bp\n")
            f.write(f"MeccDNA max length: {self.config.max_length_multi:,} bp\n")
            if self.config.num_chimeric > 0:
                f.write(
                    "Chimeric fragments: unique={}, unique_chrom_fraction>={:.2f}, same_chrom_prob={:.2f}, repeat_prob={:.2f}\n".format(
                        self.config.enforce_unique_chimeric_fragments,
                        self.config.chimeric_unique_chrom_fraction,
                        self.config.chimeric_same_chrom_probability,
                        self.config.chimeric_fragment_repeat_probability,
                    )
                )
                f.write(
                    "Chimeric fragment length: {}-{} bp\n".format(
                        self.config.chimeric_fragment_min_length,
                        self.config.chimeric_fragment_max_length,
                    )
                )
            f.write("\n")

            # ===== minimap2 参数 =====
            f.write("## Minimap2 Parameters\n")
            f.write("-" * 40 + "\n")
            f.write(f"Threads: {self.config.threads}\n")
            f.write(f"Identity threshold: {self.config.identity_threshold}%\n")
            f.write(f"Min coverage: {self.config.min_coverage}%\n")
            length_tol = max(0.0, 100.0 - self.config.length_consistency_threshold)
            f.write(
                "Length consistency: query span >= {:.1f}% and target span within +/-{:.1f}% of query length\n".format(
                    self.config.length_consistency_threshold,
                    length_tol,
                )
            )
            f.write(f"Multi coverage threshold (non-self hits): {self.config.multi_coverage_threshold}%\n")
            f.write(f"Max secondary: {self.config.max_secondary}\n")
            f.write(f"Split by length: {self.config.split_by_length} (threshold: {self.config.split_length} bp)\n")
            f.write(f"Preset short: {self.config.minimap_preset_short}\n")
            f.write(f"Preset long: {self.config.minimap_preset_long}\n")
            f.write(f"No-hit policy: {self.config.no_hit_policy}\n")
            f.write(f"Self overlap fraction: {self.config.self_overlap_fraction}\n")
            f.write("\n")

            # ===== 结果摘要 =====
            f.write("## Results Summary\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total regions: {len(regions):,}\n")
            f.write(f"  - UeccDNA (Unique):   {len(u_regions):,}\n")
            f.write(f"  - MeccDNA (Multi):    {len(m_regions):,}\n")
            f.write(f"  - CeccDNA (Chimeric): {len(c_regions):,}\n")
            f.write("\n")
            warnings = []
            if len(u_regions) < self.config.num_unique:
                warnings.append(
                    f"Unique under-generated: {len(u_regions)} / {self.config.num_unique}"
                )
            if len(m_regions) < self.config.num_multi:
                warnings.append(
                    f"Multi under-generated: {len(m_regions)} / {self.config.num_multi}"
                )
            if len(c_regions) < self.config.num_chimeric:
                warnings.append(
                    f"Chimeric under-generated: {len(c_regions)} / {self.config.num_chimeric}"
                )
            if warnings:
                f.write("## Warnings\n")
                f.write("-" * 40 + "\n")
                for w in warnings:
                    f.write(f"- {w}\n")
                f.write("\n")

            # ===== 各类型长度统计 =====
            for type_name, type_regions in [('All', regions), ('Unique', u_regions), ('Multi', m_regions), ('Chimeric', c_regions)]:
                if not type_regions:
                    continue

                lengths = np.array([r.length for r in type_regions])
                f.write(f"## {type_name} Length Statistics (n={len(type_regions):,})\n")
                f.write("-" * 40 + "\n")
                f.write(f"Min:    {np.min(lengths):>12,} bp\n")
                f.write(f"Max:    {np.max(lengths):>12,} bp\n")
                f.write(f"Mean:   {np.mean(lengths):>12,.1f} bp\n")
                f.write(f"Median: {np.median(lengths):>12,.1f} bp\n")
                f.write(f"Std:    {np.std(lengths):>12,.1f} bp\n")
                f.write("\n")

                # 分位数
                f.write("Percentiles:\n")
                for p in [1, 5, 25, 50, 75, 95, 99]:
                    f.write(f"  {p:>2}%:  {np.percentile(lengths, p):>12,.1f} bp\n")
                f.write("\n")

            # ===== 长度区间分布 =====
            if regions:
                lengths = [r.length for r in regions]
                total = len(lengths)
                f.write("## Length Bin Distribution (All)\n")
                f.write("-" * 40 + "\n")
                bins = [
                    ('< 500 bp', sum(1 for l in lengths if l < 500)),
                    ('500 bp - 1 kb', sum(1 for l in lengths if 500 <= l < 1000)),
                    ('1 kb - 5 kb', sum(1 for l in lengths if 1000 <= l < 5000)),
                    ('5 kb - 50 kb', sum(1 for l in lengths if 5000 <= l < 50000)),
                    ('50 kb - 500 kb', sum(1 for l in lengths if 50000 <= l <= 500000)),
                ]
                for label, count in bins:
                    pct = 100 * count / total if total > 0 else 0
                    f.write(f"{label:>15}: {count:>8,} ({pct:>5.1f}%)\n")
                f.write("\n")

            # ===== 组分统计 =====
            main_regions = [r for r in regions if r.component == 'main']
            tail_regions = [r for r in regions if r.component == 'tail']
            chimeric_component = [r for r in regions if r.component == 'chimeric']
            
            if main_regions or tail_regions:
                f.write("## Component Statistics (main vs tail)\n")
                f.write("-" * 40 + "\n")
                total = len(regions)
                f.write(f"Main (lognormal):   {len(main_regions):>8,} ({100*len(main_regions)/total:>5.1f}%)\n")
                f.write(f"Tail (log-uniform): {len(tail_regions):>8,} ({100*len(tail_regions)/total:>5.1f}%)\n")
                f.write(f"Chimeric:           {len(chimeric_component):>8,} ({100*len(chimeric_component)/total:>5.1f}%)\n")
                if main_regions:
                    f.write(f"Main mean length:   {np.mean([r.length for r in main_regions]):>12,.1f} bp\n")
                if tail_regions:
                    f.write(f"Tail mean length:   {np.mean([r.length for r in tail_regions]):>12,.1f} bp\n")
                f.write("\n")

            # ===== 嵌合体片段分布 =====
            if c_regions:
                f.write("## Chimeric Fragment Distribution\n")
                f.write("-" * 40 + "\n")
                frag_counts = defaultdict(int)
                for r in c_regions:
                    frag_counts[len(r.fragments)] += 1
                for n_frags, count in sorted(frag_counts.items()):
                    pct = 100 * count / len(c_regions)
                    f.write(f"{n_frags} fragments: {count:>6,} ({pct:>5.1f}%)\n")
                f.write("\n")

            # ===== 染色体分布 =====
            if regions:
                f.write("## Chromosome Distribution (Top 15)\n")
                f.write("-" * 40 + "\n")
                chrom_counts = defaultdict(int)
                for r in regions:
                    chrom_counts[r.chrom] += 1
                total = len(regions)
                for chrom, count in sorted(chrom_counts.items(), key=lambda x: -x[1])[:15]:
                    pct = 100 * count / total
                    f.write(f"{chrom}: {count:,} ({pct:.1f}%)\n")
                f.write("\n")

            # ===== 采样器统计 =====
            if self.generator is not None:
                sampler_stats = self.generator.length_sampler.get_stats()
                f.write("## Sampler Statistics\n")
                f.write("-" * 40 + "\n")
                f.write(f"Total sampled: {sampler_stats['total_sampled']:,}\n")
                f.write(f"Main count: {sampler_stats['main_count']:,}\n")
                f.write(f"Tail count: {sampler_stats['tail_count']:,}\n")
                f.write(f"Rejection count: {sampler_stats['rejection_count']:,}\n")
                f.write("\n")

            f.write("=" * 70 + "\n")
            f.write("End of Report\n")

    def _write_length_data(self, regions: List[EccRegion]):
        """写入长度数据文件"""
        output_path = f"{self.config.output_prefix}.lengths.tsv"
        logging.info(f"Writing length data: {output_path}")
        with open(output_path, 'w') as f:
            f.write("region_id\tlength\tlog10_length\ttype\tcomponent\tchrom\tn_fragments\tn_hits\n")
            for r in regions:
                n_frags = len(r.fragments) if r.fragments else 1
                n_hits = len(r.multi_hits) if r.multi_hits else 0
                f.write(f"{r.region_id}\t{r.length}\t{np.log10(r.length):.4f}\t{r.ecc_type}\t{r.component}\t{r.chrom}\t{n_frags}\t{n_hits}\n")

    def _write_multi_hits(self, m_regions: List[EccRegion]):
        """写入 MeccDNA 多处匹配位置文件"""
        output_path = f"{self.config.output_prefix}.multi_hits.tsv"
        logging.info(f"Writing multi-hit locations: {output_path}")
        
        with open(output_path, 'w') as f:
            # 写入表头
            f.write("# MeccDNA Multi-mapping Locations\n")
            f.write(
                "# Multi is called only if there is at least one non-self hit with >= {:.1f}% identity, "
                ">= {:.1f}% length consistency, and >= {:.1f}% coverage\n".format(
                    self.config.identity_threshold,
                    self.config.length_consistency_threshold,
                    self.config.multi_coverage_threshold,
                )
            )
            f.write("# hit_* columns exclude the source locus (see source_* columns)\n")
            f.write("#\n")
            f.write("region_id\tsource_chrom\tsource_start\tsource_end\tsource_length\t")
            f.write("n_hits\thit_chrom\thit_start\thit_end\thit_identity\thit_coverage\n")
            
            for r in m_regions:
                source_info = f"{r.region_id}\t{r.chrom}\t{r.start}\t{r.end}\t{r.length}"
                n_hits = len(r.multi_hits)
                
                if r.multi_hits:
                    for hit_chrom, hit_start, hit_end, hit_ident, hit_cov in r.multi_hits:
                        f.write(f"{source_info}\t{n_hits}\t{hit_chrom}\t{hit_start}\t{hit_end}\t{hit_ident:.2f}\t{hit_cov:.2f}\n")
                else:
                    # 不应该发生，但防御性处理
                    f.write(f"{source_info}\t0\t.\t.\t.\t.\t.\n")
        
        # 同时输出一个汇总文件
        summary_path = f"{self.config.output_prefix}.multi_hits_summary.tsv"
        logging.info(f"Writing multi-hit summary: {summary_path}")
        
        with open(summary_path, 'w') as f:
            f.write("region_id\tsource\tlength\tn_hits\thits_detail\n")
            for r in m_regions:
                source = f"{r.chrom}:{r.start}-{r.end}"
                n_hits = len(r.multi_hits)
                
                # 格式化匹配详情
                hits_detail = ';'.join([
                    f"{c}:{s}-{e}({ident:.1f}%)"
                    for c, s, e, ident, cov in r.multi_hits
                ])
                
                f.write(f"{r.region_id}\t{source}\t{r.length}\t{n_hits}\t{hits_detail}\n")

    def _generate_ml_training_data(self, regions: List[EccRegion]):
        """生成 ML 训练数据：二倍化序列 → minimap2 比对 → BLAST 格式输出

        这个功能模拟 CircleSeeker 的处理流程：
        1. 将 eccDNA 序列二倍化（模拟 circular DNA 的 rolling circle 特性）
        2. 用 minimap2 比对到参考基因组
        3. 将 PAF 转换为 BLAST outfmt 6 + mapq 格式
        4. 过滤低 identity 的比对

        输出文件可用于 ML 模型训练，因为每个 eccDNA 的真实类型（U/M/C）已知。
        """
        if not regions:
            logging.warning("No regions to generate ML training data")
            return

        assert self.sequences is not None

        logging.info("=" * 60)
        logging.info("Generating ML Training Data")
        logging.info("=" * 60)

        prefix = self.config.output_prefix

        # Step 1: 写入二倍化 FASTA（按 CircleSeeker 格式）
        doubled_fasta = f"{prefix}.doubled.fa"
        self._write_doubled_fasta(regions, doubled_fasta)

        # Step 2: minimap2 比对
        paf_file = f"{prefix}.ml_align.paf"
        self._run_ml_alignment(doubled_fasta, paf_file)

        # Step 3: PAF → BLAST 格式 + identity 过滤
        blast_file = f"{prefix}.ml_align.tsv"
        n_written = self._paf_to_blast_filtered(paf_file, blast_file)

        logging.info(f"ML training data generated: {n_written} alignments")
        logging.info(f"  Doubled FASTA: {doubled_fasta}")
        logging.info(f"  Alignment TSV: {blast_file}")

    def _write_doubled_fasta(self, regions: List[EccRegion], output_path: str):
        """写入二倍化的 FASTA 文件

        格式：>{region_id}|1|{length}|2.0|circular
        序列：随机旋转起点后重复两次（模拟 circular DNA rolling circle 的随机线性化）

        真实的 eccDNA 检测中，环状 DNA 的线性化起点是随机的，
        所以我们随机选择一个起点来旋转序列，使模拟更真实。
        """
        assert self.sequences is not None
        logging.info(f"Writing doubled FASTA: {output_path}")

        with open(output_path, 'w') as f:
            for r in regions:
                seq = r.get_sequence(self.sequences)
                length = len(seq)

                # 随机选择一个旋转起点（模拟环状 DNA 的随机线性化）
                rotate_pos = int(self.rng.integers(0, length))
                rotated_seq = seq[rotate_pos:] + seq[:rotate_pos]

                # CircleSeeker 格式: readName|repN|consLen|copyNum|circular
                # repN=1 表示第一个重复单元
                # copyNum=2.0 表示二倍化
                header = f">{r.region_id}|1|{length}|2.0|circular"
                doubled_seq = rotated_seq + rotated_seq  # 二倍化
                f.write(header + "\n")
                for i in range(0, len(doubled_seq), 80):
                    f.write(doubled_seq[i:i+80] + "\n")

        logging.info(f"Wrote {len(regions)} doubled sequences (with random rotation)")

    def _run_ml_alignment(self, query_fasta: str, output_paf: str):
        """运行 minimap2 比对用于 ML 训练数据"""
        assert self.aligner is not None or self.reference_for_align is not None

        # 使用已有的索引或创建新的 aligner
        if self.aligner is None:
            self.aligner = Minimap2Manager(self.config, self.config.reference)
            self.aligner.setup()

        assert self.aligner.index_short is not None
        assert self.aligner.work_dir is not None

        logging.info("Running minimap2 for ML training data...")

        # 使用 map-hifi preset，与 CircleSeeker HiFi 检测保持一致
        cmd = [
            'minimap2',
            '-x', 'map-hifi',
            '--secondary=yes',
            '-N', str(self.config.max_secondary),
            '-c',  # 输出 CIGAR
            '--cs',  # 输出 cs tag
            '-t', str(self.config.threads),
            '-o', output_paf,
            self.aligner.reference_for_align,
            query_fasta
        ]

        logging.info(f"minimap2 command: {' '.join(cmd[:10])}...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f"minimap2 failed: {result.stderr}")
            raise RuntimeError(f"minimap2 failed: {result.stderr}")

        logging.info(f"PAF written to: {output_paf}")

    def _paf_to_blast_filtered(self, paf_path: str, output_path: str) -> int:
        """将 PAF 转换为 BLAST outfmt 6 + mapq + truth_type 格式，并过滤低 identity 比对

        输出列（TSV，无表头）：
          query_id, subject_id, identity, alignment_length, mismatches, gap_opens,
          q_start, q_end, s_start, s_end, evalue, bit_score, sstrand, mapq, truth_type

        truth_type 列：U (Unique), M (Multi), C (Chimeric)，从 query_id 中提取

        Args:
            paf_path: 输入 PAF 文件路径
            output_path: 输出 TSV 文件路径

        Returns:
            写入的比对记录数
        """
        import re

        logging.info(f"Converting PAF to BLAST format: {paf_path} -> {output_path}")

        def extract_truth_type(query_id: str) -> str:
            """从 query_id 提取真实类型"""
            if query_id.startswith('UeccDNA_'):
                return 'U'
            elif query_id.startswith('MeccDNA_'):
                return 'M'
            elif query_id.startswith('CeccDNA_'):
                return 'C'
            return 'X'  # Unknown

        identity_threshold = self.config.ml_identity_threshold

        # cs tag 解析辅助函数
        CS_OPS = {":", "*", "+", "-", "~"}

        def parse_cs_tag(cs: str) -> tuple:
            """解析 minimap2 cs tag"""
            mismatches = 0
            gap_opens = 0
            gap_bases = 0
            i = 0
            while i < len(cs):
                op = cs[i]
                if op == ":":
                    i += 1
                    while i < len(cs) and cs[i].isdigit():
                        i += 1
                elif op == "*":
                    mismatches += 1
                    i += 3 if i + 2 < len(cs) else len(cs)
                elif op in "+-":
                    gap_opens += 1
                    i += 1
                    start = i
                    while i < len(cs) and cs[i] not in CS_OPS:
                        i += 1
                    gap_bases += i - start
                elif op == "~":
                    gap_opens += 1
                    i += 1
                    start = i
                    while i < len(cs) and cs[i] not in CS_OPS:
                        i += 1
                    segment = cs[start:i]
                    match = re.search(r"(\d+)", segment)
                    if match:
                        gap_bases += int(match.group(1))
                else:
                    i += 1
            return mismatches, gap_opens, gap_bases

        def extract_cs_tag(tags: List[str]) -> Optional[str]:
            for tag in tags:
                if tag.startswith("cs:Z:"):
                    return tag[5:]
            return None

        written = 0
        skipped = 0
        identity_filtered = 0

        with open(paf_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line_num, line in enumerate(infile, 1):
                if not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 12:
                    skipped += 1
                    continue

                try:
                    qname = parts[0]
                    qlen = int(parts[1])
                    qstart = int(parts[2])
                    qend = int(parts[3])
                    strand = parts[4]
                    tname = parts[5]
                    tstart = int(parts[7])
                    tend = int(parts[8])
                    nmatch = int(parts[9])
                    alen = int(parts[10])
                    mapq = int(parts[11])
                except (TypeError, ValueError, IndexError):
                    skipped += 1
                    continue

                # 解析 cs tag
                cs_tag = extract_cs_tag(parts[12:])
                if cs_tag:
                    mismatches, gap_opens, gap_bases = parse_cs_tag(cs_tag)
                else:
                    gap_bases = 0
                    mismatches = max(alen - nmatch, 0)
                    gap_opens = 0

                identity = (nmatch / alen * 100.0) if alen else 0.0

                # Identity 过滤
                if identity < identity_threshold:
                    identity_filtered += 1
                    continue

                # 转换坐标（1-based, strand-oriented，与 BLAST 一致）
                q_start = qstart + 1
                q_end = qend

                if strand == "+":
                    s_start = tstart + 1
                    s_end = tend
                    sstrand = "plus"
                else:
                    s_start = tend
                    s_end = tstart + 1
                    sstrand = "minus"

                # 从 query_id 提取真实类型（去掉二倍化后缀 |1|xxx|2.0|circular）
                base_qname = qname.split('|')[0]
                truth_type = extract_truth_type(base_qname)

                out_fields = [
                    qname,
                    tname,
                    f"{identity:.2f}",
                    str(alen),
                    str(mismatches),
                    str(gap_opens),
                    str(q_start),
                    str(q_end),
                    str(s_start),
                    str(s_end),
                    "0",  # evalue
                    "0",  # bit_score
                    sstrand,
                    str(mapq),
                    truth_type,  # 真实类型标签
                ]
                outfile.write("\t".join(out_fields) + "\n")
                written += 1

        if skipped > 0:
            logging.warning(f"Skipped {skipped} malformed lines in PAF")
        if identity_filtered > 0:
            logging.info(f"Filtered {identity_filtered} alignments with identity < {identity_threshold}%")

        return written

    def _print_summary(self, regions: List[EccRegion]):
        """打印摘要"""
        u = sum(1 for r in regions if r.ecc_type == 'U')
        m = sum(1 for r in regions if r.ecc_type == 'M')
        c = sum(1 for r in regions if r.ecc_type == 'C')
        
        logging.info("")
        logging.info("=" * 60)
        logging.info("Simulation Complete!")
        logging.info("=" * 60)
        logging.info(f"Generated: {len(regions)} total eccDNA regions")
        logging.info(f"  - UeccDNA (Unique):   {u}")
        logging.info(f"  - MeccDNA (Multi):    {m}")
        logging.info(f"  - CeccDNA (Chimeric): {c}")
        if u < self.config.num_unique:
            logging.warning(f"Unique under-generated: {u} / {self.config.num_unique}")
        if m < self.config.num_multi:
            logging.warning(f"Multi under-generated: {m} / {self.config.num_multi}")
        if c < self.config.num_chimeric:
            logging.warning(f"Chimeric under-generated: {c} / {self.config.num_chimeric}")
        logging.info("")
        logging.info("Output files:")
        logging.info(f"  {self.config.output_prefix}.all.bed")
        logging.info(f"  {self.config.output_prefix}.all.fa")
        logging.info(f"  {self.config.output_prefix}.qc.log")
        logging.info(f"  {self.config.output_prefix}.lengths.tsv")
        if u > 0:
            logging.info(f"  {self.config.output_prefix}.unique.bed/fa")
        if m > 0:
            logging.info(f"  {self.config.output_prefix}.multi.bed/fa")
            logging.info(f"  {self.config.output_prefix}.multi_hits.tsv")
            logging.info(f"  {self.config.output_prefix}.multi_hits_summary.tsv")
        if c > 0:
            logging.info(f"  {self.config.output_prefix}.chimeric.bed/fa")
        if self.config.generate_ml_data:
            logging.info(f"  {self.config.output_prefix}.doubled.fa (ML training)")
            logging.info(f"  {self.config.output_prefix}.ml_align.tsv (ML training)")

    def _write_error_log(self) -> None:
        """在异常情况下写入最小错误日志，避免无日志输出"""
        import traceback

        output_path = f"{self.config.output_prefix}.qc.log"
        try:
            with open(output_path, "w") as f:
                f.write("=" * 70 + "\n")
                f.write("eccDNA Simulation QC Report - Error\n")
                f.write("=" * 70 + "\n\n")
                f.write("An error occurred during simulation.\n\n")
                f.write(traceback.format_exc())
            logging.error(f"Simulation failed; error log written to: {output_path}")
        except Exception as exc:
            logging.error(f"Failed to write error log: {exc}")


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    """解析命令行参数"""
    p = argparse.ArgumentParser(
        description="Simulate eccDNA regions with U/M/C classification (minimap2 v4 Final)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 基本用法（自动创建 sim_ecc/ 目录）
  ecc sim-region -r hg38.fa -u 10000 -o sim_ecc -t 16
  # 输出: sim_ecc/sim_ecc.all.bed, sim_ecc/sim_ecc.all.fa, ...

  # 指定输出目录（用路径作为 prefix）
  ecc sim-region -r hg38.fa -u 10000 -o /path/to/output/sim_ecc -t 16
  # 输出: /path/to/output/sim_ecc/sim_ecc.all.bed, ...

  # 生成所有类型
  ecc sim-region -r hg38.fa -u 7000 -m 2000 -c 1000 -o sim_ecc -t 16

  # 自定义分布参数
  ecc sim-region -r hg38.fa -u 10000 -o sim_ecc \\
      --mode 400 --sigma 0.8 --tail-weight 0.05

Types:
  UeccDNA (Unique): Only one genomic location matches at >=99% identity with >=99% length consistency
  MeccDNA (Multi):  Multiple genomic locations match at >=99% identity with >=99% length consistency
  CeccDNA (Chimeric): Composed of 2-5 fragments from different locations
        """
    )

    # 输入输出
    p.add_argument('-r', '--reference', required=True,
                   help='Reference FASTA (.fa/.fasta/.fna or .gz)')
    p.add_argument('-o', '--output-prefix', required=True,
                   help='Output prefix (will also be used as directory name if --output-dir not specified)')
    p.add_argument('-d', '--output-dir', type=str, default=None,
                   help='Output directory (default: create directory named after prefix)')

    # 目标数量
    p.add_argument('-u', '--num-unique', type=int, default=1000,
                   help='Number of UeccDNA (Unique) to generate (default: 1000)')
    p.add_argument('-m', '--num-multi', type=int, default=0,
                   help='Number of MeccDNA (Multi-mapped) to generate (default: 0)')
    p.add_argument('-c', '--num-chimeric', type=int, default=0,
                   help='Number of CeccDNA (Chimeric) to generate (default: 0)')

    # minimap2 参数
    p.add_argument('-t', '--threads', type=int, default=4,
                   help='Number of threads for minimap2 (default: 4)')
    p.add_argument('--identity', type=float, default=99.0,
                   help='Identity threshold %% (default: 99.0)')
    p.add_argument('--min-coverage', type=float, default=90.0,
                   help='Minimum coverage %% (default: 90.0)')
    p.add_argument('--length-consistency', type=float, default=99.0,
                   help='Length consistency threshold %% (default: 99.0)')
    p.add_argument('--max-secondary', type=int, default=50,
                   help='Max secondary alignments (default: 50)')

    # 分流设置
    p.add_argument('--split-by-length', action='store_true',
                   help='Enable short/long split mapping (default: disabled)')
    p.add_argument('--no-split-by-length', action='store_true',
                   help='Disable short/long split mapping (kept for compatibility)')
    p.add_argument('--split-length', type=int, default=5000,
                   help='Length threshold for long preset (default: 5000)')
    p.add_argument('--preset-short', type=str, default='map-hifi',
                   help='Minimap2 preset for short reads (default: map-hifi)')
    p.add_argument('--preset-long', type=str, default='map-hifi',
                   help='Minimap2 preset for long reads (default: map-hifi)')

    # no-hit 处理
    p.add_argument('--no-hit-policy', choices=['skip', 'unique'], default='skip',
                   help='How to handle no-hit regions (default: skip)')

    # 索引设置
    p.add_argument('--index-dir', type=str, default=None,
                   help='Directory to store/reuse .mmi (default: reference directory)')
    p.add_argument('--no-reuse-index', action='store_true',
                   help='Disable index reuse (always rebuild)')
    p.add_argument('--force-rebuild-index', action='store_true',
                   help='Force rebuild index even if exists')

    # 分布参数
    p.add_argument('--seed', type=int, default=42,
                   help='Random seed (default: 42)')
    p.add_argument('--mode', type=float, default=400.0,
                   help='Lognormal mode/peak in bp (default: 400)')
    p.add_argument('--sigma', type=float, default=0.8,
                   help='Lognormal sigma (default: 0.8)')
    p.add_argument('--tail-weight', type=float, default=0.05,
                   help='Tail distribution weight (default: 0.05)')
    p.add_argument('--tail-min', type=int, default=5000,
                   help='Tail minimum length bp (default: 5000)')
    p.add_argument('--tail-max', type=int, default=500000,
                   help='Tail maximum length bp (default: 500000)')

    # 长度约束
    p.add_argument('--min-length', type=int, default=100,
                   help='Global minimum length bp (default: 100)')
    p.add_argument('--max-length', type=int, default=500000,
                   help='Global maximum length bp (default: 500000)')
    p.add_argument('--max-length-multi', type=int, default=10000,
                   help='MeccDNA maximum length bp (default: 10000)')

    # 候选倍数
    p.add_argument('--multiplier-u', type=float, default=2.0,
                   help='Candidate multiplier for Unique (default: 2.0)')
    p.add_argument('--multiplier-m', type=float, default=10.0,
                   help='Candidate multiplier for Multi (default: 10.0)')

    # 其他
    p.add_argument('--keep-tmp', action='store_true',
                   help='Keep temporary minimap2 files (for debugging)')
    p.add_argument('--exclude-chroms', nargs='+',
                   default=['chrM', 'chrUn', '_random', '_alt', '_hap'],
                   help='Chromosome patterns to exclude')
    p.add_argument('-v', '--verbose', action='store_true',
                   help='Verbose logging')

    # ML 训练数据生成
    p.add_argument('--generate-ml-data', action='store_true',
                   help='Generate ML training data (doubled sequences + alignment TSV)')
    p.add_argument('--ml-identity', type=float, default=99.0,
                   help='Identity threshold %% for ML training data (default: 99.0)')

    # 内部串联重复过滤
    p.add_argument('--no-filter-internal-repeats', action='store_true',
                   help='Disable filtering of regions with internal tandem repeats')
    p.add_argument('--internal-repeat-ratio', type=float, default=0.8,
                   help='Filter regions where TideHunter consensus / region_length < ratio (default: 0.8)')

    args = p.parse_args()
    if args.split_by_length and args.no_split_by_length:
        p.error("Cannot set both --split-by-length and --no-split-by-length.")
    return args


def main():
    """主函数"""
    args = parse_args()

    # 设置日志
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    split_by_length = args.split_by_length and not args.no_split_by_length

    # 创建配置
    config = SimulationConfig(
        reference=args.reference,
        output_prefix=args.output_prefix,
        output_dir=args.output_dir,
        num_unique=args.num_unique,
        num_multi=args.num_multi,
        num_chimeric=args.num_chimeric,
        seed=args.seed,
        mode=args.mode,
        sigma=args.sigma,
        tail_weight=args.tail_weight,
        tail_min=args.tail_min,
        tail_max=args.tail_max,
        min_length=args.min_length,
        max_length=args.max_length,
        max_length_multi=args.max_length_multi,
        threads=args.threads,
        identity_threshold=args.identity,
        min_coverage=args.min_coverage,
        length_consistency_threshold=args.length_consistency,
        max_secondary=args.max_secondary,
        minimap_preset_short=args.preset_short,
        minimap_preset_long=args.preset_long,
        split_by_length=split_by_length,
        split_length=args.split_length,
        no_hit_policy=args.no_hit_policy,
        index_dir=args.index_dir,
        reuse_index=(not args.no_reuse_index),
        force_rebuild_index=args.force_rebuild_index,
        candidate_multiplier_u=args.multiplier_u,
        candidate_multiplier_m=args.multiplier_m,
        exclude_chroms=args.exclude_chroms,
        keep_tmp=args.keep_tmp,
        generate_ml_data=args.generate_ml_data,
        ml_identity_threshold=args.ml_identity,
        filter_internal_repeats=(not args.no_filter_internal_repeats),
        internal_repeat_ratio_threshold=args.internal_repeat_ratio,
    )

    # 运行模拟
    sim = EccDNASimulator(config)
    sim.run()


if __name__ == '__main__':
    main()
