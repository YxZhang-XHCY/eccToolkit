"""
背景线性DNA模块

用于模拟真实测序中的线性基因组DNA污染/背景。
这些线性DNA不经过RCA扩增，直接进行测序。
"""

import random
import math
from typing import List, Dict, Optional, Tuple, Generator
from dataclasses import dataclass
import numpy as np

import logging

from .models import (
    LinearMolecule, Segment, SegmentType, Strand,
    EccDNA, register_ecc_length
)
from .config import BackgroundDNAParams

logger = logging.getLogger(__name__)


@dataclass
class ChromosomeInfo:
    """染色体信息"""
    name: str
    length: int
    sequence: str


class BackgroundDNAGenerator:
    """
    背景线性DNA生成器

    从参考基因组中随机采样线性DNA片段，
    用于模拟测序中的非eccDNA背景信号。
    """

    def __init__(
        self,
        fasta_path: str,
        params: BackgroundDNAParams,
        seed: Optional[int] = None
    ):
        """
        初始化背景DNA生成器

        Args:
            fasta_path: 参考基因组FASTA文件路径
            params: 背景DNA参数
            seed: 随机种子
        """
        self.fasta_path = fasta_path
        self.params = params
        self.rng = np.random.default_rng(seed)
        self.random = random.Random(seed)

        # 加载染色体信息
        self.chromosomes: List[ChromosomeInfo] = []
        self.total_genome_length = 0
        self.chrom_cumsum: List[int] = []  # 用于加权随机采样

        self._load_reference()

    def _load_reference(self):
        """加载参考基因组"""
        logger.info(f"Loading reference genome from {self.fasta_path}...")

        current_name = None
        current_seq = []

        with open(self.fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('>'):
                    # 保存上一条染色体
                    if current_name is not None and current_seq:
                        seq = ''.join(current_seq)
                        self.chromosomes.append(ChromosomeInfo(
                            name=current_name,
                            length=len(seq),
                            sequence=seq
                        ))

                    # 开始新染色体
                    current_name = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line.upper())

        # 保存最后一条染色体
        if current_name is not None and current_seq:
            seq = ''.join(current_seq)
            self.chromosomes.append(ChromosomeInfo(
                name=current_name,
                length=len(seq),
                sequence=seq
            ))

        # 计算累积长度（用于加权采样）
        self.total_genome_length = sum(c.length for c in self.chromosomes)
        cumsum = 0
        for chrom in self.chromosomes:
            cumsum += chrom.length
            self.chrom_cumsum.append(cumsum)

        logger.info(f"Loaded {len(self.chromosomes)} chromosomes, "
                    f"total length: {self.total_genome_length:,} bp")

    def _sample_length(self) -> int:
        """采样片段长度"""
        if self.params.length_distribution == "uniform":
            return self.random.randint(
                self.params.min_length,
                self.params.max_length
            )
        elif self.params.length_distribution == "lognormal":
            # 将均值和标准差转换为对数正态分布参数
            mean = self.params.length_mean
            std = self.params.length_std

            # 对数正态分布参数
            # mu = ln(mean^2 / sqrt(mean^2 + std^2))
            # sigma = sqrt(ln(1 + std^2/mean^2))
            variance = std ** 2
            mu = math.log(mean ** 2 / math.sqrt(mean ** 2 + variance))
            sigma = math.sqrt(math.log(1 + variance / (mean ** 2)))

            while True:
                length = int(self.rng.lognormal(mu, sigma))
                if self.params.min_length <= length <= self.params.max_length:
                    return length
        else:
            # 默认使用均匀分布
            return self.random.randint(
                self.params.min_length,
                self.params.max_length
            )

    def _sample_region(self, length: int) -> Tuple[ChromosomeInfo, int, int]:
        """
        采样一个基因组区域

        Returns:
            (染色体, 起始位置, 结束位置)
        """
        # 按染色体长度加权随机选择染色体
        rand_pos = self.random.randint(0, self.total_genome_length - 1)

        chrom_idx = 0
        for i, cumsum in enumerate(self.chrom_cumsum):
            if rand_pos < cumsum:
                chrom_idx = i
                break

        chrom = self.chromosomes[chrom_idx]

        # 确保片段不超出染色体边界
        max_start = chrom.length - length
        if max_start < 0:
            # 染色体太短，使用整个染色体
            return chrom, 0, chrom.length

        start = self.random.randint(0, max_start)
        end = start + length

        return chrom, start, end

    def generate_fragment(self, fragment_id: int) -> LinearMolecule:
        """
        生成一个背景DNA片段

        Args:
            fragment_id: 片段ID

        Returns:
            LinearMolecule对象
        """
        # 采样长度
        length = self._sample_length()

        # 采样区域
        chrom, start, end = self._sample_region(length)
        actual_length = end - start

        # 提取序列
        sequence = chrom.sequence[start:end]

        # 随机决定链方向
        strand = Strand.FORWARD if self.random.random() < 0.5 else Strand.REVERSE

        # 创建虚拟的 EccDNA 条目用于兼容现有框架
        # 使用特殊的ID前缀 "BG_" 表示背景DNA
        bg_id = f"BG_{chrom.name}_{start}_{end}"

        # 创建 Segment
        segment = Segment(
            ecc_id=bg_id,
            ecc_offset=0,
            length=actual_length,
            strand=strand,
            segment_type=SegmentType.BACKGROUND
        )

        # 创建 LinearMolecule
        molecule = LinearMolecule(
            molecule_id=f"background_{fragment_id}",
            segments=[segment],
            source_graph_id=f"background_{fragment_id}",
            is_from_branch=False,
            is_background=True,
            background_chrom=chrom.name,
            background_start=start,
            background_end=end
        )

        return molecule, sequence, bg_id

    def generate_molecules(
        self,
        count: int,
        ecc_db: Dict[str, EccDNA]
    ) -> Tuple[List[LinearMolecule], Dict[str, EccDNA]]:
        """
        生成指定数量的背景DNA分子

        Args:
            count: 需要生成的分子数量
            ecc_db: 现有的EccDNA数据库（会被扩展）

        Returns:
            (背景分子列表, 更新后的ecc_db)
        """
        molecules = []

        logger.info(f"Generating {count} background DNA molecules...")

        for i in range(count):
            molecule, sequence, bg_id = self.generate_fragment(i)
            molecules.append(molecule)

            # 将背景DNA序列加入ecc_db以便后续读段生成
            if bg_id not in ecc_db:
                ecc_db[bg_id] = EccDNA(
                    id=bg_id,
                    seq=sequence,
                    weight=0.0  # 权重为0，不参与eccDNA采样
                )
                register_ecc_length(bg_id, len(sequence))

        logger.info(f"Generated {len(molecules)} background DNA molecules")

        return molecules, ecc_db


def calculate_background_count(
    total_target_reads: int,
    background_ratio: float
) -> int:
    """
    计算需要生成的背景DNA分子数量

    Args:
        total_target_reads: 总目标读段数
        background_ratio: 背景DNA占比

    Returns:
        需要生成的背景分子数量
    """
    # 背景读段数 = 总读段数 × 比例
    background_reads = int(total_target_reads * background_ratio)

    # 假设每个背景分子平均产生1条读段（NGS为1对）
    # 可以根据实际情况调整这个比例
    return background_reads
