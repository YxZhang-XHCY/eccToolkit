"""
Nanopore错误模型

特点：
- 较高的错误率（~5-15%取决于版本）
- 系统性同聚物误差
- 随机替换和indel
"""

from typing import Tuple, Optional, List
import numpy as np

from .base import BaseErrorModel


class NanoporeErrorModel(BaseErrorModel):
    """Nanopore测序错误模型"""
    
    def __init__(
        self,
        substitution_rate: float = 0.03,
        insertion_rate: float = 0.02,
        deletion_rate: float = 0.03,
        homopolymer_error_rate: float = 0.15,
        rng: Optional[np.random.Generator] = None
    ):
        """
        Args:
            substitution_rate: 替换错误率
            insertion_rate: 插入错误率
            deletion_rate: 删除错误率
            homopolymer_error_rate: 同聚物区域额外错误率
        """
        super().__init__(rng)
        self.substitution_rate = substitution_rate
        self.insertion_rate = insertion_rate
        self.deletion_rate = deletion_rate
        self.homopolymer_error_rate = homopolymer_error_rate
        self._bases = ['A', 'T', 'G', 'C']
    
    def _find_homopolymers(self, sequence: str, min_length: int = 4) -> List[Tuple[int, int]]:
        """找出同聚物区域"""
        regions = []
        i = 0
        while i < len(sequence):
            base = sequence[i]
            run_start = i
            while i < len(sequence) and sequence[i] == base:
                i += 1
            run_len = i - run_start
            if run_len >= min_length:
                regions.append((run_start, i))
        return regions
    
    def apply(self, sequence: str, quality: Optional[str] = None) -> Tuple[str, str]:
        """
        应用 Nanopore 错误模型（概率归一化版本）

        使用两阶段采样避免概率溢出：
        1. 先计算各类错误的调整后概率（min(1.0, ...)）
        2. 使用 multinomial 采样决定发生哪种事件

        事件类型：
        - no_error: 保持原样
        - substitution: 替换为其他碱基
        - insertion: 在此位置后插入碱基
        - deletion: 删除此碱基
        """
        sequence = sequence.upper()

        # 找出同聚物区域
        homopolymer_regions = self._find_homopolymers(sequence)
        in_homopolymer = set()
        for start, end in homopolymer_regions:
            for i in range(start, end):
                in_homopolymer.add(i)

        result = []
        qual_values = []

        for i, base in enumerate(sequence):
            is_homo = i in in_homopolymer

            # 计算调整后的错误率，确保不超过1.0
            if is_homo:
                sub_rate = min(1.0, self.substitution_rate + self.homopolymer_error_rate)
                ins_rate = min(1.0, self.insertion_rate + self.homopolymer_error_rate)
                del_rate = min(1.0, self.deletion_rate + self.homopolymer_error_rate * 2)
            else:
                sub_rate = self.substitution_rate
                ins_rate = self.insertion_rate
                del_rate = self.deletion_rate

            # 归一化：确保总错误率不超过1.0
            total_error_rate = sub_rate + ins_rate + del_rate
            if total_error_rate > 1.0:
                # 按比例缩放
                scale = 0.95 / total_error_rate  # 留5%给无错误
                sub_rate *= scale
                ins_rate *= scale
                del_rate *= scale

            no_error_rate = 1.0 - sub_rate - ins_rate - del_rate

            # multinomial 采样
            r = self.rng.random()
            if r < del_rate:
                # 删除：跳过此碱基
                continue
            elif r < del_rate + ins_rate:
                # 插入：先添加一个插入碱基
                if is_homo:
                    inserted = base
                else:
                    inserted = self.rng.choice(self._bases)
                result.append(inserted)
                qual_values.append(8)
                # 然后继续添加原碱基（可能被替换）
                if self.rng.random() < sub_rate / (sub_rate + no_error_rate + 0.001):
                    if base in self._bases:
                        alternatives = [b for b in self._bases if b != base]
                        base = self.rng.choice(alternatives)
                    qual_values.append(10)
                else:
                    qual_values.append(15 if is_homo else 18)
                result.append(base)
            elif r < del_rate + ins_rate + sub_rate:
                # 替换
                if base in self._bases:
                    alternatives = [b for b in self._bases if b != base]
                    base = self.rng.choice(alternatives)
                result.append(base)
                qual_values.append(10)
            else:
                # 无错误
                result.append(base)
                qual_values.append(15 if is_homo else 18)

        result_seq = "".join(result)

        # 生成质量字符串
        if qual_values:
            qual_values = np.array(qual_values)
            qual_values = qual_values + self.rng.integers(-3, 4, len(qual_values))
            qual_values = np.clip(qual_values, 2, 41)
            quality = "".join(chr(q + 33) for q in qual_values)
        else:
            quality = ""

        return result_seq, quality
    
    @property
    def name(self) -> str:
        return "nanopore"


class NanoporeR10ErrorModel(NanoporeErrorModel):
    """Nanopore R10（更新版本，更低错误率）"""
    
    def __init__(self, rng: Optional[np.random.Generator] = None):
        super().__init__(
            substitution_rate=0.01,
            insertion_rate=0.01,
            deletion_rate=0.01,
            homopolymer_error_rate=0.05,
            rng=rng
        )
    
    @property
    def name(self) -> str:
        return "nanopore_r10"
