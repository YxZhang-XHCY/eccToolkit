"""
PacBio HiFi错误模型

特点：
- 非常低的错误率（~0.1%）
- 主要是随机替换
- 偶尔有小indel
"""

from typing import Tuple, Optional
import numpy as np

from .base import BaseErrorModel


class PacBioHiFiErrorModel(BaseErrorModel):
    """PacBio HiFi错误模型"""
    
    def __init__(
        self,
        substitution_rate: float = 0.001,
        insertion_rate: float = 0.0001,
        deletion_rate: float = 0.0001,
        rng: Optional[np.random.Generator] = None
    ):
        """
        Args:
            substitution_rate: 替换错误率
            insertion_rate: 插入错误率
            deletion_rate: 删除错误率
        """
        super().__init__(rng)
        self.substitution_rate = substitution_rate
        self.insertion_rate = insertion_rate
        self.deletion_rate = deletion_rate
        self._bases = ['A', 'T', 'G', 'C']
    
    def apply(self, sequence: str, quality: Optional[str] = None) -> Tuple[str, str]:
        result = []
        qual_values = []
        
        for base in sequence.upper():
            # 删除
            if self.rng.random() < self.deletion_rate:
                continue
            
            # 插入（在当前碱基前）
            if self.rng.random() < self.insertion_rate:
                inserted = self.rng.choice(self._bases)
                result.append(inserted)
                qual_values.append(20)  # 低质量标记插入
            
            # 替换
            if self.rng.random() < self.substitution_rate:
                if base in self._bases:
                    alternatives = [b for b in self._bases if b != base]
                    base = self.rng.choice(alternatives)
                    qual_values.append(25)
                else:
                    qual_values.append(35)
            else:
                qual_values.append(35)  # HiFi高质量
            
            result.append(base)
        
        result_seq = "".join(result)
        
        # 生成质量字符串
        qual_values = np.array(qual_values)
        qual_values = qual_values + self.rng.integers(-2, 3, len(qual_values))
        qual_values = np.clip(qual_values, 2, 41)
        quality = "".join(chr(q + 33) for q in qual_values)
        
        return result_seq, quality
    
    @property
    def name(self) -> str:
        return "pacbio_hifi"


class PacBioCLRErrorModel(BaseErrorModel):
    """PacBio CLR错误模型（更高错误率）"""
    
    def __init__(
        self,
        substitution_rate: float = 0.01,
        insertion_rate: float = 0.05,
        deletion_rate: float = 0.02,
        rng: Optional[np.random.Generator] = None
    ):
        super().__init__(rng)
        self.substitution_rate = substitution_rate
        self.insertion_rate = insertion_rate
        self.deletion_rate = deletion_rate
        self._bases = ['A', 'T', 'G', 'C']
    
    def apply(self, sequence: str, quality: Optional[str] = None) -> Tuple[str, str]:
        result = []
        qual_values = []
        
        for base in sequence.upper():
            if self.rng.random() < self.deletion_rate:
                continue
            
            if self.rng.random() < self.insertion_rate:
                inserted = self.rng.choice(self._bases)
                result.append(inserted)
                qual_values.append(8)
            
            if self.rng.random() < self.substitution_rate:
                if base in self._bases:
                    alternatives = [b for b in self._bases if b != base]
                    base = self.rng.choice(alternatives)
                qual_values.append(10)
            else:
                qual_values.append(15)
            
            result.append(base)
        
        result_seq = "".join(result)
        
        qual_values = np.array(qual_values)
        qual_values = qual_values + self.rng.integers(-2, 3, len(qual_values))
        qual_values = np.clip(qual_values, 2, 41)
        quality = "".join(chr(q + 33) for q in qual_values)
        
        return result_seq, quality
    
    @property
    def name(self) -> str:
        return "pacbio_clr"
