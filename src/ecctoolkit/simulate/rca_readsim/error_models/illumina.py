"""
Illumina错误模型

特点：
- 主要是替换错误
- 错误率随位置增加（尤其是3'端）
- 质量分数反映错误概率
"""

from typing import Tuple, Optional
import numpy as np

from .base import BaseErrorModel


class IlluminaErrorModel(BaseErrorModel):
    """Illumina测序错误模型"""
    
    def __init__(
        self,
        base_error_rate: float = 0.001,
        end_error_rate: float = 0.01,
        rng: Optional[np.random.Generator] = None
    ):
        """
        Args:
            base_error_rate: 基础错误率（5'端）
            end_error_rate: 末端错误率（3'端）
        """
        super().__init__(rng)
        self.base_error_rate = base_error_rate
        self.end_error_rate = end_error_rate
        self._bases = np.array(['A', 'T', 'G', 'C'])
    
    def apply(self, sequence: str, quality: Optional[str] = None) -> Tuple[str, str]:
        seq_array = np.array(list(sequence.upper()))
        length = len(seq_array)
        
        # 计算位置相关的错误率
        positions = np.arange(length)
        error_rates = self.base_error_rate + (
            (self.end_error_rate - self.base_error_rate) * 
            positions / max(length - 1, 1)
        )
        
        # 决定哪些位置出错
        errors = self.rng.random(length) < error_rates
        
        # 对出错位置随机替换
        for i in np.where(errors)[0]:
            original = seq_array[i]
            alternatives = self._bases[self._bases != original]
            seq_array[i] = self.rng.choice(alternatives)
        
        result_seq = "".join(seq_array)
        
        # 生成质量字符串
        if quality is None:
            # 根据错误率计算质量值 Q = -10 * log10(p)
            q_values = -10 * np.log10(np.maximum(error_rates, 1e-10))
            q_values = np.clip(q_values, 2, 41).astype(int)
            # 添加一些随机波动
            q_values = q_values + self.rng.integers(-3, 4, length)
            q_values = np.clip(q_values, 2, 41)
            quality = "".join(chr(q + 33) for q in q_values)
        
        return result_seq, quality
    
    @property
    def name(self) -> str:
        return "illumina"


class IlluminaQualityAwareErrorModel(BaseErrorModel):
    """基于质量分数的Illumina错误模型"""
    
    def __init__(self, rng: Optional[np.random.Generator] = None):
        super().__init__(rng)
        self._bases = np.array(['A', 'T', 'G', 'C'])
    
    def apply(self, sequence: str, quality: Optional[str] = None) -> Tuple[str, str]:
        if quality is None:
            # 如果没有质量，使用默认模型
            model = IlluminaErrorModel(rng=self.rng)
            return model.apply(sequence)
        
        seq_array = np.array(list(sequence.upper()))
        length = len(seq_array)
        
        # 从质量字符串计算错误率
        q_values = np.array([ord(c) - 33 for c in quality])
        error_rates = 10 ** (-q_values / 10)
        
        # 引入错误
        errors = self.rng.random(length) < error_rates
        
        for i in np.where(errors)[0]:
            original = seq_array[i]
            if original in self._bases:
                alternatives = self._bases[self._bases != original]
                seq_array[i] = self.rng.choice(alternatives)
        
        return "".join(seq_array), quality
    
    @property
    def name(self) -> str:
        return "illumina_qa"
