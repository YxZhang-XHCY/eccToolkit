"""
错误模型基类

定义序列错误引入的接口
"""

from abc import ABC, abstractmethod
from typing import Tuple, Optional
import numpy as np


class BaseErrorModel(ABC):
    """错误模型基类"""
    
    def __init__(self, rng: Optional[np.random.Generator] = None):
        self.rng = rng if rng is not None else np.random.default_rng()
    
    @abstractmethod
    def apply(self, sequence: str, quality: Optional[str] = None) -> Tuple[str, str]:
        """
        对序列应用错误模型
        
        Args:
            sequence: 原始序列
            quality: 原始质量字符串（可选）
        
        Returns:
            (带错误的序列, 质量字符串)
        """
        pass
    
    @property
    @abstractmethod
    def name(self) -> str:
        """错误模型名称"""
        pass


class IdentityErrorModel(BaseErrorModel):
    """无错误模型（直接返回原序列）"""
    
    def apply(self, sequence: str, quality: Optional[str] = None) -> Tuple[str, str]:
        if quality is None:
            quality = "I" * len(sequence)  # 默认高质量
        return sequence, quality
    
    @property
    def name(self) -> str:
        return "identity"
