"""
错误模型模块

提供三个平台的测序错误模拟
"""

from .base import BaseErrorModel, IdentityErrorModel
from .illumina import IlluminaErrorModel, IlluminaQualityAwareErrorModel
from .pacbio import PacBioHiFiErrorModel, PacBioCLRErrorModel
from .nanopore import NanoporeErrorModel, NanoporeR10ErrorModel


def get_error_model(name: str, **kwargs) -> BaseErrorModel:
    """
    根据名称获取错误模型
    
    Args:
        name: 模型名称 (identity, illumina, pacbio_hifi, pacbio_clr, 
                       nanopore, nanopore_r10)
        **kwargs: 传递给模型的参数
    
    Returns:
        错误模型实例
    """
    models = {
        "identity": IdentityErrorModel,
        "illumina": IlluminaErrorModel,
        "illumina_qa": IlluminaQualityAwareErrorModel,
        "pacbio_hifi": PacBioHiFiErrorModel,
        "pacbio_clr": PacBioCLRErrorModel,
        "nanopore": NanoporeErrorModel,
        "nanopore_r10": NanoporeR10ErrorModel,
    }
    
    if name not in models:
        raise ValueError(f"Unknown error model: {name}. Available: {list(models.keys())}")
    
    return models[name](**kwargs)


__all__ = [
    'BaseErrorModel',
    'IdentityErrorModel',
    'IlluminaErrorModel',
    'IlluminaQualityAwareErrorModel',
    'PacBioHiFiErrorModel',
    'PacBioCLRErrorModel',
    'NanoporeErrorModel',
    'NanoporeR10ErrorModel',
    'get_error_model'
]
