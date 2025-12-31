"""
序列工具函数
"""

from typing import List, Tuple, Optional, Dict
import numpy as np


def reverse_complement(seq: str) -> str:
    """反向互补"""
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
        'N': 'N', 'n': 'n',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D'
    }
    return "".join(complement.get(base, 'N') for base in reversed(seq))


def gc_content(seq: str) -> float:
    """计算GC含量"""
    if len(seq) == 0:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in 'GC')
    return gc / len(seq)


def count_homopolymers(seq: str, min_length: int = 4) -> Dict[str, int]:
    """统计homopolymer"""
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    seq = seq.upper()
    i = 0
    while i < len(seq):
        base = seq[i]
        run_len = 1
        while i + run_len < len(seq) and seq[i + run_len] == base:
            run_len += 1
        if run_len >= min_length and base in counts:
            counts[base] += 1
        i += run_len
    return counts


def extract_circular_region(
    seq: str,
    offset: int,
    length: int,
    allow_wrap: bool = True
) -> str:
    """
    从环状序列提取区域
    
    Args:
        seq: 环状序列
        offset: 起始偏移（0-based）
        length: 提取长度
        allow_wrap: 是否允许绕环
    
    Returns:
        提取的序列
    """
    L = len(seq)
    if L == 0:
        return ""
    
    offset = offset % L
    
    if not allow_wrap:
        return seq[offset:offset + length]
    
    if offset + length <= L:
        return seq[offset:offset + length]
    
    # 需要绕环
    result = seq[offset:]
    remaining = length - len(result)
    full_copies = remaining // L
    result += seq * full_copies
    result += seq[:remaining % L]
    return result


def compute_repeat_count(read_length: int, template_length: int) -> int:
    """
    计算read覆盖的repeat次数
    
    Args:
        read_length: read长度
        template_length: 模板长度
    
    Returns:
        覆盖的完整repeat次数
    """
    if template_length == 0:
        return 0
    return (read_length + template_length - 1) // template_length


def find_tandem_repeats(seq: str, min_unit: int = 10, max_unit: int = 5000) -> List[dict]:
    """
    简单检测串联重复（用于验证）
    
    使用简单的k-mer方法，不是完整的TR检测
    
    Args:
        seq: 序列
        min_unit: 最小重复单元长度
        max_unit: 最大重复单元长度
    
    Returns:
        检测到的重复列表
    """
    results = []
    L = len(seq)
    
    for unit_len in range(min_unit, min(max_unit, L // 2) + 1):
        # 简单检查：比较首尾是否能形成周期
        unit = seq[:unit_len]
        matches = 0
        for i in range(0, L - unit_len, unit_len):
            if seq[i:i+unit_len] == unit:
                matches += 1
            else:
                break
        
        if matches >= 2:
            results.append({
                'unit_length': unit_len,
                'repeat_count': matches,
                'total_length': matches * unit_len
            })
    
    return results


def quality_string(length: int, mean_qual: int = 30, std_qual: float = 5.0, 
                   rng: Optional[np.random.Generator] = None) -> str:
    """
    生成质量字符串
    
    Args:
        length: 长度
        mean_qual: 平均质量值
        std_qual: 质量标准差
        rng: 随机数生成器
    
    Returns:
        FASTQ质量字符串
    """
    if rng is None:
        rng = np.random.default_rng()
    
    quals = rng.normal(mean_qual, std_qual, length)
    quals = np.clip(quals, 2, 41).astype(int)
    
    return "".join(chr(q + 33) for q in quals)


def simulate_quality_decay(
    length: int,
    start_qual: int = 35,
    end_qual: int = 20,
    noise_std: float = 3.0,
    rng: Optional[np.random.Generator] = None
) -> str:
    """
    模拟Illumina风格的质量衰减
    
    Args:
        length: 长度
        start_qual: 起始质量
        end_qual: 终止质量
        noise_std: 噪声标准差
        rng: 随机数生成器
    
    Returns:
        质量字符串
    """
    if rng is None:
        rng = np.random.default_rng()
    
    # 线性衰减 + 噪声
    base_quals = np.linspace(start_qual, end_qual, length)
    noise = rng.normal(0, noise_std, length)
    quals = base_quals + noise
    quals = np.clip(quals, 2, 41).astype(int)
    
    return "".join(chr(q + 33) for q in quals)
