"""
输入输出工具模块

- FASTA读写
- FASTQ写入
- 配置文件处理
"""

import logging
from typing import List, Dict, Tuple, Optional, TextIO, Generator, Union
from pathlib import Path
import gzip

from .models import EccDNA, SequencedRead

logger = logging.getLogger(__name__)

# 有效的DNA碱基
VALID_BASES = set('ACGTN')


def validate_sequence(seq: str, seq_id: str) -> str:
    """
    验证并清理序列

    Args:
        seq: 序列字符串
        seq_id: 序列ID（用于警告消息）

    Returns:
        清理后的序列（大写，去除空白）
    """
    seq = seq.upper().strip()

    # 检查非法字符
    invalid_chars = set(seq) - VALID_BASES
    if invalid_chars:
        logger.warning(
            f"Sequence '{seq_id}' contains non-standard bases: {invalid_chars}. "
            f"These will be converted to 'N'."
        )
        # 将非法字符替换为 N
        seq = ''.join(c if c in VALID_BASES else 'N' for c in seq)

    return seq


def parse_fasta(path: Union[str, Path]) -> List[EccDNA]:
    """
    解析FASTA文件
    
    支持.fa, .fasta, .fa.gz, .fasta.gz
    
    Args:
        path: FASTA文件路径
    
    Returns:
        EccDNA列表
    """
    path = Path(path)
    
    # 检查是否是gzip压缩
    opener = gzip.open if path.suffix == '.gz' else open
    mode = 'rt' if path.suffix == '.gz' else 'r'
    
    eccdnas = []
    current_id = None
    current_seq = []
    current_weight = 1.0
    
    with opener(path, mode) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # 保存上一条序列
                if current_id is not None:
                    seq = "".join(current_seq)
                    seq = validate_sequence(seq, current_id)
                    if seq:  # 只添加非空序列
                        eccdnas.append(EccDNA(
                            id=current_id,
                            seq=seq,
                            weight=current_weight
                        ))
                    else:
                        logger.warning(f"Skipping empty sequence: {current_id}")
                
                # 解析新的header
                header = line[1:].split()
                current_id = header[0]
                current_seq = []
                current_weight = 1.0
                
                # 尝试解析权重（格式：>id weight=1.5 或 >id 1.5）
                for part in header[1:]:
                    if part.startswith('weight='):
                        try:
                            current_weight = float(part.split('=')[1])
                        except (ValueError, IndexError):
                            pass
                    else:
                        try:
                            current_weight = float(part)
                        except ValueError:
                            pass
            else:
                current_seq.append(line.upper())
    
    # 保存最后一条
    if current_id is not None:
        seq = "".join(current_seq)
        seq = validate_sequence(seq, current_id)
        if seq:  # 只添加非空序列
            eccdnas.append(EccDNA(
                id=current_id,
                seq=seq,
                weight=current_weight
            ))
        else:
            logger.warning(f"Skipping empty sequence: {current_id}")

    if not eccdnas:
        logger.warning(f"No valid sequences found in {path}")

    return eccdnas


def write_fasta(eccdnas: List[EccDNA], path: Union[str, Path], line_width: int = 80):
    """
    写入FASTA文件
    
    Args:
        eccdnas: EccDNA列表
        path: 输出路径
        line_width: 每行宽度
    """
    path = Path(path)
    
    with open(path, 'w') as f:
        for ecc in eccdnas:
            f.write(f">{ecc.id} weight={ecc.weight}\n")
            seq = ecc.seq
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + "\n")


def write_fastq(
    reads: List[SequencedRead],
    path: Union[str, Path],
    compress: bool = False
):
    """
    写入FASTQ文件
    
    Args:
        reads: SequencedRead列表
        path: 输出路径
        compress: 是否gzip压缩
    """
    path = Path(path)
    if compress and not path.suffix == '.gz':
        path = Path(str(path) + '.gz')
    
    opener = gzip.open if compress else open
    mode = 'wt' if compress else 'w'
    
    with opener(path, mode) as f:
        for read in reads:
            f.write(read.to_fastq())


def write_paired_fastq(
    reads: List[SequencedRead],
    path_r1: Union[str, Path],
    path_r2: Union[str, Path],
    compress: bool = False
):
    """
    写入paired-end FASTQ文件
    
    Args:
        reads: SequencedRead列表（R1和R2交替或混合）
        path_r1: R1输出路径
        path_r2: R2输出路径
        compress: 是否gzip压缩
    """
    # 分离R1和R2
    r1_reads = [r for r in reads if r.read_number == 1]
    r2_reads = [r for r in reads if r.read_number == 2]
    
    # 按read_id排序确保配对
    r1_reads.sort(key=lambda r: r.read_id)
    r2_reads.sort(key=lambda r: r.read_id)
    
    write_fastq(r1_reads, path_r1, compress)
    write_fastq(r2_reads, path_r2, compress)


def iter_fastq(path: Union[str, Path]) -> Generator[Tuple[str, str, str], None, None]:
    """
    迭代读取FASTQ文件
    
    Yields:
        (read_id, sequence, quality)
    """
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    mode = 'rt' if path.suffix == '.gz' else 'r'
    
    with opener(path, mode) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            qual = f.readline().strip()
            
            read_id = header[1:].split()[0]  # 去掉@和后续注释
            yield read_id, seq, qual


def build_ecc_db(eccdnas: List[EccDNA]) -> Dict[str, EccDNA]:
    """
    构建eccDNA字典
    
    Args:
        eccdnas: EccDNA列表
    
    Returns:
        {id: EccDNA}字典
    """
    return {e.id: e for e in eccdnas}


def summarize_eccdnas(eccdnas: List[EccDNA]) -> str:
    """
    生成eccDNA摘要
    
    Args:
        eccdnas: EccDNA列表
    
    Returns:
        摘要字符串
    """
    if not eccdnas:
        return "No eccDNAs"
    
    lengths = [e.length for e in eccdnas]
    weights = [e.weight for e in eccdnas]
    
    return (
        f"eccDNAs: {len(eccdnas)}, "
        f"Length: {min(lengths)}-{max(lengths)}bp "
        f"(mean {sum(lengths)/len(lengths):.0f}bp), "
        f"Total weight: {sum(weights):.2f}"
    )
