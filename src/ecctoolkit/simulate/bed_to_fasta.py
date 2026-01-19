"""
BED 文件转 eccDNA FASTA

支持两种 BED 格式：
1. Simple eccDNA: chr start end [name] [score] [strand]
2. Chimeric eccDNA: 多行相同 name，会被合并为一个嵌合体

输出 FASTA 格式：
>eccDNA_1 chr1:100-500
ATGC...
"""

import logging
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class BedRecord:
    """BED 记录"""
    chrom: str
    start: int
    end: int
    name: str = ""
    score: float = 0.0
    strand: str = "+"

    @property
    def length(self) -> int:
        return self.end - self.start


def parse_bed_file(bed_path: str) -> List[BedRecord]:
    """解析 BED 文件"""
    records = []
    with open(bed_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 3:
                logger.warning(f"Line {line_num}: insufficient fields, skipping")
                continue

            try:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3] if len(fields) > 3 else f"ecc_{line_num}"
                score = float(fields[4]) if len(fields) > 4 else 0.0
                strand = fields[5] if len(fields) > 5 else "+"

                records.append(BedRecord(
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=name,
                    score=score,
                    strand=strand
                ))
            except ValueError as e:
                logger.warning(f"Line {line_num}: parse error ({e}), skipping")
                continue

    return records


def load_reference_index(fasta_path: str) -> Dict[str, Tuple[int, int, int]]:
    """
    加载参考基因组索引（简单实现，适合小基因组）
    返回 {chrom: (file_offset, seq_length, line_length)}
    """
    index = {}
    with open(fasta_path, 'r') as f:
        current_chrom = None
        seq_start = 0
        seq_length = 0
        line_length = 0

        while True:
            pos = f.tell()
            line = f.readline()
            if not line:
                break

            if line.startswith('>'):
                if current_chrom:
                    index[current_chrom] = (seq_start, seq_length, line_length)
                current_chrom = line[1:].split()[0]
                seq_start = f.tell()
                seq_length = 0
                line_length = 0
            else:
                seq_length += len(line.strip())
                if line_length == 0:
                    line_length = len(line)  # 包含换行符

        if current_chrom:
            index[current_chrom] = (seq_start, seq_length, line_length)

    return index


def fetch_sequence(
    fasta_path: str,
    chrom: str,
    start: int,
    end: int,
    strand: str = "+",
    index: Optional[Dict] = None
) -> str:
    """从参考基因组获取序列"""
    if index is None:
        index = load_reference_index(fasta_path)

    if chrom not in index:
        raise ValueError(f"Chromosome {chrom} not found in reference")

    file_offset, seq_length, line_length = index[chrom]

    if end > seq_length:
        raise ValueError(f"Region {chrom}:{start}-{end} exceeds chromosome length {seq_length}")

    # 计算文件位置
    bases_per_line = line_length - 1  # 减去换行符
    start_line = start // bases_per_line
    start_col = start % bases_per_line
    file_start = file_offset + start_line * line_length + start_col

    # 读取序列
    seq_parts = []
    remaining = end - start

    with open(fasta_path, 'r') as f:
        f.seek(file_start)
        while remaining > 0:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            take = min(len(line), remaining)
            seq_parts.append(line[:take])
            remaining -= take

    seq = ''.join(seq_parts)

    # 处理负链
    if strand == '-':
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                      'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
        seq = ''.join(complement.get(b, b) for b in reversed(seq))

    return seq.upper()


def bed_to_fasta(
    bed_path: str,
    reference_path: str,
    output_path: str,
    merge_chimeric: bool = True
) -> Tuple[str, int]:
    """
    将 BED 文件转换为 eccDNA FASTA

    Args:
        bed_path: BED 文件路径
        reference_path: 参考基因组 FASTA 路径
        output_path: 输出 FASTA 路径
        merge_chimeric: 是否合并相同 name 的记录为嵌合体

    Returns:
        (输出文件路径, eccDNA 数量)
    """
    records = parse_bed_file(bed_path)
    if not records:
        raise ValueError(f"No valid records found in {bed_path}")

    logger.info(f"Parsed {len(records)} BED records from {bed_path}")

    # 加载参考基因组索引
    index = load_reference_index(reference_path)
    logger.info(f"Loaded reference index: {len(index)} chromosomes")

    # 按 name 分组（用于嵌合体）
    if merge_chimeric:
        grouped: Dict[str, List[BedRecord]] = {}
        for rec in records:
            if rec.name not in grouped:
                grouped[rec.name] = []
            grouped[rec.name].append(rec)
    else:
        grouped = {rec.name: [rec] for rec in records}

    # 生成 FASTA
    output_path = Path(output_path)
    ecc_count = 0

    with open(output_path, 'w') as f:
        for name, recs in grouped.items():
            try:
                # 获取序列
                seq_parts = []
                regions = []
                for rec in recs:
                    seq = fetch_sequence(
                        reference_path, rec.chrom, rec.start, rec.end,
                        rec.strand, index
                    )
                    seq_parts.append(seq)
                    regions.append(f"{rec.chrom}:{rec.start}-{rec.end}({rec.strand})")

                full_seq = ''.join(seq_parts)
                region_str = ','.join(regions)

                # 写入 FASTA
                # 使用 weight 字段（如果 BED 有 score）
                weight = recs[0].score if recs[0].score > 0 else 1.0
                header = f">{name} region={region_str} length={len(full_seq)} weight={weight}"
                f.write(f"{header}\n")

                # 每行 80 字符
                for i in range(0, len(full_seq), 80):
                    f.write(full_seq[i:i+80] + '\n')

                ecc_count += 1

            except Exception as e:
                logger.warning(f"Failed to process {name}: {e}")
                continue

    logger.info(f"Generated {ecc_count} eccDNA sequences to {output_path}")
    return str(output_path), ecc_count


def detect_input_type(input_path: str) -> str:
    """检测输入文件类型"""
    path = Path(input_path)
    suffix = path.suffix.lower()

    if suffix in ['.fa', '.fasta', '.fna']:
        return 'fasta'
    elif suffix in ['.bed']:
        return 'bed'
    else:
        # 尝试读取第一行判断
        with open(input_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('>'):
                return 'fasta'
            else:
                return 'bed'
