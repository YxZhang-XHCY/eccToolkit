import os
import sys
import glob
import shutil
import logging
from typing import Optional, Union, Callable, Dict, Any

import pandas as pd
import numpy as np
import subprocess as sp
import multiprocessing as mp
from importlib.resources import files
from os.path import exists, join, split
from tqdm import tqdm

from .lite_utils import utilities

logger = logging.getLogger(__name__)

# 内置 HiFi profile ID（基于 2000 条真实 HiFi reads）
BUILTIN_HIFI_PROFILE_ID = "HeLa_HiFi_2k"


def _check_tool_available(tool_name: str) -> bool:
    """检测外部工具是否可用"""
    return shutil.which(tool_name) is not None


def _generate_molecule_seed(base_seed: int, mol_id: str, platform: str = "") -> int:
    """
    为分子生成唯一的随机种子

    使用哈希函数确保：
    1. 相同的 base_seed + mol_id + platform 总是产生相同的种子
    2. 不同的组合产生不同的种子，避免碰撞
    3. 种子在有效范围内 [0, 2^31-1]

    Args:
        base_seed: 基础种子（用户指定）
        mol_id: 分子 ID（唯一标识符）
        platform: 平台标识（可选，用于区分同一分子在不同平台的种子）

    Returns:
        32位有符号整数范围内的种子
    """
    import hashlib
    # 组合所有信息生成唯一字符串
    seed_str = f"{base_seed}:{mol_id}:{platform}"
    # 使用 MD5 哈希（足够快且分布均匀）
    hash_bytes = hashlib.md5(seed_str.encode('utf-8')).digest()
    # 取前4字节作为32位整数
    seed_value = int.from_bytes(hash_bytes[:4], byteorder='little')
    # 确保在有效范围内 [0, 2^31-1]
    return seed_value % (2**31)


# =============================================================================
# Read ID 解析函数
# =============================================================================

def _clean_read_id(read_id: str) -> str:
    """
    清理 read ID：去除空格和后续注释。

    Args:
        read_id: 原始 read ID (可能含 @ 前缀、空格、注释)

    Returns:
        清理后的 read ID
    """
    return read_id.lstrip('@').split()[0]


def _parse_int_or_str(value: str):
    """
    尝试将字符串转换为整数，失败则返回原字符串。
    """
    return int(value) if value.isdigit() else value


def parse_art_read_id(read_id: str) -> dict:
    """
    解析 ART (art_illumina) 生成的 read ID

    ART 生成的 ID 格式: {mol_id}-{read_index}/{read_number}
    例如: mol_123-45/1, mol_123-45/2

    Args:
        read_id: FASTQ 中的 read ID (不含 @ 前缀)

    Returns:
        dict with keys: mol_id, read_index, read_number (1 or 2)
    """
    read_id = _clean_read_id(read_id)
    result = {'mol_id': None, 'read_index': None, 'read_number': 1}

    # 提取 read_number (如 /1 或 /2)
    if '/' in read_id:
        base, read_num = read_id.rsplit('/', 1)
        result['read_number'] = _parse_int_or_str(read_num) if read_num else 1
    else:
        base = read_id

    # 提取 mol_id 和 read_index
    if '-' in base:
        parts = base.rsplit('-', 1)
        result['mol_id'] = parts[0]
        result['read_index'] = _parse_int_or_str(parts[1]) if parts[1] else None
    else:
        result['mol_id'] = base

    return result


def parse_pbsim_read_id(read_id: str) -> dict:
    """
    解析 PBSIM2 生成的 read ID

    PBSIM2 生成的 ID 格式: {mol_id}.{seq_num}_{read_num}
    例如: ecc_1.1_1, ecc_1.1_2, mol_123.1_5

    旧格式兼容: {mol_id}_{read_num} (无 . 分隔符)

    Args:
        read_id: FASTQ 中的 read ID (不含 @ 前缀)

    Returns:
        dict with keys: mol_id, read_index
    """
    read_id = _clean_read_id(read_id)
    result = {'mol_id': None, 'read_index': None}

    # 新格式: {mol_id}.{seq_num}_{read_num}
    if '.' in read_id and '_' in read_id:
        base, read_num = read_id.rsplit('_', 1)
        result['read_index'] = _parse_int_or_str(read_num)
        # 去掉 seq_num 部分
        result['mol_id'] = base.rsplit('.', 1)[0] if '.' in base else base
    elif '_' in read_id:
        # 旧格式兼容: {mol_id}_{read_num}
        parts = read_id.rsplit('_', 1)
        result['mol_id'] = parts[0]
        result['read_index'] = _parse_int_or_str(parts[1]) if parts[1] else None
    else:
        result['mol_id'] = read_id

    return result


def parse_hifi_simple_read_id(read_id: str) -> dict:
    """
    解析 HiFi simple 模式生成的 read ID

    HiFi simple 模式 ID 格式: {mol_id}_hifi_{read_index}
    例如: ecc_1_hifi_0, mol_123_hifi_15

    Args:
        read_id: FASTQ 中的 read ID (不含 @ 前缀)

    Returns:
        dict with keys: mol_id, read_index
    """
    read_id = _clean_read_id(read_id)
    result = {'mol_id': None, 'read_index': None}

    # HiFi simple 格式: {mol_id}_hifi_{index}
    if '_hifi_' in read_id:
        parts = read_id.split('_hifi_', 1)
        result['mol_id'] = parts[0]
        result['read_index'] = _parse_int_or_str(parts[1]) if len(parts) > 1 and parts[1] else None
    else:
        result['mol_id'] = read_id

    return result


def generate_truth_from_fastq(
    fastq_path: str,
    pool_csv_path: str,
    output_truth_path: str,
    platform: str,
    read_id_parser: str = 'auto'
) -> int:
    """
    从 FASTQ 文件和 pool CSV 生成 truth 文件

    Args:
        fastq_path: FASTQ 文件路径
        pool_csv_path: 分子池 CSV 文件路径
        output_truth_path: 输出 truth TSV 文件路径
        platform: 平台名称 (NGS, HiFi, ONT)
        read_id_parser: read ID 解析器 ('art', 'pbsim', 'hifi_simple', 'auto')

    Returns:
        写入的 truth 记录数
    """
    import csv

    if not exists(fastq_path):
        logger.warning(f"FASTQ file not found: {fastq_path}")
        return 0

    if not exists(pool_csv_path):
        logger.warning(f"Pool CSV not found: {pool_csv_path}")
        return 0

    # 读取 pool CSV 并建立分子 ID 索引
    pool_df = pd.read_csv(pool_csv_path, sep='\t')
    mol_index = {}
    for idx, row in pool_df.iterrows():
        mol_id = str(row['id'])
        mol_index[mol_id] = row.to_dict()

    # 选择 read ID 解析器
    if read_id_parser == 'auto':
        if platform.upper() == 'NGS':
            parser_func = parse_art_read_id
        elif platform.upper() == 'HIFI':
            parser_func = parse_hifi_simple_read_id  # 默认使用 simple 模式解析
        else:  # ONT
            parser_func = parse_pbsim_read_id
    elif read_id_parser == 'art':
        parser_func = parse_art_read_id
    elif read_id_parser == 'pbsim':
        parser_func = parse_pbsim_read_id
    elif read_id_parser == 'hifi_simple':
        parser_func = parse_hifi_simple_read_id
    else:
        parser_func = parse_pbsim_read_id

    # 定义 truth 输出列
    truth_columns = [
        'read_id',
        'platform',
        'source_mol_id',
        'ecc_ids',
        'repeat_count',
        'has_chimera',
        'is_background',
        'source_ecc_length',
        'background_region'
    ]

    def _clean_value(value):
        """Normalize NaN/None to empty string for TSV output."""
        if value is None:
            return ''
        try:
            if pd.isna(value):
                return ''
        except Exception:
            pass
        return value

    count = 0
    with open(output_truth_path, 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.DictWriter(out_file, fieldnames=truth_columns, delimiter='\t')
        writer.writeheader()

        # 读取 FASTQ 并解析 read ID
        with open(fastq_path, 'r') as fq:
            line_num = 0
            for line in fq:
                line_num += 1
                if line_num % 4 == 1:  # header 行
                    read_id = line.strip().lstrip('@')
                    parsed = parser_func(read_id)
                    mol_id = parsed['mol_id']

                    # 在 pool 中查找分子信息
                    mol_info = mol_index.get(mol_id, {})

                    record = {
                        'read_id': read_id,
                        'platform': platform,
                        'source_mol_id': mol_id,
                        'ecc_ids': _clean_value(mol_info.get('ecc_ids', '')),
                        'repeat_count': _clean_value(mol_info.get('repeat_count', '')),
                        'has_chimera': _clean_value(mol_info.get('has_chimera', '')),
                        'is_background': _clean_value(mol_info.get('is_background', '')),
                        'source_ecc_length': _clean_value(mol_info.get('source_ecc_length', '')),
                        'background_region': _clean_value(mol_info.get('background_region', ''))
                    }
                    writer.writerow(record)
                    count += 1

    logger.info(f"Generated truth file: {output_truth_path} ({count} records)")
    return count


class seqsim:
    """
    序列模拟类 - 生成 eccDNA 和线性 DNA 序列

    从参考基因组随机采样区域生成:
    - eccDNA (circular): 环形 DNA，存储在 pos.bed/csv
    - 线性 DNA (linear): 背景线性 DNA，存储在 neg.bed/csv

    Args:
        sample: 样本名称
        reference: 参考基因组 FASTA 路径
        path: 输出目录路径
        circular_number: eccDNA 数量 (默认 5000)
        linear_number: 线性 DNA 数量 (默认 5000)
        seed: 随机种子
        simple_ratio: 简单 eccDNA 比例 (暂不使用)
        simple_template: 简单 eccDNA 模板 BED (暂不使用)
        chimeric_template: 嵌合 eccDNA 模板 BED (暂不使用)
    """

    def __init__(
        self,
        sample: str,
        reference: str,
        path: str,
        circular_number: int = 5000,
        linear_number: int = 5000,
        seed: int = None,
        simple_ratio: float = None,
        simple_template: str = None,
        chimeric_template: str = None,
    ):
        import os
        import numpy as np

        self.sample = sample
        self.reference = reference
        self.path = path
        self.circular_number = int(circular_number)
        self.linear_number = int(linear_number)
        self.seed = seed
        self.simple_ratio = simple_ratio
        self.simple_template = simple_template
        self.chimeric_template = chimeric_template

        # 创建输出目录
        self.sample_dir = os.path.join(path, sample)
        os.makedirs(self.sample_dir, exist_ok=True)

        # 初始化工具类
        self.utils = utilities(reference)

        # 设置随机种子
        if seed is not None:
            np.random.seed(seed)

        # 运行序列生成
        self._generate_sequences()

    def _sample_length(self, rng) -> int:
        """
        采样 eccDNA/线性 DNA 长度

        使用与 sim-region 类似的长度分布:
        - 主体: lognormal 分布，峰值约 400bp
        - 尾部: 均匀分布，5kb-500kb
        """
        # 90% 使用 lognormal 分布 (峰值 ~400bp)
        if rng.random() < 0.95:
            # lognormal: mode=400, sigma=0.8
            length = int(rng.lognormal(mean=np.log(400), sigma=0.8))
            length = max(100, min(length, 50000))  # 限制在 100bp - 50kb
        else:
            # 5% 使用尾部分布 (5kb - 500kb)
            length = int(rng.integers(5000, 500000))

        return length

    def _generate_sequences(self):
        """生成 eccDNA 和线性 DNA 序列"""
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(self.seed)

        # 获取参考基因组信息
        genome_info = self.utils.genome_length()
        chroms = genome_info.index.tolist()
        chrom_lengths = genome_info['length'].values.astype(float)

        # 按染色体长度加权
        chrom_weights = chrom_lengths / chrom_lengths.sum()

        # ===== 生成 eccDNA (pos) =====
        logger.info(f"生成 {self.circular_number} 个 eccDNA...")
        pos_bed_rows = []
        pos_csv_rows = []

        for i in tqdm(range(self.circular_number), desc="生成 eccDNA"):
            length = self._sample_length(rng)
            ecc_id = f"eccDNA_{i+1:06d}"

            # 随机选择染色体
            chrom_idx = rng.choice(len(chroms), p=chrom_weights)
            chrom = chroms[chrom_idx]

            # 确保长度不超过染色体
            max_len = int(chrom_lengths[chrom_idx]) - 1000
            if length > max_len:
                length = max(100, max_len)

            try:
                region = self.utils.random_region(chrom, length)
                chrom, start, end, frag_len, seq = region
            except Exception as e:
                # 失败时使用较短长度重试
                logger.debug(f"采样失败 ({chrom}, {length}): {e}")
                region = self.utils.random_region(chroms[0], 500)
                chrom, start, end, frag_len, seq = region

            pos_bed_rows.append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "length": int(frag_len),
                "id": ecc_id,
            })
            pos_csv_rows.append({
                "id": ecc_id,
                "fragN": 1,
                "region": f"{chrom}:{start}-{end}",
                "length": int(frag_len),
                "seq": seq,
            })

        # ===== 生成线性 DNA (neg) =====
        logger.info(f"生成 {self.linear_number} 个线性 DNA...")
        neg_bed_rows = []
        neg_csv_rows = []

        for i in tqdm(range(self.linear_number), desc="生成线性 DNA"):
            length = self._sample_length(rng)
            linear_id = f"neg_{i+1:06d}"

            # 随机选择染色体
            chrom_idx = rng.choice(len(chroms), p=chrom_weights)
            chrom = chroms[chrom_idx]

            # 确保长度不超过染色体
            max_len = int(chrom_lengths[chrom_idx]) - 1000
            if length > max_len:
                length = max(100, max_len)

            try:
                region = self.utils.random_region(chrom, length)
                chrom, start, end, frag_len, seq = region
            except Exception as e:
                logger.debug(f"采样失败 ({chrom}, {length}): {e}")
                region = self.utils.random_region(chroms[0], 500)
                chrom, start, end, frag_len, seq = region

            neg_bed_rows.append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "length": int(frag_len),
                "id": linear_id,
            })
            neg_csv_rows.append({
                "id": linear_id,
                "fragN": 1,
                "region": f"{chrom}:{start}-{end}",
                "length": int(frag_len),
                "seq": seq,
            })

        # ===== 保存文件 =====
        pos_bed = pd.DataFrame(pos_bed_rows, columns=["chrom", "start", "end", "length", "id"])
        pos_csv = pd.DataFrame(pos_csv_rows, columns=["id", "fragN", "region", "length", "seq"])
        neg_bed = pd.DataFrame(neg_bed_rows, columns=["chrom", "start", "end", "length", "id"])
        neg_csv = pd.DataFrame(neg_csv_rows, columns=["id", "fragN", "region", "length", "seq"])

        pos_bed_path = os.path.join(self.sample_dir, f"{self.sample}.pos.bed")
        pos_csv_path = os.path.join(self.sample_dir, f"{self.sample}.pos.csv")
        neg_bed_path = os.path.join(self.sample_dir, f"{self.sample}.neg.bed")
        neg_csv_path = os.path.join(self.sample_dir, f"{self.sample}.neg.csv")

        pos_bed.to_csv(pos_bed_path, sep="\t", header=False, index=False)
        pos_csv.to_csv(pos_csv_path, sep="\t", index=False)
        neg_bed.to_csv(neg_bed_path, sep="\t", header=False, index=False)
        neg_csv.to_csv(neg_csv_path, sep="\t", index=False)

        logger.info(f"序列生成完成:")
        logger.info(f"  eccDNA: {len(pos_csv_rows)} 个 -> {pos_csv_path}")
        logger.info(f"  线性 DNA: {len(neg_csv_rows)} 个 -> {neg_csv_path}")


def _rca_worker(args):
    """
    RCA 扩增的 worker 函数（用于多进程）

    Args:
        args: (ecc_id, ecc_data, min_rca_length, min_repeats, worker_seed)

    Returns:
        (table_row, bed_rows): 处理结果
    """
    ecc_id, ecc_data, min_rca_length, min_repeats, worker_seed = args

    # 设置此 worker 的随机种子
    rng = np.random.default_rng(worker_seed)

    # 转换为 DataFrame
    region_table = pd.DataFrame(ecc_data)

    # === sim_breakpoint 逻辑 ===
    region_table = region_table.reset_index(drop=True)
    total_len = region_table.length.sum()

    if total_len > 0:
        breakpoint = int(rng.integers(0, total_len))
        cumsum = 0
        for i in range(len(region_table)):
            cumsum += region_table.loc[i, 'length']
            if cumsum > breakpoint:
                offset = breakpoint - (cumsum - region_table.loc[i, 'length'])
                current_row = region_table.loc[i]

                part1_len = offset
                part2_len = current_row['length'] - offset

                if part1_len > 0 and part2_len > 0:
                    part1 = {
                        'chrom': current_row['chrom'],
                        'start': current_row['start'],
                        'end': int(current_row['start']) + part1_len,
                        'length': part1_len,
                        'seq': current_row['seq'][:part1_len],
                        'id': current_row['id']
                    }
                    part2 = {
                        'chrom': current_row['chrom'],
                        'start': int(current_row['start']) + part1_len,
                        'end': current_row['end'],
                        'length': part2_len,
                        'seq': current_row['seq'][part1_len:],
                        'id': current_row['id']
                    }

                    after_current = region_table.iloc[i+1:].copy() if i+1 < len(region_table) else pd.DataFrame()
                    before_current = region_table.iloc[:i].copy() if i > 0 else pd.DataFrame()

                    parts = [pd.DataFrame([part1])]
                    if not after_current.empty:
                        parts.append(after_current)
                    if not before_current.empty:
                        parts.append(before_current)
                    parts.append(pd.DataFrame([part2]))

                    region_table = pd.concat(parts, ignore_index=True)
                break

    # === sim_RCA 逻辑 ===
    region_table = region_table.reset_index(drop=True)
    original_length = region_table.length.sum()

    target_by_length = min_rca_length
    target_by_repeats = original_length * min_repeats
    target_amp = max(target_by_length, target_by_repeats)

    roundN = int(target_amp / original_length) if original_length > 0 else 1
    residual = target_amp - roundN * original_length if original_length > 0 else 0

    amplified_parts = []
    for _ in range(roundN):
        amplified_parts.append(region_table.copy())

    if residual > 0:
        cumsum = 0
        for i in range(len(region_table)):
            cumsum += region_table.loc[i, 'length']
            if cumsum >= residual:
                residual_table = region_table.iloc[:i+1].copy()
                frag_residual = residual - (cumsum - region_table.loc[i, 'length'])
                residual_table.iloc[-1, residual_table.columns.get_loc('length')] = frag_residual
                residual_table.iloc[-1, residual_table.columns.get_loc('end')] = int(residual_table.iloc[-1]['start']) + frag_residual
                residual_table.iloc[-1, residual_table.columns.get_loc('seq')] = residual_table.iloc[-1]['seq'][:frag_residual]
                amplified_parts.append(residual_table)
                break

    if amplified_parts:
        amplified_table = pd.concat(amplified_parts, ignore_index=True)
    else:
        amplified_table = region_table.copy()

    # 构建结果
    amplified_table['region'] = amplified_table['chrom'] + ':' + amplified_table['start'].astype(str) + '-' + amplified_table['end'].astype(str)
    rca_len = amplified_table['length'].sum()

    table_row = {
        'id': ecc_id,
        'fragN': len(amplified_table),
        'region': '|'.join(amplified_table['region']),
        'length': rca_len,
        'seq': ''.join(amplified_table['seq']),
        'original_length': original_length
    }

    bed_rows = amplified_table[['chrom', 'start', 'end', 'length', 'id']].to_dict('records')

    return table_row, bed_rows


class libsim:
    """
    分子库模拟类

    Args:
        sample: 样本名称
        reference: 参考基因组 FASTA 路径
        path: 输出目录路径
        seed: 随机种子（用于结果复现）
        meancov: 平均覆盖度
        amp: 最小 RCA 扩增长度（默认 50000 bp），作为 min_rca_length 的别名
        min_rca_length: 最小 RCA 扩增长度
        min_repeats: 每个 eccDNA 的最少重复次数（默认 5）
        threads: 线程数（用于并行 RCA 扩增，默认 1）
    """
    def __init__(
        self,
        sample: str,
        reference: str,
        path: str,
        seed: Optional[int] = None,
        meancov: float = 25,
        amp: int = 50000,
        min_rca_length: Optional[int] = None,
        min_repeats: int = 5,
        threads: int = 1,
    ) -> None:

        self.sample = sample

        self.seed = seed
        self.meancov = meancov
        # min_rca_length 优先，否则使用 amp
        self.min_rca_length = min_rca_length if min_rca_length is not None else amp
        self.min_repeats = max(1, int(min_repeats))
        self.threads = max(1, int(threads))

        self.path = path  # 直接使用 path，不创建 sample 子目录
        self.prefix = [join(self.path, self.sample + '.neg'), join(self.path, self.sample + '.pos'), join(self.path,  self.sample + '.lib')]
        self.utils = utilities(reference)

        self.sim_library()
             
    def sim_breakpoint(self, region_table: pd.DataFrame) -> pd.DataFrame:
        """
        Simulate random breakpoint on circular DNA and return shifted region table.

        对于环形 DNA，随机选择一个断裂点位置，将环从该位置断开，
        重新排列片段顺序（断裂点后的部分移到前面）。

        Args:
            region_table: 包含 chrom, start, end, length, seq, id 列的 DataFrame

        Returns:
            重新排列后的 region_table
        """
        region_table = region_table.reset_index(drop=True)
        total_length = region_table['length'].sum()

        # 边界检查：如果只有一个碱基或更少，直接返回
        if total_length <= 1:
            return region_table.copy()

        shift = np.random.randint(1, total_length)

        # 特殊处理：单行 region 的情况
        if len(region_table) == 1:
            row = region_table.iloc[0]
            if shift >= row['length']:
                # 断裂点在末尾或超出，返回原表
                return region_table.copy()

            # 在单个片段内部断裂，分成两部分
            part1 = {
                'chrom': row['chrom'],
                'start': row['start'] + shift,
                'end': row['end'],
                'length': row['length'] - shift,
                'seq': row['seq'][shift:],
                'id': row['id']
            }
            part2 = {
                'chrom': row['chrom'],
                'start': row['start'],
                'end': row['start'] + shift,
                'length': shift,
                'seq': row['seq'][:shift],
                'id': row['id']
            }
            return pd.DataFrame([part1, part2])

        # 多行 region 的处理
        shifted_table = pd.DataFrame(columns=['chrom', 'start', 'end', 'length', 'seq', 'id'])
        cumulative_length = 0

        for i in range(len(region_table)):
            cumulative_length += region_table.loc[i, 'length']

            if shift == cumulative_length:
                # 断裂点正好在第 i 行末尾，无需分割该行
                after_break = region_table.iloc[i+1:].copy() if i+1 < len(region_table) else pd.DataFrame()
                before_break = region_table.iloc[:i+1].copy()
                shifted_table = pd.concat([after_break, before_break], ignore_index=True)
                break

            elif shift < cumulative_length:
                # 断裂点在第 i 行内部
                prev_cumulative = cumulative_length - region_table.loc[i, 'length']
                shift_in_row = shift - prev_cumulative

                # 当前行分割成两部分
                current_row = region_table.loc[i]
                part1 = {
                    'chrom': current_row['chrom'],
                    'start': current_row['start'] + shift_in_row,
                    'end': current_row['end'],
                    'length': current_row['length'] - shift_in_row,
                    'seq': current_row['seq'][shift_in_row:],
                    'id': current_row['id']
                }
                part2 = {
                    'chrom': current_row['chrom'],
                    'start': current_row['start'],
                    'end': current_row['start'] + shift_in_row,
                    'length': shift_in_row,
                    'seq': current_row['seq'][:shift_in_row],
                    'id': current_row['id']
                }

                # 组装结果：[part1] + [i+1 到末尾] + [0 到 i-1] + [part2]
                after_current = region_table.iloc[i+1:].copy() if i+1 < len(region_table) else pd.DataFrame()
                before_current = region_table.iloc[:i].copy() if i > 0 else pd.DataFrame()

                parts = [pd.DataFrame([part1])]
                if not after_current.empty:
                    parts.append(after_current)
                if not before_current.empty:
                    parts.append(before_current)
                parts.append(pd.DataFrame([part2]))

                shifted_table = pd.concat(parts, ignore_index=True)
                break

        return shifted_table
    
    def sim_RCA(self, region_table: pd.DataFrame) -> pd.DataFrame:
        """
        Simulate RCA (Rolling Circle Amplification).

        扩增逻辑：
        - 计算实际扩增目标长度 = max(min_rca_length, original_length * min_repeats)
        - 这确保每个 eccDNA:
          1. 扩增后长度至少为 min_rca_length (默认 50kb)
          2. 至少被完整重复 min_repeats 次 (默认 5 次)

        Args:
            region_table: 包含 eccDNA 片段信息的 DataFrame

        Returns:
            RCA 扩增后的 DataFrame
        """
        region_table = region_table.reset_index().drop('index', axis=1)
        original_length = region_table.length.sum()

        # 计算实际扩增目标：取两个条件中较大的值
        target_by_length = self.min_rca_length
        target_by_repeats = original_length * self.min_repeats
        target_amp = max(target_by_length, target_by_repeats)

        # 计算完整重复次数和残余部分
        roundN = int(target_amp / original_length)
        residual = target_amp - roundN * original_length

        # 构建扩增后的表
        amplified_table = pd.DataFrame(columns=['chrom', 'start', 'end', 'length', 'seq', 'id'])

        # 添加完整的重复
        for _ in range(roundN):
            amplified_table = pd.concat([amplified_table, region_table], axis=0)

        # 添加残余部分（如果有）
        if residual > 0:
            for i in range(1, len(region_table) + 1):
                if region_table[0:i].length.sum() >= residual:
                    residual_table = region_table[0:i].copy()
                    # 调整最后一个 fragment 的长度
                    last_idx = i - 1
                    frag_residual = residual - region_table[0:i-1].length.sum()
                    residual_table.loc[last_idx, 'length'] = frag_residual
                    residual_table.loc[last_idx, 'end'] = int(residual_table.loc[last_idx, 'start']) + frag_residual
                    residual_table.loc[last_idx, 'seq'] = residual_table.loc[last_idx, 'seq'][:frag_residual]
                    amplified_table = pd.concat([amplified_table, residual_table], axis=0)
                    break

        amplified_table = amplified_table.reset_index(drop=True)
        return amplified_table
    
    def sim_library(self) -> None:
        """
        模拟分子库：RCA 扩增 + 覆盖度计算

        覆盖度计算逻辑：
        ====================

        设定：
        - L_orig: 原始 eccDNA 长度
        - L_rca:  RCA 扩增后长度 (L_rca = L_orig * repeat_count)
        - realcov: 用户期望的目标覆盖度（相对于原始 eccDNA）
        - tempcov: 传给 ART/PBSIM 的覆盖度参数

        问题：
        - ART/PBSIM 会按照 tempcov 对 RCA 产物 (L_rca) 进行测序
        - 产生的 reads 数量 ∝ tempcov * L_rca
        - 但这些 reads 映射回原始 eccDNA 时，覆盖度会被放大 repeat_count 倍

        解决方案：
        - 设置 tempcov = realcov * L_orig / L_rca
        - 则实际产生的 reads 数量 ∝ tempcov * L_rca = realcov * L_orig
        - 映射回原始 eccDNA 后，覆盖度 = realcov（符合预期）

        RCA 扩增保证：
        - 每个 eccDNA 扩增后长度 >= min_rca_length (默认 50kb)
        - 每个 eccDNA 至少重复 min_repeats 次 (默认 5 次)
        """
        if self.seed is not None:
            np.random.seed(self.seed)

        # 读取 positive library (eccDNA)
        pos_bed = pd.read_csv(self.prefix[1] + '.bed', sep='\t', names=['chrom', 'start', 'end', 'length', 'id'])

        # 获取序列（显示进度）
        logger.info(f"获取 {len(pos_bed)} 个区域的序列...")
        tqdm.pandas(desc="获取序列")
        pos_bed.insert(pos_bed.shape[1], 'seq', pos_bed.progress_apply(lambda x: self.utils.get_seq(x.chrom, x.start, x.end), axis=1))

        # 对每个 eccDNA 进行 RCA 扩增
        unique_ids = pos_bed.id.unique()
        logger.info(f"对 {len(unique_ids)} 个 eccDNA 进行 RCA 扩增 (threads={self.threads})...")

        # 准备并行任务参数
        tasks = []
        base_seed = self.seed if self.seed is not None else 42
        for idx, ecc_id in enumerate(unique_ids):
            ecc_data = pos_bed[pos_bed.id == ecc_id].to_dict('records')
            worker_seed = base_seed + idx  # 每个 eccDNA 使用不同但可复现的种子
            tasks.append((ecc_id, ecc_data, self.min_rca_length, self.min_repeats, worker_seed))

        # 并行或串行处理
        table_rows = []
        bed_rows_all = []

        if self.threads > 1:
            # 多进程并行处理
            from concurrent.futures import ProcessPoolExecutor, as_completed

            with ProcessPoolExecutor(max_workers=self.threads) as executor:
                futures = {executor.submit(_rca_worker, task): task[0] for task in tasks}
                for future in tqdm(as_completed(futures), total=len(futures), desc="RCA 扩增"):
                    table_row, bed_rows = future.result()
                    table_rows.append(table_row)
                    bed_rows_all.extend(bed_rows)
        else:
            # 单线程串行处理
            for task in tqdm(tasks, desc="RCA 扩增"):
                table_row, bed_rows = _rca_worker(task)
                table_rows.append(table_row)
                bed_rows_all.extend(bed_rows)

        # 构建输出 DataFrame
        output_table = pd.DataFrame(table_rows, columns=['id', 'fragN', 'region', 'length', 'seq', 'original_length'])
        output_bed = pd.DataFrame(bed_rows_all, columns=['chrom', 'start', 'end', 'length', 'id'])
        ## merge negative library & positive library
        neg_table = pd.read_csv(self.prefix[0] + '.csv', sep='\t')
        # 为 neg_table 添加 original_length 列（线性 DNA 原始长度等于自身长度）
        if 'original_length' not in neg_table.columns:
            neg_table['original_length'] = neg_table['length'] if 'length' in neg_table.columns else 0
        output_table = pd.concat([output_table, neg_table])

        # 覆盖度计算修正
        # realcov: 目标覆盖度（相对于原始 eccDNA）
        output_table.insert(output_table.shape[1], 'realcov', np.random.gamma(2.5, self.meancov / 2.5, len(output_table)))

        # tempcov: ART 使用的覆盖度
        # 公式: tempcov = realcov * original_length / rca_length
        # 这确保最终覆盖度 = tempcov * rca_length / original_length = realcov
        output_table['tempcov'] = output_table.apply(
            lambda row: row['realcov'] * row['original_length'] / row['length']
                        if row['length'] > 0 else row['realcov'],
            axis=1
        )

        # 添加 truth 追踪所需的字段
        # ecc_ids: 等于分子 ID（Lite 模式中每个分子对应一个 eccDNA）
        output_table['ecc_ids'] = output_table['id']
        # repeat_count: RCA 扩增倍数（约等于 length / original_length）
        output_table['repeat_count'] = output_table.apply(
            lambda row: int(row['length'] / row['original_length'])
                        if row['original_length'] > 0 else 1,
            axis=1
        )
        # 其他 truth 字段
        output_table['has_chimera'] = False  # Lite 模式不模拟 chimera
        output_table['is_background'] = output_table['id'].str.startswith('neg_')
        output_table['source_ecc_length'] = output_table['original_length']
        output_table['background_region'] = ''

        output_table.to_csv(self.prefix[2] + '.csv', index=None, sep='\t')
        self.utils.write_fasta(output_table, self.prefix[2] + '.fasta')
        output_bed.to_csv(self.prefix[2] + '.bed', header=None, index=None, sep='\t')
        self.utils.transfer_files(self.prefix[0] + '.bed', self.prefix[2] + '.bed')
        return


# =============================================================================
# ART 并行化 Worker 函数
# =============================================================================

def _run_art_batch_worker(args: tuple) -> tuple:
    """
    独立的 ART 批次执行 worker 函数（用于多进程并行）。

    Args:
        args: 元组包含 (batch_idx, fasta_path, output_prefix, cov,
                       sr_platform, sr_readlen, sr_mean, sr_std, batch_seed,
                       cwd, art_env)

    Returns:
        (batch_idx, success, r1_path, r2_path, error_msg)
    """
    (batch_idx, fasta_path, output_prefix, cov,
     sr_platform, sr_readlen, sr_mean, sr_std, batch_seed,
     cwd, art_env) = args

    r1_path = f'{output_prefix}1.fq'
    r2_path = f'{output_prefix}2.fq'

    try:
        cmd = [
            'art_illumina',
            '-na',
            '-q',
            '-nf', '0',
            '-p',
            '-i', str(fasta_path),
            '-o', str(output_prefix),
            '-ss', str(sr_platform),
            '-l', str(sr_readlen),
            '-f', str(cov),
            '-m', str(sr_mean),
            '-s', str(sr_std),
        ]
        if batch_seed is not None:
            cmd.extend(['-rs', str(batch_seed)])

        sp.run(cmd, check=True, cwd=cwd, capture_output=True, env=art_env)

        # 检查输出文件是否存在
        if not os.path.exists(r1_path) or not os.path.exists(r2_path):
            return (batch_idx, False, None, None, 'ART output files not found')

        return (batch_idx, True, r1_path, r2_path, None)
    except Exception as e:
        return (batch_idx, False, None, None, str(e))


class fqsim:
    """
    从 eccDNA 分子库生成 FASTQ reads

    Args:
        sample: 样本名称
        csv: libsim 生成的模板 CSV 文件路径
        path: 输出目录路径
        seed: 随机种子
        thread: 最大线程数
        skip_sr: 跳过短读 (NGS) 模拟
        skip_hifi: 跳过 PacBio HiFi 模拟
        skip_ont: 跳过 ONT 模拟
        ont_mean: ONT 读长均值
        ont_std: ONT 读长标准差
        ont_model: PBSIM2 ONT 模型
        hifi_sample_fastq: HiFi 采样 FASTQ 文件
        hifi_mode: HiFi 模式 ('auto', 'sampling', 'simple')
        hifi_profile_id: HiFi profile ID
        hifi_profile_root: HiFi profile 目录
        hifi_len_min/peak_min/peak_max/max: HiFi 长度分布参数
        hifi_qmin/qmean/qsd: HiFi 质量参数
        hifi_total_reads: HiFi 总读数（simple 模式）
        sr_mean/std: NGS insert size 参数
        sr_readlen: NGS 读长
        sr_platform: ART 平台
        generate_truth: 是否生成 truth 文件
    """
    def __init__(
        self,
        sample: str,
        csv: str,
        path: str,
        seed: Optional[int] = None,
        thread: int = 8,
        skip_sr: bool = False,
        skip_hifi: bool = False,
        skip_ont: bool = False,
        ont_mean: float = 3000,
        ont_std: float = 2500,
        ont_model: str = 'R94',
        hifi_sample_fastq: Optional[str] = None,
        hifi_mode: str = 'auto',
        hifi_profile_id: Optional[str] = None,
        hifi_profile_root: Optional[str] = None,
        hifi_len_min: int = 5000,
        hifi_len_peak_min: int = 10000,
        hifi_len_peak_max: int = 25000,
        hifi_len_max: int = 60000,
        hifi_qmin: int = 20,
        hifi_qmean: int = 30,
        hifi_qsd: float = 0.0,
        hifi_total_reads: Optional[int] = None,
        sr_mean: float = 400,
        sr_std: float = 125,
        sr_readlen: int = 150,
        sr_platform: str = 'HS25',
        generate_truth: bool = False,
    ) -> None:

        self.sample = sample
        self.path = path  # 直接使用 path，不创建 sample 子目录
        self.seed = seed
        self.thread = int(thread)
        self.skip_sr = bool(skip_sr)
        self.skip_hifi = bool(skip_hifi)
        self.skip_ont = bool(skip_ont)
        self.ont_mean = float(ont_mean)
        self.ont_std = float(ont_std)
        self.ont_model = ont_model
        self.hifi_sample_fastq = hifi_sample_fastq
        self.hifi_mode = hifi_mode
        self.hifi_profile_id = hifi_profile_id
        self.hifi_profile_root = hifi_profile_root
        self.hifi_len_min = int(hifi_len_min)
        self.hifi_len_peak_min = int(hifi_len_peak_min)
        self.hifi_len_peak_max = int(hifi_len_peak_max)
        self.hifi_len_max = int(hifi_len_max)
        self.hifi_qmin = int(hifi_qmin)
        self.hifi_qmean = int(hifi_qmean)
        self.hifi_qsd = float(hifi_qsd)
        self.hifi_total_reads = None if hifi_total_reads is None else int(hifi_total_reads)
        self.sr_mean = float(sr_mean)
        self.sr_std = float(sr_std)
        self.sr_readlen = int(sr_readlen)
        self.sr_platform = sr_platform
        self.generate_truth = bool(generate_truth)
        self.pool_csv_path = os.path.abspath(csv)  # 保存 pool CSV 路径供 truth 生成使用

        self.utils = utilities('')
        self.ec = pd.read_csv(csv, sep='\t')
        if 'id' not in self.ec.columns or 'seq' not in self.ec.columns:
            raise ValueError("pool CSV/TSV must contain columns: 'id' and 'seq'")
        if 'tempcov' not in self.ec.columns:
            self.ec['tempcov'] = 1.0
        else:
            self.ec['tempcov'] = pd.to_numeric(self.ec['tempcov'], errors='coerce').fillna(0.0)
        # Convert to absolute paths to avoid issues with os.chdir
        self.path = os.path.abspath(self.path)
        self.unifa = join(self.path, f'.unifa_{self.sample}')  # 隐藏临时目录
        self.tmp = join(self.path, f'.tmp_{self.sample}')      # 隐藏临时目录
        self.tool_path = split(__file__)[0]

        if not exists(self.unifa):
            os.makedirs(self.unifa)
        if not exists(self.tmp):
            os.makedirs(self.tmp)

        # Only create per-molecule FASTA files when required by external tools.
        # - NGS: we use batched ART on multi-FASTA (no per-molecule files needed).
        # - HiFi sampling (PBSIM2) and ONT (PBSIM2) still require per-molecule FASTA.
        need_unifasta = (not self.skip_ont) or (
            (not self.skip_hifi) and (self.hifi_effective_mode() == 'sampling')
        )
        if need_unifasta:
            self.unifasta()
        if not self.skip_sr:
            self.sim_fastq_sr()
        if not self.skip_ont:
            self.sim_fastq_ont()
        if not self.skip_hifi:
            self.sim_fastq_hifi()
        self.sort()
        if self.generate_truth:
            self._generate_truth_files()
        self.generate_readme()

    def unifasta(self):
        '''
        Function to write each single eccDNA into a fasta file
        '''
        for i in self.ec.index:
            tmp_fasta = join(self.unifa, self.ec.loc[i,'id'] + '.fasta')
            self.utils.write_fasta(self.ec.loc[i:i,:], tmp_fasta)
        return

    def _art_env(self) -> dict:
        """
        Build environment for running ART.

        On macOS + conda installs, art_illumina can fail to locate libcblas.3.dylib
        due to broken symlinks. We create a local shim in the temp directory and
        prepend it to DYLD_LIBRARY_PATH when needed.
        """
        env = os.environ.copy()
        if sys.platform != 'darwin':
            return env

        # If libcblas.3.dylib resolves in the current Python prefix, nothing to do.
        lib_dir = join(sys.prefix, 'lib')
        cblas_path = join(lib_dir, 'libcblas.3.dylib')
        if exists(cblas_path):
            return env

        try:
            import glob

            shim_dir = join(self.tmp, '.dyld')
            os.makedirs(shim_dir, exist_ok=True)

            # Find a usable OpenBLAS dylib and link it as libcblas.3.dylib.
            candidates = []
            if os.path.isdir(lib_dir):
                candidates.extend(glob.glob(join(lib_dir, 'libopenblasp*.dylib')))
                candidates.extend(glob.glob(join(lib_dir, 'libopenblas*.dylib')))

            target = next((p for p in candidates if exists(p)), None)
            if target:
                link_path = join(shim_dir, 'libcblas.3.dylib')
                if os.path.lexists(link_path):
                    os.remove(link_path)
                os.symlink(target, link_path)

                prev = env.get('DYLD_LIBRARY_PATH', '')
                env['DYLD_LIBRARY_PATH'] = shim_dir if not prev else f"{shim_dir}:{prev}"
            else:
                logger.warning('ART dependency shim: no libopenblas*.dylib found under %s', lib_dir)
        except Exception as e:
            logger.warning('ART dependency shim setup failed: %s', e)
        return env

    def _write_multi_fasta(self, df: 'pd.DataFrame', fasta_path: str) -> int:
        """Write a multi-FASTA from a pool DataFrame subset. Returns number of records written."""
        records = 0
        with open(fasta_path, 'w') as out:
            for _, row in df.iterrows():
                mol_id = str(row.get('id', '')).strip()
                if not mol_id:
                    continue
                seq_value = row.get('seq', None)
                if pd.isna(seq_value):
                    continue
                seq = str(seq_value).strip()
                if not seq:
                    continue
                out.write(f">{mol_id}\n{seq}\n")
                records += 1
        return records
    
    def multi_para_art(self):
        para_file = []
        for idx, i in enumerate(self.ec.index):
            mol_id = str(self.ec.loc[i, 'id'])
            input_ = join(self.unifa, mol_id + '.fasta')
            output_ = join(self.tmp, mol_id + '.R')
            # Generate per-molecule seed using hash-based method for better uniqueness
            mol_seed = None if self.seed is None else _generate_molecule_seed(self.seed, mol_id, "NGS")
            para_file.append((input_, output_, self.sr_platform, self.sr_readlen, self.ec.loc[i,'tempcov'], self.sr_mean, self.sr_std, mol_seed, self.tmp))
        return para_file

    def art(self, input_, output_, sr_platform, length, cov, fral, frastd, mol_seed, cwd):
        mol_id = split(input_)[1].replace('.fasta', '')
        try:
            cmd = [
                'art_illumina',
                '-na',
                '-q',
                '-nf',
                '0',
                '-p',
                '-i',
                str(input_),
                '-o',
                str(output_),
                '-ss',
                str(sr_platform),
                '-l',
                str(length),
                '-f',
                str(cov),
                '-m',
                str(fral),
                '-s',
                str(frastd),
            ]
            if mol_seed is not None:
                cmd.extend(['-rs', str(mol_seed)])
            sp.run(cmd, check=True, cwd=cwd, capture_output=True, env=self._art_env())
            self.utils.transfer_files('{0}1.fq'.format(output_), '{0}/{1}.NGS.R1.fastq'.format(self.path, self.sample))
            self.utils.transfer_files('{0}2.fq'.format(output_), '{0}/{1}.NGS.R2.fastq'.format(self.path, self.sample))
            self.utils.remove_glob(output_ + '*')
            return (mol_id, True, None)
        except Exception as e:
            self.utils.remove_glob(output_ + '*')
            return (mol_id, False, str(e))

    def sim_fastq_sr(self):
        '''
        Generate NGS fastq reads using ART.

        并行化策略：
        1. 按 coverage 分组成多个 chunk
        2. 每个 chunk 并行执行 ART（单线程）
        3. 最后按顺序合并结果，避免冲突
        '''
        if not _check_tool_available('art_illumina'):
            logger.warning(
                'art_illumina not found or not runnable. '
                'NGS simulation requires ART. Install with: conda install -c bioconda art'
            )
            return

        # Clean/validate pool
        df = self.ec[['id', 'seq', 'tempcov']].copy()
        df['id'] = df['id'].astype(str)
        df['tempcov'] = pd.to_numeric(df['tempcov'], errors='coerce').fillna(0.0)
        df = df[df['tempcov'] > 0]
        df = df[df['seq'].notna()]
        if df.empty:
            logger.warning('NGS simulation skipped: no molecules with tempcov > 0.')
            return

        # Remove existing outputs to avoid accidental appends.
        out_r1 = '{0}/{1}.NGS.R1.fastq'.format(self.path, self.sample)
        out_r2 = '{0}/{1}.NGS.R2.fastq'.format(self.path, self.sample)
        for out_path in (out_r1, out_r2):
            if exists(out_path):
                os.remove(out_path)

        # Batch ART calls by (rounded) coverage to avoid spawning thousands of processes.
        max_groups = max(1, min(64, int(self.thread) * 4))
        round_digits = 6
        for d in (6, 4, 3, 2, 1, 0):
            if df['tempcov'].round(d).nunique() <= max_groups:
                round_digits = d
                break
        df['_cov_key'] = df['tempcov'].round(round_digits)

        # =================================================================
        # 阶段 1: 准备所有批次任务（串行写入 fasta 文件）
        # =================================================================
        batch_tasks = []
        art_env = self._art_env()  # 预先获取环境变量（避免在子进程中重复创建）

        for idx, (cov_key, group) in enumerate(df.groupby('_cov_key', sort=True), start=1):
            cov = float(cov_key)
            if cov <= 0:
                continue

            fasta_path = join(self.tmp, f'art_batch_{idx}.fasta')
            output_prefix = join(self.tmp, f'art_batch_{idx}')
            written = self._write_multi_fasta(group, fasta_path)
            if written <= 0:
                self.utils.remove_glob(fasta_path)
                continue

            batch_seed = None
            if self.seed is not None:
                batch_seed = _generate_molecule_seed(int(self.seed), f"batch_{idx}", "NGS")

            # 构建任务参数元组
            task = (
                idx,                  # batch_idx
                fasta_path,           # fasta_path
                output_prefix,        # output_prefix
                cov,                  # coverage
                self.sr_platform,     # sr_platform
                self.sr_readlen,      # sr_readlen
                self.sr_mean,         # sr_mean
                self.sr_std,          # sr_std
                batch_seed,           # batch_seed
                self.tmp,             # cwd
                art_env,              # art_env
            )
            batch_tasks.append(task)

        if not batch_tasks:
            logger.warning('NGS simulation: no valid batches to process.')
            return

        # =================================================================
        # 阶段 2: 并行执行 ART 批次
        # =================================================================
        results = []
        n_workers = max(1, int(self.thread))

        if n_workers > 1 and len(batch_tasks) > 1:
            # 多进程并行执行
            logger.info('NGS simulation: running %d batches in parallel with %d workers',
                        len(batch_tasks), n_workers)
            with mp.Pool(n_workers) as pool:
                results = list(tqdm(
                    pool.imap(_run_art_batch_worker, batch_tasks),
                    total=len(batch_tasks),
                    desc="ART simulation"
                ))
        else:
            # 单线程串行执行
            logger.info('NGS simulation: running %d batches sequentially', len(batch_tasks))
            for task in tqdm(batch_tasks, desc="ART simulation"):
                results.append(_run_art_batch_worker(task))

        # =================================================================
        # 阶段 3: 按顺序合并结果（避免合并冲突）
        # =================================================================
        # 按 batch_idx 排序，确保合并顺序一致
        results.sort(key=lambda x: x[0])

        failed_batches = 0
        successful_r1_files = []
        successful_r2_files = []

        for batch_idx, success, r1_path, r2_path, error_msg in results:
            if success and r1_path and r2_path:
                successful_r1_files.append(r1_path)
                successful_r2_files.append(r2_path)
            else:
                failed_batches += 1
                logger.warning('NGS batch %d failed: %s', batch_idx, error_msg)

        # 串行合并所有成功的输出文件
        if successful_r1_files:
            self.utils.concat_files(successful_r1_files, out_r1, append=False)
        if successful_r2_files:
            self.utils.concat_files(successful_r2_files, out_r2, append=False)

        # =================================================================
        # 阶段 4: 清理临时文件
        # =================================================================
        for task in batch_tasks:
            batch_idx = task[0]
            fasta_path = task[1]
            output_prefix = task[2]
            self.utils.remove_glob(output_prefix + '*')
            self.utils.remove_glob(fasta_path)

        if failed_batches:
            logger.warning('NGS simulation had %d failed ART batches out of %d total.',
                           failed_batches, len(batch_tasks))
        else:
            logger.info('NGS simulation completed: %d batches processed successfully.',
                        len(batch_tasks))
        return

    def pbsim_hifi_env(self):
        if not self.hifi_profile_root:
            return None
        root = os.path.abspath(os.path.expanduser(self.hifi_profile_root))
        os.makedirs(root, exist_ok=True)
        env = os.environ.copy()
        env['HOME'] = root
        env['XDG_CACHE_HOME'] = join(root, '.cache')
        env['XDG_CONFIG_HOME'] = join(root, '.config')
        env['XDG_DATA_HOME'] = join(root, '.local', 'share')
        os.makedirs(env['XDG_CACHE_HOME'], exist_ok=True)
        os.makedirs(env['XDG_CONFIG_HOME'], exist_ok=True)
        os.makedirs(env['XDG_DATA_HOME'], exist_ok=True)
        return env

    def hifi_effective_mode(self):
        """
        确定 HiFi 模拟的有效模式。

        auto 模式逻辑:
        - 如果提供了 sample_fastq 或 profile_id -> sampling
        - 否则使用内置 profile -> sampling (默认)
        - 只有显式指定 --hifi-mode simple 才会使用 simple 模式
        """
        if self.hifi_mode == 'auto':
            # auto 模式默认使用 sampling（内置或自定义 profile）
            return 'sampling'
        return self.hifi_mode

    def _get_hifi_profile_cache_dir(self) -> str:
        """获取 HiFi profile 缓存目录 (~/.ecctoolkit/hifi_profiles/)"""
        cache_dir = os.path.expanduser('~/.ecctoolkit/hifi_profiles')
        os.makedirs(cache_dir, exist_ok=True)
        return cache_dir

    def _compute_fastq_hash(self, fastq_path: str) -> str:
        """计算 FASTQ 文件的简短 hash 作为 profile ID"""
        import hashlib
        # 使用文件路径和大小计算 hash（避免读取整个文件）
        file_stat = os.stat(fastq_path)
        hash_input = f"{os.path.abspath(fastq_path)}:{file_stat.st_size}:{file_stat.st_mtime}"
        return hashlib.md5(hash_input.encode()).hexdigest()[:12]

    def ensure_hifi_profile(self):
        """
        确保 PBSIM2 sampling profile 已准备好。

        策略:
        1. 如果指定了 hifi_sample_fastq，自动生成并缓存 profile
        2. 如果指定了 hifi_profile_id，检查是否存在
        3. 如果都没有指定，使用内置的默认 profile
        """
        cache_dir = self._get_hifi_profile_cache_dir()
        builtin_profiles = join(self.tool_path, 'resource', 'profiles')

        # 确定 profile_id
        if self.hifi_sample_fastq:
            # 用户提供了自己的 HiFi 数据，生成自定义 profile
            if self.hifi_profile_id:
                profile_id = self.hifi_profile_id
            else:
                profile_id = 'auto_' + self._compute_fastq_hash(self.hifi_sample_fastq)
            self.hifi_profile_id = profile_id
            logger.info(f'Using custom HiFi data: {self.hifi_sample_fastq}')
            logger.info(f'Profile ID: {profile_id}')
        elif self.hifi_profile_id:
            # 用户指定了 profile_id
            profile_id = self.hifi_profile_id
        else:
            # 使用内置默认 profile
            profile_id = BUILTIN_HIFI_PROFILE_ID
            self.hifi_profile_id = profile_id
            logger.info(f'Using built-in HiFi profile: {profile_id} (based on 2000 real HiFi reads)')
            logger.info('Tip: For more realistic simulation, provide your own HiFi data with --hifi-sample-fastq')

        # 检查 profile 是否已存在（缓存或内置）
        profile_locations = [
            join(cache_dir, f'sample_profile_{profile_id}'),
            join(builtin_profiles, f'sample_profile_{profile_id}'),
        ]

        for loc in profile_locations:
            if exists(loc + '.fastq') and exists(loc + '.stats'):
                logger.debug(f'Found HiFi profile at: {loc}')
                return

        # Profile 不存在，需要生成
        if not self.hifi_sample_fastq:
            logger.warning(f'HiFi profile "{profile_id}" not found and no sample fastq provided.')
            logger.warning('Falling back to simple mode.')
            self.hifi_mode = 'simple'  # 强制使用 simple 模式
            return

        logger.info(f'Generating HiFi profile from: {self.hifi_sample_fastq}')
        logger.info(f'Profile will be cached at: {cache_dir}')

        # 创建 dummy fasta 用于 profile 生成
        dummy_fasta = join(self.tmp, '__ecsim_hifi_profile__.fasta')
        dummy_prefix = join(self.tmp, '__ecsim_hifi_profile__')

        # 使用随机序列避免 PBSIM 在短序列上崩溃
        import random
        random.seed(42)
        dummy_seq = ''.join(random.choices('ACGT', k=50000))
        with open(dummy_fasta, 'w') as f:
            f.write('>ecsim_hifi_profile\n')
            f.write(dummy_seq)
            f.write('\n')

        cmd_parts = [
            'pbsim',
            '--depth', '0.1',
            '--prefix', str(dummy_prefix),
            '--id-prefix', 'ecsim_hifi_profile',
            '--sample-fastq', str(self.hifi_sample_fastq),
            '--sample-profile-id', str(profile_id),
            str(dummy_fasta),
        ]

        try:
            # 在缓存目录中执行，这样 profile 文件直接生成在那里
            sp.run(cmd_parts, check=True, env=self.pbsim_hifi_env(), cwd=cache_dir, capture_output=True)
            logger.info(f'HiFi profile cached: {profile_id}')
        except Exception as e:
            logger.warning(f'Failed to generate HiFi profile: {e}')
            logger.warning('Falling back to simple mode.')
            self.hifi_mode = 'simple'
        finally:
            # 清理临时文件
            self.utils.remove_glob(dummy_prefix + '*')
            self.utils.remove_glob(dummy_fasta)
            # 清理缓存目录中的临时输出
            for pattern in ['__ecsim_hifi_profile__*', 'ecsim_hifi_profile*']:
                for f in glob.glob(join(cache_dir, pattern)):
                    try:
                        os.remove(f)
                    except Exception:
                        pass

    def _copy_hifi_profile_to_tmp(self):
        """
        复制 PBSIM2 HiFi profile 文件到 tmp 目录供多进程使用。

        搜索顺序:
        1. 用户缓存目录 (~/.ecctoolkit/hifi_profiles/)
        2. 当前工作目录
        3. 内置 profiles 目录 (resource/profiles/)
        4. 用户指定的 hifi_profile_root 目录
        """
        if not self.hifi_profile_id:
            return

        profile_base = f'sample_profile_{self.hifi_profile_id}'
        cache_dir = self._get_hifi_profile_cache_dir()
        builtin_profiles = join(self.tool_path, 'resource', 'profiles')

        search_paths = [
            cache_dir,         # 用户缓存目录（优先）
            '.',               # 当前目录
            builtin_profiles,  # 内置 profiles
        ]
        if self.hifi_profile_root:
            search_paths.append(self.hifi_profile_root)

        for ext in ['.fastq', '.stats']:
            dst = join(self.tmp, profile_base + ext)
            if exists(dst):
                continue
            # 搜索 profile 文件
            for search_dir in search_paths:
                src = join(search_dir, profile_base + ext)
                if exists(src):
                    shutil.copy2(src, dst)
                    logger.debug(f'Copied HiFi profile: {src} -> {dst}')
                    break
        return

    def multi_para_pbsim2(self):
        para_file = []
        model_file = join(self.tool_path, 'resource', 'pbsim2', self.ont_model + '.model')
        for idx, i in enumerate(self.ec.index):
            mol_id = str(self.ec.loc[i, 'id'])
            input_ = join(self.unifa, mol_id + '.fasta')
            output_ = join(self.tmp, mol_id)
            # Generate per-molecule seed using hash-based method for better uniqueness
            mol_seed = None if self.seed is None else _generate_molecule_seed(self.seed, mol_id, "ONT")
            para_file.append(
                (
                    self.ec.loc[i, 'tempcov'],
                    output_,
                    mol_id,
                    model_file,
                    self.ont_mean,
                    self.ont_std,
                    input_,
                    mol_seed,
                    self.tmp,
                )
            )
        return para_file

    def pbsim2(self, cov, output_, id_prefix, model, ont_mean, ont_std, input_, mol_seed, cwd):
        try:
            # 在 id_prefix 后加 '.' 分隔符，避免 PBSIM2 直接拼接序列号导致 ID 混淆
            # 例如: ecc_1 -> ecc_1. -> ecc_1.1_1 (而不是 ecc_11_1)
            cmd = [
                'pbsim',
                '--depth',
                str(cov),
                '--prefix',
                str(output_),
                '--id-prefix',
                str(id_prefix) + '.',
                '--hmm_model',
                str(model),
                '--length-mean',
                str(ont_mean),
                '--length-sd',
                str(ont_std),
            ]
            if mol_seed is not None:
                cmd.extend(['--seed', str(mol_seed)])
            cmd.append(str(input_))
            sp.run(cmd, check=True, cwd=cwd, capture_output=True)
            self.utils.transfer_files('{0}_0001.fastq'.format(output_), '{0}/{1}.ONT.fastq'.format(self.path, self.sample))
            self.utils.remove_glob(output_ + '*')
            return (id_prefix, True, None)
        except Exception as e:
            self.utils.remove_glob(output_ + '*')
            return (id_prefix, False, str(e))

    def sim_fastq_ont(self):
        '''
        Generate ONT fastq reads using PBSIM2.
        '''
        # 检测 pbsim2 是否可用（需要支持 --hmm_model 参数）
        if not self._check_pbsim2_available():
            logger.warning(
                'PBSIM2 not found or incompatible version. '
                'ONT simulation requires PBSIM2 with --hmm_model support. '
                'Install with: conda install -c bioconda pbsim2'
            )
            logger.warning('Skipping ONT read simulation.')
            return

        self.pbsim2_para_file = self.multi_para_pbsim2()
        with mp.Pool(self.thread) as pool:
            results = pool.starmap(self.pbsim2, self.pbsim2_para_file)
        failed = [(mol_id, err) for mol_id, success, err in results if not success]
        if failed:
            logger.warning('ONT simulation failed for %d molecules: %s', len(failed), [m for m, _ in failed[:5]])
        return

    def _check_pbsim2_available(self) -> bool:
        """检测 PBSIM2 是否可用（支持 --hmm_model 参数，用于 ONT 模拟）"""
        # 先检查 pbsim 命令是否存在
        if not _check_tool_available('pbsim'):
            return False
        # 检查是否支持 --hmm_model 参数（PBSIM2 特有，用于 ONT）
        try:
            result = sp.run(['pbsim', '--help'], capture_output=True, text=True)
            help_text = result.stdout + result.stderr
            return '--hmm_model' in help_text or 'hmm_model' in help_text
        except Exception:
            return False

    def _check_pbsim_hifi_sampling_available(self) -> bool:
        """检测 PBSIM2/3 HiFi sampling 模式是否可用（支持 --sample-fastq 参数）"""
        if not _check_tool_available('pbsim'):
            return False
        try:
            result = sp.run(['pbsim', '--help'], capture_output=True, text=True)
            help_text = result.stdout + result.stderr
            # PBSIM2/3 的 sampling 模式需要 --sample-fastq 或 --sample-profile-id
            return '--sample-fastq' in help_text or '--sample-profile-id' in help_text
        except Exception:
            return False

    def _get_pbsim_version_info(self) -> str:
        """获取 PBSIM 版本信息用于调试"""
        if not _check_tool_available('pbsim'):
            return "pbsim not found"
        try:
            result = sp.run(['pbsim', '--version'], capture_output=True, text=True)
            version_text = (result.stdout + result.stderr).strip()
            if version_text:
                return version_text.split('\n')[0]
            # 如果 --version 不工作，尝试从 --help 获取
            result = sp.run(['pbsim', '--help'], capture_output=True, text=True)
            help_text = result.stdout + result.stderr
            for line in help_text.split('\n')[:5]:
                if 'version' in line.lower() or 'pbsim' in line.lower():
                    return line.strip()
            return "version unknown"
        except Exception as e:
            return f"error: {e}"

    def multi_para_pbsim2_hifi(self):
        para_file = []
        sample_fastq_for_jobs = None if self.hifi_profile_id else self.hifi_sample_fastq
        for idx, i in enumerate(self.ec.index):
            mol_id = str(self.ec.loc[i, 'id'])
            input_ = join(self.unifa, mol_id + '.fasta')
            output_ = join(self.tmp, mol_id)
            # Generate per-molecule seed using hash-based method for better uniqueness
            mol_seed = None if self.seed is None else _generate_molecule_seed(self.seed, mol_id, "HiFi")
            para_file.append(
                (
                    self.ec.loc[i, 'tempcov'],
                    output_,
                    mol_id,
                    sample_fastq_for_jobs,
                    self.hifi_profile_id,
                    input_,
                    mol_seed,
                    self.tmp,
                )
            )
        return para_file

    def pbsim2_hifi(self, cov, output_, id_prefix, sample_fastq, profile_id, input_, mol_seed, cwd):
        try:
            if not sample_fastq and not profile_id:
                raise ValueError(
                    'HiFi simulation requires either hifi_sample_fastq (PBSIM2 --sample-fastq) or hifi_profile_id (PBSIM2 --sample-profile-id).'
                )
            # 在 id_prefix 后加 '.' 分隔符，避免 PBSIM2 直接拼接序列号导致 ID 混淆
            cmd_parts = [
                'pbsim',
                '--depth',
                str(cov),
                '--prefix',
                str(output_),
                '--id-prefix',
                str(id_prefix) + '.',
            ]
            if sample_fastq:
                cmd_parts += ['--sample-fastq', str(sample_fastq)]
            if profile_id:
                cmd_parts += ['--sample-profile-id', str(profile_id)]
            if mol_seed is not None:
                cmd_parts += ['--seed', str(mol_seed)]
            cmd_parts.append(str(input_))
            sp.run(cmd_parts, check=True, env=self.pbsim_hifi_env(), cwd=cwd, capture_output=True)
            # Only transfer if output file exists (PBSIM2 may not produce reads for very short molecules)
            output_fastq = '{0}_0001.fastq'.format(output_)
            if exists(output_fastq):
                self.utils.transfer_files(
                    output_fastq,
                    '{0}/{1}.HiFi.fastq'.format(self.path, self.sample),
                )
            self.utils.remove_glob(output_ + '*')
            return (id_prefix, True, None)
        except Exception as e:
            self.utils.remove_glob(output_ + '*')
            return (id_prefix, False, str(e))

    def sim_fastq_hifi(self):
        '''
        Generate HiFi fastq reads using PBSIM2 or simple mode.
        '''
        mode = self.hifi_effective_mode()
        if mode == 'sampling':
            if self.hifi_total_reads is not None:
                raise ValueError('hifi_total_reads is supported in HiFi simple-mode only.')
            # 检测 PBSIM2/3 是否支持 sampling 模式
            if not self._check_pbsim_hifi_sampling_available():
                pbsim_info = self._get_pbsim_version_info()
                logger.warning(
                    f'PBSIM HiFi sampling mode not available ({pbsim_info}). '
                    'HiFi sampling requires PBSIM2/3 with --sample-fastq support. '
                    'Falling back to simple mode.'
                )
                self.sim_fastq_hifi_simple()
                self._generate_hifi_fasta()
                return
            # 确保 profile 已准备好（会设置内置 profile 或用户自定义 profile）
            self.ensure_hifi_profile()
            # Copy profile files to tmp directory for PBSIM2 to find
            self._copy_hifi_profile_to_tmp()
            self.pbsim2_hifi_para_file = self.multi_para_pbsim2_hifi()
            with mp.Pool(self.thread) as pool:
                results = pool.starmap(self.pbsim2_hifi, self.pbsim2_hifi_para_file)
            failed = [(mol_id, err) for mol_id, success, err in results if not success]
            if failed:
                logger.warning('HiFi simulation failed for %d molecules: %s', len(failed), [m for m, _ in failed[:5]])
            self._generate_hifi_fasta()
            return
        if mode == 'simple':
            logger.info('HiFi simple-mode enabled (no sample fastq/profile provided); results are approximate.')
            self.sim_fastq_hifi_simple()
            self._generate_hifi_fasta()
            return
        raise ValueError("unknown hifi_mode: {0}".format(mode))

    def _generate_hifi_fasta(self):
        '''
        Convert HiFi fastq to fasta format
        '''
        fastq_path = '{0}/{1}.HiFi.fastq'.format(self.path, self.sample)
        fasta_path = '{0}/{1}.HiFi.fasta'.format(self.path, self.sample)
        if not exists(fastq_path):
            return

        # Try seqkit first
        if _check_tool_available('seqkit'):
            try:
                with open(fasta_path, 'wb') as out:
                    sp.run(['seqkit', 'fq2fa', fastq_path], check=True, stdout=out, stderr=sp.DEVNULL)
                return
            except sp.CalledProcessError as e:
                logger.warning(f'seqkit failed (exit code {e.returncode}); falling back to pure Python conversion.')
            except Exception as e:
                logger.warning(f'seqkit error: {e}; falling back to pure Python conversion.')

        # Fallback: simple FASTQ -> FASTA conversion (no quality filtering).
        logger.info('Converting HiFi FASTQ to FASTA in pure Python...')
        with open(fastq_path, 'r', encoding='utf-8', errors='replace') as fq, open(
            fasta_path, 'w', encoding='utf-8'
        ) as fa:
            while True:
                header = fq.readline()
                if not header:
                    break
                seq = fq.readline()
                plus = fq.readline()
                qual = fq.readline()
                if not qual:
                    break
                rid = header.strip().lstrip('@').split()[0]
                fa.write(f">{rid}\n{seq.strip()}\n")

    def sample_hifi_read_length(self, rng):
        if self.hifi_len_min <= 0:
            raise ValueError('hifi_len_min must be > 0')
        if not (self.hifi_len_min <= self.hifi_len_peak_min <= self.hifi_len_peak_max <= self.hifi_len_max):
            raise ValueError('HiFi length bounds must satisfy min <= peak-min <= peak-max <= max')
        r = rng.random()
        if r < 0.85:
            return int(rng.integers(self.hifi_len_peak_min, self.hifi_len_peak_max + 1))
        if r < 0.95:
            return int(rng.integers(self.hifi_len_min, self.hifi_len_peak_min + 1))
        return int(rng.integers(self.hifi_len_peak_max, self.hifi_len_max + 1))

    def sample_hifi_qualities(self, rng, length):
        qmin = max(0, int(self.hifi_qmin))
        qmean = max(qmin, int(self.hifi_qmean))
        qsd = float(self.hifi_qsd)
        if qsd <= 0:
            q = qmean
            if q > 93:
                q = 93
            return bytes([q + 33]) * length
        quals = rng.normal(loc=qmean, scale=qsd, size=length)
        quals = np.clip(np.rint(quals), qmin, 93).astype(np.int16)
        return bytes((int(q) + 33 for q in quals))

    def introduce_substitution_errors(self, rng, seq_bytes, qmean):
        if not seq_bytes:
            return seq_bytes
        qmean = max(0, int(qmean))
        p_error = 10 ** (-qmean / 10.0) if qmean > 0 else 1.0
        expected = p_error * len(seq_bytes)
        n_err = int(rng.poisson(expected))
        if n_err <= 0:
            return seq_bytes
        if n_err >= len(seq_bytes):
            n_err = len(seq_bytes) - 1
        positions = rng.choice(len(seq_bytes), size=n_err, replace=False)
        alt = {
            ord('A'): b'CGT',
            ord('C'): b'AGT',
            ord('G'): b'ACT',
            ord('T'): b'ACG',
        }
        for pos in positions:
            b = seq_bytes[pos]
            choices = alt.get(b)
            if not choices:
                continue
            seq_bytes[pos] = choices[int(rng.integers(0, 3))]
        return seq_bytes

    def sim_fastq_hifi_simple(self):
        out_fastq = '{0}/{1}.HiFi.fastq'.format(self.path, self.sample)
        rng = np.random.default_rng(self.seed)
        with open(out_fastq, 'w') as out:
            molecules = []
            weights = []
            for i in self.ec.index:
                ecc_id = str(self.ec.loc[i, 'id'])
                seq_value = self.ec.loc[i, 'seq']
                if pd.isna(seq_value):
                    continue
                seq = str(seq_value).upper()
                if not seq:
                    continue
                seq_len = len(seq)
                cov = float(self.ec.loc[i, 'tempcov'])
                if not np.isfinite(cov):
                    cov = 0.0
                molecules.append((ecc_id, seq, cov))
                weights.append(max(0.0, cov) * seq_len)

            if not molecules:
                raise ValueError('HiFi simple-mode: no valid molecules found in pool.')

            if self.hifi_total_reads is not None:
                if self.hifi_total_reads <= 0:
                    raise ValueError('hifi_total_reads must be > 0')
                wsum = float(np.sum(weights))
                if wsum <= 0:
                    weights = [len(seq) for _, seq, _ in molecules]
                    wsum = float(np.sum(weights))
                probs = np.array(weights, dtype=float) / wsum
                expected = probs * self.hifi_total_reads
                counts = np.floor(expected).astype(int)
                remainder = int(self.hifi_total_reads - int(counts.sum()))
                if remainder > 0:
                    extras = rng.choice(len(molecules), size=remainder, replace=True, p=probs)
                    for idx in extras:
                        counts[int(idx)] += 1

                for (ecc_id, seq, _), n_reads in zip(molecules, counts):
                    if n_reads <= 0:
                        continue
                    seq_len = len(seq)
                    for read_i in range(int(n_reads)):
                        readlen = self.sample_hifi_read_length(rng)
                        if readlen > seq_len:
                            readlen = seq_len
                        if readlen <= 0:
                            continue
                        if seq_len > readlen:
                            start = int(rng.integers(0, seq_len - readlen + 1))
                        else:
                            start = 0
                        read_seq = seq[start : start + readlen]
                        seq_bytes = bytearray(read_seq.encode('ascii'))
                        seq_bytes = self.introduce_substitution_errors(rng, seq_bytes, self.hifi_qmean)
                        qual_bytes = self.sample_hifi_qualities(rng, readlen)
                        out.write('@{0}_hifi_{1}\n'.format(ecc_id, read_i))
                        out.write(seq_bytes.decode('ascii'))
                        out.write('\n+\n')
                        out.write(qual_bytes.decode('ascii'))
                        out.write('\n')
                return

            for ecc_id, seq, cov in molecules:
                if (not np.isfinite(cov)) or cov <= 0:
                    continue
                seq_len = len(seq)
                target_bases = int(np.ceil(cov * seq_len))
                produced = 0
                read_i = 0
                while produced < target_bases:
                    readlen = self.sample_hifi_read_length(rng)
                    if readlen > seq_len:
                        readlen = seq_len
                    if readlen <= 0:
                        break
                    if seq_len > readlen:
                        start = int(rng.integers(0, seq_len - readlen + 1))
                    else:
                        start = 0
                    read_seq = seq[start : start + readlen]
                    seq_bytes = bytearray(read_seq.encode('ascii'))
                    seq_bytes = self.introduce_substitution_errors(rng, seq_bytes, self.hifi_qmean)
                    qual_bytes = self.sample_hifi_qualities(rng, readlen)
                    out.write('@{0}_hifi_{1}\n'.format(ecc_id, read_i))
                    out.write(seq_bytes.decode('ascii'))
                    out.write('\n+\n')
                    out.write(qual_bytes.decode('ascii'))
                    out.write('\n')
                    read_i += 1
                    produced += readlen
        return
    
    def sort(self):
        '''
        functions to sort fastq files by read name
        '''
        sort_targets = []
        # NGS files
        if exists('{0}/{1}.NGS.R1.fastq'.format(self.path, self.sample)):
            sort_targets.append('{0}/{1}.NGS.R1.fastq'.format(self.path, self.sample))
        if exists('{0}/{1}.NGS.R2.fastq'.format(self.path, self.sample)):
            sort_targets.append('{0}/{1}.NGS.R2.fastq'.format(self.path, self.sample))
        # ONT files
        if exists('{0}/{1}.ONT.fastq'.format(self.path, self.sample)):
            sort_targets.append('{0}/{1}.ONT.fastq'.format(self.path, self.sample))
        # HiFi files
        if exists('{0}/{1}.HiFi.fastq'.format(self.path, self.sample)):
            sort_targets.append('{0}/{1}.HiFi.fastq'.format(self.path, self.sample))

        if _check_tool_available('seqkit'):
            # Sort each file in place using temp file
            for file_path in sort_targets:
                tmp_path = file_path + '.tmp'
                with open(tmp_path, 'wb') as out:
                    sp.run(['seqkit', 'sort', file_path], check=True, stdout=out)
                os.rename(tmp_path, file_path)
        elif sort_targets:
            # Sorting is for reproducibility only; skip when seqkit is unavailable.
            logger.warning('seqkit not found; skipping FASTQ sorting.')

        # Cleanup temporary directories
        cleanup_paths = [
            self.unifa,
            self.tmp,
        ]
        for path in cleanup_paths:
            self.utils.remove_glob(path)
        return

    def _generate_truth_files(self):
        '''
        Generate per-platform truth files by parsing FASTQ and matching to pool CSV.

        Generates:
            {sample}.NGS.truth.tsv   - NGS reads truth
            {sample}.HiFi.truth.tsv  - HiFi reads truth
            {sample}.ONT.truth.tsv   - ONT reads truth
        '''
        # NGS truth (for R1, R2 shares same truth since they're paired)
        if not self.skip_sr:
            ngs_r1_fastq = '{0}/{1}.NGS.R1.fastq'.format(self.path, self.sample)
            if exists(ngs_r1_fastq):
                ngs_truth_path = '{0}/{1}.NGS.truth.tsv'.format(self.path, self.sample)
                generate_truth_from_fastq(
                    fastq_path=ngs_r1_fastq,
                    pool_csv_path=self.pool_csv_path,
                    output_truth_path=ngs_truth_path,
                    platform='NGS',
                    read_id_parser='art'
                )

        # HiFi truth
        if not self.skip_hifi:
            hifi_fastq = '{0}/{1}.HiFi.fastq'.format(self.path, self.sample)
            if exists(hifi_fastq):
                hifi_truth_path = '{0}/{1}.HiFi.truth.tsv'.format(self.path, self.sample)
                # 检测 HiFi 模式来选择合适的解析器
                hifi_parser = 'hifi_simple' if self.hifi_effective_mode() == 'simple' else 'pbsim'
                generate_truth_from_fastq(
                    fastq_path=hifi_fastq,
                    pool_csv_path=self.pool_csv_path,
                    output_truth_path=hifi_truth_path,
                    platform='HiFi',
                    read_id_parser=hifi_parser
                )

        # ONT truth
        if not self.skip_ont:
            ont_fastq = '{0}/{1}.ONT.fastq'.format(self.path, self.sample)
            if exists(ont_fastq):
                ont_truth_path = '{0}/{1}.ONT.truth.tsv'.format(self.path, self.sample)
                generate_truth_from_fastq(
                    fastq_path=ont_fastq,
                    pool_csv_path=self.pool_csv_path,
                    output_truth_path=ont_truth_path,
                    platform='ONT',
                    read_id_parser='pbsim'
                )

        return

    def generate_readme(self):
        '''
        Generate README.md file with simulation details
        '''
        readme_path = '{0}/README.md'.format(self.path)

        # Collect statistics from lib.csv
        lib_csv_path = '{0}/{1}.lib.csv'.format(self.path, self.sample)
        stats = {}
        if exists(lib_csv_path):
            lib_df = pd.read_csv(lib_csv_path, sep='\t')
            stats['total_molecules'] = len(lib_df)
            stats['total_length'] = int(lib_df['length'].sum()) if 'length' in lib_df.columns else 'N/A'
            stats['mean_length'] = lib_df['length'].mean() if 'length' in lib_df.columns else 0
            stats['mean_realcov'] = lib_df['realcov'].mean() if 'realcov' in lib_df.columns else 0
            stats['mean_tempcov'] = lib_df['tempcov'].mean() if 'tempcov' in lib_df.columns else 0

        with open(readme_path, 'w') as f:
            f.write('# {0} - eccDNA Simulation Results\n\n'.format(self.sample))
            f.write('Generated by **ecsim**\n\n')

            # Output Files Section
            f.write('## Output Files\n\n')
            f.write('| File | Description |\n')
            f.write('|------|-------------|\n')

            # Check and list actual generated files
            if exists('{0}/{1}.NGS.R1.fastq'.format(self.path, self.sample)):
                f.write('| `{0}.NGS.R1.fastq` | NGS paired-end reads (Read 1, Illumina) |\n'.format(self.sample))
                f.write('| `{0}.NGS.R2.fastq` | NGS paired-end reads (Read 2, Illumina) |\n'.format(self.sample))
            if exists('{0}/{1}.ONT.fastq'.format(self.path, self.sample)):
                f.write('| `{0}.ONT.fastq` | Oxford Nanopore long reads |\n'.format(self.sample))
            if exists('{0}/{1}.HiFi.fastq'.format(self.path, self.sample)):
                f.write('| `{0}.HiFi.fastq` | PacBio HiFi long reads (FASTQ format) |\n'.format(self.sample))
            if exists('{0}/{1}.HiFi.fasta'.format(self.path, self.sample)):
                f.write('| `{0}.HiFi.fasta` | PacBio HiFi long reads (FASTA format) |\n'.format(self.sample))
            if exists('{0}/{1}.lib.csv'.format(self.path, self.sample)):
                f.write('| `{0}.lib.csv` | Molecule library information (Ground Truth) |\n'.format(self.sample))
                f.write('| `{0}.lib.bed` | Molecule library genomic coordinates |\n'.format(self.sample))
                f.write('| `{0}.lib.fasta` | Molecule library sequences |\n'.format(self.sample))
            if exists('{0}/{1}.neg.csv'.format(self.path, self.sample)):
                f.write('| `{0}.neg.*` | Negative strand eccDNA molecules |\n'.format(self.sample))
                f.write('| `{0}.pos.*` | Positive strand eccDNA molecules |\n'.format(self.sample))

            f.write('\n')

            # CSV Column Descriptions
            f.write('## CSV Column Descriptions\n\n')
            f.write('### lib.csv / neg.csv / pos.csv\n\n')
            f.write('| Column | Description |\n')
            f.write('|--------|-------------|\n')
            f.write('| `id` | Unique molecule identifier |\n')
            f.write('| `fragN` | Number of fragments (1 for simple eccDNA, >1 for chimeric) |\n')
            f.write('| `region` | Genomic region(s) in format `chr:start-end` (pipe-separated for chimeric) |\n')
            f.write('| `length` | Total sequence length (bp) |\n')
            f.write('| `seq` | DNA sequence |\n')
            f.write('| `realcov` | Assigned sequencing coverage for this molecule |\n')
            f.write('| `tempcov` | Template coverage weight after RCA adjustment |\n')
            f.write('\n')

            # BED Format
            f.write('### BED File Format\n\n')
            f.write('| Column | Description |\n')
            f.write('|--------|-------------|\n')
            f.write('| 1 | Chromosome |\n')
            f.write('| 2 | Start position (0-based) |\n')
            f.write('| 3 | End position |\n')
            f.write('| 4 | Fragment length (bp) |\n')
            f.write('| 5 | Molecule ID |\n')
            f.write('\n')

            # Simulation Parameters
            f.write('## Simulation Parameters\n\n')
            f.write('| Parameter | Value |\n')
            f.write('|-----------|-------|\n')
            f.write('| Sample Name | `{0}` |\n'.format(self.sample))
            f.write('| Random Seed | `{0}` |\n'.format(self.seed if self.seed else 'None'))
            f.write('| Threads | `{0}` |\n'.format(self.thread))
            # Platform flags
            platforms = []
            if not self.skip_sr:
                platforms.append('NGS')
            if not self.skip_hifi:
                platforms.append('HiFi')
            if not self.skip_ont:
                platforms.append('ONT')
            f.write('| Platforms | `{0}` |\n'.format(', '.join(platforms) if platforms else 'None'))
            if not self.skip_sr:
                f.write('| SR Platform | `{0}` |\n'.format(self.sr_platform))
                f.write('| SR Read Length | `{0}` bp |\n'.format(self.sr_readlen))
                f.write('| SR Insert Mean | `{0}` bp |\n'.format(self.sr_mean))
                f.write('| SR Insert SD | `{0}` bp |\n'.format(self.sr_std))
            if not self.skip_ont:
                f.write('| ONT Model | `{0}` |\n'.format(self.ont_model))
                f.write('| ONT Mean Length | `{0}` bp |\n'.format(self.ont_mean))
                f.write('| ONT Length SD | `{0}` bp |\n'.format(self.ont_std))
            if not self.skip_hifi:
                f.write('| HiFi Mode | `{0}` |\n'.format(self.hifi_effective_mode()))
                f.write('| HiFi Length Range | `{0}-{1}` bp |\n'.format(self.hifi_len_min, self.hifi_len_max))
                f.write('| HiFi Peak Range | `{0}-{1}` bp |\n'.format(self.hifi_len_peak_min, self.hifi_len_peak_max))
                f.write('| HiFi Quality (mean) | Q`{0}` |\n'.format(self.hifi_qmean))
            f.write('\n')

            # Coverage Statistics
            if stats:
                f.write('## Coverage Statistics\n\n')
                f.write('| Metric | Value |\n')
                f.write('|--------|-------|\n')
                f.write('| Total Molecules | `{0}` |\n'.format(stats.get('total_molecules', 'N/A')))
                f.write('| Total Sequence Length | `{0:,}` bp |\n'.format(stats.get('total_length', 'N/A')))
                f.write('| Mean Molecule Length | `{0:.1f}` bp |\n'.format(stats.get('mean_length', 0)))
                f.write('| Mean Real Coverage | `{0:.2f}x` |\n'.format(stats.get('mean_realcov', 0)))
                f.write('| Mean Template Coverage | `{0:.2f}` |\n'.format(stats.get('mean_tempcov', 0)))
                f.write('\n')

            # Usage Recommendations
            f.write('## Usage Recommendations\n\n')
            f.write('### For eccDNA Detection Pipelines\n\n')
            f.write('1. **Short-read analysis**: Use `{0}.NGS.R1.fastq` and `{0}.NGS.R2.fastq` with Circle-Map, ecc_finder, or similar tools.\n\n'.format(self.sample))
            f.write('2. **Long-read analysis**: Use `{0}.ONT.fastq` or `{0}.HiFi.fastq` with long-read eccDNA detection tools.\n\n'.format(self.sample))
            f.write('3. **Ground truth validation**: Use `{0}.lib.csv` containing all simulated eccDNA molecules for benchmarking.\n\n'.format(self.sample))

            f.write('### File Format Notes\n\n')
            f.write('- All FASTQ files are sorted by read name for reproducibility.\n')
            f.write('- HiFi reads are provided in both FASTQ and FASTA formats.\n')
            f.write('- BED files use 0-based coordinates.\n')
            f.write('- CSV files use tab as delimiter.\n')
            f.write('\n')

            f.write('---\n')
            f.write('*Generated by ecsim*\n')

        return
