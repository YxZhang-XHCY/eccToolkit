import pandas as pd
from itertools import combinations
import random
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import numpy as np
from pathlib import Path
import argparse
import logging
from typing import Tuple, List, Dict, Set
import sys

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 预定义颜色
BRIGHT_COLORS = [
    "255,0,0",      # 鲜红
    "0,255,0",      # 鲜绿
    "0,0,255",      # 鲜蓝
    "255,165,0",    # 橙色
    "255,0,255",    # 品红
    "0,255,255",    # 青色
    "255,255,0",    # 黄色
    "128,0,128",    # 紫色
]

GRAY_COLOR = "128,128,128"  # 灰色

def filter_chr_specific_groups(df: pd.DataFrame, target_chr: str = "chr1") -> pd.DataFrame:
    """筛选出指定染色体上的eName组（组内所有元素都在该染色体上）"""
    logger.info(f"开始筛选{target_chr}上的eName组...")
    
    # 按eName分组
    grouped = df.groupby('eName')
    
    # 筛选出所有元素都在目标染色体上的组
    chr_specific_groups = []
    total_groups = 0
    
    for name, group in grouped:
        total_groups += 1
        unique_chrs = group['eChr'].unique()
        
        # 检查是否所有元素都在目标染色体上
        if len(unique_chrs) == 1 and unique_chrs[0] == target_chr:
            chr_specific_groups.append(name)
    
    logger.info(f"总eName组数: {total_groups}")
    logger.info(f"{target_chr}专属组数: {len(chr_specific_groups)}")
    
    # 返回筛选后的数据
    filtered_df = df[df['eName'].isin(chr_specific_groups)]
    logger.info(f"筛选后数据行数: {len(filtered_df)}")
    
    return filtered_df

def get_largest_group_info(df: pd.DataFrame) -> Tuple[str, int]:
    """找出行数最多的eName组"""
    group_sizes = df.groupby('eName').size()
    largest_group_name = group_sizes.idxmax()
    largest_group_size = group_sizes.max()
    
    logger.info(f"最大组: {largest_group_name} (行数: {largest_group_size})")
    return largest_group_name, largest_group_size

def validate_input_data(df: pd.DataFrame) -> bool:
    """验证输入数据的完整性"""
    required_columns = ['eName', 'eChr', 'eStart', 'eEnd']
    missing_columns = [col for col in required_columns if col not in df.columns]
    
    if missing_columns:
        logger.error(f"缺少必需的列: {missing_columns}")
        return False
    
    # 检查数据类型
    if not pd.api.types.is_numeric_dtype(df['eStart']) or not pd.api.types.is_numeric_dtype(df['eEnd']):
        logger.error("eStart和eEnd必须是数值类型")
        return False
    
    # 检查坐标合理性
    invalid_coords = df[df['eStart'] >= df['eEnd']]
    if not invalid_coords.empty:
        logger.warning(f"发现{len(invalid_coords)}行坐标异常(eStart >= eEnd)")
    
    return True

def process_group_with_color_scheme(args: Tuple[str, pd.DataFrame, str, str]) -> List[List]:
    """处理每个分组的数据，使用指定的颜色方案"""
    name, group_df, largest_group_name, color_scheme = args
    
    # 跳过只有一个元素的组
    if len(group_df) < 2:
        return []
    
    # 根据颜色方案分配颜色
    if color_scheme == "highlight_largest":
        if name == largest_group_name:
            # 最大组使用鲜艳颜色
            color = random.choice(BRIGHT_COLORS)
        else:
            # 其他组使用灰色
            color = GRAY_COLOR
    else:
        # 原始随机颜色方案
        color = f"{random.randint(0,255)},{random.randint(0,255)},{random.randint(0,255)}"
    
    output = []
    # 使用numpy数组操作提高效率
    data_array = group_df[['eChr', 'eStart', 'eEnd']].values
    
    for i in range(len(data_array)):
        for j in range(i + 1, len(data_array)):
            output.append([
                data_array[i][0], int(data_array[i][1]), int(data_array[i][2]),
                data_array[j][0], int(data_array[j][1]), int(data_array[j][2]),
                color
            ])
    
    return output

def process_file_with_chr_filter(input_file: str, output_file: str, 
                                target_chr: str = "chr1", color_scheme: str = "highlight_largest",
                                n_processes: int = None):
    """处理文件，支持染色体筛选和颜色方案"""
    if n_processes is None:
        n_processes = cpu_count()
    
    logger.info(f"读取输入文件: {input_file}")
    df = pd.read_csv(input_file)
    
    # 验证数据
    if not validate_input_data(df):
        raise ValueError("输入文件格式不正确")
    
    # 输出原始统计信息
    logger.info(f"原始数据统计:")
    logger.info(f"  - 总行数: {len(df)}")
    logger.info(f"  - 唯一eName数: {df['eName'].nunique()}")
    logger.info(f"  - 染色体分布: {dict(df['eChr'].value_counts())}")
    
    # 筛选指定染色体的数据
    if target_chr:
        df_filtered = filter_chr_specific_groups(df, target_chr)
        if df_filtered.empty:
            logger.warning(f"没有找到完全在{target_chr}上的eName组")
            # 创建空的输出文件
            pd.DataFrame(columns=[
                "Chr1", "Start1", "End1", "Chr2", "Start2", "End2", "RGB"
            ]).to_csv(output_file, sep="\t", index=False)
            return
    else:
        df_filtered = df
    
    # 分组
    groups = list(df_filtered.groupby("eName"))
    logger.info(f"筛选后分组数: {len(groups)}")
    
    # 过滤掉只有一个元素的组
    groups_filtered = [(name, group) for name, group in groups if len(group) > 1]
    logger.info(f"有效分组数(元素>1): {len(groups_filtered)}")
    
    if not groups_filtered:
        logger.warning("没有有效的分组数据（所有组都只有一个元素）")
        pd.DataFrame(columns=[
            "Chr1", "Start1", "End1", "Chr2", "Start2", "End2", "RGB"
        ]).to_csv(output_file, sep="\t", index=False)
        return
    
    # 找出最大的组
    largest_group_name, largest_group_size = get_largest_group_info(df_filtered)
    
    # 输出分组统计
    group_sizes = {name: len(group) for name, group in groups_filtered}
    sorted_groups = sorted(group_sizes.items(), key=lambda x: x[1], reverse=True)
    logger.info(f"前5大组:")
    for i, (name, size) in enumerate(sorted_groups[:5]):
        marker = " ⭐" if name == largest_group_name else ""
        logger.info(f"  {i+1}. {name}: {size}行{marker}")
    
    # 准备参数
    process_args = [(name, group_df, largest_group_name, color_scheme) 
                   for name, group_df in groups_filtered]
    
    # 使用进程池并行处理
    logger.info(f"开始并行处理 ({n_processes} 进程)...")
    with Pool(processes=n_processes) as pool:
        results = list(tqdm(
            pool.imap_unordered(process_group_with_color_scheme, process_args), 
            total=len(groups_filtered),
            desc="Processing groups"
        ))
    
    # 展平结果
    flat_results = [row for group_result in results for row in group_result]
    
    # 创建输出DataFrame
    output_df = pd.DataFrame(flat_results, columns=[
        "Chr1", "Start1", "End1", "Chr2", "Start2", "End2", "RGB"
    ])
    
    # 统计输出信息
    logger.info(f"输出统计:")
    logger.info(f"  - 总链接数: {len(output_df)}")
    logger.info(f"  - 唯一颜色数: {output_df['RGB'].nunique()}")
    
    # 颜色使用统计
    color_counts = output_df['RGB'].value_counts()
    logger.info(f"颜色使用情况:")
    for color, count in color_counts.head().items():
        if color in BRIGHT_COLORS:
            logger.info(f"  - {color} (鲜艳色): {count}条链接")
        elif color == GRAY_COLOR:
            logger.info(f"  - {color} (灰色): {count}条链接")
        else:
            logger.info(f"  - {color}: {count}条链接")
    
    # 保存结果
    output_df.to_csv(output_file, sep="\t", index=False)
    logger.info(f"✅ 完成！输出保存到: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='生成基因组链接关系图数据（支持染色体筛选和颜色优化）')
    parser.add_argument('input', help='输入CSV文件路径')
    parser.add_argument('-o', '--output', help='输出TSV文件路径', default=None)
    parser.add_argument('-p', '--processes', type=int, default=None,
                       help='使用的进程数（默认为CPU核心数）')
    parser.add_argument('--target-chr', default="chr1",
                       help='目标染色体（默认为chr1，设为空字符串则不筛选）')
    parser.add_argument('--color-scheme', choices=['highlight_largest', 'random'], 
                       default='highlight_largest',
                       help='颜色方案：highlight_largest(突出最大组), random(随机颜色)')
    parser.add_argument('--seed', type=int, default=None,
                       help='随机种子（用于可重复的颜色生成）')
    
    args = parser.parse_args()
    
    # 设置随机种子
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)
    
    # 自动生成输出文件名
    if args.output is None:
        input_path = Path(args.input)
        if args.target_chr:
            suffix = f'.{args.target_chr}.{args.color_scheme}.tsv'
        else:
            suffix = f'.{args.color_scheme}.tsv'
        args.output = str(input_path.with_suffix(suffix))
    
    # 检查输入文件是否存在
    if not Path(args.input).exists():
        logger.error(f"输入文件不存在: {args.input}")
        sys.exit(1)
    
    try:
        # 检查文件大小
        file_size = Path(args.input).stat().st_size / (1024 * 1024)  # MB
        logger.info(f"输入文件大小: {file_size:.1f}MB")
        
        # 处理文件
        process_file_with_chr_filter(
            args.input, args.output, 
            args.target_chr if args.target_chr else None,
            args.color_scheme,
            args.processes
        )
            
    except Exception as e:
        logger.error(f"处理过程中出错: {str(e)}")
        raise

if __name__ == "__main__":
    main()
