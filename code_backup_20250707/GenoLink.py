import pandas as pd
from itertools import combinations
import random
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import numpy as np
from pathlib import Path
import argparse
import logging
from typing import Tuple, List, Dict
import sys

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 预生成颜色池，避免重复计算
COLOR_POOL = [(r, g, b) for r in range(0, 256, 15) for g in range(0, 256, 15) for b in range(0, 256, 15)]

def generate_random_rgb() -> str:
    """生成随机RGB颜色字符串"""
    if COLOR_POOL:
        r, g, b = random.choice(COLOR_POOL)
    else:
        r, g, b = random.randint(0, 255), random.randint(0, 255), random.randint(0, 255)
    return f"{r},{g},{b}"

def generate_distinct_color(used_colors: set) -> str:
    """生成与已使用颜色有足够区分度的新颜色"""
    max_attempts = 100
    min_distance = 50  # 最小颜色距离
    
    for _ in range(max_attempts):
        new_color = generate_random_rgb()
        r1, g1, b1 = map(int, new_color.split(','))
        
        is_distinct = True
        for used in used_colors:
            r2, g2, b2 = map(int, used.split(','))
            # 计算欧氏距离
            distance = ((r1-r2)**2 + (g1-g2)**2 + (b1-b2)**2) ** 0.5
            if distance < min_distance:
                is_distinct = False
                break
        
        if is_distinct:
            return new_color
    
    # 如果找不到足够不同的颜色，返回随机颜色
    return generate_random_rgb()

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

def process_group(args: Tuple[str, pd.DataFrame, bool]) -> List[List]:
    """处理每个分组的数据"""
    name, group_df, use_distinct_colors = args
    
    # 跳过只有一个元素的组
    if len(group_df) < 2:
        return []
    
    # 生成颜色
    if use_distinct_colors:
        # 使用静态颜色集合（在主进程中预生成）
        color = generate_random_rgb()
    else:
        color = generate_random_rgb()
    
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

def process_large_file_by_group(input_file: str, output_file: str, 
                               use_distinct_colors: bool = False, n_processes: int = None,
                               batch_size: int = 1000):
    """处理大文件的函数，按组批量处理以保持组的完整性"""
    if n_processes is None:
        n_processes = cpu_count()
    
    # 先读取文件头以验证格式
    df_head = pd.read_csv(input_file, nrows=10)
    if not validate_input_data(df_head):
        raise ValueError("输入文件格式不正确")
    
    logger.info("第一步：扫描文件，统计每个eName的行数...")
    
    # 第一遍扫描：统计每个eName的大小
    ename_counts = {}
    total_rows = 0
    
    with open(input_file, 'r') as f:
        header = f.readline()
        header_parts = header.strip().split(',')
        ename_idx = header_parts.index('eName')
        
        for line in tqdm(f, desc="统计eName分布"):
            total_rows += 1
            parts = line.strip().split(',')
            if len(parts) > ename_idx:
                ename = parts[ename_idx]
                ename_counts[ename] = ename_counts.get(ename, 0) + 1
    
    logger.info(f"文件总行数: {total_rows}")
    logger.info(f"唯一eName数: {len(ename_counts)}")
    
    # 按批次处理eName
    enames = list(ename_counts.keys())
    valid_groups = [name for name in enames if ename_counts[name] > 1]
    logger.info(f"有效分组数(元素>1): {len(valid_groups)}")
    
    if not valid_groups:
        logger.warning("没有有效的分组数据")
        pd.DataFrame(columns=[
            "Chr1", "Start1", "End1", "Chr2", "Start2", "End2", "RGB"
        ]).to_csv(output_file, sep="\t", index=False)
        return
    
    # 分批处理
    all_results = []
    
    for i in range(0, len(valid_groups), batch_size):
        batch_enames = set(valid_groups[i:i+batch_size])
        logger.info(f"处理批次 {i//batch_size + 1}/{(len(valid_groups)-1)//batch_size + 1}")
        
        # 读取这批eName的所有数据
        batch_df = pd.read_csv(input_file)
        batch_df = batch_df[batch_df['eName'].isin(batch_enames)]
        
        if batch_df.empty:
            continue
        
        # 分组处理
        groups = list(batch_df.groupby("eName"))
        process_args = [(name, group_df, use_distinct_colors) for name, group_df in groups]
        
        # 并行处理
        with Pool(processes=n_processes) as pool:
            batch_results = list(tqdm(
                pool.imap_unordered(process_group, process_args), 
                total=len(groups),
                desc=f"Processing batch {i//batch_size + 1}"
            ))
        
        # 收集结果
        for group_result in batch_results:
            all_results.extend(group_result)
    
    # 写入结果
    output_df = pd.DataFrame(all_results, columns=[
        "Chr1", "Start1", "End1", "Chr2", "Start2", "End2", "RGB"
    ])
    output_df.to_csv(output_file, sep="\t", index=False)
    logger.info(f"✅ 处理完成，共生成 {len(output_df)} 条链接")

def process_large_file_streaming(input_file: str, output_file: str,
                                use_distinct_colors: bool = False, n_processes: int = None):
    """使用流式处理大文件，将完整的组写入临时文件后批量处理"""
    import tempfile
    import shutil
    
    if n_processes is None:
        n_processes = cpu_count()
    
    # 先读取文件头
    with open(input_file, 'r') as f:
        header = f.readline().strip()
    
    df_head = pd.read_csv(input_file, nrows=10)
    if not validate_input_data(df_head):
        raise ValueError("输入文件格式不正确")
    
    logger.info("使用流式处理模式...")
    
    # 创建临时目录
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # 第一步：按eName分割文件
        logger.info("第一步：按eName分割数据...")
        ename_files = {}
        ename_counts = {}
        
        with open(input_file, 'r') as f:
            header_line = f.readline()
            header_parts = header_line.strip().split(',')
            ename_idx = header_parts.index('eName')
            
            for line in tqdm(f, desc="分割数据"):
                parts = line.strip().split(',')
                if len(parts) > ename_idx:
                    ename = parts[ename_idx]
                    
                    if ename not in ename_files:
                        temp_file = temp_path / f"{ename}.csv"
                        ename_files[ename] = open(temp_file, 'w')
                        ename_files[ename].write(header_line)
                        ename_counts[ename] = 0
                    
                    ename_files[ename].write(line)
                    ename_counts[ename] += 1
        
        # 关闭所有文件
        for f in ename_files.values():
            f.close()
        
        # 筛选有效组
        valid_enames = [name for name, count in ename_counts.items() if count > 1]
        logger.info(f"有效分组数: {len(valid_enames)}")
        
        # 第二步：批量处理每个组
        logger.info("第二步：处理每个分组...")
        all_results = []
        
        # 准备处理参数
        process_args = []
        for ename in valid_enames:
            temp_file = temp_path / f"{ename}.csv"
            if temp_file.exists():
                group_df = pd.read_csv(temp_file)
                process_args.append((ename, group_df, use_distinct_colors))
        
        # 并行处理
        with Pool(processes=n_processes) as pool:
            results = list(tqdm(
                pool.imap_unordered(process_group, process_args),
                total=len(process_args),
                desc="Processing groups"
            ))
        
        # 收集结果
        for group_result in results:
            all_results.extend(group_result)
    
    # 写入最终结果
    output_df = pd.DataFrame(all_results, columns=[
        "Chr1", "Start1", "End1", "Chr2", "Start2", "End2", "RGB"
    ])
    output_df.to_csv(output_file, sep="\t", index=False)
    logger.info(f"✅ 处理完成，共生成 {len(output_df)} 条链接")

def main():
    parser = argparse.ArgumentParser(description='生成基因组链接关系图数据')
    parser.add_argument('input', help='输入CSV文件路径')
    parser.add_argument('-o', '--output', help='输出TSV文件路径', default=None)
    parser.add_argument('-p', '--processes', type=int, default=None,
                       help='使用的进程数（默认为CPU核心数）')
    parser.add_argument('--distinct-colors', action='store_true',
                       help='尝试生成更有区分度的颜色')
    parser.add_argument('--seed', type=int, default=None,
                       help='随机种子（用于可重复的颜色生成）')
    parser.add_argument('--mode', choices=['auto', 'memory', 'streaming'], default='auto',
                       help='处理模式：auto(自动选择), memory(全部载入内存), streaming(流式处理)')
    parser.add_argument('--batch-size', type=int, default=1000,
                       help='批处理大小（用于大文件处理）')
    
    args = parser.parse_args()
    
    # 设置随机种子
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)
    
    # 自动生成输出文件名
    if args.output is None:
        input_path = Path(args.input)
        args.output = str(input_path.with_suffix('.Links.RandomColor.tsv'))
    
    # 检查输入文件是否存在
    if not Path(args.input).exists():
        logger.error(f"输入文件不存在: {args.input}")
        sys.exit(1)
    
    try:
        # 检查文件大小
        file_size = Path(args.input).stat().st_size / (1024 * 1024)  # MB
        logger.info(f"输入文件大小: {file_size:.1f}MB")
        
        # 决定处理模式
        if args.mode == 'auto':
            if file_size > 500:  # 大于500MB使用流式处理
                mode = 'streaming'
            elif file_size > 100:  # 100-500MB使用批处理
                mode = 'batch'
            else:
                mode = 'memory'
        else:
            mode = args.mode
        
        if mode == 'streaming':
            logger.info("使用流式处理模式（适合超大文件）")
            process_large_file_streaming(args.input, args.output, 
                                       args.distinct_colors, args.processes)
        elif mode == 'batch' or (mode == 'memory' and file_size > 100):
            logger.info("使用批量处理模式")
            process_large_file_by_group(args.input, args.output, 
                                      args.distinct_colors, args.processes,
                                      args.batch_size)
        else:
            # 小文件直接处理
            logger.info("使用内存处理模式")
            df = pd.read_csv(args.input)
            
            # 验证数据
            if not validate_input_data(df):
                sys.exit(1)
            
            # 统计信息
            logger.info(f"数据统计:")
            logger.info(f"  - 总行数: {len(df)}")
            logger.info(f"  - 唯一eName数: {df['eName'].nunique()}")
            logger.info(f"  - 染色体: {sorted(df['eChr'].unique())}")
            
            # 分组
            groups = list(df.groupby("eName"))
            logger.info(f"📊 总分组数: {len(groups)} | CPU核心数: {cpu_count()}")
            
            # 过滤掉只有一个元素的组
            groups_filtered = [(name, group) for name, group in groups if len(group) > 1]
            logger.info(f"有效分组数(元素>1): {len(groups_filtered)}")
            
            if not groups_filtered:
                logger.warning("没有有效的分组数据（所有组都只有一个元素）")
                # 创建空的输出文件
                pd.DataFrame(columns=[
                    "Chr1", "Start1", "End1", "Chr2", "Start2", "End2", "RGB"
                ]).to_csv(args.output, sep="\t", index=False)
                return
            
            # 准备参数
            process_args = [(name, group_df, args.distinct_colors) 
                          for name, group_df in groups_filtered]
            
            # 使用进程池并行处理
            n_processes = args.processes or cpu_count()
            with Pool(processes=n_processes) as pool:
                results = list(tqdm(
                    pool.imap_unordered(process_group, process_args), 
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
            
            # 保存结果
            output_df.to_csv(args.output, sep="\t", index=False)
            logger.info(f"✅ 完成！输出保存到: {args.output}")
            
    except Exception as e:
        logger.error(f"处理过程中出错: {str(e)}")
        raise

if __name__ == "__main__":
    main()
