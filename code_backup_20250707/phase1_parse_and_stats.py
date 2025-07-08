#!/usr/bin/env python3
"""
Phase 1: eccDNA转座子数据解析和基础统计
"""

import os
import re
import pickle
import pandas as pd
from pathlib import Path
from collections import defaultdict
import argparse
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class GFFParser:
    """GFF文件解析器"""
    
    def __init__(self):
        # Target字段的正则表达式
        self.target_pattern = re.compile(r'Target "Motif:(\S+)" (\d+) (\d+)')
        
    def parse_file(self, filepath):
        """解析单个GFF文件"""
        records = []
        
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                # 解析基本字段
                record = {
                    'seq_id': parts[0],
                    'source': parts[1],
                    'feature': parts[2],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'score': float(parts[5]) if parts[5] != '.' else 0.0,
                    'strand': parts[6],
                    'frame': parts[7]
                }
                
                # 解析attributes字段中的Target信息
                attributes = parts[8]
                target_match = self.target_pattern.search(attributes)
                
                if target_match:
                    record['motif'] = target_match.group(1)
                    record['target_start'] = int(target_match.group(2))
                    record['target_end'] = int(target_match.group(3))
                    record['target_length'] = record['target_end'] - record['target_start'] + 1
                    record['match_length'] = record['end'] - record['start'] + 1
                else:
                    continue
                
                records.append(record)
                    
        return records


class DataProcessor:
    """数据处理和统计"""
    
    def __init__(self):
        self.parser = GFFParser()
        
    def process_all_files(self, input_dir):
        """处理所有GFF文件"""
        data = {}
        
        # 遍历所有GFF文件
        gff_files = list(Path(input_dir).glob('*.gff'))
        logger.info(f"找到 {len(gff_files)} 个GFF文件")
        
        for filepath in gff_files:
            filename = filepath.stem  # 去掉.gff后缀
            
            # 解析文件名，提取样本信息
            # 格式: AT1.Mecc, U87_1.Mecc, HeLa1.Uecc 等
            if '.Mecc' in filename:
                sample_name = filename.replace('.Mecc', '')
                sample_type = 'Mecc'
            elif '.Uecc' in filename:
                sample_name = filename.replace('.Uecc', '')
                sample_type = 'Uecc'
            else:
                logger.warning(f"无法识别文件类型: {filename}")
                continue
            
            # 提取样本组（AT, U87, HeLa, ZM）
            if sample_name.startswith('AT'):
                sample_group = 'AT'
            elif sample_name.startswith('U87'):
                sample_group = 'U87'
            elif sample_name.startswith('HeLa'):
                sample_group = 'HeLa'
            elif sample_name.startswith('ZM'):
                sample_group = 'ZM'
            else:
                sample_group = 'Unknown'
            
            logger.info(f"处理文件: {filepath.name} (样本: {sample_name}, 类型: {sample_type})")
            records = self.parser.parse_file(filepath)
            
            # 存储数据
            if sample_name not in data:
                data[sample_name] = {
                    'group': sample_group,
                    'Mecc': [],
                    'Uecc': []
                }
            
            data[sample_name][sample_type] = records
            logger.info(f"  - 提取到 {len(records)} 条记录")
        
        return data
    
    def calculate_basic_stats(self, data):
        """计算基础统计信息"""
        stats = []
        
        for sample_name, sample_data in data.items():
            for ecc_type in ['Mecc', 'Uecc']:
                records = sample_data.get(ecc_type, [])
                if not records:
                    continue
                
                df = pd.DataFrame(records)
                
                stat = {
                    'sample': sample_name,
                    'group': sample_data['group'],
                    'type': ecc_type,
                    'total_records': len(records),
                    'unique_eccdna': df['seq_id'].nunique(),
                    'unique_motifs': df['motif'].nunique(),
                    'avg_match_length': df['match_length'].mean(),
                    'median_match_length': df['match_length'].median(),
                    'avg_target_length': df['target_length'].mean(),
                    'median_target_length': df['target_length'].median(),
                    'avg_score': df['score'].mean(),
                    'total_bases_covered': df['match_length'].sum()
                }
                
                # 添加最常见的motif
                top_motifs = df['motif'].value_counts().head(5)
                stat['top_motifs'] = ','.join(top_motifs.index.tolist())
                stat['top_motif_counts'] = ','.join(top_motifs.values.astype(str).tolist())
                
                stats.append(stat)
        
        return pd.DataFrame(stats)
    
    def generate_sample_summary(self, data):
        """生成样本汇总信息"""
        summary = []
        
        for sample_name, sample_data in data.items():
            for ecc_type in ['Mecc', 'Uecc']:
                records = sample_data.get(ecc_type, [])
                if not records:
                    continue
                
                df = pd.DataFrame(records)
                
                # 统计每个eccDNA上的转座子数量
                transposon_counts = df.groupby('seq_id').size()
                
                # 计算覆盖度相关统计
                coverage_stats = df.groupby('seq_id').agg({
                    'match_length': 'sum',
                    'motif': 'count'
                })
                
                summary.append({
                    'sample': sample_name,
                    'group': sample_data['group'],
                    'type': ecc_type,
                    'total_eccdna': df['seq_id'].nunique(),
                    'total_transposons': len(df),
                    'single_transposon_eccdna': (transposon_counts == 1).sum(),
                    'multi_transposon_eccdna': (transposon_counts > 1).sum(),
                    'max_transposons_per_eccdna': transposon_counts.max(),
                    'avg_transposons_per_eccdna': transposon_counts.mean(),
                    'total_unique_motifs': df['motif'].nunique(),
                    'avg_coverage_per_eccdna': coverage_stats['match_length'].mean()
                })
        
        return pd.DataFrame(summary)
    
    def generate_motif_summary(self, data):
        """生成转座子motif汇总"""
        motif_stats = []
        
        for sample_name, sample_data in data.items():
            for ecc_type in ['Mecc', 'Uecc']:
                records = sample_data.get(ecc_type, [])
                if not records:
                    continue
                
                df = pd.DataFrame(records)
                
                # 按motif统计
                motif_groups = df.groupby('motif').agg({
                    'seq_id': 'count',  # 出现次数
                    'match_length': ['sum', 'mean', 'median'],
                    'target_length': ['mean', 'median']
                }).round(2)
                
                for motif, stats in motif_groups.iterrows():
                    motif_stats.append({
                        'sample': sample_name,
                        'group': sample_data['group'],
                        'type': ecc_type,
                        'motif': motif,
                        'count': int(stats[('seq_id', 'count')]),
                        'total_match_length': int(stats[('match_length', 'sum')]),
                        'avg_match_length': stats[('match_length', 'mean')],
                        'median_match_length': stats[('match_length', 'median')],
                        'avg_target_length': stats[('target_length', 'mean')],
                        'median_target_length': stats[('target_length', 'median')]
                    })
        
        return pd.DataFrame(motif_stats)


def main():
    parser = argparse.ArgumentParser(description='Phase 1: eccDNA转座子数据解析和基础统计')
    parser.add_argument('--input_dir', required=True, help='GFF文件目录')
    parser.add_argument('--output_dir', default='./output/phase1', help='输出目录')
    args = parser.parse_args()
    
    # 创建输出目录
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 初始化处理器
    processor = DataProcessor()
    
    # 1. 解析所有文件
    logger.info("开始解析GFF文件...")
    data = processor.process_all_files(args.input_dir)
    
    # 2. 保存处理后的完整数据
    with open(output_dir / 'processed_data.pkl', 'wb') as f:
        pickle.dump(data, f)
    logger.info("已保存处理后的数据")
    
    # 3. 计算基础统计
    logger.info("计算基础统计...")
    basic_stats = processor.calculate_basic_stats(data)
    basic_stats.to_csv(output_dir / 'basic_stats.csv', index=False)
    
    # 4. 生成样本汇总
    sample_summary = processor.generate_sample_summary(data)
    sample_summary.to_csv(output_dir / 'sample_summary.csv', index=False)
    
    # 5. 生成motif汇总
    logger.info("生成motif统计...")
    motif_summary = processor.generate_motif_summary(data)
    motif_summary.to_csv(output_dir / 'motif_summary.csv', index=False)
    
    # 6. 生成数据质量报告
    with open(output_dir / 'data_summary.txt', 'w') as f:
        f.write("eccDNA转座子数据汇总报告\n")
        f.write("=" * 50 + "\n\n")
        
        # 总体统计
        f.write(f"处理文件数: {len(basic_stats)}\n")
        f.write(f"样本数: {len(data)}\n")
        f.write(f"总记录数: {basic_stats['total_records'].sum():,}\n")
        f.write(f"Mecc记录数: {basic_stats[basic_stats['type']=='Mecc']['total_records'].sum():,}\n")
        f.write(f"Uecc记录数: {basic_stats[basic_stats['type']=='Uecc']['total_records'].sum():,}\n\n")
        
        # 按组统计
        f.write("各组样本统计:\n")
        for group in ['AT', 'U87', 'HeLa', 'ZM']:
            group_data = basic_stats[basic_stats['group'] == group]
            if not group_data.empty:
                f.write(f"\n{group}组:\n")
                f.write(f"  样本数: {group_data['sample'].nunique()}\n")
                f.write(f"  总记录数: {group_data['total_records'].sum():,}\n")
                f.write(f"  Mecc: {group_data[group_data['type']=='Mecc']['total_records'].sum():,}\n")
                f.write(f"  Uecc: {group_data[group_data['type']=='Uecc']['total_records'].sum():,}\n")
        
        # 输出文件列表
        f.write(f"\n\n输出文件:\n")
        f.write(f"- processed_data.pkl: 完整数据\n")
        f.write(f"- basic_stats.csv: 基础统计表\n")
        f.write(f"- sample_summary.csv: 样本汇总\n")
        f.write(f"- motif_summary.csv: Motif统计\n")
    
    logger.info(f"分析完成！结果保存在 {output_dir}")
    
    # 打印简要统计
    print("\n=== 数据概览 ===")
    print(f"总样本数: {len(data)}")
    print(f"总记录数: {basic_stats['total_records'].sum():,}")
    print(f"\n各组记录数:")
    group_stats = basic_stats.groupby(['group', 'type'])['total_records'].sum().unstack(fill_value=0)
    print(group_stats)


if __name__ == '__main__':
    main()
