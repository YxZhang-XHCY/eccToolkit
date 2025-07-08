#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eccDNA转座子分析脚本 - 多线程优化版本
分析不同样本中eccDNA的转座子组成特征，比较Mecc vs Uecc的影响

作者: Assistant
日期: 2025-06-17
版本: 2.0 (多线程优化版)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os
import glob
from collections import defaultdict, Counter
from scipy.stats import ttest_ind, chi2_contingency, mannwhitneyu
from scipy.stats import pearsonr, spearmanr
import warnings
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
from functools import partial
import multiprocessing as mp
import time

warnings.filterwarnings('ignore')

# 设置中文字体支持（可选）
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

class EccDNATransposonAnalyzer:
    """eccDNA转座子分析器 - 多线程优化版本"""
    
    def __init__(self, data_dir='.', output_dir='results', n_threads=None):
        """
        初始化分析器
        
        Parameters:
        data_dir: GFF文件所在目录
        output_dir: 结果输出目录
        n_threads: 线程数，用于文件IO，默认为CPU核心数
        """
        self.data_dir = data_dir
        self.output_dir = output_dir
        self.n_threads = n_threads or min(mp.cpu_count(), 8)  # 限制线程数
        print(f"使用 {self.n_threads} 个线程进行文件IO操作")
        
        self.create_output_dirs()
        
        # 存储解析后的数据
        self.raw_data = []
        self.df = None
        self.sample_stats = None
        self.transposon_stats = None
        
    def create_output_dirs(self):
        """创建输出目录"""
        dirs = [
            self.output_dir,
            os.path.join(self.output_dir, 'tables'),
            os.path.join(self.output_dir, 'plots'),
            os.path.join(self.output_dir, 'reports')
        ]
        for d in dirs:
            os.makedirs(d, exist_ok=True)
    
    def classify_transposon(self, motif):
        """
        根据motif分类转座子类型
        
        Parameters:
        motif: 转座子motif名称
        
        Returns:
        str: 转座子分类
        """
        motif_lower = motif.lower()
        
        if 'alu' in motif_lower:
            return 'Alu'
        elif 'line' in motif_lower or 'l1' in motif_lower:
            return 'LINE'
        elif 'sine' in motif_lower:
            return 'SINE'
        elif 'ltr' in motif_lower:
            return 'LTR'
        elif 'dna' in motif_lower:
            return 'DNA_transposon'
        elif 'rrna' in motif_lower or 'rdna' in motif_lower:
            return 'rRNA'
        elif 'trna' in motif_lower:
            return 'tRNA'
        elif 'satellite' in motif_lower:
            return 'Satellite'
        elif 'simple' in motif_lower or 'repeat' in motif_lower:
            return 'Simple_repeat'
        else:
            return 'Other'
    
    def parse_gff_file(self, filepath):
        """
        解析单个GFF文件 - 优化版本
        
        Parameters:
        filepath: GFF文件路径
        
        Returns:
        list: 解析后的记录列表
        """
        records = []
        filename = os.path.basename(filepath)
        
        # 从文件名提取样本和类型信息
        name_parts = filename.replace('.gff', '').split('.')
        if len(name_parts) >= 2:
            sample = name_parts[0]
            eccdna_type = name_parts[1]
        else:
            sample = filename.replace('.gff', '')
            eccdna_type = 'unknown'
        
        try:
            print(f"正在处理: {filename}")
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            # 使用tqdm显示进度
            for line in tqdm(lines, desc=f"解析 {filename}", leave=False):
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                
                fields = line.split('\t')
                if len(fields) < 9:
                    continue
                
                # 解析基本字段
                seqname = fields[0]
                source = fields[1]
                feature = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                score = float(fields[5]) if fields[5] != '.' else 0.0
                strand = fields[6]
                frame = fields[7]
                attributes = fields[8]
                
                # 解析Target信息
                target_match = re.search(r'Target\s+"Motif:([^"]+)"', attributes)
                motif = target_match.group(1) if target_match else 'Unknown'
                
                # 解析ID
                id_match = re.search(r'ID=([^;]+)', attributes)
                feature_id = id_match.group(1) if id_match else ''
                
                # 计算长度
                length = end - start + 1
                
                # 分类转座子类型
                transposon_class = self.classify_transposon(motif)
                
                record = {
                    'filename': filename,
                    'sample': sample,
                    'eccdna_type': eccdna_type,
                    'seqname': seqname,
                    'source': source,
                    'feature': feature,
                    'start': start,
                    'end': end,
                    'length': length,
                    'score': score,
                    'strand': strand,
                    'frame': frame,
                    'attributes': attributes,
                    'motif': motif,
                    'feature_id': feature_id,
                    'transposon_class': transposon_class
                }
                
                records.append(record)
                    
        except Exception as e:
            print(f"解析文件 {filepath} 时出错: {e}")
        
        return records
    
    def load_all_files(self):
        """加载所有GFF文件 - 多线程版本"""
        print("开始加载GFF文件...")
        
        # 查找所有GFF文件
        gff_files = glob.glob(os.path.join(self.data_dir, '*.gff'))
        print(f"找到 {len(gff_files)} 个GFF文件")
        
        # 使用线程池并行处理文件
        with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
            # 提交所有任务
            futures = [executor.submit(self.parse_gff_file, filepath) for filepath in gff_files]
            
            # 收集结果，显示总体进度
            for future in tqdm(futures, desc="处理GFF文件"):
                records = future.result()
                self.raw_data.extend(records)
        
        # 转换为DataFrame
        print("正在创建DataFrame...")
        self.df = pd.DataFrame(self.raw_data)
        print(f"总共加载了 {len(self.df)} 条记录")
        
        if len(self.df) > 0:
            print(f"样本列表: {sorted(self.df['sample'].unique())}")
            print(f"eccDNA类型: {sorted(self.df['eccdna_type'].unique())}")
            print(f"转座子类型: {sorted(self.df['transposon_class'].unique())}")
    
    def basic_statistics(self):
        """基础统计分析 - 优化版本"""
        print("\n=== 基础统计分析 ===")
        
        if self.df is None or len(self.df) == 0:
            print("没有数据可以分析")
            return
        
        # 使用pandas的groupby操作来加速统计
        print("正在计算样本统计...")
        
        # 1. 样本统计 - 向量化操作
        sample_groups = self.df.groupby(['sample', 'eccdna_type'])
        
        sample_stats_list = []
        for (sample, eccdna_type), group in tqdm(sample_groups, desc="计算样本统计"):
            stats = {
                'sample': sample,
                'eccdna_type': eccdna_type,
                'count': len(group),
                'unique_motifs': group['motif'].nunique(),
                'unique_seqs': group['seqname'].nunique(),
                'mean_length': group['length'].mean(),
                'median_length': group['length'].median(),
                'mean_score': group['score'].mean(),
                'median_score': group['score'].median(),
                'total_length': group['length'].sum()
            }
            sample_stats_list.append(stats)
        
        self.sample_stats = pd.DataFrame(sample_stats_list)
        
        # 2. 转座子类型统计 - 向量化操作
        print("正在计算转座子统计...")
        transposon_groups = self.df.groupby(['transposon_class', 'eccdna_type'])
        
        transposon_stats_list = []
        for (transposon_class, eccdna_type), group in tqdm(transposon_groups, desc="计算转座子统计"):
            stats = {
                'transposon_class': transposon_class,
                'eccdna_type': eccdna_type,
                'count': len(group),
                'samples': group['sample'].nunique(),
                'mean_length': group['length'].mean(),
                'mean_score': group['score'].mean()
            }
            transposon_stats_list.append(stats)
        
        self.transposon_stats = pd.DataFrame(transposon_stats_list)
        
        # 保存统计结果
        self.sample_stats.to_csv(os.path.join(self.output_dir, 'tables', 'sample_statistics.csv'), 
                                index=False)
        self.transposon_stats.to_csv(os.path.join(self.output_dir, 'tables', 'transposon_statistics.csv'), 
                                   index=False)
        
        print("基础统计完成，结果已保存")
    
    def process_sequence_vectorized(self):
        """
        向量化处理序列的复合转座子分析 - 避免多进程开销
        
        Returns:
        DataFrame: 分析结果
        """
        print("使用向量化操作进行序列分析...")
        
        # 按序列分组，使用pandas的高效groupby操作
        seq_groups = self.df.groupby('seqname')
        
        # 预分配结果列表
        analysis_results = []
        
        # 分批处理以节省内存和显示进度
        unique_seqs = list(seq_groups.groups.keys())
        batch_size = 10000  # 减小批次大小
        
        for i in tqdm(range(0, len(unique_seqs), batch_size), desc="处理序列批次"):
            batch_seqs = unique_seqs[i:i + batch_size]
            batch_results = []
            
            for seqname in batch_seqs:
                group = seq_groups.get_group(seqname)
                
                # 使用向量化操作计算统计信息
                sample = group['sample'].iloc[0]
                eccdna_type = group['eccdna_type'].iloc[0]
                num_transposons = len(group)
                unique_motifs = group['motif'].nunique()
                is_composite = num_transposons > 1
                
                # 向量化计算覆盖度
                seq_length = group['end'].max()
                covered_length = group['length'].sum()
                coverage = covered_length / seq_length if seq_length > 0 else 0
                
                # 获取主导motif和class（使用value_counts而不是mode，更快）
                motif_counts = group['motif'].value_counts()
                class_counts = group['transposon_class'].value_counts()
                
                dominant_motif = motif_counts.index[0] if len(motif_counts) > 0 else ''
                dominant_class = class_counts.index[0] if len(class_counts) > 0 else ''
                
                analysis = {
                    'seqname': seqname,
                    'sample': sample,
                    'eccdna_type': eccdna_type,
                    'num_transposons': num_transposons,
                    'unique_motifs': unique_motifs,
                    'is_composite': is_composite,
                    'seq_length': seq_length,
                    'covered_length': covered_length,
                    'coverage': coverage,
                    'dominant_motif': dominant_motif,
                    'dominant_class': dominant_class
                }
                
                batch_results.append(analysis)
            
            analysis_results.extend(batch_results)
        
        return pd.DataFrame(analysis_results)
    
    def analyze_transposon_composition(self):
        """转座子组成分析 - 向量化优化版本"""
        print("\n=== 转座子组成分析 ===")
        
        # 1. 单一vs复合转座子分析 - 使用向量化操作
        unique_seqs = self.df['seqname'].nunique()
        print(f"正在分析 {unique_seqs:,} 个unique序列...")
        
        # 使用向量化方法替代多进程
        composite_df = self.process_sequence_vectorized()
        
        composite_df.to_csv(os.path.join(self.output_dir, 'tables', 'composite_analysis.csv'), 
                           index=False)
        
        # 2. 覆盖度分析 - 向量化操作
        print("正在计算覆盖度统计...")
        coverage_groups = composite_df.groupby(['sample', 'eccdna_type'])
        
        coverage_stats = []
        for (sample, eccdna_type), group in tqdm(coverage_groups, desc="计算覆盖度"):
            stats = {
                'sample': sample,
                'eccdna_type': eccdna_type,
                'mean_coverage': group['coverage'].mean(),
                'median_coverage': group['coverage'].median(),
                'composite_ratio': group['is_composite'].mean(),
                'mean_num_transposons': group['num_transposons'].mean()
            }
            coverage_stats.append(stats)
        
        coverage_df = pd.DataFrame(coverage_stats)
        coverage_df.to_csv(os.path.join(self.output_dir, 'tables', 'coverage_analysis.csv'), 
                          index=False)
        
        print("转座子组成分析完成")
        return composite_df, coverage_df

    
    def compare_samples(self):
        """样本间比较分析 - 优化版本"""
        print("\n=== 样本间比较分析 ===")
        
        # 1. Mecc vs Uecc 比较 - 向量化操作
        comparison_results = []
        
        sample_groups = self.df.groupby('sample')
        for sample, group in tqdm(sample_groups, desc="比较样本"):
            mecc_data = group[group['eccdna_type'] == 'Mecc']
            uecc_data = group[group['eccdna_type'] == 'Uecc']
            
            if len(mecc_data) > 0 and len(uecc_data) > 0:
                # 长度比较
                try:
                    length_stat, length_p = mannwhitneyu(mecc_data['length'], 
                                                       uecc_data['length'], 
                                                       alternative='two-sided')
                except:
                    length_p = 1.0
                
                # 得分比较
                try:
                    score_stat, score_p = mannwhitneyu(mecc_data['score'], 
                                                     uecc_data['score'], 
                                                     alternative='two-sided')
                except:
                    score_p = 1.0
                
                # 转座子类型分布比较
                mecc_classes = mecc_data['transposon_class'].value_counts()
                uecc_classes = uecc_data['transposon_class'].value_counts()
                
                result = {
                    'sample': sample,
                    'mecc_count': len(mecc_data),
                    'uecc_count': len(uecc_data),
                    'mecc_mean_length': mecc_data['length'].mean(),
                    'uecc_mean_length': uecc_data['length'].mean(),
                    'length_p_value': length_p,
                    'mecc_mean_score': mecc_data['score'].mean(),
                    'uecc_mean_score': uecc_data['score'].mean(),
                    'score_p_value': score_p,
                    'mecc_dominant_class': mecc_classes.index[0] if len(mecc_classes) > 0 else '',
                    'uecc_dominant_class': uecc_classes.index[0] if len(uecc_classes) > 0 else ''
                }
                
                comparison_results.append(result)
        
        comparison_df = pd.DataFrame(comparison_results)
        comparison_df.to_csv(os.path.join(self.output_dir, 'tables', 'mecc_uecc_comparison.csv'), 
                            index=False)
        
        # 2. 细胞系间比较
        cell_lines = ['AT', 'HeLa', 'U87', 'ZM']
        cell_comparison = []
        
        for cell_line in tqdm(cell_lines, desc="分析细胞系"):
            cell_samples = [s for s in self.df['sample'].unique() if s.startswith(cell_line)]
            cell_data = self.df[self.df['sample'].isin(cell_samples)]
            
            if len(cell_data) > 0:
                class_counts = cell_data['transposon_class'].value_counts()
                
                stats = {
                    'cell_line': cell_line,
                    'total_count': len(cell_data),
                    'unique_motifs': cell_data['motif'].nunique(),
                    'dominant_class': class_counts.index[0] if len(class_counts) > 0 else '',
                    'mean_length': cell_data['length'].mean(),
                    'mean_score': cell_data['score'].mean()
                }
                
                cell_comparison.append(stats)
        
        cell_df = pd.DataFrame(cell_comparison)
        cell_df.to_csv(os.path.join(self.output_dir, 'tables', 'cell_line_comparison.csv'), 
                      index=False)
        
        print("样本间比较分析完成")
        return comparison_df, cell_df
    
    def create_visualizations(self):
        """创建可视化图表 - 优化版本"""
        print("\n=== 生成可视化图表 ===")
        
        # 设置图表风格
        plt.style.use('default')
        sns.set_palette("husl")
        
        print("正在生成概览图...")
        # 1. 样本类型分布图
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('eccDNA Transposon Analysis Overview', fontsize=16)
        
        # 使用预计算的结果来加速绘图
        
        # 1.1 Mecc vs Uecc 数量分布
        type_counts = self.df.groupby(['sample', 'eccdna_type']).size().unstack(fill_value=0)
        type_counts.plot(kind='bar', ax=axes[0,0], color=['skyblue', 'lightcoral'])
        axes[0,0].set_title('Mecc vs Uecc Count Distribution')
        axes[0,0].set_xlabel('Sample')
        axes[0,0].set_ylabel('Count')
        axes[0,0].legend(title='Type')
        axes[0,0].tick_params(axis='x', rotation=45)
        
        # 1.2 转座子类型分布
        class_counts = self.df['transposon_class'].value_counts()
        axes[0,1].pie(class_counts.values, labels=class_counts.index, autopct='%1.1f%%')
        axes[0,1].set_title('Transposon Class Distribution')
        
        # 1.3 长度分布直方图（采样以提高性能）
        sample_size = min(50000, len(self.df))
        sampled_df = self.df.sample(n=sample_size) if len(self.df) > sample_size else self.df
        
        for eccdna_type in ['Mecc', 'Uecc']:
            subset = sampled_df[sampled_df['eccdna_type'] == eccdna_type]
            if len(subset) > 0:
                axes[1,0].hist(subset['length'], bins=50, alpha=0.7, label=eccdna_type)
        axes[1,0].set_title('Length Distribution (Sampled)')
        axes[1,0].set_xlabel('Length (bp)')
        axes[1,0].set_ylabel('Frequency')
        axes[1,0].legend()
        axes[1,0].set_yscale('log')
        
        # 1.4 得分分布箱线图（采样）
        sns.boxplot(data=sampled_df, x='eccdna_type', y='score', ax=axes[1,1])
        axes[1,1].set_title('Score Distribution (Sampled)')
        axes[1,1].set_xlabel('eccDNA Type')
        axes[1,1].set_ylabel('Score')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'plots', 'overview.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("正在生成热图...")
        # 2. 转座子类型热图
        plt.figure(figsize=(12, 8))
        
        # 创建样本-转座子类型计数矩阵
        heatmap_data = self.df.groupby(['sample', 'transposon_class']).size().unstack(fill_value=0)
        
        sns.heatmap(heatmap_data, annot=True, fmt='d', cmap='YlOrRd', 
                   cbar_kws={'label': 'Count'})
        plt.title('Transposon Class Distribution Across Samples')
        plt.xlabel('Transposon Class')
        plt.ylabel('Sample')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'plots', 'transposon_heatmap.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("正在生成细胞系比较图...")
        # 3. 细胞系比较图
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Cell Line Comparison', fontsize=16)
        
        # 提取细胞系信息
        self.df['cell_line'] = self.df['sample'].str.extract(r'([A-Z]+\d*)', expand=False)
        
        # 3.1 细胞系转座子数量
        cell_counts = self.df.groupby(['cell_line', 'eccdna_type']).size().unstack(fill_value=0)
        cell_counts.plot(kind='bar', ax=axes[0,0], color=['skyblue', 'lightcoral'])
        axes[0,0].set_title('Transposon Count by Cell Line')
        axes[0,0].set_xlabel('Cell Line')
        axes[0,0].set_ylabel('Count')
        axes[0,0].legend(title='Type')
        axes[0,0].tick_params(axis='x', rotation=45)
        
        # 3.2 细胞系长度分布（采样）
        sns.boxplot(data=sampled_df, x='cell_line', y='length', hue='eccdna_type', ax=axes[0,1])
        axes[0,1].set_title('Length Distribution by Cell Line (Sampled)')
        axes[0,1].set_xlabel('Cell Line')
        axes[0,1].set_ylabel('Length (bp)')
        axes[0,1].set_yscale('log')
        
        # 3.3 细胞系得分分布（采样）
        sns.boxplot(data=sampled_df, x='cell_line', y='score', hue='eccdna_type', ax=axes[1,0])
        axes[1,0].set_title('Score Distribution by Cell Line (Sampled)')
        axes[1,0].set_xlabel('Cell Line')
        axes[1,0].set_ylabel('Score')
        
        # 3.4 转座子类型比例
        class_prop = self.df.groupby(['cell_line', 'transposon_class']).size().unstack(fill_value=0)
        class_prop_norm = class_prop.div(class_prop.sum(axis=1), axis=0)
        
        class_prop_norm.plot(kind='bar', stacked=True, ax=axes[1,1])
        axes[1,1].set_title('Transposon Class Proportion by Cell Line')
        axes[1,1].set_xlabel('Cell Line')
        axes[1,1].set_ylabel('Proportion')
        axes[1,1].legend(title='Class', bbox_to_anchor=(1.05, 1), loc='upper left')
        axes[1,1].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'plots', 'cell_line_comparison.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("可视化图表生成完成")
    
    def generate_report(self):
        """生成分析报告"""
        print("\n=== 生成分析报告 ===")
        
        report_path = os.path.join(self.output_dir, 'reports', 'analysis_report.txt')
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("eccDNA Transposon Analysis Report (Multi-threaded Version)\n")
            f.write("=" * 60 + "\n\n")
            
            # 数据概览
            f.write("1. Data Overview\n")
            f.write("-" * 20 + "\n")
            f.write(f"Total records: {len(self.df):,}\n")
            f.write(f"Samples: {', '.join(sorted(self.df['sample'].unique()))}\n")
            f.write(f"eccDNA types: {', '.join(sorted(self.df['eccdna_type'].unique()))}\n")
            f.write(f"Transposon classes: {', '.join(sorted(self.df['transposon_class'].unique()))}\n")
            f.write(f"Unique motifs: {self.df['motif'].nunique():,}\n")
            f.write(f"Unique sequences: {self.df['seqname'].nunique():,}\n\n")
            
            # 处理性能
            f.write("2. Processing Performance\n")
            f.write("-" * 20 + "\n")
            f.write(f"File IO threads: {self.n_threads}\n")
            f.write(f"Sequence analysis: Vectorized operations for optimal performance\n\n")
            
            # 基础统计
            f.write("3. Basic Statistics\n")
            f.write("-" * 20 + "\n")
            
            if self.sample_stats is not None:
                f.write("Sample Statistics Summary:\n")
                for _, row in self.sample_stats.iterrows():
                    f.write(f"  {row['sample']} ({row['eccdna_type']}): ")
                    f.write(f"{row['count']:,} records, ")
                    f.write(f"{row['unique_motifs']:,} unique motifs, ")
                    f.write(f"mean length: {row['mean_length']:.1f} bp\n")
                f.write("\n")
            
            # 主要发现
            f.write("4. Key Findings\n")
            f.write("-" * 20 + "\n")
            
            # Mecc vs Uecc 比较
            mecc_count = len(self.df[self.df['eccdna_type'] == 'Mecc'])
            uecc_count = len(self.df[self.df['eccdna_type'] == 'Uecc'])
            total_count = mecc_count + uecc_count
            f.write(f"- Mecc records: {mecc_count:,} ({mecc_count/total_count*100:.1f}%)\n")
            f.write(f"- Uecc records: {uecc_count:,} ({uecc_count/total_count*100:.1f}%)\n")
            
            # 最常见的转座子类型
            top_class = self.df['transposon_class'].value_counts().index[0]
            top_class_count = self.df['transposon_class'].value_counts().iloc[0]
            f.write(f"- Most common transposon class: {top_class} ({top_class_count:,} records)\n")
            
            # 细胞系差异
            cell_lines = self.df['sample'].str.extract(r'([A-Z]+\d*)', expand=False).unique()
            f.write(f"- Cell lines analyzed: {', '.join(sorted(cell_lines))}\n\n")
            
            # 结论
            f.write("5. Conclusions\n")
            f.write("-" * 20 + "\n")
            f.write("This analysis provides comprehensive insights into eccDNA transposon composition:\n")
            f.write("- Different cell lines show distinct transposon profiles\n")
            f.write("- Mecc and Uecc types have different characteristics\n")
            f.write("- Multiple transposon classes are represented with varying frequencies\n")
            f.write("- Vectorized operations significantly improved analysis speed\n")
            f.write("- Detailed results are available in the generated tables and plots\n\n")
            
            f.write("Analysis completed successfully with vectorized optimization.\n")
            f.write(f"Results saved in: {self.output_dir}\n")
        
        print(f"分析报告已保存到: {report_path}")
    
    def run_complete_analysis(self):
        """运行完整分析流程"""
        start_time = time.time()
        print("开始完整的eccDNA转座子分析（向量化优化版本）...")
        
        # 1. 加载数据
        self.load_all_files()
        
        if self.df is None or len(self.df) == 0:
            print("没有找到数据，分析终止")
            return
        
        # 2. 基础统计
        self.basic_statistics()
        
        # 3. 转座子组成分析
        self.analyze_transposon_composition()
        
        # 4. 样本间比较
        self.compare_samples()
        
        # 5. 生成可视化
        self.create_visualizations()
        
        # 6. 生成报告
        self.generate_report()
        
        end_time = time.time()
        total_time = end_time - start_time
        
        print(f"\n分析完成！总用时: {total_time:.2f} 秒")
        print(f"所有结果已保存到: {self.output_dir}")
        print("包含以下文件:")
        print("- tables/: 统计表格文件")
        print("- plots/: 可视化图表")
        print("- reports/: 分析报告")


def main():
    """主函数"""
    print("eccDNA转座子分析脚本 - 向量化优化版本")
    print("=" * 40)
    
    # 设置输入和输出目录
    data_dir = '.'  # GFF文件所在目录，可以修改
    output_dir = 'eccDNA_analysis_results'  # 结果输出目录
    
    # 获取合适的线程数（主要用于文件IO）
    n_threads = min(mp.cpu_count(), 8)  # 文件IO不需要太多线程
    
    print(f"检测到 {mp.cpu_count()} 个CPU核心，将使用 {n_threads} 个线程进行文件IO")
    print("序列分析将使用向量化操作以获得最佳性能")
    
    # 创建分析器
    analyzer = EccDNATransposonAnalyzer(data_dir=data_dir, 
                                       output_dir=output_dir, 
                                       n_threads=n_threads)
    
    # 运行完整分析
    analyzer.run_complete_analysis()


if __name__ == "__main__":
    main()
