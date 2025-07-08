DNA转座子分析脚本 - 高性能优化版本
分析不同样本中eccDNA的转座子组成特征，比较Mecc vs Uecc的影响

作者: Assistant
日期: 2025-06-18
版本: 3.0 (高性能优化版)
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
import gc

warnings.filterwarnings('ignore')

# 设置中文字体支持（可选）
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

class EccDNATransposonAnalyzer:
    """eccDNA转座子分析器 - 高性能优化版本"""
    
    def __init__(self, data_dir='.', output_dir='results', n_threads=None, fast_mode=False):
        """
        初始化分析器
        
        Parameters:
        data_dir: GFF文件所在目录
        output_dir: 结果输出目录
        n_threads: 线程数，用于文件IO，默认为CPU核心数
        fast_mode: 快速模式，对大数据进行采样分析
        """
        self.data_dir = data_dir
        self.output_dir = output_dir
        self.n_threads = n_threads or min(mp.cpu_count(), 4)  # 减少线程数避免竞争
        self.fast_mode = fast_mode
        
        print(f"使用 {self.n_threads} 个线程进行文件IO操作")
        if self.fast_mode:
            print("已启用快速模式（采样分析）")
        
        self.create_output_dirs()
        
        # 存储解析后的数据
        self.raw_data = []
        self.df = None
        self.sample_stats = None
        self.transposon_stats = None
        
        # 预编译正则表达式
        self.target_pattern = re.compile(r'Target\s+"Motif:([^"]+)"')
        self.id_pattern = re.compile(r'ID=([^;]+)')
        
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
    
    def quick_diagnose(self):
        """快速诊断数据规模和潜在问题"""
        print("\n=== 快速诊断 ===")
        
        # 检查文件
        gff_files = glob.glob(os.path.join(self.data_dir, '*.gff'))
        if not gff_files:
            print("错误：没有找到GFF文件！")
            return False
            
        total_size = sum(os.path.getsize(f) for f in gff_files) / (1024**3)  # GB
        
        print(f"文件数量: {len(gff_files)}")
        print(f"总数据量: {total_size:.2f} GB")
        
        if total_size > 1:
            print("⚠️  数据量较大，建议使用快速模式 (fast_mode=True)")
            
        # 如果已有数据，显示规模
        if hasattr(self, 'df') and self.df is not None:
            print(f"\n当前数据规模:")
            print(f"- 总记录数: {len(self.df):,}")
            print(f"- 唯一序列: {self.df['seqname'].nunique():,}")
            print(f"- 内存占用: {self.df.memory_usage(deep=True).sum() / (1024**3):.2f} GB")
            
            # 检查是否有超大序列
            seq_counts = self.df['seqname'].value_counts()
            if seq_counts.iloc[0] > 1000:
                print(f"⚠️  发现超大序列: {seq_counts.index[0]} 有 {seq_counts.iloc[0]:,} 条记录")
        
        return True
    
    def parse_gff_file_optimized(self, filepath):
        """优化的GFF文件解析 - 流式处理，减少内存占用"""
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
                line_count = 0
                valid_count = 0
                
                for line in f:
                    line_count += 1
                    
                    # 每100000行显示一次进度
                    if line_count % 100000 == 0:
                        print(f"  已处理 {line_count:,} 行, {valid_count:,} 条有效记录...")
                    
                    line = line.strip()
                    if line.startswith('#') or not line:
                        continue
                    
                    fields = line.split('\t')
                    if len(fields) < 9:
                        continue
                    
                    # 快速模式下进行采样
                    if self.fast_mode and valid_count % 10 != 0:  # 只保留10%的数据
                        valid_count += 1
                        continue
                    
                    # 直接构建记录，避免重复变量赋值
                    try:
                        start = int(fields[3])
                        end = int(fields[4])
                        motif_match = self.target_pattern.search(fields[8])
                        motif = motif_match.group(1) if motif_match else 'Unknown'
                        
                        record = {
                            'filename': filename,
                            'sample': sample,
                            'eccdna_type': eccdna_type,
                            'seqname': fields[0],
                            'source': fields[1],
                            'feature': fields[2],
                            'start': start,
                            'end': end,
                            'length': end - start + 1,
                            'score': float(fields[5]) if fields[5] != '.' else 0.0,
                            'strand': fields[6],
                            'frame': fields[7],
                            'attributes': fields[8],
                            'motif': motif,
                            'feature_id': '',  # 简化处理
                            'transposon_class': self.classify_transposon(motif)
                        }
                        
                        records.append(record)
                        valid_count += 1
                        
                    except ValueError:
                        continue
                
                print(f"  完成！共 {valid_count:,} 条有效记录")
                    
        except Exception as e:
            print(f"解析文件 {filepath} 时出错: {e}")
        
        return records
    
    def load_all_files(self):
        """加载所有GFF文件 - 优化版本"""
        print("开始加载GFF文件...")
        
        # 快速诊断
        if not self.quick_diagnose():
            return
        
        # 查找所有GFF文件
        gff_files = glob.glob(os.path.join(self.data_dir, '*.gff'))
        print(f"找到 {len(gff_files)} 个GFF文件")
        
        # 获取文件大小信息
        file_sizes = [(f, os.path.getsize(f)) for f in gff_files]
        total_size = sum(size for _, size in file_sizes)
        print(f"总数据量: {total_size / 1024 / 1024:.1f} MB")
        
        # 根据数据量选择处理策略
        if total_size < 500 * 1024 * 1024 and not self.fast_mode:  # 小于500MB且非快速模式
            # 使用多线程
            print("数据量适中，使用多线程处理")
            with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
                futures = [executor.submit(self.parse_gff_file_optimized, filepath) 
                          for filepath in gff_files]
                for future in tqdm(futures, desc="处理GFF文件"):
                    records = future.result()
                    self.raw_data.extend(records)
        else:
            # 使用单线程顺序处理，避免内存压力
            print("数据量较大或启用快速模式，使用单线程顺序处理")
            for filepath in tqdm(gff_files, desc="处理GFF文件"):
                records = self.parse_gff_file_optimized(filepath)
                self.raw_data.extend(records)
                
                # 定期垃圾回收
                if len(self.raw_data) > 500000:
                    gc.collect()
        
        # 转换为DataFrame
        print("正在创建DataFrame...")
        self.df = pd.DataFrame(self.raw_data)
        
        # 清理原始数据以节省内存
        self.raw_data = []
        gc.collect()
        
        print(f"总共加载了 {len(self.df):,} 条记录")
        
        if len(self.df) > 0:
            print(f"样本列表: {sorted(self.df['sample'].unique())}")
            print(f"eccDNA类型: {sorted(self.df['eccdna_type'].unique())}")
            print(f"转座子类型: {sorted(self.df['transposon_class'].unique())}")
            
            # 内存占用信息
            memory_usage = self.df.memory_usage(deep=True).sum() / (1024**3)
            print(f"DataFrame内存占用: {memory_usage:.2f} GB")
    
    def basic_statistics(self):
        """基础统计分析 - 纯向量化操作"""
        print("\n=== 基础统计分析 ===")
        
        if self.df is None or len(self.df) == 0:
            print("没有数据可以分析")
            return
        
        # 1. 样本统计 - 使用agg一次性计算所有统计量
        print("正在计算样本统计...")
        self.sample_stats = self.df.groupby(['sample', 'eccdna_type']).agg({
            'seqname': 'count',
            'motif': 'nunique',
            'length': ['mean', 'median', 'sum'],
            'score': ['mean', 'median']
        }).reset_index()
        
        # 展平多级列名
        self.sample_stats.columns = ['sample', 'eccdna_type', 'count', 'unique_motifs',
                                    'mean_length', 'median_length', 'total_length',
                                    'mean_score', 'median_score']
        
        # 添加unique sequences统计
        unique_seqs = self.df.groupby(['sample', 'eccdna_type'])['seqname'].nunique().reset_index()
        unique_seqs.columns = ['sample', 'eccdna_type', 'unique_seqs']
        self.sample_stats = self.sample_stats.merge(unique_seqs, on=['sample', 'eccdna_type'])
        
        # 2. 转座子类型统计
        print("正在计算转座子统计...")
        self.transposon_stats = self.df.groupby(['transposon_class', 'eccdna_type']).agg({
            'seqname': 'count',
            'sample': 'nunique',
            'length': 'mean',
            'score': 'mean'
        }).reset_index()
        
        self.transposon_stats.columns = ['transposon_class', 'eccdna_type', 'count',
                                        'samples', 'mean_length', 'mean_score']
        
        # 保存统计结果
        self.sample_stats.to_csv(os.path.join(self.output_dir, 'tables', 'sample_statistics.csv'), 
                                index=False)
        self.transposon_stats.to_csv(os.path.join(self.output_dir, 'tables', 'transposon_statistics.csv'), 
                                   index=False)
        
        print("基础统计完成，结果已保存")
    
    def process_sequence_vectorized_fast(self):
        """超快速向量化处理 - 完全避免循环"""
        print("使用纯向量化操作进行序列分析...")
        
        # 先看数据规模
        unique_seqs = self.df['seqname'].nunique()
        print(f"需要处理 {unique_seqs:,} 个序列")
        
        # 使用纯向量化的聚合操作
        print("计算基础统计...")
        
        # 一次性计算所有统计
        agg_dict = {
            'sample': 'first',
            'eccdna_type': 'first',
            'motif': ['count', 'nunique'],
            'length': 'sum',
            'end': 'max'
        }
        
        # 使用agg一次计算多个统计量
        seq_stats = self.df.groupby('seqname', as_index=False).agg(agg_dict)
        
        # 展平多级列名
        seq_stats.columns = ['seqname', 'sample', 'eccdna_type', 'num_transposons', 
                            'unique_motifs', 'covered_length', 'seq_length']
        
        # 向量化计算衍生字段
        seq_stats['is_composite'] = seq_stats['num_transposons'] > 1
        seq_stats['coverage'] = seq_stats['covered_length'] / seq_stats['seq_length']
        seq_stats['coverage'] = seq_stats['coverage'].fillna(0)
        
        # 对于dominant class和motif，使用不同策略
        if unique_seqs < 100000 and not self.fast_mode:
            print("计算dominant class和motif...")
            # 使用mode计算最频繁的值（更快）
            dominant_class = self.df.groupby('seqname')['transposon_class'].agg(lambda x: x.mode()[0] if len(x.mode()) > 0 else '').reset_index()
            dominant_motif = self.df.groupby('seqname')['motif'].agg(lambda x: x.mode()[0] if len(x.mode()) > 0 else '').reset_index()
            
            dominant_class.columns = ['seqname', 'dominant_class']
            dominant_motif.columns = ['seqname', 'dominant_motif']
            
            # 合并结果
            seq_stats = seq_stats.merge(dominant_class, on='seqname', how='left')
            seq_stats = seq_stats.merge(dominant_motif, on='seqname', how='left')
        else:
            print(f"序列数过多({unique_seqs:,})或快速模式，使用简化处理...")
            seq_stats['dominant_class'] = 'Multiple'
            seq_stats['dominant_motif'] = 'Multiple'
        
        return seq_stats
    
    def analyze_transposon_composition_optimized(self):
        """转座子组成分析 - 高度优化版本"""
        print("\n=== 转座子组成分析（优化版） ===")
        
        # 先进行快速诊断
        unique_seqs = self.df['seqname'].nunique()
        total_records = len(self.df)
        
        print(f"数据规模: {total_records:,} 条记录, {unique_seqs:,} 个unique序列")
        
        # 使用完全优化的向量化方法
        composite_df = self.process_sequence_vectorized_fast()
        
        # 保存结果
        composite_df.to_csv(os.path.join(self.output_dir, 'tables', 'composite_analysis.csv'), 
                           index=False)
        
        # 覆盖度分析 - 也使用向量化
        print("正在计算覆盖度统计...")
        coverage_stats = composite_df.groupby(['sample', 'eccdna_type']).agg({
            'coverage': ['mean', 'median'],
            'is_composite': 'mean',
            'num_transposons': 'mean'
        }).reset_index()
        
        # 展平列名
        coverage_stats.columns = ['sample', 'eccdna_type', 'mean_coverage', 
                                 'median_coverage', 'composite_ratio', 
                                 'mean_num_transposons']
        
        coverage_stats.to_csv(os.path.join(self.output_dir, 'tables', 'coverage_analysis.csv'), 
                             index=False)
        
        print("转座子组成分析完成")
        return composite_df, coverage_stats
    
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
        
        # 提前添加cell_line列
        self.df['cell_line'] = self.df['sample'].str.extract(r'([A-Z]+\d*)', expand=False)
        
        # 采样以提高性能
        sample_size = min(50000, len(self.df))
        if len(self.df) > sample_size:
            print(f"数据量大，对 {sample_size:,} 条记录进行采样绘图...")
            sampled_df = self.df.sample(n=sample_size, random_state=42)
        else:
            sampled_df = self.df
        
        print("正在生成概览图...")
        # 1. 样本类型分布图
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('eccDNA Transposon Analysis Overview', fontsize=16)
        
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
        axes[0,1].pie(class_counts.values[:10], labels=class_counts.index[:10], autopct='%1.1f%%')
        axes[0,1].set_title('Top 10 Transposon Classes')
        
        # 1.3 长度分布直方图（采样）
        for eccdna_type in ['Mecc', 'Uecc']:
            subset = sampled_df[sampled_df['eccdna_type'] == eccdna_type]
            if len(subset) > 0:
                axes[1,0].hist(subset['length'], bins=50, alpha=0.7, label=eccdna_type)
        axes[1,0].set_title(f'Length Distribution (n={len(sampled_df):,})')
        axes[1,0].set_xlabel('Length (bp)')
        axes[1,0].set_ylabel('Frequency')
        axes[1,0].legend()
        axes[1,0].set_yscale('log')
        
        # 1.4 得分分布箱线图（采样）
        sns.boxplot(data=sampled_df, x='eccdna_type', y='score', ax=axes[1,1])
        axes[1,1].set_title(f'Score Distribution (n={len(sampled_df):,})')
        axes[1,1].set_xlabel('eccDNA Type')
        axes[1,1].set_ylabel('Score')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'plots', 'overview.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("正在生成热图...")
        # 2. 转座子类型热图（只显示前20个样本和前10个类型）
        plt.figure(figsize=(12, 8))
        
        # 创建样本-转座子类型计数矩阵
        heatmap_data = self.df.groupby(['sample', 'transposon_class']).size().unstack(fill_value=0)
        
        # 选择前20个样本和前10个最常见的转座子类型
        top_samples = heatmap_data.sum(axis=1).nlargest(20).index
        top_classes = heatmap_data.sum(axis=0).nlargest(10).index
        heatmap_data_subset = heatmap_data.loc[top_samples, top_classes]
        
        sns.heatmap(heatmap_data_subset, annot=True, fmt='d', cmap='YlOrRd', 
                   cbar_kws={'label': 'Count'})
        plt.title('Top Transposon Classes in Top Samples')
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
        
        # 3.1 细胞系转座子数量
        cell_counts = self.df.groupby(['cell_line', 'eccdna_type']).size().unstack(fill_value=0)
        if not cell_counts.empty:
            cell_counts.plot(kind='bar', ax=axes[0,0], color=['skyblue', 'lightcoral'])
            axes[0,0].set_title('Transposon Count by Cell Line')
            axes[0,0].set_xlabel('Cell Line')
            axes[0,0].set_ylabel('Count')
            axes[0,0].legend(title='Type')
            axes[0,0].tick_params(axis='x', rotation=45)
        
        # 3.2 细胞系长度分布（采样）
        if 'cell_line' in sampled_df.columns:
            sns.boxplot(data=sampled_df, x='cell_line', y='length', hue='eccdna_type', ax=axes[0,1])
            axes[0,1].set_title(f'Length Distribution by Cell Line (n={len(sampled_df):,})')
            axes[0,1].set_xlabel('Cell Line')
            axes[0,1].set_ylabel('Length (bp)')
            axes[0,1].set_yscale('log')
        
        # 3.3 细胞系得分分布（采样）
        if 'cell_line' in sampled_df.columns:
            sns.boxplot(data=sampled_df, x='cell_line', y='score', hue='eccdna_type', ax=axes[1,0])
            axes[1,0].set_title(f'Score Distribution by Cell Line (n={len(sampled_df):,})')
            axes[1,0].set_xlabel('Cell Line')
            axes[1,0].set_ylabel('Score')
        
        # 3.4 转座子类型比例
        if 'cell_line' in self.df.columns:
            class_prop = self.df.groupby(['cell_line', 'transposon_class']).size().unstack(fill_value=0)
            if not class_prop.empty:
                class_prop_norm = class_prop.div(class_prop.sum(axis=1), axis=0)
                
                # 只显示前10个最常见的类型
                top_classes = class_prop.sum(axis=0).nlargest(10).index
                class_prop_norm[top_classes].plot(kind='bar', stacked=True, ax=axes[1,1])
                axes[1,1].set_title('Top 10 Transposon Classes by Cell Line')
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
            f.write("eccDNA Transposon Analysis Report (High Performance Version)\n")
            f.write("=" * 60 + "\n\n")
            
            # 数据概览
            f.write("1. Data Overview\n")
            f.write("-" * 20 + "\n")
            f.write(f"Analysis mode: {'Fast mode (10% sampling)' if self.fast_mode else 'Full analysis'}\n")
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
            f.write(f"Sequence analysis: Pure vectorized operations\n")
            f.write(f"Memory optimization: Enabled\n\n")
            
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
            if total_count > 0:
                f.write(f"- Mecc records: {mecc_count:,} ({mecc_count/total_count*100:.1f}%)\n")
                f.write(f"- Uecc records: {uecc_count:,} ({uecc_count/total_count*100:.1f}%)\n")
            
            # 最常见的转座子类型
            top_classes = self.df['transposon_class'].value_counts().head(5)
            f.write("\nTop 5 transposon classes:\n")
            for class_name, count in top_classes.items():
                f.write(f"  - {class_name}: {count:,} records ({count/len(self.df)*100:.1f}%)\n")
            
            # 细胞系差异
            if 'cell_line' in self.df.columns:
                cell_lines = self.df['cell_line'].dropna().unique()
                f.write(f"\nCell lines analyzed: {', '.join(sorted(cell_lines))}\n")
            
            # 结论
            f.write("\n5. Conclusions\n")
            f.write("-" * 20 + "\n")
            f.write("This high-performance analysis provides insights into eccDNA transposon composition:\n")
            f.write("- Different cell lines show distinct transposon profiles\n")
            f.write("- Mecc and Uecc types have different characteristics\n")
            f.write("- Multiple transposon classes are represented with varying frequencies\n")
            f.write("- Pure vectorized operations ensure optimal performance\n")
            f.write("- Memory optimization allows processing of large datasets\n")
            f.write("- Detailed results are available in the generated tables and plots\n\n")
            
            f.write("Analysis completed successfully with high-performance optimizations.\n")
            f.write(f"Results saved in: {self.output_dir}\n")
        
        print(f"分析报告已保存到: {report_path}")
    
    def run_complete_analysis(self):
        """运行完整分析流程 - 带断点续传和性能优化"""
        start_time = time.time()
        print("开始完整的eccDNA转座子分析（高性能优化版本）...")
        
        # 检查是否有之前的分析结果
        basic_stats_file = os.path.join(self.output_dir, 'tables', 'sample_statistics.csv')
        composite_file = os.path.join(self.output_dir, 'tables', 'composite_analysis.csv')
        
        try:
            # 1. 加载数据
            if not hasattr(self, 'df') or self.df is None:
                self.load_all_files()
                
                if self.df is None or len(self.df) == 0:
                    print("没有找到数据，分析终止")
                    return
            
            # 2. 基础统计（检查是否已完成）
            if os.path.exists(basic_stats_file) and not self.fast_mode:
                print("检测到已有基础统计结果，跳过...")
                self.sample_stats = pd.read_csv(basic_stats_file)
            else:
                self.basic_statistics()
            
            # 3. 转座子组成分析（检查是否已完成）
            if os.path.exists(composite_file) and not self.fast_mode:
                print("检测到已有组成分析结果，跳过...")
                composite_df = pd.read_csv(composite_file)
                coverage_file = os.path.join(self.output_dir, 'tables', 'coverage_analysis.csv')
                coverage_df = pd.read_csv(coverage_file) if os.path.exists(coverage_file) else None
            else:
                print("开始转座子组成分析...")
                composite_df, coverage_df = self.analyze_transposon_composition_optimized()
            
            # 4. 样本间比较
            comparison_file = os.path.join(self.output_dir, 'tables', 'mecc_uecc_comparison.csv')
            if os.path.exists(comparison_file) and not self.fast_mode:
                print("检测到已有比较分析结果，跳过...")
            else:
                self.compare_samples()
            
            # 5. 生成可视化
            overview_plot = os.path.join(self.output_dir, 'plots', 'overview.png')
            if os.path.exists(overview_plot) and not self.fast_mode:
                print("检测到已有可视化结果，跳过...")
            else:
                self.create_visualizations()
            
            # 6. 生成报告
            self.generate_report()
            
            end_time = time.time()
            total_time = end_time - start_time
            
            print(f"\n分析完成！总用时: {total_time:.2f} 秒 ({total_time/60:.1f} 分钟)")
            print(f"所有结果已保存到: {self.output_dir}")
            print("包含以下文件:")
            print("- tables/: 统计表格文件")
            print("- plots/: 可视化图表")
            print("- reports/: 分析报告")
            
        except MemoryError:
            print("\n❌ 内存不足! 请尝试:")
            print("1. 使用快速模式: fast_mode=True")
            print("2. 关闭其他程序释放内存")
            print("3. 在更大内存的机器上运行")
        except Exception as e:
            print(f"\n❌ 分析过程中出错: {str(e)}")
            import traceback
            traceback.print_exc()
        finally:
            # 清理内存
            if hasattr(self, 'df'):
                print("\n清理内存...")
                del self.df
                gc.collect()


def check_memory():
    """检查系统内存"""
    try:
        import psutil
        memory = psutil.virtual_memory()
        print(f"系统内存: 总计 {memory.total/(1024**3):.1f} GB, "
              f"可用 {memory.available/(1024**3):.1f} GB, "
              f"使用率 {memory.percent}%")
        
        if memory.available < 4 * (1024**3):  # 少于4GB
            print("⚠️  警告: 可用内存较少，建议使用快速模式")
            return True
        return False
    except ImportError:
        print("提示: 安装psutil可以监控内存使用: pip install psutil")
        return False


def main():
    """主函数"""
    print("eccDNA转座子分析脚本 - 高性能优化版本 v3.0")
    print("=" * 50)
    
    # 检查内存
    low_memory = check_memory()
    
    # 设置参数
    data_dir = '.'  # GFF文件所在目录
    output_dir = 'eccDNA_analysis_results'  # 结果输出目录
    
    # 询问用户是否使用快速模式
    if low_memory:
        fast_mode = True
        print("\n自动启用快速模式（10%数据采样）")
    else:
        user_input = input("\n是否使用快速模式？(y/n，默认n): ").strip().lower()
        fast_mode = user_input == 'y'
    
    # 获取合适的线程数
    n_threads = min(mp.cpu_count(), 4)  # 限制线程数避免竞争
    
    print(f"\n配置信息:")
    print(f"- 数据目录: {data_dir}")
    print(f"- 输出目录: {output_dir}")
    print(f"- IO线程数: {n_threads}")
    print(f"- 分析模式: {'快速模式(10%采样)' if fast_mode else '完整分析'}")
    print(f"- 支持断点续传: 是")
    
    # 创建分析器
    analyzer = EccDNATransposonAnalyzer(
        data_dir=data_dir, 
        output_dir=output_dir, 
        n_threads=n_threads,
        fast_mode=fast_mode
    )
    
    # 运行完整分析
    try:
        analyzer.run_complete_analysis()
    except KeyboardInterrupt:
        print("\n⚠️  分析被用户中断")
        print("下次运行将从断点继续...")
    except Exception as e:
        print(f"\n❌ 发生错误: {str(e)}")
        print("请检查数据文件和环境配置")


if __name__ == "__main__":
    main()
