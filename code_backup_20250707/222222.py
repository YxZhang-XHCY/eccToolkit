#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eccDNA复合转座子分析脚本
分析不同样本不同类型的eccDNA所富含的复合转座子组成
"""

import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

def load_and_process_data(file_path):
    """
    加载并处理eccDNA数据
    """
    print("正在加载数据...")
    df = pd.read_csv(file_path)
    print(f"原始数据shape: {df.shape}")
    print(f"包含样本: {df['Sample'].unique()}")
    print(f"包含类型: {df['Class'].unique()}")
    
    # 检查数据类型并处理异常值
    print("\n检查关键列的数据类型:")
    print(f"anno_Percent列类型: {df['anno_Percent'].dtype}")
    print(f"seq_length列类型: {df['seq_length'].dtype}")
    
    # 清理anno_Percent列
    print("\n清理anno_Percent数据...")
    df['anno_Percent'] = pd.to_numeric(df['anno_Percent'], errors='coerce')
    na_count = df['anno_Percent'].isna().sum()
    if na_count > 0:
        print(f"发现 {na_count} 个无效的anno_Percent值，将被删除")
        df = df.dropna(subset=['anno_Percent'])
    
    # 清理seq_length列
    print("清理seq_length数据...")
    df['seq_length'] = pd.to_numeric(df['seq_length'], errors='coerce')
    na_count = df['seq_length'].isna().sum()
    if na_count > 0:
        print(f"发现 {na_count} 个无效的seq_length值，将被删除")
        df = df.dropna(subset=['seq_length'])
    
    print(f"清理后数据shape: {df.shape}")
    return df

def calculate_eccdna_composition(df):
    """
    计算每个eccDNA的anno_Percent总和
    """
    print("\n正在计算每个eccDNA的anno_Percent总和...")
    
    # 检查分组前的数据
    print("检查分组键的唯一值数量:")
    print(f"Sample: {df['Sample'].nunique()}")
    print(f"Class: {df['Class'].nunique()}")
    print(f"seqname: {df['seqname'].nunique()}")
    
    # 清理motif列，确保没有NaN值
    df['motif'] = df['motif'].fillna('Unknown')
    
    # 按eccDNA分组计算anno_Percent总和
    try:
        # 先按start位置排序，确保motif按位置顺序排列
        df_sorted = df.sort_values(['Sample', 'Class', 'seqname', 'start'])
        
        # 使用更高效的聚合方式
        eccdna_summary = df_sorted.groupby(['Sample', 'Class', 'seqname']).agg({
            'anno_Percent': 'sum',
            'motif': lambda x: '|'.join([str(m) for m in x if pd.notna(m)]),  # 保留所有motif，已按位置排序
            'seq_length': 'first'
        }).reset_index()
        
        eccdna_summary.columns = ['Sample', 'Class', 'seqname', 'total_anno_percent', 'motif_combination', 'seq_length']
        
        print("正在生成按数量统计的motif组合...")
        # 生成按数量统计的组合（不考虑位置，方便模式识别）
        def create_count_combination(motifs_str):
            if pd.isna(motifs_str) or motifs_str == '':
                return ''
            motifs = motifs_str.split('|')
            motif_counts = {}
            for motif in motifs:
                motif_counts[motif] = motif_counts.get(motif, 0) + 1
            
            # 按motif名称排序，格式：motif(count)
            sorted_items = []
            for motif in sorted(motif_counts.keys()):
                count = motif_counts[motif]
                if count == 1:
                    sorted_items.append(motif)
                else:
                    sorted_items.append(f"{motif}({count})")
            return '|'.join(sorted_items)
        
        # 计算复合转座子数量
        def count_transposons(motifs_str):
            if pd.isna(motifs_str) or motifs_str == '':
                return 0
            return len(motifs_str.split('|'))
        
        # 计算独特转座子类型数量
        def count_unique_transposon_types(motifs_str):
            if pd.isna(motifs_str) or motifs_str == '':
                return 0
            motifs = motifs_str.split('|')
            return len(set(motifs))
        
        eccdna_summary['motif_combination_by_count'] = eccdna_summary['motif_combination'].apply(create_count_combination)
        eccdna_summary['transposon_count'] = eccdna_summary['motif_combination'].apply(count_transposons)
        eccdna_summary['unique_transposon_types'] = eccdna_summary['motif_combination'].apply(count_unique_transposon_types)
        
        print(f"总eccDNA数量: {len(eccdna_summary)}")
        print(f"anno_Percent总和分布:")
        print(eccdna_summary['total_anno_percent'].describe())
        print(f"\n复合转座子数量分布:")
        print(eccdna_summary['transposon_count'].describe())
        print(f"\n独特转座子类型数量分布:")
        print(eccdna_summary['unique_transposon_types'].describe())
        
        # 显示一些示例数据
        print("\n前5行数据示例:")
        print("位置排序组合 vs 数量统计组合 vs 转座子数量:")
        for i in range(min(5, len(eccdna_summary))):
            row = eccdna_summary.iloc[i]
            print(f"  {i+1}. 位置排序: {row['motif_combination']}")
            print(f"     数量统计: {row['motif_combination_by_count']}")
            print(f"     转座子总数: {row['transposon_count']}")
            print(f"     独特类型数: {row['unique_transposon_types']}")
            print()
        
        return eccdna_summary
        
    except Exception as e:
        print(f"分组聚合时出错: {str(e)}")
        # 调试信息
        print("\n显示前几行原始数据以便调试:")
        print(df[['Sample', 'Class', 'seqname', 'anno_Percent', 'motif', 'seq_length']].head(10))
        raise

def filter_high_coverage_eccdna(eccdna_summary, threshold=60):
    """
    筛选出anno_Percent总和大于阈值的eccDNA
    """
    print(f"\n正在筛选anno_Percent总和大于{threshold}%的eccDNA...")
    
    filtered_df = eccdna_summary[eccdna_summary['total_anno_percent'] > threshold].copy()
    
    print(f"筛选后eccDNA数量: {len(filtered_df)} ({len(filtered_df)/len(eccdna_summary)*100:.1f}%)")
    
    # 按样本和类型统计
    sample_class_counts = filtered_df.groupby(['Sample', 'Class']).size().reset_index(name='count')
    print("\n各样本各类型的高覆盖eccDNA数量:")
    print(sample_class_counts)
    
    # 统计转座子数量分布
    print("\n高覆盖eccDNA中的转座子数量分布:")
    transposon_dist = Counter(filtered_df['transposon_count'])
    for count, freq in sorted(transposon_dist.items()):
        print(f"  {count}个转座子: {freq} eccDNA ({freq/len(filtered_df)*100:.1f}%)")
    
    print(f"\n高覆盖eccDNA中转座子数量统计:")
    print(f"  平均转座子数量: {filtered_df['transposon_count'].mean():.2f}")
    print(f"  中位数转座子数量: {filtered_df['transposon_count'].median():.1f}")
    print(f"  最多转座子数量: {filtered_df['transposon_count'].max()}")
    print(f"  最少转座子数量: {filtered_df['transposon_count'].min()}")
    
    return filtered_df

def analyze_motif_combinations(filtered_df):
    """
    分析motif组合
    """
    print("\n正在分析motif组合...")
    
    # 统计整体最常见的组合（按位置）
    overall_combinations = Counter(filtered_df['motif_combination'])
    print("\n整体最常见的top10 motif组合（按位置排序）:")
    for combo, count in overall_combinations.most_common(10):
        print(f"{combo}: {count} ({count/len(filtered_df)*100:.1f}%)")
    
    # 统计整体最常见的组合（按数量）
    overall_combinations_count = Counter(filtered_df['motif_combination_by_count'])
    print("\n整体最常见的top10 motif组合（按数量统计）:")
    for combo, count in overall_combinations_count.most_common(10):
        print(f"{combo}: {count} ({count/len(filtered_df)*100:.1f}%)")
    
    # 按样本和类型分析
    results = {}
    for (sample, class_type), group in filtered_df.groupby(['Sample', 'Class']):
        combinations_pos = Counter(group['motif_combination'])
        combinations_count = Counter(group['motif_combination_by_count'])
        
        results[(sample, class_type)] = {
            'total_count': len(group),
            'top3_combinations_by_position': combinations_pos.most_common(3),
            'top3_combinations_by_count': combinations_count.most_common(3),
            'unique_combinations_by_position': len(combinations_pos),
            'unique_combinations_by_count': len(combinations_count),
            'avg_transposon_count': group['transposon_count'].mean(),
            'avg_unique_types': group['unique_transposon_types'].mean()
        }
        
        print(f"\n{sample} - {class_type} (总数: {len(group)}):")
        print(f"  独特位置组合数: {len(combinations_pos)}")
        print(f"  独特数量组合数: {len(combinations_count)}")
        print(f"  平均转座子数量: {group['transposon_count'].mean():.2f}")
        print(f"  平均独特类型数: {group['unique_transposon_types'].mean():.2f}")
        
        print("  Top3组合（按位置）:")
        for i, (combo, count) in enumerate(combinations_pos.most_common(3), 1):
            print(f"    {i}. {combo}: {count} ({count/len(group)*100:.1f}%)")
        
        print("  Top3组合（按数量）:")
        for i, (combo, count) in enumerate(combinations_count.most_common(3), 1):
            print(f"    {i}. {combo}: {count} ({count/len(group)*100:.1f}%)")
    
    return results, overall_combinations, overall_combinations_count

def generate_composition_report(filtered_df, results):
    """
    生成组成分析报告
    """
    print("\n" + "="*80)
    print("eccDNA组成分析报告")
    print("="*80)
    
    # 统计单个motif的出现频率（从位置排序组合中统计）
    all_motifs = []
    for combo in filtered_df['motif_combination']:
        all_motifs.extend(combo.split('|'))
    
    motif_counter = Counter(all_motifs)
    print(f"\n最常见的单个motif (Top10):")
    for motif, count in motif_counter.most_common(10):
        print(f"  {motif}: {count} ({count/len(all_motifs)*100:.1f}%)")
    
    # 分析组合复杂度（从位置排序组合）
    combo_complexity_pos = [len(combo.split('|')) for combo in filtered_df['motif_combination']]
    print(f"\nmotif组合复杂度分析（按位置）:")
    print(f"  平均每个eccDNA包含motif数: {np.mean(combo_complexity_pos):.2f}")
    print(f"  最少motif数: {min(combo_complexity_pos)}")
    print(f"  最多motif数: {max(combo_complexity_pos)}")
    
    complexity_dist_pos = Counter(combo_complexity_pos)
    print("  motif数量分布:")
    for num_motifs, count in sorted(complexity_dist_pos.items()):
        print(f"    {num_motifs}个motif: {count} eccDNA ({count/len(combo_complexity_pos)*100:.1f}%)")
    
    # 分析独特motif类型数量（从数量统计组合）
    unique_motif_types = []
    for combo in filtered_df['motif_combination_by_count']:
        if combo:
            motif_types = set()
            for item in combo.split('|'):
                motif_name = item.split('(')[0]  # 去掉计数部分
                motif_types.add(motif_name)
            unique_motif_types.append(len(motif_types))
    
    print(f"\n独特motif类型数量分析:")
    print(f"  平均每个eccDNA包含独特motif类型数: {np.mean(unique_motif_types):.2f}")
    print(f"  最少独特motif类型数: {min(unique_motif_types)}")
    print(f"  最多独特motif类型数: {max(unique_motif_types)}")
    
    unique_types_dist = Counter(unique_motif_types)
    print("  独特motif类型数分布:")
    for num_types, count in sorted(unique_types_dist.items()):
        print(f"    {num_types}种独特motif: {count} eccDNA ({count/len(unique_motif_types)*100:.1f}%)")
    
    # 比较两种表示方法的差异
    print(f"\n两种表示方法的比较:")
    pos_unique = len(set(filtered_df['motif_combination']))
    count_unique = len(set(filtered_df['motif_combination_by_count']))
    print(f"  按位置排序的独特组合数: {pos_unique}")
    print(f"  按数量统计的独特组合数: {count_unique}")
    print(f"  压缩比例: {count_unique/pos_unique*100:.1f}% (数量统计相对于位置排序)")
    
    # 转座子数量相关统计
    print(f"\n转座子数量统计:")
    print(f"  总转座子数量: {filtered_df['transposon_count'].sum()}")
    print(f"  平均每个eccDNA的转座子数量: {filtered_df['transposon_count'].mean():.2f}")
    print(f"  转座子数量最多的eccDNA: {filtered_df['transposon_count'].max()}")
    print(f"  包含单个转座子的eccDNA数量: {(filtered_df['transposon_count'] == 1).sum()}")
    print(f"  包含多个转座子的eccDNA数量: {(filtered_df['transposon_count'] > 1).sum()}")

def save_results(filtered_df, output_file='high_coverage_eccdna.csv'):
    """
    保存筛选后的结果
    """
    print(f"\n正在保存结果到 {output_file}...")
    
    # 重新排列列的顺序，让转座子数量列更显眼
    columns_order = [
        'Sample', 'Class', 'seqname', 'total_anno_percent', 
        'transposon_count', 'unique_transposon_types',
        'motif_combination', 'motif_combination_by_count', 'seq_length'
    ]
    
    # 确保所有列都存在
    existing_columns = [col for col in columns_order if col in filtered_df.columns]
    output_df = filtered_df[existing_columns].copy()
    
    output_df.to_csv(output_file, index=False, encoding='utf-8')
    print("保存完成!")
    print(f"输出文件包含以下列:")
    for i, col in enumerate(output_df.columns, 1):
        print(f"  {i}. {col}")

def create_visualizations(filtered_df, results):
    """
    创建可视化图表
    """
    print("\n正在创建可视化图表...")
    
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei']  # 支持中文
    plt.rcParams['axes.unicode_minus'] = False
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # 1. 样本和类型分布
    sample_class_counts = filtered_df.groupby(['Sample', 'Class']).size().reset_index(name='count')
    sns.barplot(data=sample_class_counts, x='Sample', y='count', hue='Class', ax=axes[0, 0])
    axes[0, 0].set_title('High Coverage eccDNA Count by Sample and Class')
    axes[0, 0].tick_params(axis='x', rotation=45)
    
    # 2. anno_Percent分布
    axes[0, 1].hist(filtered_df['total_anno_percent'], bins=20, alpha=0.7, edgecolor='black')
    axes[0, 1].set_xlabel('Total anno_Percent')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Distribution of Total anno_Percent')
    axes[0, 1].axvline(x=80, color='red', linestyle='--', label='Threshold (80%)')
    axes[0, 1].legend()
    
    # 3. 转座子数量分布
    transposon_counts = Counter(filtered_df['transposon_count'])
    axes[0, 2].bar(transposon_counts.keys(), transposon_counts.values())
    axes[0, 2].set_xlabel('Number of Transposons per eccDNA')
    axes[0, 2].set_ylabel('Count')
    axes[0, 2].set_title('Distribution of Transposon Count')
    
    # 4. motif组合复杂度
    combo_complexity = [len(combo.split('|')) for combo in filtered_df['motif_combination']]
    complexity_counts = Counter(combo_complexity)
    axes[1, 0].bar(complexity_counts.keys(), complexity_counts.values())
    axes[1, 0].set_xlabel('Number of Motifs per eccDNA')
    axes[1, 0].set_ylabel('Count')
    axes[1, 0].set_title('Motif Combination Complexity')
    
    # 5. 序列长度分布
    axes[1, 1].scatter(filtered_df['seq_length'], filtered_df['total_anno_percent'], 
                      c=filtered_df['transposon_count'], cmap='viridis', alpha=0.6)
    axes[1, 1].set_xlabel('Sequence Length (bp)')
    axes[1, 1].set_ylabel('Total anno_Percent')
    axes[1, 1].set_title('Length vs Anno_Percent (colored by transposon count)')
    cbar = plt.colorbar(axes[1, 1].collections[0], ax=axes[1, 1])
    cbar.set_label('Transposon Count')
    
    # 6. 转座子数量 vs anno_Percent
    axes[1, 2].scatter(filtered_df['transposon_count'], filtered_df['total_anno_percent'], alpha=0.6)
    axes[1, 2].set_xlabel('Transposon Count')
    axes[1, 2].set_ylabel('Total anno_Percent')
    axes[1, 2].set_title('Transposon Count vs Total anno_Percent')
    
    plt.tight_layout()
    plt.savefig('eccdna_analysis_plots.png', dpi=300, bbox_inches='tight')
    print("图表已保存为 eccdna_analysis_plots.png")

def main():
    """
    主函数
    """
    # 文件路径 - 请修改为你的实际文件路径
    file_path = 'TE.multiple_compound.csv'
    
    try:
        # 1. 加载数据
        df = load_and_process_data(file_path)
        
        # 2. 计算eccDNA组成
        eccdna_summary = calculate_eccdna_composition(df)
        
        # 3. 筛选高覆盖率eccDNA
        filtered_df = filter_high_coverage_eccdna(eccdna_summary, threshold=80)
        
        # 4. 分析motif组合
        results, overall_combinations, overall_combinations_count = analyze_motif_combinations(filtered_df)
        
        # 5. 生成详细报告
        generate_composition_report(filtered_df, results)
        
        # 6. 保存结果
        save_results(filtered_df)
        
        # 7. 创建可视化
        create_visualizations(filtered_df, results)
        
        print("\n分析完成! 请查看生成的文件:")
        print("- high_coverage_eccdna.csv: 筛选后的eccDNA数据（包含转座子数量列）")
        print("- eccdna_analysis_plots.png: 可视化图表")
        
    except FileNotFoundError:
        print(f"错误: 找不到文件 {file_path}")
        print("请确保文件路径正确，或将文件放在当前目录下")
    except Exception as e:
        print(f"分析过程中出现错误: {str(e)}")

if __name__ == "__main__":
    main()
