#!/usr/bin/env python3
"""
单一复合转座子占比分析脚本
分析motif_percent和anno_Percent在不同样本和类型中的分布
使用5个bins: 0-20%, 20-40%, 40-60%, 60-80%, 80-100% (包括>100%)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
import time

# 设置中文字体支持
rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
rcParams['axes.unicode_minus'] = False

def create_percentage_bins(values):
    """
    创建百分比分组
    0-20%, 20-40%, 40-60%, 60-80%, 80-100% (包括>100%)
    """
    bins = [0, 20, 40, 60, 80, float('inf')]
    labels = ['0-20%', '20-40%', '40-60%', '60-80%', '80-100%']
    
    # 处理大于100的值
    values_capped = np.where(values > 100, 100, values)
    
    return pd.cut(values_capped, bins=bins, labels=labels, include_lowest=True, right=False)

def analyze_single_compound_te(input_file):
    """
    分析单一复合转座子的占比分布
    """
    start_time = time.time()
    
    # 读取数据
    print("正在读取单一复合转座子数据...")
    df = pd.read_csv(input_file)
    print(f"单一复合转座子记录数: {len(df):,}")
    
    # 基本统计信息
    print("\n=== 基本统计信息 ===")
    print(f"样本数量: {df['Sample'].nunique()}")
    print(f"转座子类型数量: {df['Class'].nunique()}")
    print(f"样本列表: {', '.join(sorted(df['Sample'].unique()))}")
    print(f"转座子类型列表: {', '.join(sorted(df['Class'].unique()))}")
    
    # 处理缺失值
    print("\n=== 数据质量检查 ===")
    motif_na = df['motif_percent'].isna().sum()
    anno_na = df['anno_Percent'].isna().sum()
    print(f"motif_percent缺失值: {motif_na:,} ({motif_na/len(df)*100:.1f}%)")
    print(f"anno_Percent缺失值: {anno_na:,} ({anno_na/len(df)*100:.1f}%)")
    
    # 移除缺失值进行分析
    df_motif = df.dropna(subset=['motif_percent']).copy()
    df_anno = df.dropna(subset=['anno_Percent']).copy()
    
    print(f"用于motif_percent分析的记录数: {len(df_motif):,}")
    print(f"用于anno_Percent分析的记录数: {len(df_anno):,}")
    
    # 创建百分比分组
    print("\n=== 创建百分比分组 ===")
    df_motif['motif_bin'] = create_percentage_bins(df_motif['motif_percent'])
    df_anno['anno_bin'] = create_percentage_bins(df_anno['anno_Percent'])
    
    # 统计超过100%的情况
    motif_over100 = (df_motif['motif_percent'] > 100).sum()
    anno_over100 = (df_anno['anno_Percent'] > 100).sum()
    print(f"motif_percent > 100%的记录数: {motif_over100:,}")
    print(f"anno_Percent > 100%的记录数: {anno_over100:,}")
    
    # 分析motif_percent分布
    print("\n" + "="*50)
    print("MOTIF_PERCENT 分布分析")
    print("="*50)
    
    analyze_percentage_distribution(df_motif, 'motif_percent', 'motif_bin', 'motif')
    
    # 分析anno_Percent分布  
    print("\n" + "="*50)
    print("ANNO_PERCENT 分布分析")
    print("="*50)
    
    analyze_percentage_distribution(df_anno, 'anno_Percent', 'anno_bin', 'anno')
    
    # 生成可视化
    create_visualizations(df_motif, df_anno)
    
    # 生成详细统计表
    generate_detailed_tables(df_motif, df_anno)
    
    end_time = time.time()
    print(f"\n✅ 分析完成！总耗时: {end_time - start_time:.2f} 秒")

def analyze_percentage_distribution(df, percent_col, bin_col, analysis_type):
    """
    分析百分比分布
    """
    print(f"\n--- {percent_col.upper()} 描述性统计 ---")
    stats = df[percent_col].describe()
    print(f"平均值: {stats['mean']:.2f}%")
    print(f"中位数: {stats['50%']:.2f}%")
    print(f"标准差: {stats['std']:.2f}%")
    print(f"最小值: {stats['min']:.2f}%")
    print(f"最大值: {stats['max']:.2f}%")
    
    print(f"\n--- {percent_col.upper()} 分组分布 ---")
    bin_counts = df[bin_col].value_counts().sort_index()
    bin_percent = (bin_counts / len(df) * 100).round(1)
    
    for bin_name, count in bin_counts.items():
        pct = bin_percent[bin_name]
        print(f"{bin_name}: {count:,} 条记录 ({pct}%)")
    
    print(f"\n--- 按样本和类型的 {percent_col.upper()} 分布 ---")
    
    # 按Sample分组统计
    sample_stats = df.groupby(['Sample', bin_col]).size().unstack(fill_value=0)
    sample_percent = sample_stats.div(sample_stats.sum(axis=1), axis=0) * 100
    
    print("各样本的分组分布 (百分比):")
    print(sample_percent.round(1))
    
    # 按Class分组统计
    class_stats = df.groupby(['Class', bin_col]).size().unstack(fill_value=0)
    class_percent = class_stats.div(class_stats.sum(axis=1), axis=0) * 100
    
    print(f"\n各转座子类型的分组分布 (百分比):")
    print(class_percent.round(1))
    
    # 按Sample和Class组合统计
    print(f"\n--- 按样本×类型的 {percent_col.upper()} 详细分布 ---")
    combined_stats = df.groupby(['Sample', 'Class', bin_col]).size().unstack(fill_value=0)
    combined_percent = combined_stats.div(combined_stats.sum(axis=1), axis=0) * 100
    
    print("样本×类型组合的分组分布 (百分比):")
    print(combined_percent.round(1))
    
    # 保存详细表格
    sample_percent.to_csv(f'{analysis_type}_percent_by_sample.csv')
    class_percent.to_csv(f'{analysis_type}_percent_by_class.csv')
    combined_percent.to_csv(f'{analysis_type}_percent_by_sample_class.csv')
    
    print(f"\n📊 详细统计表已保存:")
    print(f"  - {analysis_type}_percent_by_sample.csv")
    print(f"  - {analysis_type}_percent_by_class.csv") 
    print(f"  - {analysis_type}_percent_by_sample_class.csv")

def create_visualizations(df_motif, df_anno):
    """
    创建可视化图表
    """
    print("\n=== 生成可视化图表 ===")
    
    # 设置图表样式
    plt.style.use('default')
    sns.set_palette("husl")
    
    # 创建图表
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('单一复合转座子占比分布分析', fontsize=16, fontweight='bold')
    
    # 1. motif_percent总体分布
    axes[0,0].hist(df_motif['motif_percent'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0,0].set_title('Motif Percent 分布直方图')
    axes[0,0].set_xlabel('Motif Percent (%)')
    axes[0,0].set_ylabel('频次')
    axes[0,0].axvline(df_motif['motif_percent'].mean(), color='red', linestyle='--', label=f'平均值: {df_motif["motif_percent"].mean():.1f}%')
    axes[0,0].legend()
    
    # 2. motif_percent分组条形图
    motif_bin_counts = df_motif['motif_bin'].value_counts().sort_index()
    axes[0,1].bar(range(len(motif_bin_counts)), motif_bin_counts.values, color='lightgreen', edgecolor='black')
    axes[0,1].set_title('Motif Percent 分组分布')
    axes[0,1].set_xlabel('分组')
    axes[0,1].set_ylabel('记录数')
    axes[0,1].set_xticks(range(len(motif_bin_counts)))
    axes[0,1].set_xticklabels(motif_bin_counts.index, rotation=45)
    
    # 3. motif_percent按样本分布
    sample_motif = df_motif.groupby(['Sample', 'motif_bin']).size().unstack(fill_value=0)
    sample_motif_pct = sample_motif.div(sample_motif.sum(axis=1), axis=0) * 100
    sample_motif_pct.plot(kind='bar', stacked=True, ax=axes[0,2], colormap='viridis')
    axes[0,2].set_title('Motif Percent 按样本分布')
    axes[0,2].set_xlabel('样本')
    axes[0,2].set_ylabel('百分比 (%)')
    axes[0,2].legend(title='分组', bbox_to_anchor=(1.05, 1), loc='upper left')
    axes[0,2].tick_params(axis='x', rotation=45)
    
    # 4. anno_Percent总体分布
    axes[1,0].hist(df_anno['anno_Percent'], bins=50, alpha=0.7, color='orange', edgecolor='black')
    axes[1,0].set_title('Anno Percent 分布直方图')
    axes[1,0].set_xlabel('Anno Percent (%)')
    axes[1,0].set_ylabel('频次')
    axes[1,0].axvline(df_anno['anno_Percent'].mean(), color='red', linestyle='--', label=f'平均值: {df_anno["anno_Percent"].mean():.1f}%')
    axes[1,0].legend()
    
    # 5. anno_Percent分组条形图
    anno_bin_counts = df_anno['anno_bin'].value_counts().sort_index()
    axes[1,1].bar(range(len(anno_bin_counts)), anno_bin_counts.values, color='salmon', edgecolor='black')
    axes[1,1].set_title('Anno Percent 分组分布')
    axes[1,1].set_xlabel('分组')
    axes[1,1].set_ylabel('记录数')
    axes[1,1].set_xticks(range(len(anno_bin_counts)))
    axes[1,1].set_xticklabels(anno_bin_counts.index, rotation=45)
    
    # 6. anno_Percent按样本分布
    sample_anno = df_anno.groupby(['Sample', 'anno_bin']).size().unstack(fill_value=0)
    sample_anno_pct = sample_anno.div(sample_anno.sum(axis=1), axis=0) * 100
    sample_anno_pct.plot(kind='bar', stacked=True, ax=axes[1,2], colormap='plasma')
    axes[1,2].set_title('Anno Percent 按样本分布')
    axes[1,2].set_xlabel('样本')
    axes[1,2].set_ylabel('百分比 (%)')
    axes[1,2].legend(title='分组', bbox_to_anchor=(1.05, 1), loc='upper left')
    axes[1,2].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig('single_compound_te_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig('single_compound_te_analysis.pdf', bbox_inches='tight')
    print("📈 可视化图表已保存: single_compound_te_analysis.png/pdf")
    plt.show()

def generate_detailed_tables(df_motif, df_anno):
    """
    生成详细的统计表格
    """
    print("\n=== 生成汇总统计表 ===")
    
    # 创建汇总表
    summary_stats = []
    
    # 整体统计
    summary_stats.append({
        '分析类型': 'motif_percent',
        '样本': 'ALL',
        '转座子类型': 'ALL',
        '记录数': len(df_motif),
        '平均值': df_motif['motif_percent'].mean(),
        '中位数': df_motif['motif_percent'].median(),
        '标准差': df_motif['motif_percent'].std(),
        '0-20%': (df_motif['motif_bin'] == '0-20%').sum(),
        '20-40%': (df_motif['motif_bin'] == '20-40%').sum(),
        '40-60%': (df_motif['motif_bin'] == '40-60%').sum(),
        '60-80%': (df_motif['motif_bin'] == '60-80%').sum(),
        '80-100%': (df_motif['motif_bin'] == '80-100%').sum()
    })
    
    summary_stats.append({
        '分析类型': 'anno_Percent',
        '样本': 'ALL',
        '转座子类型': 'ALL', 
        '记录数': len(df_anno),
        '平均值': df_anno['anno_Percent'].mean(),
        '中位数': df_anno['anno_Percent'].median(),
        '标准差': df_anno['anno_Percent'].std(),
        '0-20%': (df_anno['anno_bin'] == '0-20%').sum(),
        '20-40%': (df_anno['anno_bin'] == '20-40%').sum(),
        '40-60%': (df_anno['anno_bin'] == '40-60%').sum(),
        '60-80%': (df_anno['anno_bin'] == '60-80%').sum(),
        '80-100%': (df_anno['anno_bin'] == '80-100%').sum()
    })
    
    # 按样本统计
    for sample in df_motif['Sample'].unique():
        sample_motif = df_motif[df_motif['Sample'] == sample]
        summary_stats.append({
            '分析类型': 'motif_percent',
            '样本': sample,
            '转座子类型': 'ALL',
            '记录数': len(sample_motif),
            '平均值': sample_motif['motif_percent'].mean(),
            '中位数': sample_motif['motif_percent'].median(),
            '标准差': sample_motif['motif_percent'].std(),
            '0-20%': (sample_motif['motif_bin'] == '0-20%').sum(),
            '20-40%': (sample_motif['motif_bin'] == '20-40%').sum(),
            '40-60%': (sample_motif['motif_bin'] == '40-60%').sum(),
            '60-80%': (sample_motif['motif_bin'] == '60-80%').sum(),
            '80-100%': (sample_motif['motif_bin'] == '80-100%').sum()
        })
        
        sample_anno = df_anno[df_anno['Sample'] == sample]
        if len(sample_anno) > 0:
            summary_stats.append({
                '分析类型': 'anno_Percent',
                '样本': sample,
                '转座子类型': 'ALL',
                '记录数': len(sample_anno),
                '平均值': sample_anno['anno_Percent'].mean(),
                '中位数': sample_anno['anno_Percent'].median(),
                '标准差': sample_anno['anno_Percent'].std(),
                '0-20%': (sample_anno['anno_bin'] == '0-20%').sum(),
                '20-40%': (sample_anno['anno_bin'] == '20-40%').sum(),
                '40-60%': (sample_anno['anno_bin'] == '40-60%').sum(),
                '60-80%': (sample_anno['anno_bin'] == '60-80%').sum(),
                '80-100%': (sample_anno['anno_bin'] == '80-100%').sum()
            })
    
    # 保存汇总表
    summary_df = pd.DataFrame(summary_stats)
    summary_df = summary_df.round(2)
    summary_df.to_csv('single_compound_te_summary.csv', index=False)
    
    print("📋 汇总统计表已保存: single_compound_te_summary.csv")
    print("\n汇总表预览:")
    print(summary_df.head(10))

def main():
    """
    主函数
    """
    input_file = "TE.single_compound.csv"
    
    try:
        analyze_single_compound_te(input_file)
        
        print(f"\n✅ 单一复合转座子分析完成！")
        print("📁 生成的文件:")
        print("  - motif_percent_by_sample.csv")
        print("  - motif_percent_by_class.csv")
        print("  - motif_percent_by_sample_class.csv")
        print("  - anno_percent_by_sample.csv")
        print("  - anno_percent_by_class.csv")
        print("  - anno_percent_by_sample_class.csv")
        print("  - single_compound_te_summary.csv")
        print("  - single_compound_te_analysis.png/pdf")
        
    except FileNotFoundError as e:
        print(f"❌ 错误：找不到文件 {e}")
        print("请确保 TE.single_compound.csv 文件存在于当前目录")
        
    except Exception as e:
        print(f"❌ 处理过程中发生错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
