#!/usr/bin/env python3
"""
eccDNA富集/耗竭与RNA-seq差异表达基因关联分析 - 逐样本分析版本
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import hypergeom, fisher_exact, chi2_contingency, spearmanr
from matplotlib_venn import venn2, venn3
from statsmodels.stats.multitest import multipletests
import os
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def load_eccdna_data(filepath):
    """加载eccDNA富集数据"""
    df = pd.read_csv(filepath)
    return df

def load_degs_data(filepath):
    """加载DEGs数据"""
    df = pd.read_csv(filepath, sep='\t', header=None, 
                     names=['gene', 'ensembl_id', 'tumor_exp', 'normal_exp', 'log2fc', 'qvalue'])
    return df

def analyze_eccdna_degs_overlap(eccdna_df, degs_df, fc_threshold=1):
    """分析eccDNA与DEGs重叠"""
    
    # 获取显著的eccDNA基因
    eccdna_enriched = set(eccdna_df[(eccdna_df['direction'] == 'enrichment') & 
                                    (eccdna_df['significant'] == True)]['gene'])
    eccdna_depleted = set(eccdna_df[(eccdna_df['direction'] == 'depletion') & 
                                    (eccdna_df['significant'] == True)]['gene'])
    
    # 获取DEGs（使用指定的FC阈值）
    degs_up = set(degs_df[degs_df['log2fc'] > fc_threshold]['gene'])
    degs_down = set(degs_df[degs_df['log2fc'] < -fc_threshold]['gene'])
    
    # 计算交集
    results = {
        'enriched_up': eccdna_enriched & degs_up,
        'enriched_down': eccdna_enriched & degs_down,
        'depleted_up': eccdna_depleted & degs_up,
        'depleted_down': eccdna_depleted & degs_down
    }
    
    # 统计信息
    stats = {
        'n_eccdna_enriched': len(eccdna_enriched),
        'n_eccdna_depleted': len(eccdna_depleted),
        'n_degs_up': len(degs_up),
        'n_degs_down': len(degs_down),
        'n_enriched_up': len(results['enriched_up']),
        'n_enriched_down': len(results['enriched_down']),
        'n_depleted_up': len(results['depleted_up']),
        'n_depleted_down': len(results['depleted_down'])
    }
    
    return results, stats, eccdna_enriched, eccdna_depleted, degs_up, degs_down

def fisher_test(overlap, set1_size, set2_size, total_genes):
    """Fisher精确检验"""
    table = [[overlap, set1_size - overlap],
             [set2_size - overlap, total_genes - set1_size - set2_size + overlap]]
    odds_ratio, p_value = fisher_exact(table, alternative='greater')
    return odds_ratio, p_value

def perform_gradient_fc_analysis_single_sample(eccdna_df, degs_df, sample_name, fc_thresholds=[1, 2, 4, 6]):
    """对单个样本执行梯度FC阈值分析"""
    results_summary = []
    
    # 只使用该样本的数据
    sample_eccdna = eccdna_df[eccdna_df['Sample'] == sample_name]
    
    # 计算背景基因数
    all_genes = set(sample_eccdna['gene']) | set(degs_df['gene'])
    total_genes = len(all_genes)
    
    # 获取该样本的eccDNA基因集
    eccdna_enriched_all = set(sample_eccdna[(sample_eccdna['direction'] == 'enrichment') & 
                                            (sample_eccdna['significant'] == True)]['gene'])
    eccdna_depleted_all = set(sample_eccdna[(sample_eccdna['direction'] == 'depletion') & 
                                            (sample_eccdna['significant'] == True)]['gene'])
    
    print(f"\n--- Sample: {sample_name} ---")
    print(f"Total background genes: {total_genes}")
    print(f"eccDNA enriched genes: {len(eccdna_enriched_all)}")
    print(f"eccDNA depleted genes: {len(eccdna_depleted_all)}")
    
    for fc in fc_thresholds:
        # 分析当前FC阈值
        results, stats, _, _, _, _ = analyze_eccdna_degs_overlap(sample_eccdna, degs_df, fc_threshold=fc)
        
        # 计算统计显著性
        analysis_results = {
            'sample': sample_name,
            'fc_threshold': fc,
            'n_degs_up': stats['n_degs_up'],
            'n_degs_down': stats['n_degs_down'],
            'n_enriched_up': stats['n_enriched_up'],
            'n_enriched_down': stats['n_enriched_down'],
            'n_depleted_up': stats['n_depleted_up'],
            'n_depleted_down': stats['n_depleted_down']
        }
        
        # Fisher检验
        if stats['n_eccdna_enriched'] > 0 and stats['n_degs_up'] > 0:
            or_eu, p_eu = fisher_test(stats['n_enriched_up'], stats['n_eccdna_enriched'], 
                                     stats['n_degs_up'], total_genes)
            analysis_results['enriched_up_OR'] = or_eu
            analysis_results['enriched_up_pvalue'] = p_eu
        
        if stats['n_eccdna_depleted'] > 0 and stats['n_degs_down'] > 0:
            or_dd, p_dd = fisher_test(stats['n_depleted_down'], stats['n_eccdna_depleted'], 
                                     stats['n_degs_down'], total_genes)
            analysis_results['depleted_down_OR'] = or_dd
            analysis_results['depleted_down_pvalue'] = p_dd
        
        results_summary.append(analysis_results)
    
    return pd.DataFrame(results_summary)

def plot_sample_specific_analysis(sample_name, gradient_df, eccdna_df, degs_df):
    """为单个样本创建分析图表"""
    # 创建样本特定的输出目录
    output_dir = f"sample_{sample_name}_analysis"
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 梯度分析图
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Sample {sample_name}: eccDNA-DEG Association Across FC Thresholds', fontsize=16)
    
    # 重叠基因数量
    ax = axes[0, 0]
    ax.plot(gradient_df['fc_threshold'], gradient_df['n_enriched_up'], 'o-', label='Enriched ∩ Up', color='red')
    ax.plot(gradient_df['fc_threshold'], gradient_df['n_depleted_down'], 's-', label='Depleted ∩ Down', color='blue')
    ax.plot(gradient_df['fc_threshold'], gradient_df['n_enriched_down'], '^-', label='Enriched ∩ Down', color='orange')
    ax.plot(gradient_df['fc_threshold'], gradient_df['n_depleted_up'], 'v-', label='Depleted ∩ Up', color='green')
    ax.set_xlabel('FC Threshold')
    ax.set_ylabel('Number of Overlapping Genes')
    ax.set_title('Overlapping Genes vs FC Threshold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Odds Ratio
    ax = axes[0, 1]
    if 'enriched_up_OR' in gradient_df.columns:
        mask = gradient_df['enriched_up_OR'].notna()
        ax.plot(gradient_df.loc[mask, 'fc_threshold'], gradient_df.loc[mask, 'enriched_up_OR'], 
                'o-', label='Enriched-Up', color='red')
    if 'depleted_down_OR' in gradient_df.columns:
        mask = gradient_df['depleted_down_OR'].notna()
        ax.plot(gradient_df.loc[mask, 'fc_threshold'], gradient_df.loc[mask, 'depleted_down_OR'], 
                's-', label='Depleted-Down', color='blue')
    ax.axhline(y=1, color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel('FC Threshold')
    ax.set_ylabel('Odds Ratio')
    ax.set_title('Association Strength (OR) vs FC Threshold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    
    # P值
    ax = axes[1, 0]
    if 'enriched_up_pvalue' in gradient_df.columns:
        mask = gradient_df['enriched_up_pvalue'].notna()
        ax.plot(gradient_df.loc[mask, 'fc_threshold'], 
                -np.log10(gradient_df.loc[mask, 'enriched_up_pvalue']), 
                'o-', label='Enriched-Up', color='red')
    if 'depleted_down_pvalue' in gradient_df.columns:
        mask = gradient_df['depleted_down_pvalue'].notna()
        ax.plot(gradient_df.loc[mask, 'fc_threshold'], 
                -np.log10(gradient_df.loc[mask, 'depleted_down_pvalue']), 
                's-', label='Depleted-Down', color='blue')
    ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5, label='p=0.05')
    ax.set_xlabel('FC Threshold')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('Statistical Significance vs FC Threshold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # DEG数量
    ax = axes[1, 1]
    ax.plot(gradient_df['fc_threshold'], gradient_df['n_degs_up'], 'o-', label='Upregulated DEGs', color='red')
    ax.plot(gradient_df['fc_threshold'], gradient_df['n_degs_down'], 's-', label='Downregulated DEGs', color='blue')
    ax.set_xlabel('FC Threshold')
    ax.set_ylabel('Number of DEGs')
    ax.set_title('DEG Count vs FC Threshold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_gradient_analysis.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. 韦恩图 (FC=2)
    sample_eccdna = eccdna_df[eccdna_df['Sample'] == sample_name]
    results, stats, eccdna_enriched, eccdna_depleted, degs_up, degs_down = analyze_eccdna_degs_overlap(
        sample_eccdna, degs_df, fc_threshold=2)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'Sample {sample_name}: eccDNA vs DEGs Overlap (FC=2)', fontsize=16)
    
    # 韦恩图
    ax = axes[0, 0]
    venn2([eccdna_enriched, degs_up], set_labels=['eccDNA Enriched', 'RNA-seq Up'], ax=ax)
    ax.set_title('eccDNA Enriched vs RNA-seq Upregulated')
    
    ax = axes[0, 1]
    venn2([eccdna_depleted, degs_down], set_labels=['eccDNA Depleted', 'RNA-seq Down'], ax=ax)
    ax.set_title('eccDNA Depleted vs RNA-seq Downregulated')
    
    ax = axes[1, 0]
    venn2([eccdna_enriched, degs_down], set_labels=['eccDNA Enriched', 'RNA-seq Down'], ax=ax)
    ax.set_title('eccDNA Enriched vs RNA-seq Downregulated (Inverse)')
    
    ax = axes[1, 1]
    venn2([eccdna_depleted, degs_up], set_labels=['eccDNA Depleted', 'RNA-seq Up'], ax=ax)
    ax.set_title('eccDNA Depleted vs RNA-seq Upregulated (Inverse)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_venn_diagrams.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_dir

def create_sample_comparison_heatmap(all_sample_results):
    """创建样本间比较的热图"""
    # 提取FC=2时的数据进行样本间比较
    fc2_data = []
    for sample, df in all_sample_results.items():
        fc2_row = df[df['fc_threshold'] == 2].iloc[0]
        fc2_data.append({
            'Sample': sample,
            'Enriched∩Up': fc2_row['n_enriched_up'],
            'Enriched∩Down': fc2_row['n_enriched_down'],
            'Depleted∩Up': fc2_row['n_depleted_up'],
            'Depleted∩Down': fc2_row['n_depleted_down'],
            'Enriched-Up OR': fc2_row.get('enriched_up_OR', np.nan),
            'Depleted-Down OR': fc2_row.get('depleted_down_OR', np.nan)
        })
    
    comparison_df = pd.DataFrame(fc2_data)
    comparison_df.set_index('Sample', inplace=True)
    
    # 创建热图
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # 基因数量热图
    sns.heatmap(comparison_df[['Enriched∩Up', 'Enriched∩Down', 'Depleted∩Up', 'Depleted∩Down']], 
                annot=True, fmt='d', cmap='YlOrRd', ax=ax1)
    ax1.set_title('Number of Overlapping Genes by Sample (FC=2)')
    ax1.set_xlabel('Gene Categories')
    
    # Odds Ratio热图
    sns.heatmap(comparison_df[['Enriched-Up OR', 'Depleted-Down OR']], 
                annot=True, fmt='.2f', cmap='RdBu_r', center=1, ax=ax2)
    ax2.set_title('Odds Ratios by Sample (FC=2)')
    ax2.set_xlabel('Association Types')
    
    plt.tight_layout()
    plt.savefig('sample_comparison_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return comparison_df

def analyze_cross_sample_consistency(eccdna_df, degs_df, samples):
    """分析跨样本的一致性"""
    # 收集每个样本的基因集
    sample_genes = {
        'enriched_up': {},
        'enriched_down': {},
        'depleted_up': {},
        'depleted_down': {}
    }
    
    for sample in samples:
        sample_data = eccdna_df[eccdna_df['Sample'] == sample]
        results, _, _, _, _, _ = analyze_eccdna_degs_overlap(sample_data, degs_df, fc_threshold=2)
        
        for category in sample_genes.keys():
            sample_genes[category][sample] = results[category]
    
    # 计算一致性统计
    consistency_stats = {}
    for category in sample_genes.keys():
        all_genes = set()
        for genes in sample_genes[category].values():
            all_genes.update(genes)
        
        # 计算每个基因在多少个样本中出现
        gene_frequency = {}
        for gene in all_genes:
            count = sum(1 for sample_set in sample_genes[category].values() if gene in sample_set)
            gene_frequency[gene] = count
        
        consistency_stats[category] = {
            'total_unique_genes': len(all_genes),
            'genes_in_all_samples': sum(1 for count in gene_frequency.values() if count == len(samples)),
            'genes_in_majority': sum(1 for count in gene_frequency.values() if count > len(samples)/2),
            'gene_frequency': gene_frequency
        }
    
    return consistency_stats

def generate_sample_report(sample_name, gradient_df, stats_fc2, output_dir):
    """为单个样本生成详细报告"""
    report_path = os.path.join(output_dir, f'{sample_name}_analysis_report.txt')
    
    with open(report_path, 'w') as f:
        f.write(f"=== Sample {sample_name} Analysis Report ===\n\n")
        
        # 基本统计（FC=2）
        f.write("1. BASIC STATISTICS (FC=2)\n")
        f.write("-" * 50 + "\n")
        f.write(f"eccDNA enriched genes: {stats_fc2['n_eccdna_enriched']}\n")
        f.write(f"eccDNA depleted genes: {stats_fc2['n_eccdna_depleted']}\n")
        f.write(f"Upregulated DEGs: {stats_fc2['n_degs_up']}\n")
        f.write(f"Downregulated DEGs: {stats_fc2['n_degs_down']}\n\n")
        
        f.write("2. OVERLAP STATISTICS (FC=2)\n")
        f.write("-" * 50 + "\n")
        f.write(f"Enriched ∩ Upregulated: {stats_fc2['n_enriched_up']}\n")
        f.write(f"Enriched ∩ Downregulated: {stats_fc2['n_enriched_down']}\n")
        f.write(f"Depleted ∩ Upregulated: {stats_fc2['n_depleted_up']}\n")
        f.write(f"Depleted ∩ Downregulated: {stats_fc2['n_depleted_down']}\n\n")
        
        # 梯度分析总结
        f.write("3. GRADIENT FC ANALYSIS\n")
        f.write("-" * 50 + "\n")
        f.write(gradient_df.to_string(index=False))
        f.write("\n")

def main():
    # 文件路径
    eccdna_file = "merged_gene_eccdna_enrichment_results.csv"
    degs_file = "DEGs_GBM.GEPIA2.Log2FC_1.qValue_0.01.txt"
    
    # 加载数据
    print("Loading data...")
    eccdna_df = load_eccdna_data(eccdna_file)
    degs_df = load_degs_data(degs_file)
    
    # 获取所有样本
    samples = eccdna_df['Sample'].unique()
    print(f"Found {len(samples)} samples: {', '.join(samples)}")
    
    # 创建主输出目录
    os.makedirs("eccdna_deg_analysis_results", exist_ok=True)
    os.chdir("eccdna_deg_analysis_results")
    
    # 存储所有样本的结果
    all_sample_results = {}
    all_sample_stats = {}
    
    # 1. 对每个样本进行独立分析
    print("\n=== Individual Sample Analysis ===")
    for sample in samples:
        print(f"\nAnalyzing sample: {sample}")
        
        # 梯度FC分析
        gradient_results = perform_gradient_fc_analysis_single_sample(
            eccdna_df, degs_df, sample, fc_thresholds=[1, 2, 4, 6])
        all_sample_results[sample] = gradient_results
        
        # 获取FC=2的详细统计
        sample_data = eccdna_df[eccdna_df['Sample'] == sample]
        _, stats_fc2, _, _, _, _ = analyze_eccdna_degs_overlap(sample_data, degs_df, fc_threshold=2)
        all_sample_stats[sample] = stats_fc2
        
        # 创建样本特定的图表和报告
        output_dir = plot_sample_specific_analysis(sample, gradient_results, eccdna_df, degs_df)
        generate_sample_report(sample, gradient_results, stats_fc2, output_dir)
        
        # 保存样本特定的梯度分析结果
        gradient_results.to_csv(os.path.join(output_dir, f'{sample}_gradient_results.csv'), index=False)
    
    # 2. 创建样本间比较分析
    print("\n=== Cross-Sample Comparison ===")
    comparison_df = create_sample_comparison_heatmap(all_sample_results)
    comparison_df.to_csv('sample_comparison_fc2.csv')
    
    # 3. 分析跨样本一致性
    print("\n=== Cross-Sample Consistency Analysis ===")
    consistency_stats = analyze_cross_sample_consistency(eccdna_df, degs_df, samples)
    
    # 4. 生成综合报告
    with open('comprehensive_multi_sample_report.txt', 'w') as f:
        f.write("=== Comprehensive Multi-Sample eccDNA-DEG Analysis Report ===\n\n")
        
        # 样本概览
        f.write("1. SAMPLE OVERVIEW\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total samples analyzed: {len(samples)}\n")
        f.write(f"Samples: {', '.join(samples)}\n\n")
        
        # 样本间比较（FC=2）
        f.write("2. CROSS-SAMPLE COMPARISON (FC=2)\n")
        f.write("-" * 50 + "\n")
        f.write(comparison_df.to_string())
        f.write("\n\n")
        
        # 一致性分析
        f.write("3. CROSS-SAMPLE CONSISTENCY\n")
        f.write("-" * 50 + "\n")
        for category, stats in consistency_stats.items():
            f.write(f"\n{category}:\n")
            f.write(f"  Total unique genes across all samples: {stats['total_unique_genes']}\n")
            f.write(f"  Genes found in all samples: {stats['genes_in_all_samples']}\n")
            f.write(f"  Genes found in majority of samples: {stats['genes_in_majority']}\n")
            
            # 列出在所有样本中都出现的基因
            if stats['genes_in_all_samples'] > 0:
                consistent_genes = [g for g, count in stats['gene_frequency'].items() 
                                  if count == len(samples)]
                f.write(f"  Consistent genes: {', '.join(sorted(consistent_genes)[:10])}")
                if len(consistent_genes) > 10:
                    f.write(f"... and {len(consistent_genes)-10} more")
                f.write("\n")
        
        # 样本特异性分析
        f.write("\n4. SAMPLE-SPECIFIC PATTERNS\n")
        f.write("-" * 50 + "\n")
        for sample in samples:
            fc2_data = all_sample_results[sample][all_sample_results[sample]['fc_threshold'] == 2].iloc[0]
            f.write(f"\n{sample}:\n")
            f.write(f"  Enriched-Up association: n={fc2_data['n_enriched_up']}")
            if 'enriched_up_OR' in fc2_data and not pd.isna(fc2_data['enriched_up_OR']):
                f.write(f", OR={fc2_data['enriched_up_OR']:.2f}, p={fc2_data['enriched_up_pvalue']:.3e}")
            f.write("\n")
            f.write(f"  Depleted-Down association: n={fc2_data['n_depleted_down']}")
            if 'depleted_down_OR' in fc2_data and not pd.isna(fc2_data['depleted_down_OR']):
                f.write(f", OR={fc2_data['depleted_down_OR']:.2f}, p={fc2_data['depleted_down_pvalue']:.3e}")
            f.write("\n")
    
    # 5. 创建样本间韦恩图（展示共享的基因）
    if len(samples) >= 2:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Cross-Sample Gene Overlap (FC=2)', fontsize=16)
        
        categories = ['enriched_up', 'enriched_down', 'depleted_up', 'depleted_down']
        titles = ['Enriched∩Up', 'Enriched∩Down', 'Depleted∩Up', 'Depleted∩Down']
        
        for idx, (cat, title) in enumerate(zip(categories, titles)):
            ax = axes[idx//2, idx%2]
            
            # 获取前两个样本的基因集
            sample1, sample2 = samples[0], samples[1]
            sample1_data = eccdna_df[eccdna_df['Sample'] == sample1]
            sample2_data = eccdna_df[eccdna_df['Sample'] == sample2]
            
            results1, _, _, _, _, _ = analyze_eccdna_degs_overlap(sample1_data, degs_df, fc_threshold=2)
            results2, _, _, _, _, _ = analyze_eccdna_degs_overlap(sample2_data, degs_df, fc_threshold=2)
            
            venn2([results1[cat], results2[cat]], 
                  set_labels=[sample1, sample2], ax=ax)
            ax.set_title(f'{title} Gene Overlap')
        
        plt.tight_layout()
        plt.savefig('cross_sample_gene_overlap_venn.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 6. 汇总所有样本的梯度分析结果
    all_gradients = pd.concat(all_sample_results.values(), ignore_index=True)
    all_gradients.to_csv('all_samples_gradient_analysis.csv', index=False)
    
    print("\n=== Analysis Complete! ===")
    print("\nOutput structure:")
    print("eccdna_deg_analysis_results/")
    print("├── sample_[name]_analysis/  (for each sample)")
    print("│   ├── [sample]_gradient_analysis.png")
    print("│   ├── [sample]_venn_diagrams.png")
    print("│   ├── [sample]_gradient_results.csv")
    print("│   └── [sample]_analysis_report.txt")
    print("├── sample_comparison_heatmap.png")
    print("├── cross_sample_gene_overlap_venn.png")
    print("├── sample_comparison_fc2.csv")
    print("├── all_samples_gradient_analysis.csv")
    print("└── comprehensive_multi_sample_report.txt")

if __name__ == "__main__":
    main()
