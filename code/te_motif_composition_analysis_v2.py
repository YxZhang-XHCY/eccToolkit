import pandas as pd
from collections import Counter
import argparse
import matplotlib.pyplot as plt

def analyze_motif_composition(csvfile):
    df = pd.read_csv(csvfile)
    results = []

    # 按Sample, Class, seqname分组统计motif
    group_cols = ['Sample', 'Class', 'seqname']
    motif_per_seq = df.groupby(group_cols)['motif'].apply(list).reset_index()

    # 标记“单一”与“复合”
    motif_per_seq['motif_count'] = motif_per_seq['motif'].apply(len)
    motif_per_seq['type'] = motif_per_seq['motif_count'].apply(lambda x: 'single' if x == 1 else 'composite')

    # 合并标记结果回原表
    merged = df.merge(motif_per_seq[group_cols + ['motif_count', 'type']], on=group_cols)

    # 单一/复合转座子的所有行
    single_rows = merged[merged['type'] == 'single'].copy()
    composite_rows = merged[merged['type'] == 'composite'].copy()

    single_rows.to_csv('single_motif_eccDNA.csv', index=False)
    composite_rows.to_csv('composite_motif_eccDNA.csv', index=False)
    print('已输出单一转座子文件 single_motif_eccDNA.csv 和复合转座子文件 composite_motif_eccDNA.csv')

    # ---- motif_percent 分析 ----
    if 'motif_percent' in single_rows.columns:
        motif_percent = pd.to_numeric(single_rows['motif_percent'], errors='coerce').dropna()
        desc = motif_percent.describe(percentiles=[0.05,0.25,0.5,0.75,0.95])
        desc.to_csv('single_motif_percent_stats.csv')
        print('已输出单一转座子的motif_percent描述性统计 single_motif_percent_stats.csv')
        # 可选：画分布图
        plt.figure(figsize=(6,4))
        plt.hist(motif_percent, bins=30, alpha=0.7)
        plt.xlabel('motif_percent (%)')
        plt.ylabel('Number of eccDNAs')
        plt.title('Distribution of motif_percent in single-motif eccDNA')
        plt.tight_layout()
        plt.savefig('single_motif_percent_hist.png', dpi=150)
        print('已输出单一转座子的motif_percent直方图 single_motif_percent_hist.png')

    # ---- 总体比例和top10分析 ----
    for (sample, class_), group in motif_per_seq.groupby(['Sample', 'Class']):
        total = len(group)
        n_single = (group['type'] == 'single').sum()
        n_composite = (group['type'] == 'composite').sum()
        single_ratio = n_single / total if total else 0
        composite_ratio = n_composite / total if total else 0

        # 单一转座子motif统计
        single_motifs = [motif_list[0] for motif_list in group.loc[group['type'] == 'single', 'motif']]
        top10_single = Counter(single_motifs).most_common(10)

        # 复合转座子motif组合统计
        composite_lists = group.loc[group['type'] == 'composite', 'motif']
        composite_counts = composite_lists.apply(len)
        mean_motif_per_composite = composite_counts.mean() if len(composite_counts) else 0
        # 组合字符串排序，便于统计top组合（防止有float/None干扰）
        composite_combos = composite_lists.apply(
            lambda x: '|'.join(sorted([str(m) for m in x if pd.notnull(m)]))
        )
        top10_combos = Counter(composite_combos).most_common(10)

        results.append({
            'Sample': sample,
            'Class': class_,
            'total_seq': total,
            'single_count': n_single,
            'composite_count': n_composite,
            'single_ratio': f'{single_ratio:.2%}',
            'composite_ratio': f'{composite_ratio:.2%}',
            'top10_single_motif': top10_single,
            'mean_motif_per_composite': round(mean_motif_per_composite, 2),
            'top10_composite_combo': top10_combos
        })

    # 输出分析结果
    resdf = pd.DataFrame(results)
    resdf.to_csv('motif_composition_summary.csv', index=False)
    print('已输出 motif_composition_summary.csv，包含各Sample/Class的比例与top信息。')

    # 附加详细top10信息输出，txt文本更直观
    with open('motif_composition_detail.txt', 'w') as f:
        for r in results:
            f.write(f"Sample: {r['Sample']}\tClass: {r['Class']}\n")
            f.write(f"  Total eccDNA: {r['total_seq']}\n")
            f.write(f"  单一转座子数量: {r['single_count']} ({r['single_ratio']})\n")
            f.write(f"  复合转座子数量: {r['composite_count']} ({r['composite_ratio']})\n")
            f.write(f"  单一转座子Top10: {r['top10_single_motif']}\n")
            f.write(f"  复合转座子平均Motif数量: {r['mean_motif_per_composite']}\n")
            f.write(f"  复合转座子Top10组合: {r['top10_composite_combo']}\n")
            f.write('-'*60+'\n')
    print('已输出 motif_composition_detail.txt，包含详细分组和top10信息。')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="分析eccDNA中单一/复合转座子的比例及motif组合")
    parser.add_argument('--input', '-i', required=True, help='输入csv文件')
    args = parser.parse_args()
    analyze_motif_composition(args.input)
