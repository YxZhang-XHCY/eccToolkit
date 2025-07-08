import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# === 读取数据 ===
csv_file = 'TE.ALl.Merged.csv'   # 如果文件名不同请修改这里
df = pd.read_csv(csv_file)
sample_col = df.columns[0]       # 自动识别第一列为Sample

# === 结果表保存 ===
gap_table = []

# === 循环每个样本 ===
for sample, group in df.groupby(sample_col):
    motif_counts = group['motif'].value_counts()
    if len(motif_counts) < 2:
        continue

    # 计算相邻motif数量差值，找最大gap（断层）
    diffs = motif_counts.values[:-1] - motif_counts.values[1:]
    max_gap_idx = np.argmax(diffs)
    gap_size = diffs[max_gap_idx]

    # 输出断层相关信息
    gap_table.append([
        sample,
        max_gap_idx+1,  # 断层发生在rank max_gap_idx+1 到 max_gap_idx+2之间
        gap_size,
        motif_counts.index[max_gap_idx], motif_counts.iloc[max_gap_idx],
        motif_counts.index[max_gap_idx+1], motif_counts.iloc[max_gap_idx+1]
    ])

    # === 画rank plot图 ===
    fig, ax = plt.subplots(figsize=(10,5))
    ax.plot(np.arange(1, len(motif_counts)+1), motif_counts.values, marker='o', lw=1)
    ax.set_yscale('log')
    ax.set_xlabel('Motif Rank')
    ax.set_ylabel('Motif Count (log scale)')
    ax.set_title(f'Motif Abundance Rank Plot: {sample}')
    ax.grid(True, which="both", ls="--")
    # 标记断层点
    ax.axvline(max_gap_idx+1, color='red', ls='--', label='Gap')
    # 标记点文本
    ax.annotate(f'{motif_counts.index[max_gap_idx]}\n{motif_counts.iloc[max_gap_idx]}',
                xy=(max_gap_idx+1, motif_counts.iloc[max_gap_idx]), 
                xytext=(max_gap_idx+5, motif_counts.iloc[max_gap_idx]), 
                arrowprops=dict(arrowstyle='->',color='red'), color='red', fontsize=10)
    ax.annotate(f'{motif_counts.index[max_gap_idx+1]}\n{motif_counts.iloc[max_gap_idx+1]}',
                xy=(max_gap_idx+2, motif_counts.iloc[max_gap_idx+1]), 
                xytext=(max_gap_idx+7, motif_counts.iloc[max_gap_idx+1]), 
                arrowprops=dict(arrowstyle='->',color='red'), color='red', fontsize=10)
    plt.legend()
    # 保存图片
    fig.savefig(f'{sample}_rankplot.png', dpi=200)
    plt.close()

# === 输出断层表 ===
gap_df = pd.DataFrame(gap_table, columns=['Sample','Gap Rank','Gap Size',
                                          'Before Motif','Before Count','After Motif','After Count'])
gap_df.to_csv('motif_gap_summary.csv', index=False)
print('分析完成！每个样本的rank plot图已保存为PNG，断层信息已保存为 motif_gap_summary.csv')
