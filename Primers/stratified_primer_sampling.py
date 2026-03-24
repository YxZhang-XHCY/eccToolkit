import pandas as pd
import numpy as np

# 设置随机种子
np.random.seed(42)

# 读取数据
df = pd.read_csv('merged_primers_500bp_32t_2a_specific_unique_DR.tsv', sep='\t')

# 定义样本名称和重复类型的简化映射
sample_map = {
    '3SEP_LR': '3LR',
    'Circel-Seq_LR': 'CLR',
    'Circel-Seq_SR': 'CSR',
    'MMC-seq': 'MMC'
}

repeat_map = {
    'Other': 'O',
    'Direct-Repeat': 'DR',
    'Both': 'B',
    'Inverted-Repeat': 'IR'
}

# 定义每个重复类型需要抽取的样本数
sample_sizes = {
    'Other': 72,
    'Direct-Repeat': 80,
    'Both': 20,
    'Inverted-Repeat': 20
}

# 按sample分组进行抽样
sampled_dfs = []
for sample in df['sample'].unique():
    sample_df = df[df['sample'] == sample]
    
    # 对每个重复类型进行抽样
    sampled_by_repeat = []
    for repeat_class, size in sample_sizes.items():
        repeat_df = sample_df[sample_df['Repeat_Class'] == repeat_class]
        if len(repeat_df) > 0:
            # 如果数据量不足，则全部保留
            if len(repeat_df) < size:
                sampled = repeat_df
            else:
                sampled = repeat_df.sample(n=size, random_state=42)
            sampled_by_repeat.append(sampled)
    
    # 合并该样本的所有抽样结果
    if sampled_by_repeat:
        sampled_dfs.append(pd.concat(sampled_by_repeat))

# 合并所有结果
result_df = pd.concat(sampled_dfs, ignore_index=True)

# 创建NewID列
result_df['NewID'] = result_df.apply(lambda x: f"{sample_map[x['sample']]}_{repeat_map[x['Repeat_Class']]}", axis=1)

# 为每个组合添加自增编号
for sample in sample_map.values():
    for repeat in repeat_map.values():
        mask = result_df['NewID'] == f"{sample}_{repeat}"
        count = mask.sum()
        if count > 0:
            result_df.loc[mask, 'NewID'] = [f"{sample}_{repeat}{i+1:03d}" for i in range(count)]

# 保存主要结果
result_df.to_csv('sampled_primers_with_newid.tsv', sep='\t', index=False)

# 创建引物对文件
primer_pairs = []
for _, row in result_df.iterrows():
    primer_pairs.extend([
        (f"{row['NewID']}_F", row['FORWARD_PRIMER']),
        (f"{row['NewID']}_R", row['REVERSE_PRIMER'])
    ])

primer_df = pd.DataFrame(primer_pairs, columns=['PrimerID', 'Sequence'])
primer_df.to_csv('primer_pairs.csv', index=False)

# 打印统计信息
print("\n抽样统计:")
for sample in df['sample'].unique():
    print(f"\n{sample}:")
    sample_data = result_df[result_df['sample'] == sample]
    print(sample_data['Repeat_Class'].value_counts())

print("\n总抽样数:", len(result_df))
