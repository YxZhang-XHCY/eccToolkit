import pandas as pd

# 读取所有TSV文件并添加sample列
dfs = []

# 3SEP_LR
df1 = pd.read_csv('3SEP_LR_primers_500bp_32t_2a.tsv', sep='\t')
df1['sample'] = '3SEP_LR'
dfs.append(df1)

# Circel-Seq_LR
df2 = pd.read_csv('Circel-Seq_LR_primers_500bp_32t_2a.tsv', sep='\t')
df2['sample'] = 'Circel-Seq_LR'
dfs.append(df2)

# Circel-Seq_SR
df3 = pd.read_csv('eccDNA_Circel-Seq_NGS_primers_500bp_32t_2a.tsv', sep='\t')
df3['sample'] = 'Circel-Seq_SR'
dfs.append(df3)

# MMC-seq
df4 = pd.read_csv('MMC-seq_primers_500bp_32t_2a.tsv', sep='\t')
df4['sample'] = 'MMC-seq'
dfs.append(df4)

# 合并所有数据框
merged_df = pd.concat(dfs, ignore_index=True)

# 保存完整的合并文件
merged_df.to_csv('merged_primers_500bp_32t_2a.tsv', sep='\t', index=False)

# 筛选特异性引物并去重
specific_df = merged_df[
    (merged_df['FORWARD_SPECIFIC'] == True) & 
    (merged_df['REVERSE_SPECIFIC'] == True)
].drop_duplicates(subset=['Name'], keep='first')

# 保存筛选后的文件
specific_df.to_csv('merged_primers_500bp_32t_2a_specific_unique.tsv', sep='\t', index=False)

# 打印统计信息
print(f"原始合并数据行数: {len(merged_df)}")
print(f"筛选后的数据行数: {len(specific_df)}")
print("\n每个样本中特异性引物的数量:")
print(specific_df['sample'].value_counts())
