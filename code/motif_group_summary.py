import pandas as pd

# 读取生成的文件
df = pd.read_csv("TE.ALl.Merged.csv")

# 1. 第一列的非冗余元素（唯一值）
col1_name = df.columns[0]  # 第一列名，自动识别
unique_col1 = df[col1_name].unique()
print(f"第一列（{col1_name}）非冗余元素：")
for x in unique_col1:
    print(x)
print("共有元素数:", len(unique_col1))
print("-" * 40)

# 2. motif 列非冗余元素
unique_motif = df["motif"].dropna().unique()
print("motif 列非冗余元素：")
for x in unique_motif:
    print(x)
print("共有motif种类:", len(unique_motif))
print("-" * 40)

# 3. 按第一列分组统计motif种类和数量
for sample in unique_col1:
    subset = df[df[col1_name] == sample]
    motif_counts = subset["motif"].value_counts(dropna=False)
    print(f"Sample: {sample}")
    print(motif_counts)
    print("-" * 30)
