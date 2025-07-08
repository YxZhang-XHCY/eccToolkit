import pandas as pd
import glob
import os

# 获取所有csv文件
csv_files = glob.glob("*Merged.Uecc.csv")

dfs = []
for f in csv_files:
    # 提取样本名作为新列，比如文件名为 U87_Bulk_Ctrl_rep01_Final.Merged.Uecc.csv，sample就是 U87_Bulk_Ctrl_rep01
    basename = os.path.basename(f)
    sample = basename.split("_Final")[0]
    # 读取csv并加上 sample 列
    df = pd.read_csv(f)
    df['sample'] = sample
    dfs.append(df)

# 合并所有DataFrame
df_merged = pd.concat(dfs, ignore_index=True)

# 输出到新文件
df_merged.to_csv("AllSamples.Merged.Uecc.csv", index=False)
