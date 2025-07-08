import pandas as pd
import glob
import os

# 指定你的csv文件目录
csv_dir = "Uecccsv"

# 获取所有csv文件路径
csv_files = glob.glob(os.path.join(csv_dir, "*.csv"))

dfs = []
for file in csv_files:
    # 取文件名作为sample名（去掉扩展名即可）
    sample_name = os.path.splitext(os.path.basename(file))[0]
    # 读取csv文件
    df = pd.read_csv(file)
    # 添加sample列
    df['sample'] = sample_name
    dfs.append(df)

# 合并所有DataFrame
merged_df = pd.concat(dfs, ignore_index=True)

# 保存到新文件
merged_df.to_csv("Cell_All_Merged_Uecc_with_sample.csv", index=False)

print("Done! 合并后的文件为 All_Merged_Uecc_with_sample.csv")
