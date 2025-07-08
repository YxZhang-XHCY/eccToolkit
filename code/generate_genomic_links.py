import pandas as pd
from itertools import combinations
import random
from tqdm import tqdm


def generate_random_rgb():
    """生成一个随机 RGB 字符串"""
    return f"{random.randint(0,255)},{random.randint(0,255)},{random.randint(0,255)}"


# 读取输入文件
df = pd.read_csv("HeLa_rep1.Mecc.Clean.csv")

# 初始化结果列表
output_rows = []

# tqdm 包裹 groupby，显示总进度
for name, group in tqdm(df.groupby("eName"), desc="Processing groups"):
    color = generate_random_rgb()
    # 成对组合加 tqdm（总数为 nC2）
    pairs = list(combinations(group.iterrows(), 2))
    for (i1, row1), (i2, row2) in pairs:
        output_rows.append(
            [
                row1["eChr"],
                row1["eStart"],
                row1["eEnd"],
                row2["eChr"],
                row2["eStart"],
                row2["eEnd"],
                color,
            ]
        )

# 转换为 DataFrame 并输出
output_df = pd.DataFrame(
    output_rows, columns=["Chr1", "Start1", "End1", "Chr2", "Start2", "End2", "RGB"]
)
output_df.to_csv("HeLa_rep1.Mecc.Links.RandomColor.tsv", sep="\t", index=False)

print("✅ Saved: HeLa_rep1.Mecc.Links.RandomColor.tsv")
