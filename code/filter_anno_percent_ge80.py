#!/usr/bin/env python3
"""
filter_anno_percent_ge80.py

筛选出输入CSV中anno_Percent大于等于80的行，输出到新文件。

用法:
    python filter_anno_percent_ge80.py input.csv output.csv
"""

import sys
import pandas as pd

def main():
    if len(sys.argv) != 3:
        print("用法: python filter_anno_percent_ge80.py input.csv output.csv")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # 读取数据
    df = pd.read_csv(input_file)
    # 过滤
    filtered = df[df['anno_Percent'] >= 80]
    # 保存
    filtered.to_csv(output_file, index=False)
    print(f"已筛选出anno_Percent >= 80的{len(filtered)}条记录，保存为{output_file}")

if __name__ == "__main__":
    main()
