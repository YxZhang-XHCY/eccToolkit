#!/usr/bin/env python3
"""
TE数据处理脚本
功能：
1. 使用Mecc.seq_length.table.simplified_modified.csv补充TE.ALl.Merged.csv中的空值
2. 重新计算motif_percent
3. 计算anno_Percent = (end-start+1) / seq_length * 100
4. 填充空的chr, seq_start, seq_end为"-"
5. 统计处理结果
"""

import pandas as pd
import numpy as np
import sys


def process_te_data(te_file, ref_file, output_file):
    """
    处理TE数据
    """
    # 读取数据
    print("正在读取数据文件...")
    te_df = pd.read_csv(te_file)
    ref_df = pd.read_csv(ref_file)

    print(f"TE文件原始行数: {len(te_df)}")
    print(f"参考文件行数: {len(ref_df)}")

    # 统计原始空值情况
    original_null_count = te_df["seq_length"].isnull().sum()
    print(f"seq_length字段原始空值数量: {original_null_count}")

    # 创建参考数据的索引字典，基于Sample, Class, seqname
    ref_dict = {}
    for _, row in ref_df.iterrows():
        key = (row["Sample"], row["Class"], row["seqname"])
        ref_dict[key] = row["seq_length"]

    print(f"参考字典包含 {len(ref_dict)} 个条目")

    # 填充空值
    filled_count = 0
    for idx, row in te_df.iterrows():
        if pd.isnull(row["seq_length"]) or row["seq_length"] == "":
            key = (row["Sample"], row["Class"], row["seqname"])
            if key in ref_dict:
                te_df.at[idx, "seq_length"] = ref_dict[key]
                filled_count += 1

    print(f"成功填充 {filled_count} 行的seq_length数据")

    # 检查填充后的空值情况
    remaining_null_count = te_df["seq_length"].isnull().sum()
    print(f"填充后seq_length字段剩余空值数量: {remaining_null_count}")

    # 确保数值类型转换
    print("正在处理数据类型转换...")
    numeric_columns = ["start", "end", "seq_length", "motif_length"]
    for col in numeric_columns:
        te_df[col] = pd.to_numeric(te_df[col], errors="coerce")

    # 重新计算motif_percent
    print("正在重新计算motif_percent...")
    te_df["motif_percent"] = np.where(
        (pd.notnull(te_df["motif_length"]))
        & (pd.notnull(te_df["seq_length"]))
        & (te_df["seq_length"] > 0),
        (te_df["motif_length"] / te_df["seq_length"] * 100).round(2),
        te_df["motif_percent"],
    )

    # 计算anno_Percent = (end-start+1) / seq_length * 100
    print("正在计算anno_Percent...")
    te_df["anno_Percent"] = np.where(
        (pd.notnull(te_df["start"]))
        & (pd.notnull(te_df["end"]))
        & (pd.notnull(te_df["seq_length"]))
        & (te_df["seq_length"] > 0),
        ((te_df["end"] - te_df["start"] + 1) / te_df["seq_length"] * 100).round(2),
        np.nan,
    )

    # 填充chr, seq_start, seq_end的空值为"-"
    print("正在填充chr, seq_start, seq_end的空值...")
    chr_cols = ["chr", "seq_start", "seq_end"]
    chr_filled_count = 0

    for col in chr_cols:
        if col in te_df.columns:
            null_mask = te_df[col].isnull() | (te_df[col] == "")
            chr_filled_count += null_mask.sum()
            te_df[col] = te_df[col].fillna("-")
            te_df[col] = te_df[col].replace("", "-")

    print(f"chr/seq_start/seq_end字段共填充 {chr_filled_count} 个空值为'-'")

    # 保存结果
    print(f"正在保存结果到 {output_file}...")
    te_df.to_csv(output_file, index=False)

    # 输出统计信息
    print("\n=== 处理结果统计 ===")
    print(f"总行数: {len(te_df)}")
    print(f"seq_length原始空值: {original_null_count}")
    print(f"seq_length填充行数: {filled_count}")
    print(f"seq_length剩余空值: {remaining_null_count}")
    print(f"chr/seq_start/seq_end填充空值数: {chr_filled_count}")
    print(f"新增anno_Percent列，计算公式: (end-start+1)/seq_length*100")
    print(f"重新计算了motif_percent列")

    # 显示前几行结果
    print("\n=== 处理后数据预览 ===")
    print(te_df.head())

    return te_df


def main():
    """
    主函数
    """
    # 文件路径
    te_file = "TE.ALl.Merged.csv"
    ref_file = "Mecc.seq_length.table.simplified_modified.csv"
    output_file = "TE.ALl.Merged.processed.csv"

    try:
        result_df = process_te_data(te_file, ref_file, output_file)
        print(f"\n处理完成！结果已保存到: {output_file}")

    except FileNotFoundError as e:
        print(f"错误：找不到文件 {e}")
        print("请确保以下文件存在于当前目录：")
        print("- TE.ALl.Merged.csv")
        print("- Mecc.seq_length.table.simplified_modified.csv")

    except Exception as e:
        print(f"处理过程中发生错误: {e}")


if __name__ == "__main__":
    main()
