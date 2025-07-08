#!/usr/bin/env python3
"""
merge_csv_with_sample.py
合并目录中的CSV文件并添加样本信息
"""

import pandas as pd
import glob
import os
import argparse


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="合并目录中的CSV文件并添加样本信息",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
    python merge_csv_with_sample.py -i Uecccsv -o merged_results.csv
    python merge_csv_with_sample.py -i data/ -o output.csv -p "*.csv"
        """,
    )

    parser.add_argument("-i", "--input-dir", required=True, help="输入CSV文件目录")
    parser.add_argument("-o", "--output", required=True, help="输出合并后的CSV文件路径")
    parser.add_argument(
        "-p", "--pattern", default="*.csv", help="文件匹配模式 (default: *.csv)"
    )
    parser.add_argument(
        "--sample-column", default="sample", help="样本列名称 (default: sample)"
    )

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

    # 检查输入目录是否存在
    if not os.path.exists(args.input_dir):
        print(f"错误: 输入目录不存在: {args.input_dir}")
        return

    # 获取所有CSV文件路径
    csv_files = glob.glob(os.path.join(args.input_dir, args.pattern))

    if not csv_files:
        print(f"错误: 在目录 {args.input_dir} 中未找到匹配 {args.pattern} 的文件")
        return

    print(f"找到 {len(csv_files)} 个CSV文件")

    dfs = []
    for file in csv_files:
        try:
            # 取文件名作为样本名（去掉扩展名）
            sample_name = os.path.splitext(os.path.basename(file))[0]

            # 读取CSV文件
            df = pd.read_csv(file)

            # 添加样本列
            df[args.sample_column] = sample_name

            dfs.append(df)
            print(f"已处理: {file} ({len(df)} 行)")

        except Exception as e:
            print(f"警告: 无法读取文件 {file}: {e}")
            continue

    if not dfs:
        print("错误: 没有成功读取任何CSV文件")
        return

    # 合并所有DataFrame
    merged_df = pd.concat(dfs, ignore_index=True)

    # 创建输出目录
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 保存到新文件
    merged_df.to_csv(args.output, index=False)

    print(f"成功合并 {len(dfs)} 个文件")
    print(f"合并后数据: {len(merged_df)} 行, {len(merged_df.columns)} 列")
    print(f"结果已保存到: {args.output}")


if __name__ == "__main__":
    main()
