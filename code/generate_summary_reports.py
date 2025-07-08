import argparse
import logging
import os
import pandas as pd

# --- 配置 (Configuration) ---
logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(asctime)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def analyze_sample_counts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Counts occurrences of each sample in the DataFrame.
    中文：统计 DataFrame 中每个样本的出现次数。

    Args:
        df (pd.DataFrame): The input DataFrame containing a 'sample' column. / 包含 'sample' 列的输入 DataFrame。

    Returns:
        pd.DataFrame: A DataFrame with 'Sample' and 'Count' columns. / 包含 'Sample' 和 'Count' 列的 DataFrame。
    """
    if "sample" not in df.columns:
        logging.error("Input DataFrame is missing the 'sample' column.")
        return pd.DataFrame()

    sample_counts = df["sample"].value_counts().reset_index()
    sample_counts.columns = ["Sample", "Count"]
    return sample_counts


def analyze_length_distribution(df: pd.DataFrame) -> pd.DataFrame:
    """
    Analyzes the distribution of lengths for each sample, categorizing them into bins.
    中文：分析每个样本的长度分布，并将其分类到不同的区间。

    Args:
        df (pd.DataFrame): The input DataFrame with 'sample' and 'Length' columns. / 包含 'sample' 和 'Length' 列的输入 DataFrame。

    Returns:
        pd.DataFrame: A DataFrame detailing the length distribution and percentages. / 详细描述长度分布和百分比的 DataFrame。
    """
    if "Length" not in df.columns or "sample" not in df.columns:
        logging.error("Input DataFrame is missing 'Length' or 'sample' column.")
        return pd.DataFrame()

    # Ensure 'Length' column is numeric, coercing errors to NaN / 确保 'Length' 列是数值类型，无法转换的设为 NaN
    df["Length"] = pd.to_numeric(df["Length"], errors="coerce")
    df.dropna(
        subset=["Length"], inplace=True
    )  # Drop rows where length could not be converted / 删除长度无法转换的行
    df["Length"] = df["Length"].astype(int)

    # Define bins and labels for length categories / 定义长度区间的边界和标签
    bins = [0, 1000, 10000, float("inf")]
    labels = ["0-1000", "1001-10000", ">10000"]

    # Use pd.cut for efficient binning / 使用 pd.cut 高效地进行分箱
    df["length_category"] = pd.cut(df["Length"], bins=bins, labels=labels, right=True)

    # Group by sample and length category / 按样本和长度类别进行分组
    distribution = (
        df.groupby(["sample", "length_category"])
        .size()
        .unstack(fill_value=0)
        .stack()
        .reset_index(name="Count")
    )
    distribution.columns = ["Sample", "Length Range", "Count"]

    # Calculate percentages / 计算百分比
    sample_totals = df.groupby("sample").size().rename("Total")
    distribution = distribution.merge(sample_totals, left_on="Sample", right_index=True)
    distribution["Percentage"] = (
        distribution["Count"] / distribution["Total"] * 100
    ).round(2)

    return distribution.sort_values(["Sample", "Length Range"])


def main():
    """
    Main function to parse arguments and run the analysis.
    中文: 解析参数并运行分析的主函数。
    """
    parser = argparse.ArgumentParser(
        description="Generate summary statistics from a combined results file. / 从合并的结果文件生成摘要统计信息。",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        help="Path to the input TSV/CSV file (e.g., Combined_Basic_Results.tsv). / 输入的 TSV/CSV 文件路径。",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help="Prefix for output files (e.g., 'MyProject'). / 输出文件的前缀（例如 'MyProject'）。",
    )
    parser.add_argument(
        "--sep",
        default="\t",
        help="Separator for the input file. Default is tab. Use ',' for CSV. / 输入文件的分隔符。默认为制表符，CSV文件请使用','。",
    )

    args = parser.parse_args()

    # --- Read Input File / 读取输入文件 ---
    if not os.path.exists(args.input_file):
        logging.error(f"Input file not found: {args.input_file}")
        return

    logging.info(f"Reading input file: {args.input_file}")
    try:
        df = pd.read_csv(args.input_file, sep=args.sep)
    except Exception as e:
        logging.error(f"Failed to read file: {e}")
        return

    # --- Run Analyses / 运行分析 ---
    logging.info("Analyzing sample counts...")
    sample_counts_df = analyze_sample_counts(df)

    logging.info("Analyzing length distribution...")
    length_distribution_df = analyze_length_distribution(df)

    # --- Save and Print Results / 保存并打印结果 ---
    if not sample_counts_df.empty:
        output_path = f"{args.output_prefix}_sample_counts.csv"
        sample_counts_df.to_csv(output_path, index=False)
        logging.info(f"Sample counts saved to {output_path}")
        print("\n--- Sample Counts ---")
        print(sample_counts_df.to_string(index=False))

    if not length_distribution_df.empty:
        output_path = f"{args.output_prefix}_length_distribution.csv"
        length_distribution_df.to_csv(output_path, index=False)
        logging.info(f"Length distribution saved to {output_path}")
        print("\n--- Length Distribution by Sample ---")
        # Format for printing / 格式化输出
        for sample, group in length_distribution_df.groupby("Sample"):
            print(f"\nSample: {sample}")
            print(group[["Length Range", "Count", "Percentage"]].to_string(index=False))

    logging.info("✔ Analysis complete!")


if __name__ == "__main__":
    main()
