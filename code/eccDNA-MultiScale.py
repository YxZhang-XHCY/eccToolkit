import argparse
import pandas as pd
import subprocess
import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import warnings

warnings.filterwarnings("ignore")


def convert_csv_to_bed(csv_file, bed_file):
    """将 eccDNA CSV 文件转换为 BED 格式"""
    print(f"  Converting {csv_file} to BED format...")
    try:
        df = pd.read_csv(csv_file)
        if (
            "eChr" not in df.columns
            or "eStart" not in df.columns
            or "eEnd" not in df.columns
        ):
            raise ValueError(
                f"CSV file must contain eChr, eStart, and eEnd columns. Found: {list(df.columns)}"
            )
        df_bed = df[["eChr", "eStart", "eEnd"]]
        df_bed.to_csv(bed_file, sep="\t", header=False, index=False)
        print(f"  Converted CSV to BED: {bed_file} ({len(df_bed)} regions)")
    except Exception as e:
        print(f"Error converting CSV to BED: {e}")
        raise


def run_bedtools_makewindows(fai_file, window_size, windows_file):
    """使用 bedtools makewindows 生成窗口"""
    print(f"  Generating windows of size {window_size}bp...")
    command = ["bedtools", "makewindows", "-g", fai_file, "-w", str(window_size)]
    try:
        with open(windows_file, "w") as out_file:
            result = subprocess.run(
                command, stdout=out_file, stderr=subprocess.PIPE, check=True
            )
        print(f"  Generated windows: {windows_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running bedtools makewindows: {e.stderr.decode()}")
        raise


def run_bedtools_intersect(windows_file, eccdna_bed_file, result_file):
    """使用 bedtools intersect 计算窗口与 eccDNA 的交集"""
    print(f"  Running intersection analysis...")
    command = ["bedtools", "intersect", "-a", windows_file, "-b", eccdna_bed_file, "-c"]
    try:
        with open(result_file, "w") as out_file:
            result = subprocess.run(
                command, stdout=out_file, stderr=subprocess.PIPE, check=True
            )
        print(f"  Calculated intersections: {result_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running bedtools intersect: {e.stderr.decode()}")
        raise


def format_result(result_file, sample_name):
    """格式化结果文件为 CSV 格式，并添加样本名称"""
    result_df = pd.read_csv(result_file, sep="\t", header=None)
    result_df.columns = ["Chrom", "WindowStart", "WindowEnd", sample_name]
    return result_df


def merge_binned_data(binned_data_list):
    """将所有样本的binned数据合并到一起"""
    merged_data = binned_data_list[0]
    for data in binned_data_list[1:]:
        merged_data = pd.merge(
            merged_data, data, on=["Chrom", "WindowStart", "WindowEnd"], how="outer"
        )
    # 填充NaN为0
    merged_data = merged_data.fillna(0)
    return merged_data


def process_samples_for_window_size(
    sample_files, sample_names, fai_file, window_size, output_prefix
):
    """处理所有样本并生成特定窗口大小的矩阵"""
    binned_data_list = []

    # 只生成一次窗口文件（所有样本共用）
    windows_file = f"temp_windows_{window_size}.bed"
    print(f"  Generating genomic windows...")
    run_bedtools_makewindows(fai_file, window_size, windows_file)

    for csv_file, sample_name in zip(sample_files, sample_names):
        print(f"\n  Processing sample: {sample_name}")
        # 创建临时文件名
        eccdna_bed_file = f"temp_{sample_name}_eccDNA.bed"
        result_file = f"temp_{sample_name}_result_{window_size}.bed"

        try:
            # 转换为BED格式
            convert_csv_to_bed(csv_file, eccdna_bed_file)

            # 计算交集
            run_bedtools_intersect(windows_file, eccdna_bed_file, result_file)

            # 格式化结果
            formatted_result = format_result(result_file, sample_name)
            binned_data_list.append(formatted_result)

        finally:
            # 清理临时文件
            for temp_file in [eccdna_bed_file, result_file]:
                if os.path.exists(temp_file):
                    os.remove(temp_file)

    # 清理窗口文件
    if os.path.exists(windows_file):
        os.remove(windows_file)

    # 合并所有样本的数据
    print(f"\n  Merging data from all samples...")
    merged_data = merge_binned_data(binned_data_list)

    # 保存窗口化数据
    output_file = f"{output_prefix}_window_{window_size//1000}kb_matrix.csv"
    merged_data.to_csv(output_file, index=False)
    print(f"  Saved matrix for window size {window_size//1000}kb to {output_file}")

    return merged_data


def group_samples_by_replicate(merged_data, sample_names):
    """将样本按重复分组"""
    # 获取样本列（排除前三列：Chrom, WindowStart, WindowEnd）
    groups = {1: [], 2: [], 3: []}

    # 根据样本名称中的后缀分组
    for name in sample_names:
        if name.endswith("-1"):
            groups[1].append(name)
        elif name.endswith("-2"):
            groups[2].append(name)
        elif name.endswith("-3"):
            groups[3].append(name)

    return groups


def filter_zero_rows(merged_data, groups):
    """剔除三个组样本全部为0的行"""
    # 检查每个组是否有至少一个非零值
    mask = pd.Series([False] * len(merged_data))

    for group_samples in groups.values():
        if group_samples:  # 如果组不为空
            group_mask = (merged_data[group_samples] != 0).any(axis=1)
            mask = mask | group_mask

    filtered_data = merged_data[mask].copy()
    print(f"Filtered data: {len(merged_data)} rows -> {len(filtered_data)} rows")

    return filtered_data


def calculate_correlations(filtered_data, groups):
    """计算两两样本组之间的相关性"""
    correlations = {}

    # 获取所有样本组合
    group_pairs = [(1, 2), (1, 3), (2, 3)]

    for group1, group2 in group_pairs:
        if groups[group1] and groups[group2]:
            # 计算每个组的平均值
            mean1 = filtered_data[groups[group1]].mean(axis=1)
            mean2 = filtered_data[groups[group2]].mean(axis=1)

            # 计算Pearson相关性
            pearson_r, pearson_p = stats.pearsonr(mean1, mean2)

            # 计算Spearman相关性
            spearman_r, spearman_p = stats.spearmanr(mean1, mean2)

            correlations[f"Group{group1}_vs_Group{group2}"] = {
                "pearson_r": pearson_r,
                "pearson_p": pearson_p,
                "spearman_r": spearman_r,
                "spearman_p": spearman_p,
            }

    return correlations


def create_correlation_heatmaps(all_correlations, output_prefix):
    """创建相关性热图"""
    window_sizes = sorted(all_correlations.keys())
    group_pairs = ["Group1_vs_Group2", "Group1_vs_Group3", "Group2_vs_Group3"]

    # 准备数据矩阵
    pearson_matrix = np.zeros((len(window_sizes), len(group_pairs)))
    spearman_matrix = np.zeros((len(window_sizes), len(group_pairs)))

    for i, window_size in enumerate(window_sizes):
        for j, pair in enumerate(group_pairs):
            if pair in all_correlations[window_size]:
                pearson_matrix[i, j] = all_correlations[window_size][pair]["pearson_r"]
                spearman_matrix[i, j] = all_correlations[window_size][pair][
                    "spearman_r"
                ]

    # 创建热图
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Pearson相关性热图
    sns.heatmap(
        pearson_matrix,
        xticklabels=group_pairs,
        yticklabels=[f"{size//1000}kb" for size in window_sizes],
        annot=True,
        fmt=".3f",
        cmap="coolwarm",
        center=0,
        ax=ax1,
        vmin=-1,
        vmax=1,
    )
    ax1.set_title("Pearson Correlation Coefficients")
    ax1.set_xlabel("Sample Group Pairs")
    ax1.set_ylabel("Window Size")

    # Spearman相关性热图
    sns.heatmap(
        spearman_matrix,
        xticklabels=group_pairs,
        yticklabels=[f"{size//1000}kb" for size in window_sizes],
        annot=True,
        fmt=".3f",
        cmap="coolwarm",
        center=0,
        ax=ax2,
        vmin=-1,
        vmax=1,
    )
    ax2.set_title("Spearman Correlation Coefficients")
    ax2.set_xlabel("Sample Group Pairs")
    ax2.set_ylabel("Window Size")

    plt.tight_layout()
    plt.savefig(
        f"{output_prefix}_correlation_heatmaps.png", dpi=300, bbox_inches="tight"
    )
    plt.savefig(f"{output_prefix}_correlation_heatmaps.pdf", bbox_inches="tight")
    plt.close()

    # 保存热图数据
    pearson_df = pd.DataFrame(
        pearson_matrix,
        index=[f"{size//1000}kb" for size in window_sizes],
        columns=group_pairs,
    )
    spearman_df = pd.DataFrame(
        spearman_matrix,
        index=[f"{size//1000}kb" for size in window_sizes],
        columns=group_pairs,
    )

    pearson_df.to_csv(f"{output_prefix}_pearson_heatmap_data.csv")
    spearman_df.to_csv(f"{output_prefix}_spearman_heatmap_data.csv")

    print(f"Saved correlation heatmaps to {output_prefix}_correlation_heatmaps.png/pdf")


def create_correlation_line_plots(all_correlations, output_prefix):
    """创建相关性随窗口大小变化的折线图"""
    window_sizes = sorted(all_correlations.keys())
    window_labels = [f"{size//1000}" for size in window_sizes]

    # 准备数据
    data_dict = {
        "Group1_vs_Group2": {"pearson": [], "spearman": []},
        "Group1_vs_Group3": {"pearson": [], "spearman": []},
        "Group2_vs_Group3": {"pearson": [], "spearman": []},
    }

    for window_size in window_sizes:
        for pair in data_dict.keys():
            if pair in all_correlations[window_size]:
                data_dict[pair]["pearson"].append(
                    all_correlations[window_size][pair]["pearson_r"]
                )
                data_dict[pair]["spearman"].append(
                    all_correlations[window_size][pair]["spearman_r"]
                )
            else:
                data_dict[pair]["pearson"].append(np.nan)
                data_dict[pair]["spearman"].append(np.nan)

    # 创建折线图
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    # Pearson相关性折线图
    for pair, data in data_dict.items():
        ax1.plot(window_labels, data["pearson"], marker="o", label=pair, linewidth=2)

    ax1.set_xlabel("Window Size (kb)")
    ax1.set_ylabel("Pearson Correlation Coefficient")
    ax1.set_title("Pearson Correlation vs Window Size")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(-0.1, 1.1)

    # Spearman相关性折线图
    for pair, data in data_dict.items():
        ax2.plot(window_labels, data["spearman"], marker="s", label=pair, linewidth=2)

    ax2.set_xlabel("Window Size (kb)")
    ax2.set_ylabel("Spearman Correlation Coefficient")
    ax2.set_title("Spearman Correlation vs Window Size")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(-0.1, 1.1)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_correlation_trends.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{output_prefix}_correlation_trends.pdf", bbox_inches="tight")
    plt.close()

    # 保存折线图数据
    trend_data = pd.DataFrame(index=window_labels)
    for pair, data in data_dict.items():
        trend_data[f"{pair}_pearson"] = data["pearson"]
        trend_data[f"{pair}_spearman"] = data["spearman"]

    trend_data.to_csv(f"{output_prefix}_correlation_trends_data.csv")

    print(f"Saved correlation trends to {output_prefix}_correlation_trends.png/pdf")


def save_correlation_summary(all_correlations, output_prefix):
    """保存相关性分析的详细汇总"""
    summary_data = []

    for window_size in sorted(all_correlations.keys()):
        for pair, stats in all_correlations[window_size].items():
            summary_data.append(
                {
                    "Window_Size_kb": window_size // 1000,
                    "Comparison": pair,
                    "Pearson_R": stats["pearson_r"],
                    "Pearson_P": stats["pearson_p"],
                    "Spearman_R": stats["spearman_r"],
                    "Spearman_P": stats["spearman_p"],
                }
            )

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(f"{output_prefix}_correlation_summary.csv", index=False)
    print(f"Saved correlation summary to {output_prefix}_correlation_summary.csv")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="Multi-scale window analysis of eccDNA data with correlation analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python script.py -1 sample1-1.csv sample2-1.csv sample3-1.csv \\
                     -2 sample1-2.csv sample2-2.csv sample3-2.csv \\
                     -3 sample1-3.csv sample2-3.csv sample3-3.csv \\
                     -f genome.fai -o analysis_output

With custom window parameters:
    python script.py -1 s1-1.csv s2-1.csv -2 s1-2.csv s2-2.csv -3 s1-3.csv s2-3.csv \\
                     -f genome.fai -o output --min_window 5000 --max_window 200000 --window_step 5000
        """,
    )

    # 样本输入参数
    parser.add_argument(
        "-1",
        "--group1",
        nargs="+",
        required=True,
        help="CSV files for replicate group 1 (sample names will end with -1)",
    )
    parser.add_argument(
        "-2",
        "--group2",
        nargs="+",
        required=True,
        help="CSV files for replicate group 2 (sample names will end with -2)",
    )
    parser.add_argument(
        "-3",
        "--group3",
        nargs="+",
        required=True,
        help="CSV files for replicate group 3 (sample names will end with -3)",
    )

    # 其他必需参数
    parser.add_argument(
        "-f", "--fai", required=True, help="FAI file of the reference genome"
    )
    parser.add_argument(
        "-o",
        "--output_prefix",
        required=True,
        help="Output prefix for all generated files",
    )

    # 窗口参数
    parser.add_argument(
        "--min_window",
        type=int,
        default=10000,
        help="Minimum window size in bp (default: 10000)",
    )
    parser.add_argument(
        "--max_window",
        type=int,
        default=100000,
        help="Maximum window size in bp (default: 100000)",
    )
    parser.add_argument(
        "--window_step",
        type=int,
        default=10000,
        help="Window size step in bp (default: 10000)",
    )

    args = parser.parse_args()

    # 验证输入
    if args.min_window <= 0 or args.max_window <= 0 or args.window_step <= 0:
        print("Error: Window sizes must be positive")
        return

    if args.min_window > args.max_window:
        print("Error: Minimum window size must be less than maximum window size")
        return

    # 检查FAI文件是否存在
    if not os.path.exists(args.fai):
        print(f"Error: FAI file not found: {args.fai}")
        return

    # 检查所有输入文件是否存在
    all_files = args.group1 + args.group2 + args.group3
    for file in all_files:
        if not os.path.exists(file):
            print(f"Error: Input file not found: {file}")
            return

    # 生成样本名称
    sample_names = []
    sample_files = []

    # 处理每个组的文件
    for i, (group_files, group_num) in enumerate(
        [(args.group1, 1), (args.group2, 2), (args.group3, 3)]
    ):
        for j, file in enumerate(group_files):
            # 从文件名生成样本名（去除路径和扩展名）
            base_name = os.path.splitext(os.path.basename(file))[0]
            sample_name = f"{base_name}-{group_num}"
            sample_names.append(sample_name)
            sample_files.append(file)

    print(f"Processing {len(sample_files)} samples:")
    for name, file in zip(sample_names, sample_files):
        print(f"  {name}: {file}")

    # 生成窗口大小列表
    window_sizes = list(range(args.min_window, args.max_window + 1, args.window_step))
    print(f"\nWill analyze window sizes: {[f'{w//1000}kb' for w in window_sizes]}")

    # 存储所有窗口大小的相关性结果
    all_correlations = {}

    # 处理每个窗口大小
    for window_size in window_sizes:
        print(f"\n--- Processing window size: {window_size//1000}kb ---")

        # 生成窗口化矩阵
        merged_data = process_samples_for_window_size(
            sample_files, sample_names, args.fai, window_size, args.output_prefix
        )

        # 按重复分组
        groups = group_samples_by_replicate(merged_data, sample_names)
        print(
            f"Sample groups: Group1={len(groups[1])}, Group2={len(groups[2])}, Group3={len(groups[3])}"
        )

        # 过滤全零行
        filtered_data = filter_zero_rows(merged_data, groups)

        # 计算相关性
        correlations = calculate_correlations(filtered_data, groups)
        all_correlations[window_size] = correlations

        # 打印当前窗口的相关性结果
        print(f"Correlations for {window_size//1000}kb window:")
        for pair, stats in correlations.items():
            print(
                f"  {pair}: Pearson={stats['pearson_r']:.3f} (p={stats['pearson_p']:.3e}), "
                f"Spearman={stats['spearman_r']:.3f} (p={stats['spearman_p']:.3e})"
            )

    # 创建可视化和保存结果
    print("\n--- Creating visualizations and saving results ---")

    # 创建相关性热图
    create_correlation_heatmaps(all_correlations, args.output_prefix)

    # 创建相关性趋势折线图
    create_correlation_line_plots(all_correlations, args.output_prefix)

    # 保存详细的相关性汇总
    save_correlation_summary(all_correlations, args.output_prefix)

    print(f"\nAnalysis completed! All results saved with prefix: {args.output_prefix}")
    print("Generated files:")
    print(f"  - Matrix files: {args.output_prefix}_window_*kb_matrix.csv")
    print(
        f"  - Correlation heatmaps: {args.output_prefix}_correlation_heatmaps.png/pdf"
    )
    print(f"  - Correlation trends: {args.output_prefix}_correlation_trends.png/pdf")
    print(f"  - Heatmap data: {args.output_prefix}_pearson/spearman_heatmap_data.csv")
    print(f"  - Trend data: {args.output_prefix}_correlation_trends_data.csv")
    print(f"  - Summary: {args.output_prefix}_correlation_summary.csv")


if __name__ == "__main__":
    main()
