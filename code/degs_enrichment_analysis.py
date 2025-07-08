#!/usr/bin/env python3
"""
差异基因与富集/贫化基因关系分析脚本
分析不同log2FC阈值下的差异基因与富集/贫化基因的重叠关系
"""

import pandas as pd
import numpy as np
from collections import defaultdict
import argparse
import os


def read_degs(filepath):
    """读取差异基因文件"""
    degs = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=["gene", "ensembl", "tumor_expr", "normal_expr", "log2FC", "qvalue"],
    )
    return degs


def read_enrichment(filepath):
    """读取富集分析结果文件"""
    enrichment = pd.read_csv(filepath)
    return enrichment


def analyze_overlap(degs_df, enrichment_df, log2fc_thresholds=[1, 2, 4, 6]):
    """分析不同log2FC阈值下的基因重叠"""
    results = []

    # 获取唯一的样本组合
    samples = enrichment_df[["Sample", "Sample2"]].drop_duplicates()

    for _, sample_row in samples.iterrows():
        sample, sample2 = sample_row["Sample"], sample_row["Sample2"]

        # 筛选该样本的显著富集/贫化基因
        sample_data = enrichment_df[
            (enrichment_df["Sample"] == sample)
            & (enrichment_df["Sample2"] == sample2)
            & (enrichment_df["significant"] == True)
        ]

        enriched_genes = set(
            sample_data[sample_data["direction"] == "enrichment"]["gene"]
        )
        depleted_genes = set(
            sample_data[sample_data["direction"] == "depletion"]["gene"]
        )

        # 对每个log2FC阈值进行分析
        for threshold in log2fc_thresholds:
            # 筛选差异基因
            up_degs = set(degs_df[degs_df["log2FC"] >= threshold]["gene"])
            down_degs = set(degs_df[degs_df["log2FC"] <= -threshold]["gene"])

            # 计算交集
            up_enriched = up_degs & enriched_genes
            up_depleted = up_degs & depleted_genes
            down_enriched = down_degs & enriched_genes
            down_depleted = down_degs & depleted_genes

            result = {
                "Sample": f"{sample}_{sample2}",
                "log2FC_threshold": threshold,
                "n_up_DEGs": len(up_degs),
                "n_down_DEGs": len(down_degs),
                "n_enriched": len(enriched_genes),
                "n_depleted": len(depleted_genes),
                "n_up_enriched": len(up_enriched),
                "n_up_depleted": len(up_depleted),
                "n_down_enriched": len(down_enriched),
                "n_down_depleted": len(down_depleted),
                "up_enriched_genes": list(up_enriched),
                "up_depleted_genes": list(up_depleted),
                "down_enriched_genes": list(down_enriched),
                "down_depleted_genes": list(down_depleted),
            }
            results.append(result)

    return pd.DataFrame(results)


def print_summary(results_df):
    """打印汇总结果"""
    print("=" * 80)
    print("差异基因与富集/贫化基因关系分析结果")
    print("=" * 80)

    # 打印统计表格
    summary_cols = [
        "Sample",
        "log2FC_threshold",
        "n_up_DEGs",
        "n_down_DEGs",
        "n_enriched",
        "n_depleted",
        "n_up_enriched",
        "n_up_depleted",
        "n_down_enriched",
        "n_down_depleted",
    ]
    print("\n统计汇总:")
    print(results_df[summary_cols].to_string(index=False))

    # 打印详细基因列表（以log2FC >= 2为例）
    print("\n" + "=" * 80)
    print("详细基因列表 (log2FC >= 2):")
    print("=" * 80)

    threshold_2_results = results_df[results_df["log2FC_threshold"] == 2]

    for _, row in threshold_2_results.iterrows():
        print(f"\n样本: {row['Sample']}")
        print(
            f"上调且富集的基因 ({row['n_up_enriched']}个): {', '.join(row['up_enriched_genes']) if row['up_enriched_genes'] else '无'}"
        )
        print(
            f"上调且贫化的基因 ({row['n_up_depleted']}个): {', '.join(row['up_depleted_genes']) if row['up_depleted_genes'] else '无'}"
        )
        print(
            f"下调且富集的基因 ({row['n_down_enriched']}个): {', '.join(row['down_enriched_genes']) if row['down_enriched_genes'] else '无'}"
        )
        print(
            f"下调且贫化的基因 ({row['n_down_depleted']}个): {', '.join(row['down_depleted_genes']) if row['down_depleted_genes'] else '无'}"
        )


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="差异基因与富集/贫化基因关系分析脚本",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
    python degs_enrichment_analysis.py -d DEGs_GBM.GEPIA2.Log2FC_1.qValue_0.01.txt -e merged_gene_eccdna_enrichment_results.csv -o results.csv
    python degs_enrichment_analysis.py -d DEGs.txt -e enrichment.csv -o output.csv --thresholds 1 2 4 6
        """,
    )

    parser.add_argument(
        "-d",
        "--degs",
        required=True,
        help="差异基因文件路径 (TSV format with columns: gene, ensembl, tumor_expr, normal_expr, log2FC, qvalue)",
    )
    parser.add_argument(
        "-e", "--enrichment", required=True, help="富集分析结果文件路径 (CSV format)"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="输出文件路径 (CSV format)"
    )
    parser.add_argument(
        "--thresholds",
        nargs="+",
        type=float,
        default=[1, 2, 4, 6],
        help="Log2FC阈值列表 (default: 1 2 4 6)",
    )

    return parser.parse_args()


def main():
    # 解析命令行参数
    args = parse_arguments()

    # 检查输入文件是否存在
    if not os.path.exists(args.degs):
        print(f"错误: 差异基因文件不存在: {args.degs}")
        return

    if not os.path.exists(args.enrichment):
        print(f"错误: 富集分析文件不存在: {args.enrichment}")
        return

    # 读取数据
    print("读取数据文件...")
    degs_df = read_degs(args.degs)
    enrichment_df = read_enrichment(args.enrichment)

    print(f"差异基因数量: {len(degs_df)}")
    print(f"富集分析记录数: {len(enrichment_df)}")

    # 分析重叠
    print("\n开始分析...")
    results_df = analyze_overlap(degs_df, enrichment_df, args.thresholds)

    # 打印结果
    print_summary(results_df)

    # 保存结果
    results_df.to_csv(args.output, index=False)
    print(f"\n结果已保存到: {args.output}")


if __name__ == "__main__":
    main()
