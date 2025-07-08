#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
validate_eccdna_amplicons.py
----------------------------
A complete pipeline to validate candidate eccDNAs from targeted amplicon sequencing data.
This script implements a novel method for eccDNA validation:
1. It constructs a custom "head-to-tail" junction reference sequence for each candidate eccDNA.
2. It aligns sequencing reads to this custom reference.
3. It counts reads supporting the junction breakpoint, including those spanning the junction.
4. It makes a comprehensive determination of which eccDNAs were successfully amplified and detected.

中文简介：
一个完整的流程，用于通过靶向扩增子测序数据验证候选的 eccDNA。
该脚本实现了一种新的 eccDNA 验证方法：
1. 为每个候选 eccDNA 构建一个定制的“头尾拼接”连接点参考序列。
2. 将测序 reads 比对到该定制参考序列上。
3. 统计支持断点的 reads 数量（包括跨越拼接点的 reads）。
4. 综合判断哪些 eccDNA 被成功扩增和检测。

Dependencies / 依赖:
  - python 3.7+ (pip install pysam pandas)
  - bwa
  - samtools >= 1.10

Example / 示例 (Paired-end):
  python validate_eccdna_amplicons.py \\
      -i primer.tsv \\
      -1 eccPCR_R1.fq.gz -2 eccPCR_R2.fq.gz \\
      -r /path/to/hg38.fa \\
      -o HeLa_eccPCR_validation \\
      --inside 300 --flank 200 --threads 12
"""

import argparse
import os
import subprocess
import sys
import csv
import re
import pysam
import pandas as pd


# --- Command-line Interface / 命令行接口 ---
def get_args():
    """Parses command-line arguments."""
    p = argparse.ArgumentParser(
        description="A pipeline to validate candidate eccDNAs from amplicon sequencing.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Input/Output Arguments / 输入输出参数
    p.add_argument(
        "-i",
        "--tsv",
        required=True,
        help="Input TSV file with eccDNA candidates (Chr, Start, End, Name). / 包含候选eccDNA的TSV输入文件。",
    )
    p.add_argument(
        "-1",
        "--fastq1",
        required=True,
        help="Input FASTQ file, read 1 (gzipped OK). / 输入的FASTQ文件，Read 1。",
    )
    p.add_argument(
        "-2",
        "--fastq2",
        help="Input FASTQ file, read 2 (optional, for paired-end). / 输入的FASTQ文件，Read 2 (可选，用于双端测序)。",
    )
    p.add_argument(
        "-r",
        "--ref",
        required=True,
        help="Reference genome FASTA file (must be indexed). / 参考基因组FASTA文件(必须已建立索引)。",
    )
    p.add_argument(
        "-o",
        "--prefix",
        default="eccdna_validation",
        help="Prefix for all output files. / 所有输出文件的前缀。",
    )
    p.add_argument(
        "-t",
        "--threads",
        type=int,
        default=8,
        help="Number of threads for BWA and samtools. / BWA和samtools使用的线程数。",
    )

    # Method-specific parameters / 方法相关参数
    p.add_argument(
        "--inside",
        type=int,
        default=300,
        help="Length (bp) to keep from the inside of each eccDNA end. / 从eccDNA每个末端内侧保留的序列长度(bp)。",
    )
    p.add_argument(
        "--flank",
        type=int,
        default=200,
        help="Length (bp) to keep from the flanking region of each eccDNA end. / 从eccDNA每个末端外侧(旁翼)保留的序列长度(bp)。",
    )

    # Analysis thresholds / 分析阈值参数
    p.add_argument(
        "--window",
        type=int,
        default=20,
        help="Half-window size around the junction center to define a crossing read. / 定义跨越reads时，环绕连接点中心点的半窗口大小。",
    )
    p.add_argument(
        "--min-crossing",
        type=int,
        default=3,
        help="Minimum number of crossing reads to support a junction. / 支持一个连接点所需的最小跨越reads数。",
    )
    p.add_argument(
        "--min-unique",
        type=int,
        default=20,
        help="Minimum number of unique reads mapped to the amplicon. / 比对到扩增子上的最小唯一reads数。",
    )
    p.add_argument(
        "--min-split",
        type=int,
        default=3,
        help="Minimum number of split reads mapped to the amplicon. / 比对到扩增子上的最小剪切reads数。",
    )
    return p.parse_args()


# --- Shell Command Execution / Shell命令执行 ---
def run_command(cmd: str, log_prefix: str = None):
    """
    Executes a shell command and logs its output.
    中文：执行一个shell命令并记录其输出。
    """
    log_file = f"{log_prefix}.log" if log_prefix else os.devnull
    print(f">> Executing: {cmd}", file=sys.stderr)
    try:
        with open(log_file, "w") as f:
            subprocess.check_call(cmd, shell=True, stdout=f, stderr=f)
    except subprocess.CalledProcessError as e:
        print(f"!! Command failed with exit code {e.returncode}.", file=sys.stderr)
        print(f"!! See log for details: {log_file}", file=sys.stderr)
        raise  # Re-raise the exception to stop the pipeline / 重新抛出异常以终止流程


# --- Core Logic: Building Custom Reference / 核心逻辑：构建定制参考序列 ---
def build_amplicon_fasta(
    tsv_path: str, ref_fa_path: str, out_fa_path: str, inside_len: int, flank_len: int
):
    """
    Builds a custom FASTA file with head-to-tail concatenated sequences for each eccDNA candidate.
    The structure for each reference is: [TAIL_INSIDE][TAIL_OUTSIDE][HEAD_OUTSIDE][HEAD_INSIDE]
    This places the eccDNA junction point exactly in the center of the new reference sequence.

    中文：为每个候选eccDNA构建一个头尾拼接的定制FASTA参考文件。
    每个参考序列的结构是：[环内尾部][环外尾部][环外头部][环内头部]
    这种结构将eccDNA的连接点精确地置于新参考序列的中心。
    """
    ref = pysam.FastaFile(ref_fa_path)
    with open(tsv_path) as fh, open(out_fa_path, "w") as out:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            chrom = row["Chr"]
            start = int(row["Start"])
            end = int(row["End"])
            name = row.get("Name") or f"ecc_{row['No']}"
            chr_len = ref.get_reference_length(chrom)

            # Dynamically adjust 'inside' length for very short circles to prevent fetching negative coordinates.
            # 中文：对于非常短的环，动态调整'inside'长度，以防止坐标为负。
            effective_inside = min(inside_len, max(1, (end - start + 1) // 2))

            # Define coordinates for the four segments / 定义四个片段的坐标
            tail_in_s, tail_in_e = max(1, end - effective_inside + 1), end
            tail_out_s, tail_out_e = end + 1, min(chr_len, end + flank_len)
            head_out_s, head_out_e = max(1, start - flank_len), start - 1
            head_in_s, head_in_e = start, min(chr_len, start + effective_inside - 1)

            # Fetch and concatenate sequences / 获取并拼接序列
            seq = (
                ref.fetch(chrom, tail_in_s - 1, tail_in_e)  # Tail-Inside / 环内尾部
                + ref.fetch(
                    chrom, tail_out_s - 1, tail_out_e
                )  # Tail-Outside / 环外尾部
                + ref.fetch(
                    chrom, head_out_s - 1, head_out_e
                )  # Head-Outside / 环外头部
                + ref.fetch(chrom, head_in_s - 1, head_in_e)
            )  # Head-Inside / 环内头部

            # Write to FASTA file / 写入FASTA文件
            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i : i + 60] + "\n")
    ref.close()


# --- Core Logic: Alignment and Analysis / 核心逻辑：比对与分析 ---
def bwa_align(
    fa_path: str, fq1_path: str, fq2_path: str, threads: int, prefix: str
) -> str:
    """
    Performs BWA-MEM alignment and samtools sort/index.
    中文：执行BWA-MEM比对以及samtools排序/索引。
    """
    # Index the custom FASTA if not already done / 如果定制FASTA尚未索引，则建立索引
    if not os.path.exists(fa_path + ".bwt"):
        run_command(f"bwa index {fa_path}", f"{prefix}_bwa_index")

    bam_path = f"{prefix}.sorted.bam"
    bwa_cmd = f"bwa mem -t {threads} {fa_path} {fq1_path}" + (
        f" {fq2_path}" if fq2_path else ""
    )
    full_cmd = f"{bwa_cmd} | samtools sort -@{threads} -o {bam_path}"

    run_command(full_cmd, f"{prefix}_bwa_sort")
    run_command(f"samtools index -@ {threads} {bam_path}", f"{prefix}_samtools_index")
    return bam_path


def crosses_center(aln: pysam.AlignedSegment, center: int, window: int) -> bool:
    """Checks if a single alignment segment directly spans the junction center."""
    # 中文：检查单个比对片段是否直接跨越了连接点中心。
    return aln.reference_start < center - window and aln.reference_end > center + window


def sa_tag_crosses(aln: pysam.AlignedSegment, center: int, window: int) -> bool:
    """Checks the SA tag to see if a read is split across the junction center."""
    # 中文：检查SA标签，判断一个read是否被分割在连接点中心的两侧。
    if not aln.has_tag("SA"):
        return False
    for block in aln.get_tag("SA").split(";"):
        if not block:
            continue  # Handle trailing semicolon / 处理末尾的分号

        sa_pos = (
            int(block.split(",")[1]) - 1
        )  # SA tag position is 1-based / SA标签的位置是1-based
        # Check if primary and supplementary alignments are on opposite sides of the junction
        # 中文：检查主比对和辅助比对是否位于连接点的两侧
        if (aln.reference_end < center - window and sa_pos > center + window) or (
            aln.reference_start > center + window and sa_pos < center - window
        ):
            return True
    return False


def determine_status(unique_reads, split_reads, crossing_reads, thresholds):
    """
    Determines the detection status and confidence based on a two-tiered strategy.
    中文：基于两级判定策略，确定检测状态和置信度。
    """
    # Tier 1: Basic evidence from mapping quality / 第1级：来自比对质量的基本证据
    basic_pass = (
        unique_reads >= thresholds.min_unique and split_reads >= thresholds.min_split
    )

    # Tier 2: Strong evidence from junction-crossing reads / 第2级：来自跨连接点reads的强证据
    crossing_pass = crossing_reads >= thresholds.min_crossing

    if basic_pass and crossing_pass:
        return "Detected", "High"
    elif basic_pass or crossing_pass:
        return "Detected", "Medium"
    else:
        return "Not_detected", "Low"


def analyze_bam(bam_path: str, args: argparse.Namespace) -> pd.DataFrame:
    """
    Performs a comprehensive analysis of the BAM file to validate eccDNAs.
    中文：对BAM文件进行综合分析以验证eccDNA。
    """
    bam = pysam.AlignmentFile(bam_path)
    results = []

    for ref_name in bam.references:
        unique_reads, crossing_reads = set(), set()
        split_read_count = 0

        ref_len = bam.get_reference_length(ref_name)
        # The junction is engineered to be at the exact center of the custom reference
        # 中文：通过设计，连接点精确地位于定制参考序列的中心
        center = ref_len // 2

        for aln in bam.fetch(ref_name):
            if (
                aln.is_unmapped
                or aln.is_duplicate
                or aln.is_secondary
                or aln.is_supplementary
            ):
                continue

            # Count unique reads (primary alignments only) / 统计唯一reads（仅限主比对）
            unique_reads.add(aln.query_name)

            # Count split reads (from SA tag) / 统计剪切reads（来自SA标签）
            if aln.has_tag("SA"):
                split_read_count += 1

            # Check for junction-crossing evidence / 检查跨越连接点的证据
            if crosses_center(aln, center, args.window) or sa_tag_crosses(
                aln, center, args.window
            ):
                crossing_reads.add(aln.query_name)

        # Determine status and confidence / 判定状态和置信度
        status, confidence = determine_status(
            len(unique_reads), split_read_count, len(crossing_reads), args
        )

        results.append(
            {
                "ecc_id": ref_name,
                "unique_reads": len(unique_reads),
                "split_reads": split_read_count,
                "crossing_reads": len(crossing_reads),
                "ref_length": ref_len,
                "status": status,
                "confidence": confidence,
            }
        )

        # Log progress to stderr / 将进度实时输出到标准错误
        print(
            f"  {ref_name:<20}\tUnique: {len(unique_reads):<5}\tSplit: {split_read_count:<5}\tCrossing: {len(crossing_reads):<5}\tStatus: {status} ({confidence})",
            file=sys.stderr,
        )

    bam.close()
    return pd.DataFrame(results)


# --- Main Execution Block / 主执行模块 ---
def main():
    """Main function to orchestrate the pipeline."""
    args = get_args()

    # Define output file paths / 定义输出文件路径
    fa_file = f"{args.prefix}.amplicon.fa"
    bam_file = f"{args.prefix}.sorted.bam"
    summary_file = f"{args.prefix}.summary.tsv"
    detail_file = f"{args.prefix}.detail.tsv"

    # --- Start Pipeline ---
    print(f"## eccDNA Amplicon Validation Pipeline Started", file=sys.stderr)
    print(f"## Parameters:", file=sys.stderr)
    print(
        f"   - Custom Reference: inside={args.inside}bp, flank={args.flank}bp",
        file=sys.stderr,
    )
    print(
        f"   - Detection Thresholds: unique≥{args.min_unique}, split≥{args.min_split}, crossing≥{args.min_crossing}",
        file=sys.stderr,
    )

    try:
        print(f"\n## Step 1: Building custom amplicon reference...", file=sys.stderr)
        build_amplicon_fasta(args.tsv, args.ref, fa_file, args.inside, args.flank)

        print(f"\n## Step 2: Aligning reads with BWA-MEM...", file=sys.stderr)
        bam_file = bwa_align(
            fa_file, args.fastq1, args.fastq2, args.threads, args.prefix
        )

        print(f"\n## Step 3: Analyzing junction-spanning reads...", file=sys.stderr)
        df_results = analyze_bam(bam_file, args)

        # Save detailed and summary results / 保存详细和汇总结果
        df_results.to_csv(detail_file, sep="\t", index=False)
        summary = df_results[
            ["ecc_id", "status", "confidence", "unique_reads", "crossing_reads"]
        ].copy()
        summary.to_csv(summary_file, sep="\t", index=False)

        # Print final statistics / 打印最终统计数据
        n_total = len(df_results)
        if n_total > 0:
            n_detected = (df_results["status"] == "Detected").sum()
            n_high = (df_results["confidence"] == "High").sum()
            n_medium = (df_results["confidence"] == "Medium").sum()

            print(f"\n## Validation Complete!", file=sys.stderr)
            print(f"   Total candidates analyzed: {n_total}", file=sys.stderr)
            print(
                f"   Detected: {n_detected} ({n_detected/n_total*100:.1f}%)",
                file=sys.stderr,
            )
            print(f"     - High Confidence: {n_high}", file=sys.stderr)
            print(f"     - Medium Confidence: {n_medium}", file=sys.stderr)
        else:
            print(
                f"\n## Validation Complete! No candidates were processed.",
                file=sys.stderr,
            )

        print(f"\n## Output Files:", file=sys.stderr)
        print(f"   - Summary: {summary_file}", file=sys.stderr)
        print(f"   - Details: {detail_file}", file=sys.stderr)
        print(f"   - Alignment: {bam_file}", file=sys.stderr)
        print(f"   - Custom Reference: {fa_file}", file=sys.stderr)

    except Exception as e:
        print(
            f"\n!! An error occurred during the pipeline execution: {e}",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
