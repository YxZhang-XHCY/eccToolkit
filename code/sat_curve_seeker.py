#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
satCurveSeeker.py
-----------------
A workflow script to generate data for a saturation curve by subsampling a FASTA file
and running CircleSeeker on each subset sequentially.

中文简介：
一个工作流脚本，通过对FASTA文件进行下采样，并对每个子集依次运行CircleSeeker，
来生成用于绘制饱和度曲线的数据。
"""

import argparse
import logging
import subprocess
import time
import os
import sys


# --- 辅助函数 (Helper Function) ---
def run_command(cmd: str, log_file: str = None):
    """
    Executes a shell command, logs its output, and handles errors.
    中文：执行一个shell命令，记录其输出，并处理错误。
    """
    print(f">> Executing: {cmd}", file=sys.stderr)
    try:
        # If no log file is specified, pipe output to stderr
        # 中文：如果没有指定日志文件，则将输出管道输送到标准错误
        if log_file:
            with open(log_file, "w") as f:
                subprocess.run(
                    cmd, shell=True, check=True, stdout=f, stderr=subprocess.STDOUT
                )
        else:
            subprocess.run(cmd, shell=True, check=True, stderr=sys.stderr)

    except subprocess.CalledProcessError as e:
        print(f"!! Command failed with exit code {e.returncode}.", file=sys.stderr)
        if log_file:
            print(f"!! Check log for details: {log_file}", file=sys.stderr)
        raise


# --- 主流程类 (Main Pipeline Class) ---
class SatCurveSeeker:
    """
    A class to encapsulate the downsampling and CircleSeeker execution workflow.
    中文：一个封装了下采样和调用CircleSeeker工作流的类。
    """

    def __init__(self, args):
        """
        Initializes the workflow with command-line arguments.
        中文：使用命令行参数初始化工作流。
        """
        self.input_fasta = os.path.abspath(args.input_fasta)
        self.reference = os.path.abspath(args.reference)
        self.prefix = args.prefix
        self.outdir = os.path.abspath(args.outdir)
        self.threads = args.threads
        self.seed = args.seed

        # Generate sampling percentages based on user input / 根据用户输入生成采样百分比
        self.percents = list(range(args.start, args.end + 1, args.step))

        # Setup logging / 设置日志记录器
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            handlers=[
                logging.FileHandler(
                    os.path.join(self.outdir, f"{self.prefix}_saturation_run.log")
                ),
                logging.StreamHandler(sys.stdout),
            ],
        )
        self.logger = logging.getLogger("SatCurveSeeker")

        # Create output directory if it doesn't exist / 如果输出目录不存在，则创建它
        os.makedirs(self.outdir, exist_ok=True)

    def run(self):
        """
        Main workflow: iterates through percentages, downsamples, and calls CircleSeeker sequentially.
        中文：主流程：遍历所有百分比，进行下采样，然后依次调用CircleSeeker。
        """
        self.logger.info("=== Saturation Curve Seeker (Sequential) Started ===")
        self.logger.info(f"Input FASTA: {self.input_fasta}")
        self.logger.info(f"Reference: {self.reference}")
        self.logger.info(f"Output directory: {self.outdir}")
        self.logger.info(f"Sampling percentages: {self.percents}%")

        total_start_time = time.time()
        for p in self.percents:
            try:
                # Downsample the data / 对数据进行下采样
                sampled_fasta = self.downsample(p)
                # Run CircleSeeker on the subsampled data / 对下采样后的数据运行CircleSeeker
                self.run_circleseeker(sampled_fasta, p)
            except Exception as e:
                self.logger.error(
                    f"Pipeline failed at {p}% sampling level. Aborting. Error: {e}"
                )
                sys.exit(1)

        total_time = time.time() - total_start_time
        self.logger.info(
            f"All downsampling and CircleSeeker jobs have completed in {total_time:.2f} seconds."
        )
        self.logger.info("=== Saturation Curve Seeker Finished ===")

    def downsample(self, percent: int) -> str:
        """
        Downsamples the input FASTA using seqkit.
        中文：使用seqkit对输入的FASTA文件进行下采样。
        """
        fraction = percent / 100.0
        self.logger.info(f"--- Starting downsampling to {percent}% ---")
        start_time = time.time()

        sampled_fasta = os.path.join(self.outdir, f"{self.prefix}_{percent}p.fasta")
        cmd = f"seqkit sample -p {fraction} -s {self.seed} {self.input_fasta} -o {sampled_fasta}"

        run_command(
            cmd,
            log_file=os.path.join(self.outdir, f"{self.prefix}_{percent}p_seqkit.log"),
        )

        cost = time.time() - start_time
        self.logger.info(f"Sampling to {percent}% completed in {cost:.2f} seconds.")
        return sampled_fasta

    def run_circleseeker(self, sampled_fasta: str, percent: int):
        """
        Runs CircleSeeker on a given FASTA file and waits for it to complete.
        中文：对给定的FASTA文件运行CircleSeeker，并等待其完成。
        """
        self.logger.info(f"--- Running CircleSeeker for {percent}% sample ---")
        start_time = time.time()

        # Create a dedicated output subdirectory for each CircleSeeker run
        # 中文：为每次CircleSeeker运行创建一个专用的输出子目录
        cs_outdir = os.path.join(self.outdir, f"circleseeker_{percent}p_output")
        cs_prefix = f"{self.prefix}_{percent}p"

        cs_cmd = (
            f"CircleSeeker "
            f"-i {sampled_fasta} "
            f"-p {cs_prefix} "
            f"-r {self.reference} "
            f"-o {cs_outdir} "
            f"-t {self.threads} "
            "--enable_X"
        )

        run_command(
            cs_cmd,
            log_file=os.path.join(
                self.outdir, f"{self.prefix}_{percent}p_circleseeker.log"
            ),
        )

        cost = time.time() - start_time
        self.logger.info(
            f"CircleSeeker for {percent}% sample completed in {cost:.2f} seconds."
        )


def main():
    """Main function to parse arguments and kick off the workflow."""
    parser = argparse.ArgumentParser(
        description="Generate data for a saturation curve using CircleSeeker (sequential execution).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # --- I/O Arguments / 输入输出参数 ---
    parser.add_argument(
        "-i",
        "--input_fasta",
        required=True,
        help="Input FASTA file (e.g., from an assembler). / 输入的FASTA文件（例如，来自组装软件）。",
    )
    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="Path to the reference genome. / 参考基因组的路径。",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        required=True,
        help="Sample prefix for output files. / 用于输出文件的样本前缀。",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Main output directory to store all results. / 用于存储所有结果的主输出目录。",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=True,
        help="Number of threads for CircleSeeker. / CircleSeeker使用的线程数。",
    )

    # --- Sampling Control Arguments / 下采样控制参数 ---
    parser.add_argument(
        "--start",
        type=int,
        default=10,
        help="Starting percentage for subsampling. / 下采样的起始百分比。",
    )
    parser.add_argument(
        "--end",
        type=int,
        default=90,
        help="Ending percentage for subsampling. / 下采样的结束百分比。",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=10,
        help="Step percentage for subsampling. / 下采样的步长百分比。",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for seqkit to ensure reproducibility. / 用于seqkit的随机种子，以确保结果可重复。",
    )

    args = parser.parse_args()

    # Create and run the SatCurveSeeker instance / 创建并运行SatCurveSeeker实例
    runner = SatCurveSeeker(args)
    runner.run()


if __name__ == "__main__":
    main()
