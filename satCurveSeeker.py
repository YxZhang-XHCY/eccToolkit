#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
satCurveSeeker.py (Sequential)

功能：
  1. 按 10% 至 90%（每增10%）进行下采样（使用 seqkit sample）
  2. 对抽样后的每个文件执行 CircleSeeker，严格串行：等待前一个结束后再执行下一个
  3. 记录执行过程和耗时（logging 模块）
  4. 采用面向对象结构：SatCurveSeeker类 + main入口
"""

import argparse
import logging
import subprocess
import time
import os
import sys

class SatCurveSeeker:
    """
    一个封装下采样和调用 CircleSeeker 的类（严格串行版本）
    """
    def __init__(self, input_fasta, reference, prefix, outdir, threads):
        self.input_fasta = input_fasta
        self.reference = reference
        self.prefix = prefix
        self.outdir = outdir
        self.threads = threads

        # 需要抽样的比例：10%, 20%, ..., 90%
        self.percents = [p for p in range(10, 100, 10)]
        # 固定随机种子，便于结果可重复
        self.seed = 42

        # 设置日志记录器
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            handlers=[logging.StreamHandler(sys.stdout)]
        )
        self.logger = logging.getLogger("SatCurveSeeker")

        # 如果输出目录不存在，就创建
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
    
    def run(self):
        """
        主流程：对每个百分比做下采样，并调用CircleSeeker（严格串行）
        """
        self.logger.info("=== Saturation Curve Seeker (Sequential) Started ===")
        self.logger.info("Input FASTA: %s", self.input_fasta)
        self.logger.info("Reference: %s", self.reference)
        self.logger.info("Sample prefix: %s", self.prefix)
        self.logger.info("Output directory: %s", self.outdir)
        self.logger.info("Threads: %s", self.threads)

        for p in self.percents:
            fraction = p / 100.0
            # 下采样
            sampled_fasta = self.downsample(fraction)
            # 调用 CircleSeeker，等待完成
            self.run_circleseeker(sampled_fasta, p)
        
        self.logger.info("All downsampling & CircleSeeker jobs have completed.")
        self.logger.info("=== Saturation Curve Seeker Finished ===")

    def downsample(self, fraction):
        """
        使用 seqkit sample -p fraction -s seed 进行下采样，返回下采样文件路径
        """
        start_time = time.time()

        label = int(fraction * 100)
        sampled_fasta = os.path.join(self.outdir, f"{self.prefix}_{label}Per.fasta")

        cmd = f"seqkit sample -p {fraction} -s {self.seed} {self.input_fasta} > {sampled_fasta}"
        self.logger.info("Executing seqkit command: %s", cmd)
        subprocess.run(cmd, shell=True, check=True)

        cost = time.time() - start_time
        self.logger.info(
            "Sampling %d%% completed. Time used: %.2f seconds.",
            label, cost
        )

        return sampled_fasta

    def run_circleseeker(self, sampled_fasta, label):
        """
        调用 CircleSeeker, 并【严格串行】等待它完成
        """
        start_time = time.time()
        cs_cmd = (
            f"CircleSeeker "
            f"-i {sampled_fasta} "
            f"-p {self.prefix}_{label}Per "
            f"-r {self.reference} "
            f"-o {self.outdir} "
            f"-t {self.threads} "
            "--enable_X"
        )
        self.logger.info("Running CircleSeeker: %s", cs_cmd)

        # 使用 subprocess.run，使脚本在该步骤等待CircleSeeker执行完才继续
        subprocess.run(cs_cmd, shell=True, check=True)

        cost = time.time() - start_time
        self.logger.info(
            "CircleSeeker completed for %d%%. Time used: %.2f seconds.",
            label, cost
        )

def main():
    parser = argparse.ArgumentParser(description="Saturation curve approach with CircleSeeker (sequential execution).")
    parser.add_argument("-i", "--input_fasta", required=True,
                        help="FASTA格式的输入文件")
    parser.add_argument("-r", "--reference", required=True,
                        help="参考基因组文件路径")
    parser.add_argument("-p", "--prefix", required=True,
                        help="样本前缀")
    parser.add_argument("-o", "--outdir", required=True,
                        help="输出目录")
    parser.add_argument("-t", "--threads", required=True,
                        help="线程数量")
    args = parser.parse_args()

    # 创建并运行 SatCurveSeeker
    runner = SatCurveSeeker(
        input_fasta=args.input_fasta,
        reference=args.reference,
        prefix=args.prefix,
        outdir=args.outdir,
        threads=args.threads
    )
    runner.run()

if __name__ == "__main__":
    main()
