#!/usr/bin/env python3
import argparse
import logging
import os
import sys
import subprocess
from datetime import datetime
from pathlib import Path

class EccDNANanoporeRunner:
    def __init__(self):
        # 设置脚本默认路径
        self.script_path = Path("/data9/home/yxzhang/Tools/eccDNA_RCA_nanopore")
        self.setup_argument_parser()
        self.setup_logging()
        
    def setup_argument_parser(self):
        """设置命令行参数解析"""
        self.parser = argparse.ArgumentParser(description='eccDNA RCA Nanopore Pipeline Runner')
        
        # 必需参数
        self.parser.add_argument('--fastq', required=True, help='Input fastq file path')
        self.parser.add_argument('--reference', required=True, help='Reference genome file path')
        self.parser.add_argument('--output-dir', required=True, help='Output directory')
        self.parser.add_argument('--name', required=True, help='Analysis name')
        
        # 可选参数
        self.parser.add_argument('--threads', type=int, default=40,
                               help='Number of threads to use (default: 40)')
        self.parser.add_argument('--script-dir', default=str(self.script_path),
                               help=f'Directory containing eccDNA_RCA_nanopore.py (default: {self.script_path})')
        self.parser.add_argument('--verbose', action='store_true',
                               help='Enable verbose output')
    
    def setup_logging(self):
        """设置日志记录"""
        self.logger = logging.getLogger('eccDNA_nanopore')
        self.logger.setLevel(logging.DEBUG)
        
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(console_formatter)
        self.logger.addHandler(console_handler)
    
    def create_output_dirs(self, args):
        """创建输出目录结构"""
        output_dir = Path(args.output_dir)
        log_dir = output_dir / 'logs'
        result_dir = output_dir / 'results'
        
        for directory in [output_dir, log_dir, result_dir]:
            directory.mkdir(parents=True, exist_ok=True)
            self.logger.info(f"Created directory: {directory}")
            
        return log_dir, result_dir

    def run_minimap2(self, args, result_dir):
        """运行minimap2比对"""
        fastq_name = Path(args.fastq).stem
        paf_file = result_dir / f"{fastq_name}.paf"
        
        self.logger.info("Starting minimap2 alignment...")
        cmd = [
            'minimap2',
            '-t', str(args.threads),
            '-c',
            str(args.reference),
            str(args.fastq)
        ]
        
        try:
            with open(paf_file, 'w') as f:
                process = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True, check=True)
            self.logger.info("Minimap2 alignment completed successfully")
            return paf_file
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Minimap2 alignment failed: {e.stderr}")
            raise
    
    def run_eccDNA_analysis(self, args, paf_file, result_dir):
        """运行eccDNA分析"""
        # 切换到脚本目录
        os.chdir(args.script_dir)
        self.logger.info(f"Changed working directory to: {args.script_dir}")
        
        fastq_name = Path(args.fastq).stem
        info_file = result_dir / f"{fastq_name}_info.tsv"
        seq_file = result_dir / f"{fastq_name}.fa"
        var_file = result_dir / f"{fastq_name}_var.tsv"
        log_file = result_dir / 'analysis.log'
        
        self.logger.info("Starting eccDNA analysis...")
        cmd = [
            './eccDNA_RCA_nanopore.py',  # 使用相对路径
            '--fastq', str(args.fastq),
            '--paf', str(paf_file),
            '--info', str(info_file),
            '--seq', str(seq_file),
            '--var', str(var_file),
            '--reference', str(args.reference)
        ]
        
        if args.verbose:
            cmd.append('--verbose')
        
        try:
            with open(log_file, 'w') as f:
                process = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True, check=True)
            self.logger.info("eccDNA analysis completed successfully")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"eccDNA analysis failed: {e.stderr}")
            raise
    
    def run(self):
        """运行主流程"""
        args = self.parser.parse_args()
        
        # 设置文件日志
        log_dir = Path(args.output_dir) / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / f'eccDNA_nanopore_{args.name}_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
        
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(file_formatter)
        self.logger.addHandler(file_handler)
        
        # 检查输入文件和脚本
        for file_path in [args.fastq, args.reference]:
            if not os.path.exists(file_path):
                self.logger.error(f"Input file not found: {file_path}")
                sys.exit(1)
        
        script_file = Path(args.script_dir) / 'eccDNA_RCA_nanopore.py'
        if not script_file.exists():
            self.logger.error(f"eccDNA_RCA_nanopore.py not found in: {args.script_dir}")
            sys.exit(1)
        
        try:
            # 创建目录
            log_dir, result_dir = self.create_output_dirs(args)
            
            # 记录开始时间
            self.logger.info(f"Started analysis: {args.name}")
            start_time = datetime.now()
            
            # 运行minimap2
            paf_file = self.run_minimap2(args, result_dir)
            
            # 运行eccDNA分析
            self.run_eccDNA_analysis(args, paf_file, result_dir)
            
            # 记录结束时间和运行时间
            end_time = datetime.now()
            duration = end_time - start_time
            self.logger.info(f"Analysis completed. Total runtime: {duration}")
            
        except Exception as e:
            self.logger.error(f"Analysis failed: {str(e)}")
            sys.exit(1)

if __name__ == '__main__':
    runner = EccDNANanoporeRunner()
    runner.run()
