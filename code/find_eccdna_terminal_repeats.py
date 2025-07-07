#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
find_eccdna_terminal_repeats.py
-------------------------------
A sophisticated pipeline to identify terminal direct or inverted repeats at the
breakpoints of candidate eccDNAs. This is achieved by extracting flanking sequences
around each breakpoint, performing self-alignment using BLAST, and analyzing the
results to find repeats anchored at the junction.

The method is parallelized to efficiently process large lists of eccDNA candidates.

中文简介：
一个用于鉴定候选eccDNA断点处是否存在末端重复序列（正向或反向）的复杂流程。
本方法通过提取每个断点两侧的旁翼序列，使用BLAST进行自我比对，并分析比对结果
来寻找锚定在连接点上的重复序列。

该流程经过并行化处理，可高效分析大量的候选eccDNA。
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import uuid
import pandas as pd
import pysam
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from abc import ABC, abstractmethod

# --- 全局设置 (Global Setup) ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(module)s - %(message)s')

# --- 辅助模块 (Helper Modules) ---
def check_dependencies():
    """
    Checks if required command-line tools (bedtools, blastn) are in the system's PATH.
    中文：检查所需的命令行工具 (bedtools, blastn) 是否存在于系统路径中。
    """
    tools = ["bedtools", "blastn"]
    all_found = True
    for tool in tools:
        if shutil.which(tool) is None:
            logging.error(f"FATAL: Required tool '{tool}' not found in PATH. Please install it.")
            all_found = False
    if not all_found:
        sys.exit(1)
    logging.info("All required tools (bedtools, blastn) found in PATH.")

def run_command(cmd, log_file=None):
    """
    Executes a shell command robustly, with logging.
    中文：以健壮的方式执行shell命令，并记录日志。
    """
    logging.debug(f"Executing command: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, check=False, capture_output=True, text=True, encoding='utf-8')
        if result.returncode != 0:
            # Check for common "no hits" messages which are not true errors
            # 中文：检查常见的“未找到匹配”信息，这些并非真正的错误
            no_hits_messages = ["No hits found", "Query sequence not found"]
            if any(msg.lower() in (result.stderr + result.stdout).lower() for msg in no_hits_messages):
                logging.info(f"Command finished with no hits: {' '.join(cmd)}")
                return True, result
            
            logging.error(f"Command failed: {' '.join(cmd)}")
            logging.error(f"STDOUT: {result.stdout.strip()}")
            logging.error(f"STDERR: {result.stderr.strip()}")
            if log_file:
                with open(log_file, "w") as f:
                    f.write("COMMAND:\n" + ' '.join(cmd) + "\n\n")
                    f.write("STDOUT:\n" + result.stdout + "\n")
                    f.write("STDERR:\n" + result.stderr + "\n")
            return False, result
        return True, result
    except FileNotFoundError:
        logging.error(f"Command not found: {cmd[0]}. Please ensure it is installed and in your PATH.")
        raise
    except Exception as e:
        logging.error(f"An unexpected error occurred while running command {' '.join(cmd)}: {e}", exc_info=True)
        raise

# --- 核心类定义 (Core Class Definitions) ---

class BaseProcessor(ABC):
    """
    Abstract base class for processors.
    中文：处理器的抽象基类。
    """
    @abstractmethod
    def process(self):
        pass

class BlastAnalyzer(BaseProcessor):
    """
    Handles FASTA extraction and BLAST analysis for a single eccDNA.
    中文：为单个eccDNA处理FASTA序列提取和BLAST分析。
    """
    def __init__(self, genome_file, fai_file, extend_outer, extend_inner, blast_params, parent_temp_dir):
        self.genome_file = Path(genome_file)
        self.fai_file = Path(fai_file) if fai_file else None
        self.extend_outer = extend_outer
        self.extend_inner = extend_inner
        self.blast_params = blast_params
        self.chrom_lengths = self._load_fai()
        self.parent_temp_dir = Path(parent_temp_dir)

    def _load_fai(self):
        """Loads chromosome lengths from a .fai index file."""
        chrom_lengths = {}
        if not self.fai_file or not self.fai_file.exists():
            logging.warning(f"FAI file not provided or does not exist: {self.fai_file}. Chromosome boundary checks will be limited.")
            return chrom_lengths
        try:
            with open(self.fai_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    chrom_lengths[parts[0]] = int(parts[1])
            logging.info(f"Loaded {len(chrom_lengths)} chromosome lengths from {self.fai_file}")
        except Exception as e:
            logging.error(f"Error loading FAI file {self.fai_file}: {e}")
        return chrom_lengths
        
    def _clamp_coords(self, chrom, start, end):
        """Ensures coordinates are within chromosome boundaries."""
        chrom_len = self.chrom_lengths.get(chrom) or self.chrom_lengths.get(f"chr{chrom}") or self.chrom_lengths.get(chrom.replace("chr", ""))
        if chrom_len:
            clamped_start = max(0, start)
            clamped_end = min(chrom_len, end)
            is_valid = clamped_start < clamped_end
            return clamped_start, clamped_end, is_valid
        return max(0, start), end, start < end

    def _get_fasta(self, chrom, start, end, ecc_name, suffix, temp_dir):
        """Extracts a FASTA sequence for a given genomic region using bedtools."""
        clamped_start, clamped_end, is_valid = self._clamp_coords(chrom, start, end)
        if not is_valid:
            logging.warning(f"Skipping FASTA extraction for {ecc_name} ({suffix}): region {chrom}:{start}-{end} is invalid after clamping.")
            return None

        uid = uuid.uuid4().hex[:8]
        temp_bed_path = temp_dir / f"{ecc_name}_{uid}_{suffix}.bed"
        temp_fasta_path = temp_dir / f"{ecc_name}_{uid}_{suffix}.fasta"
        
        with open(temp_bed_path, "w") as f:
            f.write(f"{chrom}\t{clamped_start}\t{clamped_end}\n")

        cmd = ["bedtools", "getfasta", "-fi", str(self.genome_file), "-bed", str(temp_bed_path), "-fo", str(temp_fasta_path)]
        success, _ = run_command(cmd)
        
        if success and temp_fasta_path.exists() and temp_fasta_path.stat().st_size > 0:
            return temp_fasta_path
        return None

    def _reverse_complement_fasta(self, infile, outfile):
        """Creates a reverse-complemented version of a FASTA file."""
        try:
            # A simple way to reverse complement using seqkit if available, otherwise fallback
            # 中文：如果seqkit可用，则用其进行反向互补，否则使用备用方案
            if shutil.which("seqkit"):
                cmd = ["seqkit", "seq", "-r", "-p", str(infile), "-o", str(outfile)]
                success, _ = run_command(cmd)
                if success: return True
            
            # Fallback python implementation / Python备用实现
            complement = str.maketrans("ATCGN", "TAGCN")
            with open(infile, 'r') as fin, open(outfile, 'w') as fout:
                header = fin.readline()
                seq = "".join(line.strip() for line in fin)
                rev_comp_seq = seq.upper().translate(complement)[::-1]
                fout.write(header)
                for i in range(0, len(rev_comp_seq), 60):
                    fout.write(rev_comp_seq[i:i+60] + "\n")
            return True
        except Exception as e:
            logging.error(f"Failed to reverse complement {infile}: {e}")
            return False

    def _run_blast(self, query, subject, ecc_name, temp_dir):
        """Runs blastn and returns results as a DataFrame."""
        if not all([query, subject]): return pd.DataFrame()

        blast_out = temp_dir / f"blast_{ecc_name}_{uuid.uuid4().hex[:8]}.out"
        cmd = [
            "blastn", "-query", str(query), "-subject", str(subject), "-out", str(blast_out),
            "-outfmt", "6 std qseq qlen sseq slen", "-task", self.blast_params['task'],
            "-word_size", str(self.blast_params['word_size']), "-evalue", str(self.blast_params['evalue'])
        ]
        success, _ = run_command(cmd)

        if not success or not blast_out.exists() or blast_out.stat().st_size == 0:
            return pd.DataFrame()
        
        try:
            df = pd.read_csv(blast_out, sep="\t", header=None, names=[
                "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                "sstart", "send", "evalue", "bitscore", "qseq", "qlen", "sseq", "slen"
            ])
            # Apply filtering / 应用筛选
            df_filtered = df[(df['pident'] >= self.blast_params['min_identity']) & (df['length'] >= self.blast_params['min_length'])].copy()
            if not df_filtered.empty:
                df_filtered['eccDNA_name_internal'] = ecc_name
            return df_filtered
        except pd.errors.EmptyDataError:
            return pd.DataFrame()

    def process(self, ecc_data_tuple): # Renamed from process_single_eccdna for clarity
        """
        Processes a single eccDNA: extracts flanks, runs BLAST for direct/inverted repeats.
        中文：处理单个eccDNA：提取旁翼序列，为正向/反向重复运行BLAST。
        """
        ecc_chr, ecc_start, ecc_end, ecc_name = ecc_data_tuple
        process_temp_dir = Path(tempfile.mkdtemp(prefix=f"ecc_{ecc_name}_", dir=self.parent_temp_dir))
        
        try:
            # Get flank sequences / 获取旁翼序列
            start_flank_fasta = self._get_fasta(ecc_chr, ecc_start - self.extend_outer, ecc_start + self.extend_inner, ecc_name, "start_flank", process_temp_dir)
            end_flank_fasta = self._get_fasta(ecc_chr, ecc_end - self.extend_inner, ecc_end + self.extend_outer, ecc_name, "end_flank", process_temp_dir)

            if not start_flank_fasta or not end_flank_fasta:
                return pd.DataFrame(), pd.DataFrame()

            # --- Direct Repeats Analysis / 正向重复分析 ---
            direct_df = self._run_blast(start_flank_fasta, end_flank_fasta, ecc_name, process_temp_dir)

            # --- Inverted Repeats Analysis / 反向重复分析 ---
            reversed_end_fasta = process_temp_dir / f"{end_flank_fasta.stem}_revcomp.fasta"
            if self._reverse_complement_fasta(end_flank_fasta, reversed_end_fasta):
                inverted_df = self._run_blast(start_flank_fasta, reversed_end_fasta, ecc_name, process_temp_dir)
            else:
                inverted_df = pd.DataFrame()
            
            return direct_df, inverted_df
        finally:
            shutil.rmtree(process_temp_dir)

class RepeatResultsProcessor(BaseProcessor):
    """
    Processes a combined DataFrame of BLAST hits to filter for true terminal repeats.
    中文：处理一个合并后的BLAST结果DataFrame，以筛选出真正的末端重复。
    """
    def __init__(self, df_blast, output_path, is_inverted=False, tolerance=5):
        self.df = df_blast
        self.output_path = Path(output_path)
        self.is_inverted = is_inverted
        self.tolerance = tolerance

    def process(self):
        """
        Main processing method for filtering BLAST results.
        中文：用于筛选BLAST结果的主处理方法。
        """
        if self.df.empty:
            self.output_path.touch()
            return
        
        df = self.df.copy()

        # Remap alignment coordinates to genomic coordinates / 将比对坐标重映射至基因组坐标
        q_coords = df['qseqid'].str.extract(r'([^:]*):(\d+)-(\d+)', expand=True)
        s_coords = df['sseqid'].str.extract(r'([^:]*):(\d+)-(\d+)', expand=True)
        df['q_genome_start'] = pd.to_numeric(q_coords[1], errors='coerce') + df['qstart']
        df['q_genome_end'] = pd.to_numeric(q_coords[1], errors='coerce') + df['qend']
        
        if self.is_inverted:
            s_flank_start = pd.to_numeric(s_coords[1], errors='coerce')
            df['s_genome_start'] = s_flank_start + (df['slen'] - df['send'])
            df['s_genome_end'] = s_flank_start + (df['slen'] - df['sstart'])
        else:
            s_flank_start = pd.to_numeric(s_coords[1], errors='coerce')
            df['s_genome_start'] = s_flank_start + df['sstart']
            df['s_genome_end'] = s_flank_start + df['send']

        # Extract eccDNA breakpoint info from the name / 从名称中提取eccDNA断点信息
        ecc_info = df['eccDNA_name_internal'].str.split('_', expand=True)
        df['ecc_chr'] = ecc_info.get(0)
        df['ecc_start'] = pd.to_numeric(ecc_info.get(1), errors='coerce')
        df['ecc_end'] = pd.to_numeric(ecc_info.get(2), errors='coerce')

        # Drop rows where coordinate conversion failed / 删除坐标转换失败的行
        df.dropna(subset=['q_genome_start', 's_genome_start', 'ecc_start', 'ecc_end'], inplace=True)
        if df.empty:
            self.output_path.touch()
            return

        # --- Core Filtering Logic / 核心筛选逻辑 ---
        # A terminal repeat must have its ends anchored near the eccDNA breakpoints.
        # 中文：一个末端重复序列的末端必须锚定在eccDNA断点附近。
        BP1 = df['ecc_start'] + 1 # 1-based breakpoint 1 / 1-based断点1
        BP2 = df['ecc_end']       # 1-based breakpoint 2 / 1-based断点2

        # Check if alignment ends are within tolerance of breakpoints / 检查比对末端是否在断点的容忍范围内
        cond1 = (df['q_genome_end'] - BP1).abs() <= self.tolerance
        cond2 = (df['s_genome_start'] - BP2).abs() <= self.tolerance
        
        df_final = df[cond1 & cond2].copy()
        
        logging.info(f"Processed {len(self.df)} BLAST hits -> {len(df_final)} terminal repeats (is_inverted={self.is_inverted}).")
        df_final.to_csv(self.output_path, index=False)

class MainProcessor(BaseProcessor):
    """
    Orchestrates the entire end-to-end pipeline.
    中文：协调整个端到端流程。
    """
    def __init__(self, args):
        super().__init__(args.genome, args.input)
        self.args = args
        self.output_file = Path(args.output or Path(self.bed_file).parent / f"{Path(self.bed_file).stem}_terminal_repeats.csv")
        self.threads = args.threads if args.threads else max(1, os.cpu_count() // 2)
        self.name_map = {}
        self.run_temp_dir = Path(tempfile.mkdtemp(prefix="ecc_repeats_main_", dir=args.custom_temp_dir))
        logging.info(f"Main temporary directory: {self.run_temp_dir}")
        
    def _parse_input_line(self, line, line_num):
        """Parses and validates a single line from the input TSV."""
        fields = line.strip().split('\t')
        if len(fields) < 4: return None

        original_name, ecc_chr, start_str, end_str = fields[0], fields[1], fields[2], fields[3]
        try:
            # --- CRITICAL: Coordinate System Handling ---
            # Assumes input file uses 0-based, inclusive coordinates for start and end.
            # Converts to standard BED format (0-based start, 1-based exclusive end) for internal use.
            # 中文：--- 关键：坐标系统处理 ---
            # 假定输入文件使用0-based、包含型坐标 (start和end)。
            # 内部将其转换为标准的BED格式 (0-based start, 1-based exclusive end) 以便处理。
            start_0_incl = int(start_str)
            end_0_incl = int(end_str)
            
            # Internal representation / 内部表示
            ecc_start_0based = start_0_incl
            ecc_end_1based_exclusive = end_0_incl + 1

            if ecc_start_0based >= ecc_end_1based_exclusive: return None

            # Standardized name for internal tracking / 用于内部追踪的标准化名称
            standardized_name = f"{ecc_chr}_{ecc_start_0based}_{ecc_end_1based_exclusive}"
            self.name_map[standardized_name] = original_name
            return (ecc_chr, ecc_start_0based, ecc_end_1based_exclusive, standardized_name)
        except (ValueError, IndexError):
            logging.warning(f"Skipping malformed line {line_num}: {line.strip()}")
            return None

    def _preprocess_input(self):
        """Reads and parses the entire input file."""
        input_path = self.validate_input_file(self.bed_file)
        eccdnas = []
        with open(input_path, 'r') as f:
            header = f.readline().lower()
            if not ('chr' in header and 'start' in header and 'end' in header):
                f.seek(0) # No header detected, rewind / 未检测到表头，返回文件开头
            
            for i, line in enumerate(f, 1):
                if parsed := self._parse_input_line(line, i):
                    eccdnas.append(parsed)
        logging.info(f"Preprocessed {len(eccdnas)} eccDNA records from {input_path}")
        return eccdnas
    
    def process(self):
        """Main execution method."""
        all_eccdnas = self._preprocess_input()
        if not all_eccdnas: return
        
        blast_params = {
            'min_identity': self.args.min_blast_identity, 'min_length': self.args.min_blast_length,
            'task': self.args.blast_task, 'word_size': self.args.blast_word_size, 'evalue': self.args.blast_evalue
        }
        analyzer = BlastAnalyzer(self.args.genome, self.args.fai, self.args.extend_outer, self.args.extend_inner, blast_params, self.run_temp_dir)

        all_direct, all_inverted = [], []
        logging.info(f"Starting parallel BLAST analysis for {len(all_eccdnas)} eccDNAs using {self.threads} workers.")
        
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            for i, (direct_df, inverted_df) in enumerate(executor.map(analyzer.process, all_eccdnas)):
                if not direct_df.empty: all_direct.append(direct_df)
                if not inverted_df.empty: all_inverted.append(inverted_df)
                if (i + 1) % 100 == 0: logging.info(f"  ... {i + 1}/{len(all_eccdnas)} futures collected.")
        
        logging.info("BLAST analysis complete. Filtering results...")
        
        # Process direct repeats / 处理正向重复
        master_direct_df = pd.concat(all_direct, ignore_index=True) if all_direct else pd.DataFrame()
        direct_processor = RepeatResultsProcessor(master_direct_df, self.run_temp_dir / "direct.csv", is_inverted=False, tolerance=self.args.breakpoint_tolerance)
        direct_processor.process()

        # Process inverted repeats / 处理反向重复
        master_inverted_df = pd.concat(all_inverted, ignore_index=True) if all_inverted else pd.DataFrame()
        inverted_processor = RepeatResultsProcessor(master_inverted_df, self.run_temp_dir / "inverted.csv", is_inverted=True, tolerance=self.args.breakpoint_tolerance)
        inverted_processor.process()

        # Combine final results / 合并最终结果
        final_dfs = []
        try:
            df_direct = pd.read_csv(self.run_temp_dir / "direct.csv")
            if not df_direct.empty:
                df_direct['repeat_type'] = 'Direct'
                final_dfs.append(df_direct)
        except (FileNotFoundError, pd.errors.EmptyDataError): pass
        
        try:
            df_inverted = pd.read_csv(self.run_temp_dir / "inverted.csv")
            if not df_inverted.empty:
                df_inverted['repeat_type'] = 'Inverted'
                final_dfs.append(df_inverted)
        except (FileNotFoundError, pd.errors.EmptyDataError): pass

        if not final_dfs:
            logging.warning("No terminal repeats found after filtering.")
            self.output_file.touch()
        else:
            df_combined = pd.concat(final_dfs, ignore_index=True)
            df_combined['eccDNA_orig_name'] = df_combined['eccDNA_name_internal'].map(self.name_map)
            df_combined.rename(columns={'eccDNA_name_internal': 'eccDNA_std_name'}, inplace=True)
            # Select and reorder columns for final output / 为最终输出选择并重排序列
            final_cols = [
                'eccDNA_orig_name', 'eccDNA_std_name', 'repeat_type', 'pident', 'length', 'evalue',
                'q_genome_start', 'q_genome_end', 's_genome_start', 's_genome_end', 'ecc_start', 'ecc_end'
            ]
            df_final = df_combined[[c for c in final_cols if c in df_combined.columns]]
            df_final.to_csv(self.output_file, index=False)
            logging.info(f"✔ Analysis complete. Found {len(df_final)} terminal repeats. Results saved to {self.output_file}")
            
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        shutil.rmtree(self.run_temp_dir)
        logging.info(f"Main temporary directory removed: {self.run_temp_dir}")


def main():
    """Main entry point of the script."""
    check_dependencies()
    parser = argparse.ArgumentParser(
        description='Find terminal direct or inverted repeats at eccDNA breakpoints.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # --- Add arguments ---
    io_group = parser.add_argument_group('Input/Output')
    io_group.add_argument('-i', '--input', required=True, help="Input TSV file. Format: name, chr, start (0-based, inclusive), end (0-based, inclusive).")
    io_group.add_argument('-g', '--genome', required=True, help='Path to the reference genome FASTA file.')
    io_group.add_argument('-o', '--output', default=None, help='Output CSV file path. [Default: <input_stem>_terminal_repeats.csv]')
    io_group.add_argument('-t', '--threads', type=int, default=None, help='Number of parallel processes. [Default: half of CPU cores]')
    
    region_group = parser.add_argument_group('Region Definition')
    region_group.add_argument('--extend_outer', type=int, default=100, help="Base pairs to extend outward from eccDNA ends.")
    region_group.add_argument('--extend_inner', type=int, default=50, help="Base pairs to extend inward from eccDNA ends.")
    
    filter_group = parser.add_argument_group('Filtering Options')
    filter_group.add_argument('--breakpoint_tolerance', type=int, default=5, help="Tolerance (bp) for matching repeat ends to eccDNA breakpoints.")
    filter_group.add_argument('--min_blast_identity', type=float, default=90.0, help="Minimum BLAST percent identity.")
    filter_group.add_argument('--min_blast_length', type=int, default=10, help="Minimum BLAST alignment length.")

    blast_group = parser.add_argument_group('Advanced BLAST Options')
    blast_group.add_argument('--blast_task', default='blastn-short', help="Task for blastn (e.g., 'blastn', 'blastn-short').")
    blast_group.add_argument('--blast_word_size', type=int, default=7, help="Word size for blastn.")
    blast_group.add_argument('--blast_evalue', type=float, default=10.0, help="E-value threshold for blastn.")

    misc_group = parser.add_argument_group('Miscellaneous')
    misc_group.add_argument('--fai', default=None, help='Path to genome .fai index (recommended for accuracy).')
    misc_group.add_argument('--custom_temp_dir', default=None, help='Specify a directory for temporary files.')

    args = parser.parse_args()

    # --- Run Processor ---
    try:
        with MainProcessor(args) as processor:
            processor.process()
    except Exception:
        logging.error("The pipeline failed due to an unhandled exception.", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()