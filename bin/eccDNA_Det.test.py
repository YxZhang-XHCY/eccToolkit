import subprocess
import os
import sys
import datetime
import argparse
from pathlib import Path

def log_with_timestamp(message, log_file):
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(log_file, 'a') as log:
        log.write(f"[{timestamp}] {message}\n")

def check_file_exists(file_path, log_file):
    if not Path(file_path).is_file():
        log_with_timestamp(f"错误：文件不存在 - {file_path}", log_file)
        sys.exit(1)

def run_command(cmd, log_file):
    log_with_timestamp(f"运行命令：{cmd}", log_file)
    process = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    with open(log_file, 'a') as log:
        log.write(process.stdout)
    if process.returncode != 0:
        log_with_timestamp(f"命令失败：{cmd}", log_file)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='EccDNA Detection and Processing Script.')
    parser.add_argument('-s', '--sample', required=True, help='Sample name')
    parser.add_argument('-t', '--threads', required=True, type=int, help='Number of threads')
    parser.add_argument('-g', '--genome_ref', required=True, help='Genome reference file')
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    args = parser.parse_args()

    fastq1 = Path(args.input_dir) / f"{args.sample}.R1.fastq.gz"
    fastq2 = Path(args.input_dir) / f"{args.sample}.R2.fastq.gz"
    align_dir = Path(args.output_dir) / args.sample / "align_files"
    log_file = f"process_eccDNA_{args.sample}.log"

    # Create or clear the log file
    open(log_file, 'w').close()

    # Check if input files exist
    check_file_exists(fastq1, log_file)
    check_file_exists(fastq2, log_file)
    check_file_exists(args.genome_ref, log_file)

    # Check and clone ecc_finder if needed
    eccToolkit_dir = Path(__file__).parent.parent
    ecc_finder_dir = eccToolkit_dir / "ecc_finder"
    if not ecc_finder_dir.is_dir():
        log_with_timestamp("ecc_finder 目录不存在，正在从 GitHub 克隆...", log_file)
        run_command(f"git clone https://github.com/njaupan/ecc_finder.git {ecc_finder_dir}", log_file)

    # Change to the ecc_finder directory and run ecc_finder.py
    os.chdir(ecc_finder_dir)
    ecc_finder_cmd = f"python ecc_finder.py map-sr {args.genome_ref} {fastq1} {fastq2} -t {args.threads} -r {args.genome_ref} -o {args.output_dir}/{args.sample} -x {args.sample}"
    run_command(ecc_finder_cmd, log_file)

    # Process the output using samtools and Circle-Map
    os.chdir(align_dir)
    samtools_sort_qname_cmd = f"samtools sort -@ {args.threads} -n -o qname_eccDNA_{args.sample}.bam {args.sample}.sam"
    run_command(samtools_sort_qname_cmd, log_file)

    samtools_sort_cmd = f"samtools sort -@ {args.threads} -o sorted_eccDNA_{args.sample}.bam {args.sample}.sam"
    run_command(samtools_sort_cmd, log_file)

    samtools_index_cmd = f"samtools index -@ {args.threads} sorted_eccDNA_{args.sample}.bam"
    run_command(samtools_index_cmd, log_file)

    circle_map_read_extractor_cmd = f"Circle-Map ReadExtractor -i qname_eccDNA_{args.sample}.bam -o eccDNA_{args.sample}_candidates.bam"
    run_command(circle_map_read_extractor_cmd, log_file)

    samtools_sort_candidates_cmd = f"samtools sort -@ {args.threads} eccDNA_{args.sample}_candidates.bam -o sort_eccDNA_{args.sample}_candidates.bam"
    run_command(samtools_sort_candidates_cmd, log_file)

    rm_unsorted_candidates_cmd = f"rm eccDNA_{args.sample}_candidates.bam"
    run_command(rm_unsorted_candidates_cmd, log_file)

    samtools_index_candidates_cmd = f"samtools index -@ {args.threads} sort_eccDNA_{args.sample}_candidates.bam"
    run_command(samtools_index_candidates_cmd, log_file)

    circle_map_realign_cmd = f"Circle-Map Realign -t {args.threads} -i sort_eccDNA_{args.sample}_candidates.bam -qbam qname_eccDNA_{args.sample}.bam -sbam sorted_eccDNA_{args.sample}.bam -fasta {args.genome_ref} -o eccDNA_{args.sample}_CM.bed"
    run_command(circle_map_realign_cmd, log_file)

    log_with_timestamp(f"脚本完成处理样本：{args.sample}", log_file)

if __name__ == "__main__":
    main()

