#!/bin/bash

log_with_timestamp() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a $log_file
}

check_file_exists() {
    if [ ! -f "$1" ]; then
        log_with_timestamp "错误：文件不存在 - $1"
        exit 1
    fi
}

run_command() {
    log_with_timestamp "运行命令：$1"
    if ! eval $1 2>&1 | tee -a $log_file; then
        log_with_timestamp "命令失败：$1"
        exit 1
    fi
}

# Check if the correct number of arguments is provided
if [ "$#" -ne 10 ]; then
    echo "Usage: $0 -s sample_name -t threads -g genome_ref -i input_dir -o output_dir"
    exit 1
fi

while getopts ":s:t:g:i:o:" opt; do
  case $opt in
    s) sample="$OPTARG"
    ;;
    t) threads="$OPTARG"
    ;;
    g) genome_ref="$OPTARG"
    ;;
    i) input_dir="$OPTARG"
    ;;
    o) output_dir="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

fastq1="${input_dir}/${sample}.R1.fastq.gz"
fastq2="${input_dir}/${sample}.R2.fastq.gz"
align_dir="${output_dir}/${sample}/align_files"
log_file="process_eccDNA_${sample}.log"

# Create or clear the log file
> $log_file

# Check if input files exist
check_file_exists "$fastq1"
check_file_exists "$fastq2"
check_file_exists "$genome_ref"

# Change to the directory of ecc_finder
cd /data5/home/yxzhang2/tools/ecc_finder

# Run ecc_finder.py
cmd="python ecc_finder.py map-sr $genome_ref $fastq1 $fastq2 -t $threads -r $genome_ref -o ${output_dir}/${sample} -x $sample"
run_command "$cmd"

# Wait for ecc_finder.py to finish
wait

# Process the output using samtools and Circle-Map
cd "$align_dir" || { log_with_timestamp "无法切换到目录：$align_dir"; exit 1; }

# Sort SAM file by read name
cmd="samtools sort -@ $threads -n -o qname_eccDNA_${sample}.bam ${sample}.sam"
run_command "$cmd"

# Sort SAM file by coordinate
cmd="samtools sort -@ $threads -o sorted_eccDNA_${sample}.bam ${sample}.sam"
run_command "$cmd"

# Index the sorted BAM file
cmd="samtools index -@ $threads sorted_eccDNA_${sample}.bam"
run_command "$cmd"

# Extract reads mapping to eccDNA regions
cmd="Circle-Map ReadExtractor -i qname_eccDNA_${sample}.bam -o eccDNA_${sample}_candidates.bam"
run_command "$cmd"

# Wait for the above processes to finish
wait

# Sort the candidate BAM file
cmd="samtools sort -@ $threads eccDNA_${sample}_candidates.bam -o sort_eccDNA_${sample}_candidates.bam"
run_command "$cmd"

# Wait for the sort to finish
wait

# Remove unsorted candidates BAM
cmd="rm eccDNA_${sample}_candidates.bam"
run_command "$cmd"

# Index the sorted candidate BAM file
cmd="samtools index -@ $threads sort_eccDNA_${sample}_candidates.bam"
run_command "$cmd"

# Realign reads to improve eccDNA detection
cmd="Circle-Map Realign -t $threads -i sort_eccDNA_${sample}_candidates.bam -qbam qname_eccDNA_${sample}.bam -sbam sorted_eccDNA_${sample}.bam -fasta $genome_ref -o eccDNA_${sample}_CM.bed"
run_command "$cmd"

log_with_timestamp "脚本完成处理样本：$sample"

