import pandas as pd
import subprocess
import argparse
from concurrent.futures import ProcessPoolExecutor
import os

def get_fasta_using_bedtools(genome_file, bed_line, extend_outer, extend_inner):
    """使用bedtools从基因组文件中提取序列"""
    chr, start, end, name = bed_line
    start_seq_bed = f"temp_{name}_start.bed"
    end_seq_bed = f"temp_{name}_end.bed"

    # 创建起始和终止位置的 BED 文件
    with open(start_seq_bed, "w") as f:
        f.write(f"{chr}\t{max(0, start - extend_outer)}\t{start + extend_inner}\n")
    with open(end_seq_bed, "w") as f:
        f.write(f"{chr}\t{end - extend_inner}\t{end + extend_outer}\n")

    # 生成FASTA文件
    start_seq_fasta = f"temp_{name}_start.fasta"
    end_seq_fasta = f"temp_{name}_end.fasta"

    subprocess.run(["bedtools", "getfasta", "-fi", genome_file, "-bed", start_seq_bed, "-fo", start_seq_fasta])
    subprocess.run(["bedtools", "getfasta", "-fi", genome_file, "-bed", end_seq_bed, "-fo", end_seq_fasta])


    # 删除临时BED文件
    os.remove(start_seq_bed)
    os.remove(end_seq_bed)

    return start_seq_fasta, end_seq_fasta


def reverse_fasta(infile_path, outfile_path):
    """反转fasta序列"""
    def rev_seq(seq):
        return seq[::-1]

    with open(infile_path, 'r') as infile, open(outfile_path, 'w') as outfile:
        for line in infile:
            if line.startswith(">"):
                outfile.write(line)
            else:
                outfile.write(rev_seq(line.strip()) + '\n')

def blast_to_csv(start_fasta, end_fasta, eccDNA_name, csv_file):
    """使用BLAST比对两个序列并保存输出到CSV"""
    output_file = f"blast_output_{eccDNA_name}.txt"

    # 执行BLAST比对
    blastn_cline = ["blastn", "-query", start_fasta, "-subject", end_fasta, "-out", output_file, "-outfmt", "6", "-task", "blastn-short", "-word_size", "4", "-evalue", "1"]
    subprocess.run(blastn_cline)

    # 检查输出文件是否为空
    if os.stat(output_file).st_size > 0:
        # 读取BLAST输出并转换为DataFrame
        df = pd.read_csv(output_file, sep="\t", header=None, names=["query_id", "subject_id", "identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"])
        df['eccDNA_name'] = eccDNA_name  # 添加eccDNA名称为一列

        # 追加到CSV文件
        with open(csv_file, 'a') as f:
            df.to_csv(f, header=f.tell()==0, index=False)
    else:
        print(f"There are no repeat sequences at the ends of {eccDNA_name}.")

    # 删除临时BLAST输出文件，无论是否有输出
    os.remove(output_file)
    

def process_eccDNA(row, genome_file, direct_csv_file, inverted_csv_file, extend_outer, extend_inner):
    """处理单个eccDNA的BLAST比对并保存输出"""
    start_fasta, end_fasta = get_fasta_using_bedtools(genome_file, (row['chr'], row['start'], row['end'], row['name']), extend_outer, extend_inner)
    
    # 直接比对
    blast_to_csv(start_fasta, end_fasta, row['name'], direct_csv_file)

    # 反转end序列后进行比对
    reversed_end_fasta = f"temp_{row['name']}_reversed_end.fasta"
    reverse_fasta(end_fasta, reversed_end_fasta)
    blast_to_csv(start_fasta, reversed_end_fasta, row['name'], inverted_csv_file)

    # 删除临时FASTA文件
    os.remove(start_fasta)
    os.remove(end_fasta)
    os.remove(reversed_end_fasta)


def main(args):
    # 清空或创建CSV文件，并添加标题行
    open(args.direct_csv, 'w').close()
    open(args.inverted_csv, 'w').close()

    # 读取BED文件
    df = pd.read_csv(args.bed_file, sep="\t", header=None, names=["chr", "start", "end", "name"])

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_eccDNA, row, args.genome_file, args.direct_csv, args.inverted_csv, args.extend_outer, args.extend_inner) for index, row in df.iterrows()]
        for future in futures:
            future.result()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BLAST eccDNA sequences")
    parser.add_argument("genome_file", help="Path to the genome fasta file")
    parser.add_argument("bed_file", help="Path to the bed file")
    parser.add_argument("direct_csv", help="Output CSV file for direct BLAST results")
    parser.add_argument("inverted_csv", help="Output CSV file for inverted BLAST results")
    parser.add_argument("--extend_outer", type=int, default=100, help="Number of base pairs to extend outward from the ends (default: 100)")
    parser.add_argument("--extend_inner", type=int, default=50, help="Number of base pairs to extend inward from the ends (default: 50)")

    args = parser.parse_args()
    main(args)


