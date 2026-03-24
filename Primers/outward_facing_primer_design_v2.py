import primer3
import pandas as pd
from pyfaidx import Fasta
import argparse
import subprocess
import os
from Bio import SeqIO

def read_bed(bed_file):
    """读取BED文件，包括name列"""
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                raise ValueError(f"BED文件的行 '{line.strip()}' 少于4列")
            chrom, start, end, name = parts[:4]
            regions.append((chrom, int(start), int(end), name))
    return regions

def read_fasta(fasta_file):
    """读取FASTA文件"""
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append((record.id, str(record.seq)))
    return sequences

def get_sequence(fasta, chrom, start, end):
    """从参考基因组获取序列"""
    return str(fasta[chrom][start:end])

def prepare_sequence_for_outward_primers(sequence):
    """准备用于outward-facing引物设计的序列"""
    length = len(sequence)
    if length <= 600:
        return sequence + sequence
    else:
        return sequence[-300:] + sequence + sequence[:300]

def design_primers(sequence, seq_id, global_args):
    """设计引物"""
    seq_args = {
        'SEQUENCE_ID': seq_id,
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_INCLUDED_REGION': [0, len(sequence) - 1],
    }
    return primer3.bindings.design_primers(seq_args, global_args)

def create_blast_db(reference_genome):
    """为参考基因组创建BLAST数据库"""
    db_name = os.path.splitext(reference_genome)[0]
    cmd = f"makeblastdb -in {reference_genome} -dbtype nucl -out {db_name}"
    subprocess.run(cmd, shell=True, check=True)
    return db_name

def run_blast(primer_seq, db_name, threads=1):
    """对单个引物序列运行BLAST搜索"""
    cmd = f"echo '{primer_seq}' | blastn -task blastn-short -db {db_name} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident' -num_threads {threads}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip().split('\n')

def evaluate_blast_results(blast_results, primer_len):
    """评估BLAST结果以确定引物特异性"""
    perfect_matches = 0
    for result in blast_results:
        if result:
            fields = result.split('\t')
            identity = float(fields[2])
            alignment_length = int(fields[3])
            if identity == 100 and alignment_length == primer_len:
                perfect_matches += 1
    return perfect_matches == 1

def process_primer_results_with_blast(sequence, seq_id, name, db_name, global_args, threads, max_attempts=3):
    """处理引物设计结果并进行BLAST验证，如果没有找到特异性引物则重试"""
    for attempt in range(max_attempts):
        primer3_result = design_primers(sequence, seq_id, global_args)
        primer_results = []
        
        for i in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):
            forward_primer = primer3_result[f"PRIMER_LEFT_{i}_SEQUENCE"]
            reverse_primer = primer3_result[f"PRIMER_RIGHT_{i}_SEQUENCE"]
            
            forward_blast = run_blast(forward_primer, db_name, threads)
            reverse_blast = run_blast(reverse_primer, db_name, threads)
            
            forward_specific = evaluate_blast_results(forward_blast, len(forward_primer))
            reverse_specific = evaluate_blast_results(reverse_blast, len(reverse_primer))
            
            primer_pair = {
                'ID': f"{seq_id}_{i+1}",
                'Name': name,
                'FORWARD_PRIMER': forward_primer,
                'REVERSE_PRIMER': reverse_primer,
                'PRODUCT_SIZE': primer3_result[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                'FORWARD_SPECIFIC': forward_specific,
                'REVERSE_SPECIFIC': reverse_specific
            }
            primer_results.append(primer_pair)
        
        # 检查是否有特异性引物对
        if any(pr['FORWARD_SPECIFIC'] and pr['REVERSE_SPECIFIC'] for pr in primer_results):
            return primer_results
        
        # 如果没有找到特异性引物对，放宽参数再试
        global_args['PRIMER_MAX_SELF_ANY'] += 1
        global_args['PRIMER_MAX_SELF_END'] += 1
        global_args['PRIMER_PAIR_MAX_COMPL_ANY'] += 1
        global_args['PRIMER_PAIR_MAX_COMPL_END'] += 1
        print(f"No specific primer pair found for {seq_id}. Retrying with relaxed parameters. Attempt {attempt + 1}/{max_attempts}")
    
    print(f"Failed to design specific primers for {seq_id} after {max_attempts} attempts.")
    return primer_results  # 返回最后一次尝试的结果

def main():
    parser = argparse.ArgumentParser(description='Design outward-facing primers for circular DNA with BLAST verification')
    parser.add_argument('-i', '--input', required=True, help='Input file (BED or FASTA format)')
    parser.add_argument('-f', '--format', choices=['bed', 'fasta'], required=True, help='Input file format (bed or fasta)')
    parser.add_argument('-r', '--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output file for primers')
    parser.add_argument('-d', '--blast_db', help='Path to existing BLAST database. If not provided, a new one will be created from the reference genome.')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use for BLAST searches (default: 1)')
    parser.add_argument('-p', '--product_size', type=int, default=200, help='Target product size (default: 200)')
    args = parser.parse_args()

    # 读取参考基因组
    reference_fasta = Fasta(args.reference)

    # 处理BLAST数据库
    if args.blast_db:
        db_name = args.blast_db
    else:
        print("Creating BLAST database...")
        db_name = create_blast_db(args.reference)
        print(f"BLAST database created: {db_name}")

    # 设置引物设计参数
    product_size_min = max(100, args.product_size - 100)
    product_size_max = args.product_size + 100
    
    global_args = {
        'PRIMER_NUM_RETURN': 5,
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 58.0,
        'PRIMER_MAX_TM': 62.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_MAX_POLY_X': 4,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 8,
        'PRIMER_MAX_SELF_END': 3,
        'PRIMER_PAIR_MAX_COMPL_ANY': 8,
        'PRIMER_PAIR_MAX_COMPL_END': 3,
        'PRIMER_PRODUCT_SIZE_RANGE': [product_size_min, product_size_max],
        'PRIMER_GC_CLAMP': 2
    }

    all_primers = []

    if args.format == 'bed':
        regions = read_bed(args.input)
        for chrom, start, end, name in regions:
            seq_id = f"{chrom}:{start}-{end}"
            sequence = get_sequence(reference_fasta, chrom, start, end)
            prepared_sequence = prepare_sequence_for_outward_primers(sequence)
            primers = process_primer_results_with_blast(prepared_sequence, seq_id, name, db_name, global_args.copy(), args.threads)
            all_primers.extend(primers)
    elif args.format == 'fasta':
        sequences = read_fasta(args.input)
        for seq_id, sequence in sequences:
            prepared_sequence = prepare_sequence_for_outward_primers(sequence)
            primers = process_primer_results_with_blast(prepared_sequence, seq_id, seq_id, db_name, global_args.copy(), args.threads)
            all_primers.extend(primers)

    # 将结果写入文件
    df = pd.DataFrame(all_primers)
    df.to_csv(args.output, index=False, sep='\t')
    print(f"Primers designed and written to {args.output}")

    # 报告特异性引物对的数量
    specific_pairs = df[(df['FORWARD_SPECIFIC'] == True) & (df['REVERSE_SPECIFIC'] == True)]
    print(f"Total sequences processed: {len(df['Name'].unique())}")
    print(f"Sequences with specific primer pairs: {len(specific_pairs['Name'].unique())}")

if __name__ == "__main__":
    main()
