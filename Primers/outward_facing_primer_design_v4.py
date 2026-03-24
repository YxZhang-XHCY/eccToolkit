import primer3
import pandas as pd
from pyfaidx import Fasta
import argparse
import subprocess
import os
from Bio import SeqIO
from typing import Dict, List, Tuple, Optional
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class PrimerDesignRetry:
    """引物设计重试策略类"""
    def __init__(self, max_attempts: int = 3):
        self.max_attempts = max_attempts
        self.current_attempt = 0
        
    def should_retry(self) -> bool:
        return self.current_attempt < self.max_attempts
    
    def next_attempt(self, global_args: Dict) -> Dict:
        self.current_attempt += 1
        if self.current_attempt >= self.max_attempts:
            return global_args
            
        # 逐步放宽参数
        relaxed_args = global_args.copy()
        relaxed_args['PRIMER_MAX_SELF_ANY'] += 1
        relaxed_args['PRIMER_MAX_SELF_END'] += 1
        relaxed_args['PRIMER_PAIR_MAX_COMPL_ANY'] += 1
        relaxed_args['PRIMER_PAIR_MAX_COMPL_END'] += 1
        
        if self.current_attempt == 2:
            relaxed_args['PRIMER_MIN_TM'] -= 0.5
            relaxed_args['PRIMER_MAX_TM'] += 0.5
            relaxed_args['PRIMER_MIN_GC'] -= 2.0
            relaxed_args['PRIMER_MAX_GC'] += 2.0
            
        logger.info(f"Attempt {self.current_attempt + 1}/{self.max_attempts} with relaxed parameters")
        return relaxed_args

def read_bed(bed_file: str) -> List[Tuple]:
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                raise ValueError(f"BED文件的行 '{line.strip()}' 少于4列")
            chrom, start, end, name = parts[:4]
            regions.append((chrom, int(start), int(end), name))
    return regions

def read_fasta(fasta_file: str) -> List[Tuple]:
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append((record.id, str(record.seq)))
    return sequences

def get_sequence(fasta: Fasta, chrom: str, start: int, end: int) -> str:
    return str(fasta[chrom][start:end])

def prepare_sequence_for_outward_primers(sequence: str) -> str:
    length = len(sequence)
    # 为了获得outward-facing引物，将序列前后拼接以确保可以设计跨越起点与终点的引物
    if length <= 600:
        return sequence + sequence
    else:
        return sequence[-300:] + sequence + sequence[:300]

def analyze_sequence(sequence: str) -> Dict:
    length = len(sequence)
    gc_content = (sequence.count('G') + sequence.count('C')) / length * 100
    return {
        'length': length,
        'gc_content': gc_content,
        'gc_content_formatted': f"{gc_content:.2f}%"
    }

def design_primers(sequence: str, seq_id: str, global_args: Dict) -> Optional[Dict]:
    seq_args = {
        'SEQUENCE_ID': seq_id,
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_INCLUDED_REGION': [0, len(sequence) - 1],
    }
    try:
        result = primer3.bindings.design_primers(seq_args, global_args)
        if result is None or result.get("PRIMER_PAIR_NUM_RETURNED", 0) == 0:
            seq_info = analyze_sequence(sequence)
            logger.warning(f"No primers found for {seq_id}. "
                         f"Sequence length: {seq_info['length']}, "
                         f"GC content: {seq_info['gc_content_formatted']}")
            return None
        return result
    except Exception as e:
        seq_info = analyze_sequence(sequence)
        logger.error(f"Error designing primers for {seq_id}: {str(e)}. "
                    f"Sequence length: {seq_info['length']}, "
                    f"GC content: {seq_info['gc_content_formatted']}")
        return None

def create_blast_db(reference_genome: str) -> str:
    db_name = os.path.splitext(reference_genome)[0]
    cmd = f"makeblastdb -in {reference_genome} -dbtype nucl -out {db_name}"
    subprocess.run(cmd, shell=True, check=True)
    return db_name

def run_blast(primer_seq: str, db_name: str, threads: int = 1) -> List[str]:
    cmd = f"echo '{primer_seq}' | blastn -task blastn-short -db {db_name} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident' -num_threads {threads}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip().split('\n')

def evaluate_blast_results(blast_results: List[str], primer_len: int) -> bool:
    perfect_matches = 0
    for result in blast_results:
        if result:
            fields = result.split('\t')
            identity = float(fields[2])
            alignment_length = int(fields[3])
            if identity == 100 and alignment_length == primer_len:
                perfect_matches += 1
    return perfect_matches == 1

def process_primer_results_with_blast(sequence: str, seq_id: str, name: str, 
                                      db_name: str, global_args: Dict, 
                                      threads: int, retry_strategy: PrimerDesignRetry) -> List[Dict]:
    all_primer_results = []
    current_args = global_args.copy()
    
    while retry_strategy.should_retry():
        primer3_result = design_primers(sequence, seq_id, current_args)
        if primer3_result is None:
            current_args = retry_strategy.next_attempt(current_args)
            continue
            
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
                'REVERSE_SPECIFIC': reverse_specific,
                'ATTEMPT': retry_strategy.current_attempt + 1
            }
            all_primer_results.append(primer_pair)
        
        if any(pr['FORWARD_SPECIFIC'] and pr['REVERSE_SPECIFIC'] for pr in all_primer_results):
            logger.info(f"Found specific primers for {seq_id} on attempt {retry_strategy.current_attempt + 1}")
            return all_primer_results
            
        current_args = retry_strategy.next_attempt(current_args)
    
    # 如果所有尝试都失败
    if not all_primer_results:
        logger.warning(f"Failed to design specific primers for {seq_id} after {retry_strategy.max_attempts} attempts")
        return [{
            'ID': seq_id, 'Name': name, 'FORWARD_PRIMER': '', 'REVERSE_PRIMER': '', 
            'PRODUCT_SIZE': 0, 'FORWARD_SPECIFIC': False, 'REVERSE_SPECIFIC': False,
            'ATTEMPT': retry_strategy.max_attempts
        }]
    
    return all_primer_results

def chunk_list(data: List, chunk_size: int) -> List[List]:
    """将列表数据分成若干个chunk"""
    for i in range(0, len(data), chunk_size):
        yield data[i:i+chunk_size]

def process_chunk(chunk, reference_fasta, db_name, global_args, threads, attempts, format_type):
    """处理一个chunk的序列/区域数据"""
    results = []
    if format_type == 'bed':
        for chrom, start, end, name in chunk:
            seq_id = f"{chrom}:{start}-{end}"
            sequence = get_sequence(reference_fasta, chrom, start, end)
            prepared_sequence = prepare_sequence_for_outward_primers(sequence)
            primers = process_primer_results_with_blast(
                prepared_sequence, seq_id, name, db_name, 
                global_args.copy(), threads, 
                PrimerDesignRetry(max_attempts=attempts)
            )
            results.extend(primers)
    else:  # fasta
        for seq_id, sequence in chunk:
            prepared_sequence = prepare_sequence_for_outward_primers(sequence)
            primers = process_primer_results_with_blast(
                prepared_sequence, seq_id, seq_id, db_name, 
                global_args.copy(), threads,
                PrimerDesignRetry(max_attempts=attempts)
            )
            results.extend(primers)
    return results

def main():
    parser = argparse.ArgumentParser(description='Design outward-facing primers for circular DNA with BLAST verification')
    parser.add_argument('-i', '--input', required=True, help='Input file (BED or FASTA format)')
    parser.add_argument('-f', '--format', choices=['bed', 'fasta'], required=True, help='Input file format (bed or fasta)')
    parser.add_argument('-r', '--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output file for primers')
    parser.add_argument('-d', '--blast_db', help='Path to existing BLAST database. If not provided, a new one will be created from the reference genome.')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use for BLAST searches (default: 1)')
    parser.add_argument('-p', '--product_size', type=int, default=200, help='Target product size (default: 200)')
    parser.add_argument('-a', '--attempts', type=int, default=3, help='Maximum number of design attempts (default: 3)')
    parser.add_argument('-l', '--log', help='Log file path')
    parser.add_argument('-c', '--chunk_size', type=int, default=10, help='Number of sequences/regions per chunk (default: 10)')
    parser.add_argument('-w', '--workers', type=int, default=1, help='Number of worker threads for parallel processing of chunks (default: 1)')
    args = parser.parse_args()

    # 设置文件日志
    if args.log:
        file_handler = logging.FileHandler(args.log)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)

    reference_fasta = Fasta(args.reference)

    # 创建或使用已有BLAST数据库
    if args.blast_db:
        db_name = args.blast_db
    else:
        logger.info("Creating BLAST database...")
        db_name = create_blast_db(args.reference)
        logger.info(f"BLAST database created: {db_name}")

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

    # 读取输入数据
    if args.format == 'bed':
        data = read_bed(args.input)
    else:
        data = read_fasta(args.input)

    # 分块处理
    chunks = list(chunk_list(data, args.chunk_size))
    all_primers = []

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = []
        for chk in chunks:
            fut = executor.submit(
                process_chunk, chk, reference_fasta, db_name, global_args, 
                args.threads, args.attempts, args.format
            )
            futures.append(fut)

        for fut in as_completed(futures):
            res = fut.result()
            all_primers.extend(res)

    df = pd.DataFrame(all_primers)
    df.to_csv(args.output, index=False, sep='\t')
    logger.info(f"Primers designed and written to {args.output}")

    # 统计与报告
    total_sequences = len(df['Name'].unique())
    specific_pairs = df[(df['FORWARD_SPECIFIC'] == True) & (df['REVERSE_SPECIFIC'] == True)]
    successful_sequences = len(specific_pairs['Name'].unique())
    
    logger.info(f"=== 设计结果统计 ===")
    logger.info(f"总序列数: {total_sequences}")
    logger.info(f"成功设计特异性引物的序列数: {successful_sequences}")
    if total_sequences > 0:
        logger.info(f"成功率: {(successful_sequences/total_sequences*100):.2f}%")
    else:
        logger.info("无有效的输入序列，无法计算成功率。")
    
    if 'ATTEMPT' in df.columns:
        attempt_stats = specific_pairs['ATTEMPT'].value_counts().sort_index()
        logger.info("\n每次尝试的成功数量：")
        for attempt, count in attempt_stats.items():
            logger.info(f"第 {attempt} 次尝试成功数量: {count}")

if __name__ == "__main__":
    main()
