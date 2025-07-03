#!/usr/bin/env python3
"""
detect_eccdna_complete.py
-------------------------
完整的eccDNA检测流程：
1. 根据候选 eccDNA 区段生成"头/尾 inside + outside"拼接参考
2. 将测序 reads 比对到该参考
3. 统计支持断点的 reads 数量（包括跨越拼接点的reads）
4. 综合判断哪些 eccDNA 被 PCR-NGS 成功扩增

依赖：
  python 3.7+   pip install pysam pandas
  bwa-mem, samtools >=1.10
  hg38.fa 需已有 .fai (samtools faidx)

示例（双端测序）：
  python detect_eccdna_complete.py \
      -i primer.tsv \
      -1 eccPCR_R1.fq.gz -2 eccPCR_R2.fq.gz \
      -r /path/hg38.fa \
      -o HeLa_eccPCR \
      --inside 300 --flank 200 --threads 12
"""

import argparse, os, subprocess, sys, csv, re
import pysam, pandas as pd

# ------------------------- CLI ------------------------- #
def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-i", "--tsv", required=True, help="primer.tsv (tab-delimited)")
    p.add_argument("-1", "--fastq1", required=True, help="FASTQ read1 (gz OK)")
    p.add_argument("-2", "--fastq2", help="FASTQ read2 (omit for single-end)")
    p.add_argument("-r", "--ref", required=True, help="hg38 fasta (indexed)")
    p.add_argument("-o", "--prefix", default="eccdna", help="output file prefix")
    p.add_argument("--inside", type=int, default=300, help="bp kept inside ecc (each side)")
    p.add_argument("--flank",  type=int, default=200, help="bp kept outside ecc (each side)")
    p.add_argument("-t", "--threads", type=int, default=8, help="bwa threads")
    # 跨越分析参数
    p.add_argument("--window", type=int, default=20,
                   help="half-window around center regarded as junction")
    p.add_argument("--min-crossing", type=int, default=3,
                   help="min crossing reads to support junction")
    p.add_argument("--min-unique", type=int, default=20,
                   help="min unique reads for initial detection")
    p.add_argument("--min-split", type=int, default=3,
                   help="min split reads for initial detection")
    return p.parse_args()

# ------------------ build custom fasta ---------------- #
def build_amplicon_fasta(tsv, ref_fa, out_fa, Lin, Lflank):
    """构建头尾拼接的参考序列"""
    ref = pysam.FastaFile(ref_fa)
    with open(tsv) as fh, open(out_fa, "w") as out:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            chrom  = row['Chr']
            start  = int(row['Start'])
            end    = int(row['End'])
            name   = row.get('Name') or f"ecc_{row['No']}"
            chrlen = ref.get_reference_length(chrom)

            # 动态裁剪：环太短时缩小 inside 长度
            inside = min(Lin, max(1, (end - start + 1)//2))

            # tail inside/outside 坐标
            tail_in_s, tail_in_e   = max(1, end - inside + 1), end
            tail_out_s, tail_out_e = end + 1, min(chrlen, end + Lflank)

            # head outside/inside 坐标
            head_out_s, head_out_e = max(1, start - Lflank), start - 1
            head_in_s,  head_in_e  = start, min(chrlen, start + inside - 1)

            seq = (ref.fetch(chrom, tail_in_s-1, tail_in_e)   +  # 环内尾
                   ref.fetch(chrom, tail_out_s-1, tail_out_e) +  # 环外尾
                   ref.fetch(chrom, head_out_s-1, head_out_e) +  # 环外头
                   ref.fetch(chrom, head_in_s-1, head_in_e))     # 环内头

            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")
    ref.close()

# ---------------- shell helpers ----------------------- #
def sh(cmd):
    print(">>", cmd, file=sys.stderr)
    subprocess.check_call(cmd, shell=True)

def bwa_align(fa, fq1, fq2, threads, prefix):
    """BWA比对和排序"""
    if not os.path.exists(fa + ".bwt"):
        sh(f"bwa index {fa}")
    bam = prefix + ".sorted.bam"
    bwa_cmd = f"bwa mem -t {threads} {fa} {fq1}" + (f" {fq2}" if fq2 else "")
    sh(f"{bwa_cmd} | samtools sort -@{threads} -o {bam}")
    sh(f"samtools index {bam}")
    return bam

# -------------- crossing reads analysis --------------- #
def crosses_center(aln, center, W):
    """判断是否单条比对段直接跨越 center±W"""
    return aln.reference_start < center - W and aln.reference_end > center + W

def sa_cross(aln, center, W):
    """检查 SA 标签是否把 read 的另一段落在另一侧"""
    if not aln.has_tag("SA"):
        return False
    for block in aln.get_tag("SA").split(';'):
        if not block:  # 末尾 empty
            continue
        fields = block.split(',')
        sa_pos = int(fields[1]) - 1  # 0-based
        if (aln.reference_end < center - W and sa_pos > center + W) or \
           (aln.reference_start > center + W and sa_pos < center - W):
            return True
    return False

# -------------- comprehensive analysis ---------------- #
def analyze_eccdna(bamfile, args):
    """综合分析：基础统计 + 跨越reads分析"""
    bam = pysam.AlignmentFile(bamfile)
    results = []
    
    for ref in bam.references:
        # 基础统计
        unique_reads = set()
        split_reads = 0
        
        # 跨越分析
        ref_len = bam.get_reference_length(ref)
        center = ref_len // 2
        crossing_reads = set()
        
        for aln in bam.fetch(ref):
            if aln.is_unmapped or aln.is_duplicate:
                continue
                
            # 基础统计
            unique_reads.add(aln.query_name)
            if aln.has_tag("SA"):
                split_reads += 1
            else:
                cig = aln.cigarstring or ""
                if re.match(r'^\d+S', cig) or cig.endswith('S'):
                    split_reads += 1
            
            # 跨越分析
            if crosses_center(aln, center, args.window) or sa_cross(aln, center, args.window):
                crossing_reads.add(aln.query_name)
        
        # 综合判定
        n_unique = len(unique_reads)
        n_crossing = len(crossing_reads)
        
        # 两级判定策略
        basic_pass = n_unique >= args.min_unique and split_reads >= args.min_split
        crossing_pass = n_crossing >= args.min_crossing
        
        if basic_pass and crossing_pass:
            status = "Detected"
            confidence = "High"
        elif basic_pass or crossing_pass:
            status = "Detected"
            confidence = "Medium"
        else:
            status = "Not_detected"
            confidence = "Low"
        
        results.append({
            'ecc_id': ref,
            'unique_reads': n_unique,
            'split_reads': split_reads,
            'crossing_reads': n_crossing,
            'ref_length': ref_len,
            'status': status,
            'confidence': confidence
        })
        
        # 实时输出进度
        print(f"{ref}\t{n_unique} unique\t{split_reads} split\t"
              f"{n_crossing} crossing\t{status} ({confidence})", file=sys.stderr)
    
    bam.close()
    return pd.DataFrame(results)

# ------------------------------ main ------------------ #
def main():
    args = get_args()
    
    # 输出文件名
    fa_file = args.prefix + ".amplicon.fa"
    bam_file = args.prefix + ".sorted.bam"
    summary_file = args.prefix + ".summary.tsv"
    detail_file = args.prefix + ".detail.tsv"
    
    print(f"## eccDNA 完整检测流程", file=sys.stderr)
    print(f"## 参数设置：", file=sys.stderr)
    print(f"   - 参考序列：inside={args.inside}bp, flank={args.flank}bp", file=sys.stderr)
    print(f"   - 检测阈值：unique≥{args.min_unique}, split≥{args.min_split}, crossing≥{args.min_crossing}", file=sys.stderr)
    print(f"   - 拼接窗口：±{args.window}bp", file=sys.stderr)
    
    print(f"\n## 1) 生成拼接参考序列...", file=sys.stderr)
    build_amplicon_fasta(args.tsv, args.ref, fa_file, args.inside, args.flank)
    
    print(f"\n## 2) BWA-MEM 比对...", file=sys.stderr)
    bam_file = bwa_align(fa_file, args.fastq1, args.fastq2, args.threads, args.prefix)
    
    print(f"\n## 3) 综合分析...", file=sys.stderr)
    df = analyze_eccdna(bam_file, args)
    
    # 保存详细结果
    df.to_csv(detail_file, sep="\t", index=False)
    
    # 生成简化汇总
    summary = df[['ecc_id', 'status', 'confidence', 'unique_reads', 'crossing_reads']].copy()
    summary.to_csv(summary_file, sep="\t", index=False)
    
    # 统计输出
    n_total = len(df)
    n_detected = len(df[df['status'] == 'Detected'])
    n_high = len(df[df['confidence'] == 'High'])
    n_medium = len(df[df['confidence'] == 'Medium'])
    
    print(f"\n## 检测完成！", file=sys.stderr)
    print(f"   总计: {n_total} 个候选 eccDNA", file=sys.stderr)
    print(f"   检出: {n_detected} 个 ({n_detected/n_total*100:.1f}%)", file=sys.stderr)
    print(f"   高置信度: {n_high} 个", file=sys.stderr)
    print(f"   中置信度: {n_medium} 个", file=sys.stderr)
    print(f"\n输出文件：", file=sys.stderr)
    print(f"   - {summary_file} (简化汇总)", file=sys.stderr)
    print(f"   - {detail_file} (详细结果)", file=sys.stderr)
    print(f"   - {bam_file} (比对文件)", file=sys.stderr)
    print(f"   - {fa_file} (参考序列)", file=sys.stderr)

if __name__ == "__main__":
    main()
