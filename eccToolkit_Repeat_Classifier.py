#!/usr/bin/env python3
import os
import subprocess
import argparse
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class EccDNAProcessor:
    def __init__(self, input_file, genome_file, output_file, extend_outer, extend_inner, threads):
        self.input_file = input_file
        self.genome_file = genome_file
        self.output_file = output_file or str(Path(input_file).parent / 'combined_results.csv')
        self.extend_outer = extend_outer
        self.extend_inner = extend_inner
        self.threads = threads
        self.name_map = {}
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    def validate_file(self, f):
        p = Path(f)
        if not p.exists():
            logging.error(f"{f} not found")
            exit(1)
        return p

    def preprocess_input(self):
        self.validate_file(self.input_file)
        bed_data = []
        with open(self.input_file, 'r') as f:
            first_line = f.readline().strip()
            if not first_line.startswith('eName'):
                f.seek(0)
            for line in f:
                x = line.strip().split()
                if len(x) >= 11:
                    o, c, s, e = x[0], x[1], int(x[2]), int(x[3])
                    n = f"{c}_{s}_{e}"
                    bed_data.append(f"{c}\t{s}\t{e}\t{n}\n")
                    self.name_map[n] = o
        temp_bed = Path(self.input_file).with_suffix('.temp.bed')
        with open(temp_bed, 'w') as f:
            f.writelines(bed_data)
        self.bed_file = str(temp_bed)

    def get_fasta(self, row):
        c, s, e, n = row['chr'], row['start'], row['end'], row['name']
        b1 = f"temp_{n}_start.bed"
        b2 = f"temp_{n}_end.bed"
        with open(b1, 'w') as f:
            f.write(f"{c}\t{max(0, s - self.extend_outer)}\t{s + self.extend_inner}\n")
        with open(b2, 'w') as f:
            f.write(f"{c}\t{e - self.extend_inner}\t{e + self.extend_outer}\n")
        f1 = f"temp_{n}_start.fasta"
        f2 = f"temp_{n}_end.fasta"
        subprocess.run(["bedtools", "getfasta", "-fi", self.genome_file, "-bed", b1, "-fo", f1])
        subprocess.run(["bedtools", "getfasta", "-fi", self.genome_file, "-bed", b2, "-fo", f2])
        os.remove(b1)
        os.remove(b2)
        return f1, f2

    def reverse_fasta(self, i, o):
        with open(i, 'r') as fi, open(o, 'w') as fo:
            for l in fi:
                if l.startswith(">"):
                    fo.write(l)
                else:
                    fo.write(l.strip()[::-1] + "\n")

    def blast_to_csv(self, q, s, n, csvf):
        o = f"blast_output_{n}.txt"
        c = ["blastn", "-query", q, "-subject", s, "-out", o, "-outfmt", "6",
             "-task", "blastn-short", "-word_size", "4", "-evalue", "1"]
        subprocess.run(c)
        if os.path.exists(o) and os.stat(o).st_size > 0:
            df = pd.read_csv(o, sep="\t", header=None,
                             names=["query_id","subject_id","identity","alignment_length",
                                    "mismatches","gap_opens","q_start","q_end",
                                    "s_start","s_end","evalue","bit_score"])
            df['eccDNA_name'] = n
            with open(csvf, 'a') as f:
                df.to_csv(f, header=f.tell()==0, index=False)
        os.remove(o)

    def run_blast_analysis(self):
        open('intermediate_direct_results.csv', 'w').close()
        open('intermediate_inverted_results.csv', 'w').close()
        df = pd.read_csv(self.bed_file, sep='\t', header=None, names=["chr","start","end","name"])
        def task(row):
            f1, f2 = self.get_fasta(row)
            self.blast_to_csv(f1, f2, row['name'], 'intermediate_direct_results.csv')
            rev = f"temp_{row['name']}_rev.fasta"
            self.reverse_fasta(f2, rev)
            self.blast_to_csv(f1, rev, row['name'], 'intermediate_inverted_results.csv')
            for x in [f1, f2, rev]:
                os.remove(x)
        with ProcessPoolExecutor(max_workers=self.threads) as exe:
            futs = [exe.submit(task, r) for _, r in df.iterrows()]
            for f in futs: f.result()

    def process_repeat(self, inf, outf, inv=False):
        df = pd.read_csv(inf)
        for c in ["identity","q_end","q_start","s_end","s_start"]:
            df[c] = pd.to_numeric(df[c], errors='coerce')
        x = df[(df['identity']>90)&((df['q_end']-df['q_start'])>0)&((df['s_end']-df['s_start'])>0)].copy()
        x['q_real_start'] = x['q_start'] + pd.to_numeric(x['query_id'].str.extract(r':(\d+)-')[0], errors='coerce')
        x['q_real_end']   = x['q_end']   + pd.to_numeric(x['query_id'].str.extract(r':(\d+)-')[0], errors='coerce')
        if inv:
            x['s_real_start'] = pd.to_numeric(x['subject_id'].str.extract(r'-(\d+)')[0], errors='coerce') - x['s_end'] + 1
            x['s_real_end']   = pd.to_numeric(x['subject_id'].str.extract(r'-(\d+)')[0], errors='coerce') - x['s_start'] + 1
        else:
            x['s_real_start'] = x['s_start'] + pd.to_numeric(x['subject_id'].str.extract(r':(\d+)-')[0], errors='coerce')
            x['s_real_end']   = x['s_end']   + pd.to_numeric(x['subject_id'].str.extract(r':(\d+)-')[0], errors='coerce')
        cc = x['eccDNA_name'].str.split('_', expand=True)
        x['ecc_chr']   = cc[0]
        x['ecc_start'] = pd.to_numeric(cc[1], errors='coerce')
        x['ecc_end']   = pd.to_numeric(cc[2], errors='coerce')
        cts = (x['q_real_start']-x['ecc_start']).abs()<2
        cte = (x['s_real_end']-x['ecc_end']).abs()<2
        y = x[cts|cte].copy()
        y.to_csv(outf, index=False)

    def combine_results(self):
        df1 = pd.read_csv('processed_direct_results.csv')
        df2 = pd.read_csv('processed_inverted_results.csv')
        df1['Repeat_Class'] = 'Direct-Repeat'
        df2['Repeat_Class'] = 'Inverted-Repeat'
        df1['q_real_chr'] = df1['s_real_chr'] = df1['ecc_chr']
        df2['q_real_chr'] = df2['s_real_chr'] = df2['ecc_chr']
        d = pd.concat([df1, df2], ignore_index=True)
        d['Original_eName'] = d['eccDNA_name'].map(self.name_map)
        k = ['eccDNA_name','Original_eName','identity','alignment_length',
             'q_real_chr','q_real_start','q_real_end',
             's_real_chr','s_real_start','s_real_end','Repeat_Class']
        d.to_csv(self.output_file, index=False, columns=k)

    def cleanup(self):
        for f in [self.bed_file, 'intermediate_direct_results.csv',
                  'processed_direct_results.csv','intermediate_inverted_results.csv',
                  'processed_inverted_results.csv']:
            if os.path.exists(f): os.remove(f)

    def process(self):
        self.preprocess_input()
        self.run_blast_analysis()
        self.process_repeat('intermediate_direct_results.csv','processed_direct_results.csv')
        self.process_repeat('intermediate_inverted_results.csv','processed_inverted_results.csv',inv=True)
        self.combine_results()
        self.cleanup()

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-i','--input',required=True)
    p.add_argument('-g','--genome',required=True)
    p.add_argument('-o','--output',default=None)
    p.add_argument('-t','--threads',type=int,default=None)
    p.add_argument('--extend_outer',type=int,default=100)
    p.add_argument('--extend_inner',type=int,default=50)
    a = p.parse_args()
    e = EccDNAProcessor(a.input,a.genome,a.output,a.extend_outer,a.extend_inner,a.threads)
    e.process()

if __name__ == "__main__":
    main()