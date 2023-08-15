#!/usr/bin/env python3

import os
import sys
import shutil

import pandas as pd


def reverse_fasta(infile_path, outfile_path):
    """Reverse complement FASTA sequences."""
    
    def rev_comp(seq):
        return seq.upper()[::-1]

    with open(infile_path) as infile:
        sequences = []
        name, seq = '', ''
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if name:
                    sequences.append((name, rev_comp(seq)))
                name, seq = line, ''
            else:
                seq += line
        sequences.append((name, rev_comp(seq)))

    with open(outfile_path, 'w') as outfile:
        for name, seq in sequences:
            outfile.write(f'{name}\n{seq}\n')

            
def replace_in_file(file_path, old, new):
    """Replace text in a file."""
    
    with open(file_path) as f:
        content = f.read()
        
    content = content.replace(old, new)
    
    with open(file_path, 'w') as f:
        f.write(content)

        
def split_fasta(in_fasta, out_dir):
    """Split FASTA file into individual sequences."""
    
    os.makedirs(out_dir, exist_ok=True)
    
    with open(in_fasta) as f:
        seq_id, seq = '', ''
        for line in f:
            if line.startswith('>'):
                if seq_id:
                    with open(os.path.join(out_dir, f'{seq_id}.fa'), 'w') as out:
                        out.write(f'>{seq_id}\n{seq}')
                seq_id, seq = line[1:-1], ''
            else:
                seq += line
        with open(os.path.join(out_dir, f'{seq_id}.fa'), 'w') as out:
            out.write(f'>{seq_id}\n{seq}')

            
def generate_table(num_seqs, out_path):
    """Generate table for BLAST."""
    
    with open(out_path, 'w') as f:
        for i in range(1, num_seqs+1):
            cols = ['UpEccDNA{}'.format(i), 
                    'DoEccDNA{}'.format(i),
                    'ReEccDNA{}'.format(i),
                    'Dir{}'.format(i), 
                    'Inv{}'.format(i),
                    str(i)]
            f.write('\t'.join(cols) + '\n')

            
def run_blast(table_path, out_dir):
    """Run BLAST on sequence pairs."""
    
    os.makedirs(out_dir, exist_ok=True)
    
    with open(table_path) as f:
        for line in f:
            cols = line.strip().split()
            seq1, seq2, seq3 = [c + '.fa' for c in cols[:3]]
            
            # Make BLAST DB
            os.system(f'makeblastdb -in {seq1} -dbtype nucl -out {seq1}.db')
            
            # Run BLAST
            os.system(f'blastn -query {seq2} -db {seq1}.db -out {cols[3]}.blast -outfmt 6') 
            os.system(f'blastn -query {seq3} -db {seq1}.db -out {cols[4]}.blast -outfmt 6')
            
            # Move BLAST outputs
            for blast_file in [cols[3] + '.blast', cols[4] + '.blast']:
                shutil.move(blast_file, out_dir)
                
            # Delete intermediate files
            for fa_file in [seq1, seq2, seq3]:
                os.remove(fa_file)
            os.remove(f'{seq1}.db.*')
            
    
def merge_blast_outputs(out_dir):
    """Merge all BLAST outputs."""
    
    os.system(f'cat {out_dir}/Dir*.blast > {out_dir}/Dir.All.blast')
    os.system(f'cat {out_dir}/Inv*.blast > {out_dir}/Inv.All.blast')

    
if __name__ == '__main__':

    seq_loc = sys.argv[1]
    ref_genome = sys.argv[2]
    genome_fa = sys.argv[3]
    out_dir = sys.argv[4]

    # Get flanking regions
    os.system(f'bedtools flank -i {seq_loc} -g {ref_genome} -l 100 -r 0 > {out_dir}/100_UpUp.bed')
    os.system(f'bedtools flank -i {out_dir}/100_UpUp.bed -g {ref_genome} -l 0 -r 50 > {out_dir}/UpUp_50.bed')
    os.system(f'join -j 4 -t $'\t' -o 1.1 1.2 2.3 1.4 2.4 {out_dir}/100_UpUp.bed {out_dir}/UpUp_50.bed > {out_dir}/100_UpUp_50.bed')

    os.system(f'bedtools flank -i {seq_loc} -g {ref_genome} -l 0 -r 100 > {out_dir}/100_Down.bed')
    os.system(f'bedtools flank -i {out_dir}/100_Down.bed -g {ref_genome} -l 50 -r 0 > {out_dir}/Down_50.bed')
    os.system(f'join -j 4 -t $'\t' -o 1.1 1.2 2.3 1.4 2.4 {out_dir}/Down_50.bed {out_dir}/100_Down.bed > {out_dir}/100_Down_50.bed')

    # Get sequences
    os.system(f'bedtools getfasta -nameOnly -fi {genome_fa} -bed {out_dir}/100_UpUp_50.bed > {out_dir}/100_UpUp_50.fa')
    os.system(f'bedtools getfasta -nameOnly -fi {genome_fa} -bed {out_dir}/100_Down_50.bed > {out_dir}/100_Down_50.fa')

    # Clean up
    for f in ['100_UpUp.bed', 'UpUp_50.bed', '100_UpUp_50.bed', 
              '100_Down.bed', 'Down_50.bed', '100_Down_50.bed']:
        os.remove(os.path.join(out_dir, f))

    # Reverse complement
    reverse_fasta(f'{out_dir}/100_Down_50.fa', f'{out_dir}/100_Reve_50.fa')

    # Replace EccDNA
    replace_in_file(f'{out_dir}/100_UpUp_50.fa', 'Ecc', 'UpEcc')
    replace_in_file(f'{out_dir}/100_Down_50.fa', 'Ecc', 'DoEcc')
    replace_in_file(f'{out_dir}/100_Reve_50.fa', 'Ecc', 'ReEcc')

    # Move to BLAST dir
    os.makedirs(f'{out_dir}/blast', exist_ok=True)
    for fa_file in ['100_UpUp_50.fa', '100_Down_50.fa', '100_Reve_50.fa']:
        shutil.move(fa_file, f'{out_dir}/blast/')

    # Split FASTA files 
    os.chdir(f'{out_dir}/blast')
    split_fasta('100_UpUp_50.fa', '.')
    split_fasta('100_Down_50.fa', '.')
    split_fasta('100_Reve_50.fa', '.')

    # Generate table
    file_num = len(os.listdir('.'))
    seq_num = file_num // 3
    generate_table(seq_num, 'blast.txt')

    # BLAST
    run_blast('blast.txt', './Result')

    # Merge outputs
    merge_blast_outputs('./Result')
    
print('Pipeline completed successfully!')