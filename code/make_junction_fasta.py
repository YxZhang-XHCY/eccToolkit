# make_junction_fasta.py
from Bio import SeqIO
k = 300                                   # 末端/起始各取 300 bp
with open("all_cecc_renamed.fa") as fa_in, open("all_cecc_junc.fa","w") as fa_out:
    for rec in SeqIO.parse(fa_in, "fasta"):
        junc_seq = rec.seq[-k:] + rec.seq[:k]
        rec.id += "_JUNC"
        rec.description = "junction_of_" + rec.id
        rec.seq = junc_seq
        SeqIO.write(rec, fa_out, "fasta")
