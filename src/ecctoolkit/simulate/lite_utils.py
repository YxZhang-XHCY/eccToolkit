#!/usr/bin/env python3

import os
import sys
import glob
import shutil
import pandas as pd
import numpy as np
import pysam as ps
import subprocess as sp

try:
    import fcntl
except ImportError:  # pragma: no cover - non-Unix platforms
    fcntl = None


class utilities(object):
    def __init__(self, 
                 reference: str,
                ):
        """
        utilities
        """
        self.reference = reference
        self._fasta = None
        self._genome_length_cache = None

    def _open_fasta(self):
        if not self.reference:
            raise ValueError('reference fasta path is required for this operation')
        if self._fasta is None:
            self._fasta = ps.Fastafile(self.reference)
        return self._fasta

    def resolve_chrom(self, chrom, genome=None):
        if genome is None:
            genome = self.genome_length()
        if chrom in genome.index:
            return chrom
        alt = chrom[3:] if str(chrom).startswith('chr') else 'chr{0}'.format(chrom)
        if alt in genome.index:
            return alt
        raise KeyError("chromosome '{0}' not found in reference".format(chrom))
        
    def genome_length(self):
        """
        calculate the length of each chromosome
        """
        if self._genome_length_cache is not None:
            return self._genome_length_cache
        fasta = self._open_fasta()
        output = pd.DataFrame([fasta.lengths], columns=fasta.references, index=['length']).T
        drop_names = list(output.filter(regex='_', axis=0).index) + ['chrM']
        output = output.drop(drop_names, errors='ignore')
        self._genome_length_cache = output
        return self._genome_length_cache
    
    def get_seq(self, chrom, start, end):
        genome = self.genome_length()
        chrom = self.resolve_chrom(chrom, genome)
        fasta = self._open_fasta()
        seq = fasta.fetch(chrom, start, end)
        seq = seq.upper()
        return seq
    
    def random_region(self, chrom, length, max_attempts=1000):
        genome = self.genome_length()
        chrom = self.resolve_chrom(chrom, genome)
        length = int(length)
        if length <= 0:
            raise ValueError('length must be > 0')
        max_start = int(genome.loc[chrom, 'length']) - length
        if max_start <= 0:
            raise ValueError('length exceeds reference contig length: {0} {1}'.format(chrom, length))
        for _ in range(int(max_attempts)):
            start = int(np.random.randint(max_start))
            end = start + length
            seq = self.get_seq(chrom, start, end)
            if seq.count('N') == 0:
                return [chrom, start, end, length, seq]
        raise RuntimeError(
            'failed to sample region without N after {0} attempts: {1} {2}'.format(max_attempts, chrom, length)
        )
    
    def transfer_files(self, file1, file2):
        with open(file1, 'rb') as src, open(file2, 'a+b') as dst:
            if fcntl:
                fcntl.flock(dst.fileno(), fcntl.LOCK_EX)
            try:
                shutil.copyfileobj(src, dst, length=1024 * 1024)
                # Ensure the appended content ends with a newline to prevent
                # records from different molecules being concatenated on
                # the same line, which would corrupt FASTQ format.
                dst.seek(0, 2)  # seek to end
                if dst.tell() > 0:
                    dst.seek(-1, 2)  # seek to last byte
                    if dst.read(1) != b'\n':
                        dst.write(b'\n')
            finally:
                if fcntl:
                    fcntl.flock(dst.fileno(), fcntl.LOCK_UN)
        return

    def concat_files(self, inputs, output, append=False):
        mode = 'ab' if append else 'wb'
        with open(output, mode) as dst:
            if fcntl:
                fcntl.flock(dst.fileno(), fcntl.LOCK_EX)
            try:
                for src_path in inputs:
                    with open(src_path, 'rb') as src:
                        shutil.copyfileobj(src, dst, length=1024 * 1024)
            finally:
                if fcntl:
                    fcntl.flock(dst.fileno(), fcntl.LOCK_UN)
        return

    def remove_glob(self, pattern):
        for path in glob.glob(pattern):
            try:
                if os.path.isdir(path):
                    shutil.rmtree(path)
                else:
                    os.remove(path)
            except FileNotFoundError:
                pass
        return
    
    def write_fasta(self, data, output):
        with open(output,'w') as f:
            for _, row in data.iterrows():
                f.write('>{0}\n'.format(row['id']))
                f.write(row['seq'])
                f.write('\n')
        return
    
