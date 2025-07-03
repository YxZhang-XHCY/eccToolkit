#!/usr/bin/env python3
import os
import pandas as pd
import subprocess
import argparse
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from abc import ABC, abstractmethod
import tempfile
import shutil
import uuid

# --- Global Setup ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(module)s - %(message)s')

class BaseProcessor(ABC):
    def __init__(self, genome_file=None, bed_file=None, extend_outer=100, extend_inner=50):
        self.genome_file = Path(genome_file) if genome_file else None
        self.bed_file = Path(bed_file) if bed_file else None # Original input bed file path
        self.extend_outer = extend_outer
        self.extend_inner = extend_inner

    @abstractmethod
    def process(self):
        pass

    def validate_input_file(self, file_path):
        path = Path(file_path)
        if not path.exists():
            logging.error(f'Input file {file_path} does not exist.')
            raise FileNotFoundError(f'Input file {file_path} does not exist.')
        if path.stat().st_size == 0:
            logging.warning(f'Input file {file_path} is empty.') # Warning, might be an error depending on context
        return path

class BlastAnalyzer(BaseProcessor):
    def __init__(self, genome_file, 
                 extend_outer=100, extend_inner=50, 
                 fai_file=None, parent_temp_dir=None,
                 min_blast_identity=90.0, min_blast_length=10):
        super().__init__(genome_file=genome_file, extend_outer=extend_outer, extend_inner=extend_inner)
        self.fai_file = Path(fai_file) if fai_file else None
        self.chrom_lengths = self._load_fai() if self.fai_file and self.fai_file.exists() else {}
        self.min_blast_identity = min_blast_identity
        self.min_blast_length = min_blast_length

        if parent_temp_dir:
            Path(parent_temp_dir).mkdir(parents=True, exist_ok=True)
            self.analyzer_main_temp_dir = Path(tempfile.mkdtemp(prefix="blast_analyzer_", dir=str(parent_temp_dir)))
        else:
            self.analyzer_main_temp_dir = Path(tempfile.mkdtemp(prefix="blast_analyzer_"))
        logging.info(f"BlastAnalyzer main temporary directory: {self.analyzer_main_temp_dir}")

    def _load_fai(self):
        chrom_lengths = {}
        if not self.fai_file or not self.fai_file.exists():
            logging.warning(f"FAI file not provided or does not exist: {self.fai_file}. Chromosome boundary clamping will be limited.")
            return chrom_lengths
        try:
            with open(self.fai_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        chrom_lengths[parts[0]] = int(parts[1])
            logging.info(f"Loaded {len(chrom_lengths)} chromosome lengths from {self.fai_file}")
        except Exception as e:
            logging.error(f"Error loading FAI file {self.fai_file}: {e}")
        return chrom_lengths

    def _clamp_coords(self, chrom, start, end):
        original_start, original_end = start, end
        clamped_start, clamped_end = start, end
        is_valid = True 

        chrom_variants = [chrom]
        if chrom.startswith("chr"):
            chrom_variants.append(chrom[3:])
        else:
            chrom_variants.append(f"chr{chrom}")
        
        matched_chrom_in_fai = None
        for variant in chrom_variants:
            if variant in self.chrom_lengths:
                matched_chrom_in_fai = variant
                break
        
        if matched_chrom_in_fai:
            chrom_len = self.chrom_lengths[matched_chrom_in_fai]
            clamped_start = max(0, start)
            clamped_end = min(chrom_len, end)
        else: 
            clamped_start = max(0, start)
            logging.debug(f"Chromosome '{chrom}' not found in FAI index. Limited clamping applied.")

        if clamped_start >= clamped_end: 
             is_valid = False
        
        if not is_valid:
             logging.debug(f"Region {chrom}:{original_start}-{original_end} resulted in "
                           f"zero or negative length after clamping: {clamped_start}-{clamped_end}.")
        return chrom, clamped_start, clamped_end, is_valid


    def _get_fasta_from_bed_region(self, chrom, region_start, region_end, region_name_suffix, ecc_name, process_temp_dir):
        uid = uuid.uuid4().hex[:8]
        safe_ecc_name = "".join(c if c.isalnum() or c in ('_', '-') else '_' for c in ecc_name)
        
        temp_bed_path = process_temp_dir / f"{safe_ecc_name}_{uid}_{region_name_suffix}.bed"
        temp_fasta_path = process_temp_dir / f"{safe_ecc_name}_{uid}_{region_name_suffix}.fasta"
        bed_content = f"{chrom}\t{region_start}\t{region_end}\n"
        
        try:
            with open(temp_bed_path, "w") as f:
                f.write(bed_content)

            cmd = ["bedtools", "getfasta", "-fi", str(self.genome_file), "-bed", str(temp_bed_path), "-fo", str(temp_fasta_path)]
            logging.debug(f"Running bedtools for {ecc_name} ({region_name_suffix}): {' '.join(cmd)}")
            result = subprocess.run(cmd, check=False, capture_output=True, text=True, encoding='utf-8')
            
            if result.returncode != 0:
                logging.error(f"bedtools failed for {ecc_name} ({region_name_suffix}): STDOUT='{result.stdout.strip()}' STDERR='{result.stderr.strip()}'")
                return None
            
            if not temp_fasta_path.exists() or temp_fasta_path.stat().st_size == 0:
                logging.warning(f"FASTA file is empty or missing for {ecc_name} ({region_name_suffix}) at {temp_fasta_path}. Bedtools stderr: '{result.stderr.strip()}'")
                return None
            return temp_fasta_path
        except Exception as e:
            logging.error(f"Exception in _get_fasta_from_bed_region for {ecc_name} ({region_name_suffix}): {e}", exc_info=True)
            return None
        finally:
            if temp_bed_path.exists():
                try: 
                    temp_bed_path.unlink(missing_ok=True)
                except TypeError: 
                     if os.path.exists(temp_bed_path): os.remove(temp_bed_path) # type: ignore
                except OSError: 
                    pass 
       
    def get_flank_fastas(self, ecc_data_tuple, process_temp_dir):
        ecc_chr, ecc_0_start, ecc_1_end, ecc_name = ecc_data_tuple[:4]

        start_flank_s_calc = ecc_0_start - self.extend_outer
        start_flank_e_calc = ecc_0_start + self.extend_inner
        end_flank_s_calc = ecc_1_end - self.extend_inner
        end_flank_e_calc = ecc_1_end + self.extend_outer

        sf_chr_clamped, sf_s_clamped, sf_e_clamped, sf_valid = self._clamp_coords(ecc_chr, start_flank_s_calc, start_flank_e_calc)
        ef_chr_clamped, ef_s_clamped, ef_e_clamped, ef_valid = self._clamp_coords(ecc_chr, end_flank_s_calc, end_flank_e_calc)

        if not (sf_valid and ef_valid): 
            logging.warning(f"Skipping {ecc_name} for FASTA extraction due to invalid/zero-length flank regions after clamping. "
                            f"StartFlank: {sf_chr_clamped}:{sf_s_clamped}-{sf_e_clamped} (valid: {sf_valid}), "
                            f"EndFlank: {ef_chr_clamped}:{ef_s_clamped}-{ef_e_clamped} (valid: {ef_valid})")
            return None, None

        start_fasta = self._get_fasta_from_bed_region(sf_chr_clamped, sf_s_clamped, sf_e_clamped, "start_flank", ecc_name, process_temp_dir)
        end_fasta = self._get_fasta_from_bed_region(ef_chr_clamped, ef_s_clamped, ef_e_clamped, "end_flank", ecc_name, process_temp_dir)
        
        if not start_fasta or not end_fasta: 
            return None, None
            
        return start_fasta, end_fasta

    def reverse_fasta(self, infile_path, outfile_path):
        try:
            rev_seq = lambda seq: seq[::-1]
            with open(infile_path, 'r') as infile, open(outfile_path, 'w') as outfile:
                header = ""
                sequence_parts = []
                for line in infile:
                    stripped_line = line.strip()
                    if stripped_line.startswith(">"):
                        if header: 
                            outfile.write(header + '\n') 
                            outfile.write(rev_seq("".join(sequence_parts)) + '\n')
                        header = stripped_line 
                        sequence_parts = []
                    else:
                        sequence_parts.append(stripped_line)
                if header: 
                    outfile.write(header + '\n')
                    outfile.write(rev_seq("".join(sequence_parts)) + '\n')
            return True
        except Exception as e:
            logging.error(f"Error reversing FASTA {infile_path} to {outfile_path}: {e}", exc_info=True)
            return False

    def run_blast_and_get_df(self, query_fasta, subject_fasta, ecc_name, process_temp_dir):
        if not query_fasta or not Path(query_fasta).exists() or Path(query_fasta).stat().st_size == 0:
            logging.warning(f"Skipping BLAST for {ecc_name}: query FASTA {query_fasta} is missing or empty.")
            return pd.DataFrame()
        if not subject_fasta or not Path(subject_fasta).exists() or Path(subject_fasta).stat().st_size == 0:
            logging.warning(f"Skipping BLAST for {ecc_name}: subject FASTA {subject_fasta} is missing or empty.")
            return pd.DataFrame()

        uid = uuid.uuid4().hex[:8]
        safe_ecc_name = "".join(c if c.isalnum() or c in ('_', '-') else '_' for c in ecc_name)
        blast_output_path = process_temp_dir / f"blast_output_{safe_ecc_name}_{uid}.txt"
        
        blastn_cline = [
            "blastn",
            "-query", str(query_fasta),
            "-subject", str(subject_fasta),
            "-out", str(blast_output_path),
            "-outfmt", "6 std qseq qlen sseq slen", 
            "-task", "blastn-short", 
            "-word_size", "7",       
            "-evalue", "10",         
        ]
        logging.debug(f"Running BLAST for {ecc_name} (query: {Path(query_fasta).name}, subject: {Path(subject_fasta).name}): {' '.join(blastn_cline)}")
        
        try:
            result = subprocess.run(blastn_cline, check=False, capture_output=True, text=True, encoding='utf-8')
            if result.returncode != 0:
                no_hits_messages = ["No hits found", "Query sequence not found", "subject sequence not found"] 
                if any(msg.lower() in (result.stderr + result.stdout).lower() for msg in no_hits_messages):
                    logging.info(f"No BLAST hits found for {ecc_name} (query: {Path(query_fasta).name}, subject: {Path(subject_fasta).name}).")
                    return pd.DataFrame()
                else: 
                    logging.error(f"BLASTn failed for {ecc_name}: STDOUT='{result.stdout.strip()}' STDERR='{result.stderr.strip()}'")
                    return pd.DataFrame() 
            
            if blast_output_path.exists() and blast_output_path.stat().st_size > 0:
                col_names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                             "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                             "qseq", "qlen", "sseq", "slen"]
                try:
                    df = pd.read_csv(blast_output_path, sep="\t", header=None, names=col_names)
                except pd.errors.EmptyDataError: 
                    logging.info(f"BLAST output file {blast_output_path} was empty for {ecc_name} despite non-zero size (possibly malformed).")
                    return pd.DataFrame()

                df = df[(df['pident'] >= self.min_blast_identity) & (df['length'] >= self.min_blast_length)]

                if df.empty:
                    logging.info(f"No hits passed filter (ident>={self.min_blast_identity}, len>={self.min_blast_length}) for {ecc_name} (query: {Path(query_fasta).name}, subject: {Path(subject_fasta).name}).")
                    return pd.DataFrame()

                df['eccDNA_name_internal'] = ecc_name 
                return df 
            else: 
                logging.info(f"No repeat sequences found (BLAST output empty or not created by successful run) for {ecc_name} (query: {Path(query_fasta).name}, subject: {Path(subject_fasta).name}).")
            return pd.DataFrame()
        except FileNotFoundError: 
            logging.error("blastn command not found. Please ensure BLAST+ is installed and in your PATH.")
            return pd.DataFrame() 
        except Exception as e:
            logging.error(f"Error during BLAST for {ecc_name}: {e}", exc_info=True)
            return pd.DataFrame()

    def process_single_eccdna(self, ecc_data_tuple_for_blast):
        ecc_name = ecc_data_tuple_for_blast[3]
        safe_ecc_name = "".join(c if c.isalnum() or c in ('_', '-') else '_' for c in ecc_name)
        
        process_temp_dir = Path(tempfile.mkdtemp(prefix=f"ecc_{safe_ecc_name}_", dir=self.analyzer_main_temp_dir))
        logging.debug(f"Created process temp dir: {process_temp_dir} for {ecc_name}")

        final_direct_df = pd.DataFrame()
        final_inverted_df = pd.DataFrame()
        
        try:
            start_fasta_path, end_fasta_path = self.get_flank_fastas(ecc_data_tuple_for_blast, process_temp_dir)
            if not start_fasta_path or not end_fasta_path: 
                logging.warning(f"Could not obtain one or both FASTA flank sequences for {ecc_name}. Skipping BLAST for this eccDNA.")
                return final_direct_df, final_inverted_df 

            logging.debug(f"Processing direct repeats for {ecc_name}")
            direct_blast_df = self.run_blast_and_get_df(start_fasta_path, end_fasta_path, ecc_name, process_temp_dir)
            if not direct_blast_df.empty:
                final_direct_df = direct_blast_df
            
            logging.debug(f"Processing inverted repeats for {ecc_name}")
            reversed_end_fasta_path = process_temp_dir / f"{Path(end_fasta_path.name).stem}_reversed.fasta"
            if self.reverse_fasta(end_fasta_path, reversed_end_fasta_path):
                inverted_blast_df = self.run_blast_and_get_df(start_fasta_path, reversed_end_fasta_path, ecc_name, process_temp_dir)
                if not inverted_blast_df.empty:
                    final_inverted_df = inverted_blast_df
            else:
                logging.error(f"Failed to reverse end FASTA for {ecc_name}")
            
        except Exception as e:
            logging.error(f"Unhandled error processing eccDNA {ecc_name}: {e}", exc_info=True)
        finally:
            try:
                shutil.rmtree(process_temp_dir)
                logging.debug(f"Successfully removed process temp dir: {process_temp_dir}")
            except Exception as e:
                logging.warning(f"Failed to remove process temp dir {process_temp_dir}: {e}")
        
        return final_direct_df, final_inverted_df

    def process(self): # Required by BaseProcessor
        logging.debug("BlastAnalyzer.process() called, but primary logic is in process_single_eccdna driven by MainProcessor.")
        pass

    def cleanup(self):
        try:
            if hasattr(self, 'analyzer_main_temp_dir') and self.analyzer_main_temp_dir.exists():
                shutil.rmtree(self.analyzer_main_temp_dir)
                logging.info(f"BlastAnalyzer main temporary directory removed: {self.analyzer_main_temp_dir}")
        except Exception as e:
            logging.warning(f"Failed to remove BlastAnalyzer main temporary directory {self.analyzer_main_temp_dir}: {e}")


class RepeatResultsProcessor(BaseProcessor):
    def __init__(self, input_df_or_file, output_file_path, is_inverted=False, breakpoint_tolerance=5):
        super().__init__() 
        self.input_df_or_file = input_df_or_file 
        self.output_file_path = Path(output_file_path)
        self.is_inverted = is_inverted
        self.breakpoint_tolerance = breakpoint_tolerance

    def process(self):
        df_blast = None
        input_source_name = ""
        if isinstance(self.input_df_or_file, pd.DataFrame):
            df_blast = self.input_df_or_file.copy() 
            input_source_name = "DataFrame input"
        elif isinstance(self.input_df_or_file, (str, Path)):
            input_path = Path(self.input_df_or_file)
            input_source_name = str(input_path)
            if input_path.exists() and input_path.stat().st_size > 0:
                try:
                    df_blast = pd.read_csv(input_path)
                except pd.errors.EmptyDataError:
                    logging.warning(f"Input file {input_source_name} is empty. No results to process.")
                    self.output_file_path.touch() 
                    return
                except Exception as e:
                    logging.error(f"Error reading {input_source_name}: {e}", exc_info=True)
                    self.output_file_path.touch()
                    return
            else:
                logging.warning(f"Input file for RepeatResultsProcessor is missing or empty: {input_source_name}. Creating empty output.")
                self.output_file_path.touch()
                return
        else:
            logging.warning(f"Input for RepeatResultsProcessor is invalid type: {type(self.input_df_or_file)}. Creating empty output.")
            self.output_file_path.touch()
            return

        if df_blast is None or df_blast.empty:
            logging.info(f"No BLAST results from {input_source_name} to process for RepeatResultsProcessor.")
            self.output_file_path.touch()
            return

        required_cols = ['qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 
                         'sstart', 'send', 'eccDNA_name_internal', 'qlen', 'slen']
        if not all(col in df_blast.columns for col in required_cols):
            missing = [col for col in required_cols if col not in df_blast.columns]
            logging.error(f"Missing one or more required columns in input from {input_source_name} for RepeatResultsProcessor. "
                          f"Columns found: {df_blast.columns.tolist()}. Missing: {missing}")
            self.output_file_path.touch()
            return

        numeric_cols = ['pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen', 'evalue', 'bitscore']
        for col in numeric_cols:
            if col in df_blast.columns:
                df_blast[col] = pd.to_numeric(df_blast[col], errors='coerce')
        
        df_blast.dropna(subset=['pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen'], inplace=True)
        
        if df_blast.empty:
            logging.info(f"No data after numeric conversion/dropna from {input_source_name} for RepeatResultsProcessor.")
            self.output_file_path.touch()
            return

        try:
            q_coords = df_blast['qseqid'].str.extract(r'([^:]*):(\d+)-(\d+)', expand=True)
            s_coords = df_blast['sseqid'].str.extract(r'([^:]*):(\d+)-(\d+)', expand=True)

            df_blast['q_flank_chr'] = q_coords[0]
            q_flank_genome_start_1based = pd.to_numeric(q_coords[1], errors='coerce')
            
            df_blast['QR_start'] = q_flank_genome_start_1based + df_blast['qstart'] - 1 
            df_blast['QR_end'] = q_flank_genome_start_1based + df_blast['qend'] - 1     
            
            df_blast['s_flank_chr'] = s_coords[0]
            s_flank_genome_start_1based = pd.to_numeric(s_coords[1], errors='coerce')
            
            if self.is_inverted:
                df_blast['SR_start'] = s_flank_genome_start_1based + (df_blast['slen'] - df_blast['send']) 
                df_blast['SR_end'] = s_flank_genome_start_1based + (df_blast['slen'] - df_blast['sstart'])   
                
                mask_inverted_swap = df_blast['SR_start'] > df_blast['SR_end']
                if mask_inverted_swap.any(): 
                    df_blast.loc[mask_inverted_swap, ['SR_start', 'SR_end']] = \
                        df_blast.loc[mask_inverted_swap, ['SR_end', 'SR_start']].values
            else: 
                df_blast['SR_start'] = s_flank_genome_start_1based + df_blast['sstart'] - 1 
                df_blast['SR_end'] = s_flank_genome_start_1based + df_blast['send'] - 1     

            ecc_info = df_blast['eccDNA_name_internal'].str.split('_', expand=True)
            df_blast['ecc_chr'] = ecc_info.get(0) 
            df_blast['ecc_genome_start_0based'] = pd.to_numeric(ecc_info.get(1), errors='coerce')
            df_blast['ecc_genome_end_1based_exclusive'] = pd.to_numeric(ecc_info.get(2), errors='coerce') 
            
            coord_cols_to_check = ['QR_start', 'QR_end', 'SR_start', 'SR_end',
                                   'ecc_genome_start_0based', 'ecc_genome_end_1based_exclusive']
            df_blast.dropna(subset=coord_cols_to_check, inplace=True) 
            if df_blast.empty: 
                logging.info(f"No data after coordinate remapping/dropna from {input_source_name} for RepeatResultsProcessor.")
                self.output_file_path.touch()
                return

            BP1 = df_blast['ecc_genome_start_0based'] + 1 
            BP2 = df_blast['ecc_genome_end_1based_exclusive']
            
            cond_QR_start_near_BP1 = (df_blast['QR_start'] - BP1).abs() <= self.breakpoint_tolerance
            cond_QR_end_near_BP1   = (df_blast['QR_end'] - BP1).abs() <= self.breakpoint_tolerance
            
            cond_SR_start_near_BP2 = (df_blast['SR_start'] - BP2).abs() <= self.breakpoint_tolerance
            cond_SR_end_near_BP2   = (df_blast['SR_end'] - BP2).abs() <= self.breakpoint_tolerance

            pattern1 = cond_QR_end_near_BP1 & cond_SR_start_near_BP2
            pattern2 = cond_QR_start_near_BP1 & cond_SR_end_near_BP2
            pattern3 = cond_QR_end_near_BP1 & cond_SR_end_near_BP2
            pattern4 = cond_QR_start_near_BP1 & cond_SR_start_near_BP2
            
            combined_filter = pattern1 | pattern2 | pattern3 | pattern4
            df_final = df_blast[combined_filter].copy()
            
            if not df_final.empty: 
                df_final.loc[:, 'q_flank_region'] = df_blast.loc[df_final.index, 'qseqid'] 
                df_final.loc[:, 's_flank_region'] = df_blast.loc[df_final.index, 'sseqid']

            if df_final.empty:
                logging.info(f"No results from {input_source_name} satisfied any of the defined breakpoint proximity patterns (tolerance: {self.breakpoint_tolerance}bp).")
            else:
                logging.info(f"Processed {input_source_name}: {len(df_blast)} initial candidates -> {len(df_final)} final repeats based on flexible patterns.")
            
            df_final.to_csv(self.output_file_path, index=False)

        except Exception as e:
            logging.error(f"Unexpected error processing results from {input_source_name} in RepeatResultsProcessor: {e}", exc_info=True)
            self.output_file_path.touch()


class MainProcessor(BaseProcessor):
    def __init__(self, input_file, genome_file, output_file=None, 
                 extend_outer=100, extend_inner=50, threads=None, 
                 fai_file=None, custom_temp_dir=None,
                 breakpoint_tolerance=5, min_blast_identity=90.0, min_blast_length=10):
        
        super().__init__(genome_file, input_file, extend_outer, extend_inner)
        self.output_file = Path(output_file if output_file else Path(self.bed_file).parent / f"{Path(self.bed_file).stem}_repeat_results.csv")
        self.threads = threads if threads is not None and threads > 0 else max(1, os.cpu_count() // 2)
        self.fai_file = fai_file
        self.name_map = {} 
        self.breakpoint_tolerance = breakpoint_tolerance
        self.min_blast_identity = min_blast_identity
        self.min_blast_length = min_blast_length

        if custom_temp_dir:
            Path(custom_temp_dir).mkdir(parents=True, exist_ok=True)
            self.run_temp_dir = Path(tempfile.mkdtemp(prefix="ecc_repeats_main_", dir=str(custom_temp_dir)))
        else:
            self.run_temp_dir = Path(tempfile.mkdtemp(prefix="ecc_repeats_main_"))
        logging.info(f"Main run temporary directory: {self.run_temp_dir}")

    def _parse_input_line(self, line, line_num):
        fields = line.strip().split('\t')
        original_line_content = line.strip()

        if len(fields) < 4: 
            logging.warning(f"Skipping malformed line {line_num} (expected at least 4 columns 'name chr start end', got {len(fields)}): {original_line_content}")
            return None
        
        original_name = fields[0] 
        ecc_chr = fields[1]       
        try:
            # *** COORDINATE SYSTEM CONVERSION FOR USER'S INPUT FILE ***
            # Assumes input file's column 3 (fields[2]) is 0-based start (inclusive).
            # Assumes input file's column 4 (fields[3]) is 0-based end (inclusive).
            
            input_start_0based_inclusive = int(fields[2]) 
            input_end_0based_inclusive = int(fields[3])

            # Convert to internal representation: 0-based start (inclusive), 1-based end (exclusive).
            # This is the standard for BED files and how bedtools getfasta expects its input BED regions.
            ecc_start_0based = input_start_0based_inclusive
            ecc_end_1based_exclusive = input_end_0based_inclusive + 1

        except ValueError:
            logging.warning(f"Skipping malformed line {line_num} (non-integer coordinates for start (col 3) or end (col 4)): {original_line_content}")
            return None

        if ecc_start_0based < 0 or ecc_end_1based_exclusive <= ecc_start_0based:
            logging.warning(f"Skipping malformed/invalid line {line_num} (start < 0 or end <= start after conversion): {original_line_content}")
            return None
        
        # Standardized name uses the internal 0-based start, 1-based exclusive end for consistency
        standardized_name = f"{ecc_chr}_{ecc_start_0based}_{ecc_end_1based_exclusive}" 
        self.name_map[standardized_name] = original_name
        return (ecc_chr, ecc_start_0based, ecc_end_1based_exclusive, standardized_name)


    def preprocess_input(self):
        try:
            input_path = self.validate_input_file(self.bed_file) 
        except FileNotFoundError:
            return None 

        eccdnas_for_blast = [] 
        try:
            with open(input_path, 'r') as f_in:
                line_num = 0
                first_line_peek = f_in.readline()
                if first_line_peek.startswith('#') or \
                   any(hdr_start.lower() in first_line_peek.lower() for hdr_start in ['chrom', 'ename', 'name', 'chr\t']):
                    logging.info(f"Skipping assumed header line: {first_line_peek.strip()}")
                else:
                    f_in.seek(0) 

                for line in f_in:
                    line_num += 1
                    if line.startswith('#') or not line.strip(): 
                        continue
                    parsed_data = self._parse_input_line(line, line_num)
                    if parsed_data:
                        eccdnas_for_blast.append(parsed_data)
            
            logging.info(f"Preprocessed {len(eccdnas_for_blast)} eccDNA records from {input_path}")
            return eccdnas_for_blast
        except Exception as e:
            logging.error(f"Error preprocessing input file {input_path}: {e}", exc_info=True)
            return None

    def process(self):
        all_eccdna_tuples = self.preprocess_input() 
        if not all_eccdna_tuples:
            logging.error(f"No valid eccDNA records to process from {self.bed_file}. Exiting.")
            self.cleanup_main_temp_dir()
            return

        blast_analyzer = BlastAnalyzer(
            self.genome_file,
            extend_outer=self.extend_outer,
            extend_inner=self.extend_inner,
            fai_file=self.fai_file,
            parent_temp_dir=self.run_temp_dir,
            min_blast_identity=self.min_blast_identity,
            min_blast_length=self.min_blast_length
        )

        all_direct_blast_dfs = []
        all_inverted_blast_dfs = []
        
        logging.info(f"Starting BLAST analysis for {len(all_eccdna_tuples)} eccDNAs using up to {self.threads} threads.")
        
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            future_results = [executor.submit(blast_analyzer.process_single_eccdna, ecc_tuple) for ecc_tuple in all_eccdna_tuples]
            
            for i, future in enumerate(future_results):
                try:
                    direct_df, inverted_df = future.result() 
                    if direct_df is not None and not direct_df.empty:
                        all_direct_blast_dfs.append(direct_df)
                    if inverted_df is not None and not inverted_df.empty:
                        all_inverted_blast_dfs.append(inverted_df)
                except Exception as e: 
                    ecc_name_for_log = all_eccdna_tuples[i][3] if i < len(all_eccdna_tuples) else f"future_index_{i}"
                    logging.error(f"Error processing an eccDNA future (name: {ecc_name_for_log}): {e}", exc_info=True)

                if (i + 1) % 100 == 0 or (i + 1) == len(all_eccdna_tuples):
                    logging.info(f"BLAST processing status: {i + 1}/{len(all_eccdna_tuples)} eccDNAs futures collected.")
        
        blast_analyzer.cleanup() 
        logging.info(f"BLAST analysis complete.")

        master_direct_blast_df = pd.concat(all_direct_blast_dfs, ignore_index=True) if all_direct_blast_dfs else pd.DataFrame()
        master_inverted_blast_df = pd.concat(all_inverted_blast_dfs, ignore_index=True) if all_inverted_blast_dfs else pd.DataFrame()

        processed_direct_path = self.run_temp_dir / "filtered_direct_repeats.csv"
        processed_inverted_path = self.run_temp_dir / "filtered_inverted_repeats.csv"

        logging.info("Filtering all direct repeat BLAST hits...")
        direct_processor = RepeatResultsProcessor(
            master_direct_blast_df, 
            processed_direct_path,
            is_inverted=False,
            breakpoint_tolerance=self.breakpoint_tolerance
        )
        direct_processor.process()

        logging.info("Filtering all inverted repeat BLAST hits...")
        inverted_processor = RepeatResultsProcessor(
            master_inverted_blast_df, 
            processed_inverted_path,
            is_inverted=True,
            breakpoint_tolerance=self.breakpoint_tolerance
        )
        inverted_processor.process()

        final_dfs = []
        if processed_direct_path.exists() and processed_direct_path.stat().st_size > 0:
            try:
                df_direct = pd.read_csv(processed_direct_path)
                if not df_direct.empty:
                    df_direct['repeat_type'] = 'Direct'
                    final_dfs.append(df_direct)
            except pd.errors.EmptyDataError:
                 logging.warning(f"Filtered direct repeats file {processed_direct_path} was empty or unreadable.")
        
        if processed_inverted_path.exists() and processed_inverted_path.stat().st_size > 0:
            try:
                df_inverted = pd.read_csv(processed_inverted_path)
                if not df_inverted.empty:
                    df_inverted['repeat_type'] = 'Inverted'
                    final_dfs.append(df_inverted)
            except pd.errors.EmptyDataError:
                 logging.warning(f"Filtered inverted repeats file {processed_inverted_path} was empty or unreadable.")

        # Define comprehensive output columns, these should match columns created in RepeatResultsProcessor + added ones
        final_output_columns = [ 
            'eccDNA_std_name', 'eccDNA_orig_name', 
            'ecc_chr', 'ecc_genome_start_0based', 'ecc_genome_end_1based_exclusive', # Changed name for clarity
            'repeat_type', 
            'q_flank_chr', 'QR_start', 'QR_end', # Changed to QR_start, QR_end for consistency
            's_flank_chr', 'SR_start', 'SR_end', # Changed to SR_start, SR_end for consistency
            'pident', 'length', 'evalue', 'bitscore', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 
            'qlen', 'slen', 
            'q_flank_region', 's_flank_region', 
            'qseq', 'sseq' 
        ]

        if not final_dfs:
            logging.warning("No repeats found after all processing. Output file will be empty or just headers.")
            try:
                with open(self.output_file, 'w') as f_out:
                    f_out.write(",".join(final_output_columns) + "\n")
            except IOError as e:
                logging.error(f"Could not write empty output file {self.output_file}: {e}")
        else:
            df_combined = pd.concat(final_dfs, ignore_index=True)
            if 'eccDNA_name_internal' in df_combined.columns: 
                 df_combined['eccDNA_orig_name'] = df_combined['eccDNA_name_internal'].map(self.name_map)
                 df_combined.rename(columns={'eccDNA_name_internal': 'eccDNA_std_name',
                                             'ecc_genome_end_1based': 'ecc_genome_end_1based_exclusive', # if old name exists
                                             'q_real_start_genome_1based': 'QR_start', # if old name exists
                                             'q_real_end_genome_1based': 'QR_end',     # if old name exists
                                             's_real_start_genome_1based': 'SR_start', # if old name exists
                                             's_real_end_genome_1based': 'SR_end'},    # if old name exists
                                   inplace=True)
            else:
                logging.warning("Column 'eccDNA_name_internal' not found for mapping original names. Adding placeholder columns.")
                if 'eccDNA_orig_name' not in df_combined.columns: df_combined['eccDNA_orig_name'] = "N/A" 
                if 'eccDNA_std_name' not in df_combined.columns: df_combined['eccDNA_std_name'] = "N/A_std"
                if 'ecc_genome_end_1based_exclusive' not in df_combined.columns and 'ecc_genome_end_1based' in df_combined.columns:
                     df_combined.rename(columns={'ecc_genome_end_1based': 'ecc_genome_end_1based_exclusive'}, inplace=True)


            # Ensure all desired columns are present, fill with NA if missing before final selection
            for col in final_output_columns:
                if col not in df_combined.columns:
                    # Check for old names in case rename didn't catch all from older intermediate files
                    old_name_map = {
                        'QR_start': 'q_real_start_genome_1based', 'QR_end': 'q_real_end_genome_1based',
                        'SR_start': 's_real_start_genome_1based', 'SR_end': 's_real_end_genome_1based',
                        'ecc_genome_end_1based_exclusive':'ecc_genome_end_1based'
                    }
                    if col in old_name_map and old_name_map[col] in df_combined.columns:
                        df_combined.rename(columns={old_name_map[col]: col}, inplace=True)
                    else:
                         df_combined[col] = pd.NA 
            
            df_combined_final = df_combined[final_output_columns] 

            try:
                df_combined_final.to_csv(self.output_file, index=False)
                logging.info(f"Combined results with {len(df_combined_final)} entries saved to {self.output_file}")
            except IOError as e:
                logging.error(f"Could not write final output file {self.output_file}: {e}")


        self.cleanup_main_temp_dir()

    def cleanup_main_temp_dir(self):
        try:
            if hasattr(self, 'run_temp_dir') and self.run_temp_dir.exists():
                shutil.rmtree(self.run_temp_dir)
                logging.info(f"Main run temporary directory removed: {self.run_temp_dir}")
        except Exception as e:
            logging.warning(f"Failed to remove main run temporary directory {self.run_temp_dir}: {e}")

def check_dependencies():
    tools = ["bedtools", "blastn"]
    all_found = True
    for tool in tools:
        if shutil.which(tool) is None:
            logging.error(f"FATAL: Required tool '{tool}' not found in PATH. Please install it.")
            all_found = False
    if not all_found:
        exit(1)
    logging.info("All required tools (bedtools, blastn) found in PATH.")

def main():
    check_dependencies() 

    parser = argparse.ArgumentParser(
        description='Process eccDNA data to find terminal direct or inverted repeats.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    
    parser.add_argument('-i', '--input', required=True, 
                        help="Input tab-separated file. Expected format: name (col1), chr (col2), "
                             "start (col3, 0-based inclusive), end (col4, 0-based inclusive).")
    parser.add_argument('-g', '--genome', required=True, help='Genome FASTA file path.')
    parser.add_argument('-o', '--output', default=None, help='Output CSV file name (optional, default: <input_stem>_repeat_results.csv).')
    parser.add_argument('-t', '--threads', type=int, default=None, help='Number of threads for parallel processing (default: auto-detect, uses about half of CPU cores).')
    
    parser.add_argument('--extend_outer', type=int, default=100, help="Base pairs to extend outward from eccDNA ends.")
    parser.add_argument('--extend_inner', type=int, default=50, help="Base pairs to extend inward from eccDNA ends.")
    
    parser.add_argument('--fai', default=None, help='Path to the genome .fai index file (recommended for accurate coordinate clamping).')
    parser.add_argument('--custom_temp_dir', default=None, help='Optional: Specify a parent directory for all temporary files.')
    
    filter_group = parser.add_argument_group('Filtering options')
    filter_group.add_argument('--breakpoint_tolerance', type=int, default=5, help="Tolerance in bp for matching repeat ends to eccDNA breakpoints.")
    filter_group.add_argument('--min_blast_identity', type=float, default=90.0, help="Minimum BLAST percent identity.")
    filter_group.add_argument('--min_blast_length', type=int, default=10, help="Minimum BLAST alignment length.")

    args = parser.parse_args()

    try:
        if not Path(args.input).exists():
            logging.error(f"Input file not found: {args.input}")
            exit(1)
        if not Path(args.genome).exists():
            logging.error(f"Genome file not found: {args.genome}")
            exit(1)
        if args.fai and not Path(args.fai).exists():
            logging.warning(f"FAI file specified ({args.fai}) but not found. Proceeding with limited or no FAI-based clamping.")
    except Exception as e: 
        logging.error(f"Error with input file paths: {e}")
        exit(1)

    logging.info(f"Starting eccDNA repeat analysis with parameters: {vars(args)}")

    processor = MainProcessor(
        input_file=args.input,
        genome_file=args.genome,
        output_file=args.output,
        extend_outer=args.extend_outer,
        extend_inner=args.extend_inner,
        threads=args.threads,
        fai_file=args.fai,
        custom_temp_dir=args.custom_temp_dir,
        breakpoint_tolerance=args.breakpoint_tolerance,
        min_blast_identity=args.min_blast_identity,
        min_blast_length=args.min_blast_length
    )
    try:
        processor.process()
    except Exception as e: 
        logging.error("An unhandled exception occurred in MainProcessor.process:", exc_info=True)
        if 'processor' in locals() and hasattr(processor, 'cleanup_main_temp_dir'):
            processor.cleanup_main_temp_dir() 
    finally:
        logging.info("eccDNA Repeat Classifier (Refactored) finished.")

if __name__ == "__main__":
    main()
