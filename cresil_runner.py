#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import logging
from datetime import datetime
from pathlib import Path

class CresilPipeline:
    def __init__(self, args):
        self.threads = args.threads
        self.fastq = args.fastq
        self.reference_mmi = args.reference_mmi
        self.reference_fa = args.reference_fa
        self.reference_fai = args.reference_fai
        self.output_dir = args.output_dir
        self.rmsk_bed = args.rmsk_bed
        self.cpg_bed = args.cpg_bed
        self.gene_bed = args.gene_bed
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Setup logging
        self.setup_logging()
    
    def setup_logging(self):
        """Setup logging configuration"""
        log_file = os.path.join(self.output_dir, 'cresil_pipeline.log')
        
        # Configure logging format
        log_format = '%(asctime)s - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        
        # Setup file handler
        logging.basicConfig(
            level=logging.INFO,
            format=log_format,
            datefmt=date_format,
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        
        logging.info("Pipeline initialized")
        logging.info(f"Output directory: {self.output_dir}")
        logging.info(f"Log file created at: {log_file}")
    
    def run_command(self, command):
        """Execute command and handle potential errors"""
        try:
            logging.info(f"Executing command: {command}")
            process = subprocess.Popen(
                command,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )
            
            # Real-time output processing
            while True:
                output = process.stdout.readline()
                error = process.stderr.readline()
                
                if output == '' and error == '' and process.poll() is not None:
                    break
                
                if output:
                    logging.info(output.strip())
                if error:
                    logging.error(error.strip())
            
            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode, command)
                
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed with exit status {e.returncode}")
            logging.error(f"Command: {e.cmd}")
            sys.exit(1)
        except Exception as e:
            logging.error(f"Unexpected error: {str(e)}")
            sys.exit(1)
    
    def validate_input_files(self):
        """Validate the existence of input files"""
        files_to_check = {
            'FASTQ file': self.fastq,
            'Reference MMI': self.reference_mmi,
            'Reference FASTA': self.reference_fa,
            'Reference FAI': self.reference_fai,
            'RMSK BED': self.rmsk_bed,
            'CpG BED': self.cpg_bed,
            'Gene BED': self.gene_bed
        }
        
        for file_name, file_path in files_to_check.items():
            if not os.path.exists(file_path):
                logging.error(f"{file_name} not found: {file_path}")
                sys.exit(1)
            else:
                logging.info(f"{file_name} found: {file_path}")
    
    def trim(self):
        """Run trim step"""
        logging.info("Starting trim step")
        trim_output = os.path.join(self.output_dir, "trim_result")
        command = f"cresil trim -t {self.threads} -fq {self.fastq} -r {self.reference_mmi} -o {trim_output}"
        self.run_command(command)
        trim_file = os.path.join(trim_output, "trim.txt")
        logging.info(f"Trim step completed. Output: {trim_file}")
        return trim_file
    
    def identify(self, trim_file):
        """Run identify step"""
        logging.info("Starting identify step")
        command = (f"cresil identify -t {self.threads} -fa {self.reference_fa} "
                  f"-fai {self.reference_fai} -fq {self.fastq} -trim {trim_file}")
        self.run_command(command)
        identify_file = os.path.join(self.output_dir, "eccDNA_final.txt")
        logging.info(f"Identify step completed. Output: {identify_file}")
        return identify_file
    
    def annotate(self, identify_file):
        """Run annotate step"""
        logging.info("Starting annotate step")
        command = (f"cresil annotate -t {self.threads} -rp {self.rmsk_bed} "
                  f"-cg {self.cpg_bed} -gb {self.gene_bed} -identify {identify_file}")
        self.run_command(command)
        logging.info("Annotate step completed")
    
    def run_pipeline(self):
        """Run the complete pipeline"""
        start_time = datetime.now()
        logging.info(f"Pipeline started at: {start_time}")
        
        # Validate input files
        logging.info("Validating input files...")
        self.validate_input_files()
        
        # Execute trim step
        logging.info("=== Executing Trim Step ===")
        trim_file = self.trim()
        
        # Execute identify step
        logging.info("=== Executing Identify Step ===")
        identify_file = self.identify(trim_file)
        
        # Execute annotate step
        logging.info("=== Executing Annotate Step ===")
        self.annotate(identify_file)
        
        end_time = datetime.now()
        duration = end_time - start_time
        logging.info(f"Pipeline completed at: {end_time}")
        logging.info(f"Total duration: {duration}")
        logging.info(f"Results saved in: {self.output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Cresil eccDNA Analysis Pipeline")
    
    # Required arguments
    parser.add_argument("--threads", "-t", type=int, default=40,
                      help="Number of threads to use")
    parser.add_argument("--fastq", "-fq", required=True,
                      help="Path to input FASTQ file")
    parser.add_argument("--reference-mmi", "-r", required=True,
                      help="Path to reference MMI file")
    parser.add_argument("--reference-fa", "-fa", required=True,
                      help="Path to reference FASTA file")
    parser.add_argument("--reference-fai", "-fai", required=True,
                      help="Path to reference FAI index file")
    parser.add_argument("--output-dir", "-o", required=True,
                      help="Output directory")
    parser.add_argument("--rmsk-bed", "-rp", required=True,
                      help="Path to repeat sequences BED file")
    parser.add_argument("--cpg-bed", "-cg", required=True,
                      help="Path to CpG islands BED file")
    parser.add_argument("--gene-bed", "-gb", required=True,
                      help="Path to gene annotation BED file")
    
    args = parser.parse_args()
    
    # Run pipeline
    pipeline = CresilPipeline(args)
    pipeline.run_pipeline()

if __name__ == "__main__":
    main()
