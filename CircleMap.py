import argparse
import os
import subprocess
import sys
from pathlib import Path

class PipelineError(Exception):
    """Custom exception for pipeline errors"""
    pass

def run_command(command, log_file=None):
    """执行命令并处理可能的错误"""
    try:
        if log_file:
            with open(log_file, 'w') as f:
                result = subprocess.run(command, shell=True, check=True,
                                     stdout=f, stderr=subprocess.STDOUT)
        else:
            result = subprocess.run(command, shell=True, check=True,
                                  capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {command}")
        print(f"Error message: {e.stderr if e.stderr else 'No error message'}")
        return False

def check_output_exists(*files):
    """检查所有指定的输出文件是否存在"""
    return all(os.path.exists(f) for f in files)

class Pipeline:
    def __init__(self, args):
        self.threads = args.threads
        self.fastq1 = args.fastq1
        self.fastq2 = args.fastq2
        self.reference = args.reference
        self.output_dir = args.output_dir
        self.sample_name = args.sample_name
        
        # 创建输出目录
        os.makedirs(self.output_dir, exist_ok=True)
        os.chdir(self.output_dir)

    def run_fastp(self):
        """运行fastp步骤"""
        output_files = [
            f"{self.sample_name}.R1.fp.fastq.gz",
            f"{self.sample_name}.R2.fp.fastq.gz",
            f"{self.sample_name}.report.html",
            f"{self.sample_name}.report.json"
        ]
        
        if check_output_exists(*output_files):
            print("Fastp outputs already exist, skipping...")
            return True
            
        cmd = (f"fastp -w 16 -i {self.fastq1} -I {self.fastq2} "
               f"-o {self.sample_name}.R1.fp.fastq.gz "
               f"-O {self.sample_name}.R2.fp.fastq.gz "
               f"-h {self.sample_name}.report.html "
               f"-j {self.sample_name}.report.json")
        
        return run_command(cmd, f"fastp_{self.sample_name}.log")

    def run_bwa(self):
        """运行BWA MEM步骤"""
        output_file = f"{self.sample_name}.sam"
        
        if check_output_exists(output_file):
            print("BWA outputs already exist, skipping...")
            return True
            
        cmd = (f"bwa mem -t {self.threads} {self.reference} "
               f"{self.sample_name}.R1.fp.fastq.gz "
               f"{self.sample_name}.R2.fp.fastq.gz > {output_file}")
        
        return run_command(cmd)

    def run_samtools_sort_qname(self):
        """运行samtools按查询名排序"""
        output_file = f"qname_eccDNA_{self.sample_name}.bam"
        
        if check_output_exists(output_file):
            print("Samtools qname sort outputs already exist, skipping...")
            return True
            
        cmd = (f"samtools sort -@ {self.threads} -n "
               f"-o {output_file} {self.sample_name}.sam")
        
        return run_command(cmd)

    def run_samtools_sort_coordinate(self):
        """运行samtools按坐标排序"""
        output_files = [
            f"sorted_eccDNA_{self.sample_name}.bam",
            f"sorted_eccDNA_{self.sample_name}.bam.bai"
        ]
        
        if check_output_exists(*output_files):
            print("Samtools coordinate sort outputs already exist, skipping...")
            return True
            
        cmd1 = (f"samtools sort -@ {self.threads} "
                f"-o sorted_eccDNA_{self.sample_name}.bam {self.sample_name}.sam")
        cmd2 = (f"samtools index -@ {self.threads} "
                f"sorted_eccDNA_{self.sample_name}.bam")
        
        return run_command(cmd1) and run_command(cmd2)

    def run_circle_map_extractor(self):
        """运行Circle-Map ReadExtractor"""
        output_file = f"eccDNA_{self.sample_name}_candidates.bam"
        
        if check_output_exists(output_file):
            print("Circle-Map ReadExtractor outputs already exist, skipping...")
            return True
            
        cmd = (f"Circle-Map ReadExtractor -i qname_eccDNA_{self.sample_name}.bam "
               f"-o {output_file}")
        
        return run_command(cmd, f"circle_map_read_extractor_{self.sample_name}.log")

    def run_samtools_sort_candidates(self):
        """对candidates进行排序"""
        output_files = [
            f"sort_eccDNA_{self.sample_name}_candidates.bam",
            f"sort_eccDNA_{self.sample_name}_candidates.bam.bai"
        ]
        
        if check_output_exists(*output_files):
            print("Samtools candidates sort outputs already exist, skipping...")
            return True
            
        cmd1 = (f"samtools sort -@ {self.threads} "
                f"eccDNA_{self.sample_name}_candidates.bam "
                f"-o sort_eccDNA_{self.sample_name}_candidates.bam")
        cmd2 = (f"samtools index -@ {self.threads} "
                f"sort_eccDNA_{self.sample_name}_candidates.bam")
        
        return run_command(cmd1) and run_command(cmd2)

    def run_circle_map_realign(self):
        """运行Circle-Map Realign"""
        output_file = f"eccDNA_{self.sample_name}_CM.bed"
        
        if check_output_exists(output_file):
            print("Circle-Map Realign outputs already exist, skipping...")
            return True
            
        cmd = (f"Circle-Map Realign -t {self.threads} "
               f"-i sort_eccDNA_{self.sample_name}_candidates.bam "
               f"-qbam qname_eccDNA_{self.sample_name}.bam "
               f"-sbam sorted_eccDNA_{self.sample_name}.bam "
               f"-fasta {self.reference} "
               f"-o {output_file}")
        
        return run_command(cmd, f"circle_map_realign_{self.sample_name}.log")

    def run_pipeline(self):
        """运行完整的流程"""
        steps = [
            (self.run_fastp, "Fastp"),
            (self.run_bwa, "BWA MEM"),
            (self.run_samtools_sort_qname, "Samtools Sort by Query Name"),
            (self.run_samtools_sort_coordinate, "Samtools Sort by Coordinate"),
            (self.run_circle_map_extractor, "Circle-Map ReadExtractor"),
            (self.run_samtools_sort_candidates, "Samtools Sort Candidates"),
            (self.run_circle_map_realign, "Circle-Map Realign")
        ]

        print(f"Starting pipeline for sample {self.sample_name}")
        
        for step_func, step_name in steps:
            print(f"\nExecuting {step_name}...")
            if not step_func():
                print(f"Error in {step_name} step")
                return False
                
        print(f"\nPipeline completed successfully for sample {self.sample_name}")
        return True

def main():
    parser = argparse.ArgumentParser(description="Pipeline for processing sequencing data")
    parser.add_argument("-t", "--threads", required=True, help="Number of threads")
    parser.add_argument("-1", "--fastq1", required=True, help="Input fastq1 file")
    parser.add_argument("-2", "--fastq2", required=True, help="Input fastq2 file")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome file")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
    parser.add_argument("-s", "--sample_name", required=True, help="Sample name")

    args = parser.parse_args()

    # 验证输入文件是否存在
    for file_path in [args.fastq1, args.fastq2, args.reference]:
        if not os.path.exists(file_path):
            print(f"Error: Input file {file_path} does not exist")
            sys.exit(1)

    pipeline = Pipeline(args)
    success = pipeline.run_pipeline()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
