import argparse
import os
import subprocess
import sys
from pathlib import Path

# --- Helper Functions / 辅助函数 ---


class PipelineError(Exception):
    """
    Custom exception for pipeline errors.
    中文: 用于流程错误的自定义异常。
    """

    pass


def run_command(command, log_file=None):
    """
    Executes a shell command and handles potential errors, with optional logging.

    Args:
        command (str): The command to execute.
        log_file (str, optional): Path to a file for logging stdout/stderr. Defaults to None.

    Returns:
        bool: True if the command was successful, False otherwise.

    中文: 执行一个 shell 命令并处理潜在错误，支持可选的日志记录功能。

    参数:
        command (str): 需要执行的命令。
        log_file (str, optional): 用于记录标准输出/错误的日志文件路径。默认为 None。

    返回:
        bool: 命令成功执行则返回 True，否则返回 False。
    """
    try:
        print(
            f"  Executing: {command}"
        )  # English: Print command being executed / 中文: 打印正在执行的命令
        if log_file:
            with open(log_file, "w") as f:
                subprocess.run(
                    command, shell=True, check=True, stdout=f, stderr=subprocess.STDOUT
                )
        else:
            # When not logging, capture output for better error display if it fails
            # 中文: 不记录日志时，捕获输出以便更好地显示错误信息
            subprocess.run(
                command, shell=True, check=True, capture_output=True, text=True
            )
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: Command failed with exit code {e.returncode}", file=sys.stderr)
        print(f"  Failed Command: {command}", file=sys.stderr)
        if e.stdout:
            print(f"  Stdout: {e.stdout}", file=sys.stderr)
        if e.stderr:
            print(f"  Stderr: {e.stderr}", file=sys.stderr)
        return False


def check_output_exists(*files):
    """
    Checks if all specified output files exist.
    中文: 检查所有指定的输出文件是否存在。
    """
    return all(os.path.exists(f) for f in files)


# --- Main Pipeline Class / 主流程类 ---


class EccDNAPipeline:
    def __init__(self, args):
        """
        Initializes the pipeline with user-provided arguments.
        中文: 使用用户提供的参数初始化流程。
        """
        self.threads = args.threads
        self.fastq1 = os.path.abspath(
            args.fastq1
        )  # English: Use absolute paths for robustness / 中文: 使用绝对路径以增强稳健性
        self.fastq2 = os.path.abspath(args.fastq2)
        self.reference = os.path.abspath(args.reference)
        self.output_dir = os.path.abspath(args.output_dir)
        self.sample_name = args.sample_name

        # English: Create output directory if it doesn't exist.
        # 中文: 如果输出目录不存在，则创建它。
        os.makedirs(self.output_dir, exist_ok=True)

    def _get_path(self, filename):
        """
        Helper method to get the full path for an output file.
        中文: 获取输出文件完整路径的辅助方法。
        """
        return os.path.join(self.output_dir, filename)

    def run_fastp(self):
        """
        Runs the fastp step for quality control.
        中文: 运行 fastp 步骤进行质量控制。
        """
        output_files = [
            self._get_path(f"{self.sample_name}.R1.fp.fastq.gz"),
            self._get_path(f"{self.sample_name}.R2.fp.fastq.gz"),
        ]

        if check_output_exists(*output_files):
            print("  Fastp outputs already exist, skipping...")
            return True

        cmd = (
            f"fastp -w {self.threads} -i {self.fastq1} -I {self.fastq2} "
            f"-o {output_files[0]} "
            f"-O {output_files[1]} "
            f"-h {self._get_path(f'{self.sample_name}.report.html')} "
            f"-j {self._get_path(f'{self.sample_name}.report.json')}"
        )

        return run_command(cmd, self._get_path(f"fastp_{self.sample_name}.log"))

    def run_bwa(self):
        """
        Runs the BWA MEM step for alignment.
        中文: 运行 BWA MEM 步骤进行序列比对。
        """
        output_file = self._get_path(f"{self.sample_name}.sam")

        if check_output_exists(output_file):
            print("  BWA SAM output already exists, skipping...")
            return True

        # English: Define input files using the _get_path helper for clarity
        # 中文: 为清晰起见，使用 _get_path 辅助函数定义输入文件
        r1_fp = self._get_path(f"{self.sample_name}.R1.fp.fastq.gz")
        r2_fp = self._get_path(f"{self.sample_name}.R2.fp.fastq.gz")

        cmd = f"bwa mem -t {self.threads} {self.reference} {r1_fp} {r2_fp} > {output_file}"

        return run_command(cmd)

    def run_samtools_sort_qname(self):
        """
        Runs samtools to sort the alignment by query name.
        中文: 运行 samtools 按查询名称（query name）排序比对文件。
        """
        sam_input = self._get_path(f"{self.sample_name}.sam")
        output_file = self._get_path(f"qname_eccDNA_{self.sample_name}.bam")

        if check_output_exists(output_file):
            print("  Samtools qname-sorted BAM already exists, skipping...")
            return True

        cmd = f"samtools sort -@ {self.threads} -n -o {output_file} {sam_input}"

        return run_command(cmd)

    def run_samtools_sort_coordinate(self):
        """
        Runs samtools to sort the alignment by coordinate and then index it.
        中文: 运行 samtools 按坐标排序比对文件并建立索引。
        """
        sam_input = self._get_path(f"{self.sample_name}.sam")
        output_bam = self._get_path(f"sorted_eccDNA_{self.sample_name}.bam")

        if check_output_exists(output_bam, f"{output_bam}.bai"):
            print(
                "  Samtools coordinate-sorted BAM and index already exist, skipping..."
            )
            return True

        cmd1 = f"samtools sort -@ {self.threads} -o {output_bam} {sam_input}"
        cmd2 = f"samtools index -@ {self.threads} {output_bam}"

        # English: Both commands must succeed.
        # 中文: 两个命令都必须成功执行。
        return run_command(cmd1) and run_command(cmd2)

    def run_circle_map_extractor(self):
        """
        Runs the Circle-Map ReadExtractor to find candidate reads.
        中文: 运行 Circle-Map ReadExtractor 寻找候选 reads。
        """
        qname_bam = self._get_path(f"qname_eccDNA_{self.sample_name}.bam")
        output_file = self._get_path(f"eccDNA_{self.sample_name}_candidates.bam")

        if check_output_exists(output_file):
            print("  Circle-Map ReadExtractor output already exists, skipping...")
            return True

        cmd = f"Circle-Map ReadExtractor -i {qname_bam} -o {output_file}"

        return run_command(
            cmd, self._get_path(f"circle_map_extractor_{self.sample_name}.log")
        )

    def run_samtools_sort_candidates(self):
        """
        Sorts and indexes the candidate reads BAM file.
        中文: 排序并索引候选 reads 的 BAM 文件。
        """
        candidates_bam = self._get_path(f"eccDNA_{self.sample_name}_candidates.bam")
        output_bam = self._get_path(f"sort_eccDNA_{self.sample_name}_candidates.bam")

        if check_output_exists(output_bam, f"{output_bam}.bai"):
            print("  Sorted candidate BAM and index already exist, skipping...")
            return True

        cmd1 = f"samtools sort -@ {self.threads} {candidates_bam} -o {output_bam}"
        cmd2 = f"samtools index -@ {self.threads} {output_bam}"

        return run_command(cmd1) and run_command(cmd2)

    def run_circle_map_realign(self):
        """
        Runs the Circle-Map Realign step to identify circular DNA.
        中文: 运行 Circle-Map Realign 步骤来鉴定环状 DNA。
        """
        sorted_candidates_bam = self._get_path(
            f"sort_eccDNA_{self.sample_name}_candidates.bam"
        )
        qname_bam = self._get_path(f"qname_eccDNA_{self.sample_name}.bam")
        coord_sorted_bam = self._get_path(f"sorted_eccDNA_{self.sample_name}.bam")
        output_bed = self._get_path(f"eccDNA_{self.sample_name}_CM.bed")

        if check_output_exists(output_bed):
            print("  Circle-Map Realign output already exists, skipping...")
            return True

        cmd = (
            f"Circle-Map Realign -t {self.threads} -i {sorted_candidates_bam} "
            f"-qbam {qname_bam} -sbam {coord_sorted_bam} -fasta {self.reference} "
            f"-o {output_bed}"
        )

        return run_command(
            cmd, self._get_path(f"circle_map_realign_{self.sample_name}.log")
        )

    def execute(self):
        """
        Executes the full pipeline step-by-step.
        中文: 按顺序执行完整的流程。
        """
        steps = [
            (self.run_fastp, "Step 1: Quality Control (fastp)"),
            (self.run_bwa, "Step 2: Alignment (BWA-MEM)"),
            (self.run_samtools_sort_qname, "Step 3: Sort by Query Name (samtools)"),
            (
                self.run_samtools_sort_coordinate,
                "Step 4: Sort by Coordinate (samtools)",
            ),
            (self.run_circle_map_extractor, "Step 5: Extract Candidates (Circle-Map)"),
            (self.run_samtools_sort_candidates, "Step 6: Sort Candidates (samtools)"),
            (self.run_circle_map_realign, "Step 7: Realign and Detect (Circle-Map)"),
        ]

        print(f"🚀 Starting eccDNA detection pipeline for sample: {self.sample_name}")
        print(f"   Output directory: {self.output_dir}")

        for step_func, step_name in steps:
            print(f"\n▶️ Executing {step_name}...")
            if not step_func():
                print(f"❌ Pipeline failed at step: {step_name}", file=sys.stderr)
                return False
            print(f"✅ {step_name} completed successfully.")

        print(f"\n🎉 Pipeline completed successfully for sample: {self.sample_name}")
        return True


# --- Entry Point / 程序入口 ---


def main():
    """
    Main function to parse arguments and run the pipeline.
    中文: 解析命令行参数并运行流程的主函数。
    """
    parser = argparse.ArgumentParser(
        description="A comprehensive pipeline to detect eccDNA using Circle-Map.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=True,
        help="Number of threads for all steps.",
    )
    parser.add_argument(
        "-1", "--fastq1", required=True, help="Path to input FASTQ file (Read 1)."
    )
    parser.add_argument(
        "-2", "--fastq2", required=True, help="Path to input FASTQ file (Read 2)."
    )
    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="Path to the reference genome FASTA file.",
    )
    parser.add_argument(
        "-o", "--output_dir", required=True, help="Directory to store all output files."
    )
    parser.add_argument(
        "-s",
        "--sample_name",
        required=True,
        help="A unique name for the sample, used for output files.",
    )

    args = parser.parse_args()

    # English: Validate that input files exist before starting.
    # 中文: 在开始前，验证所有输入文件都存在。
    for file_path in [args.fastq1, args.fastq2, args.reference]:
        if not os.path.exists(file_path):
            print(f"Error: Input file does not exist at '{file_path}'", file=sys.stderr)
            sys.exit(1)

    try:
        pipeline = EccDNAPipeline(args)
        success = pipeline.execute()
        if success:
            sys.exit(0)
        else:
            sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
