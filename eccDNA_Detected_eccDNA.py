import argparse
import subprocess
import os
from pathlib import Path

# 初始化参数解析器
parser = argparse.ArgumentParser(description='eccDNA Detection and Processing Script.')
parser.add_argument('-s', '--sample', required=True, help='Sample name')
parser.add_argument('-t', '--threads', type=int, default=12, help='Number of threads (default: 12)')
parser.add_argument('-g', '--genome', required=True, help='Genome reference file')
parser.add_argument('-i', '--input_dir', required=True, help='Input directory')
parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
parser.add_argument('-f', '--overlap_ratio', type=float, default=0.9, help='Overlap ratio for BedTools processing (default: 0.9)')
args = parser.parse_args()

# 定义脚本路径
bin_path = Path(__file__).parent / 'bin'
eccDNA_toolkit_path = bin_path / 'eccDNA_Detect_eccDNA.sh'
bedtools_processor_path = bin_path / 'eccDNA_BedTools_Processor.sh'

# 运行eccDNA检测脚本
subprocess.run([
    "bash", eccDNA_toolkit_path,
    "-s", args.sample,
    "-t", str(args.threads),
    "-g", args.genome,
    "-i", args.input_dir,
    "-o", args.output_dir
], check=True)

# 目录和文件操作
os.chdir(args.output_dir)
os.makedirs("temp_Res_" + args.sample, exist_ok=True)
os.chdir(os.path.join(args.output_dir, args.sample, "align_files"))
for file in os.listdir("."):
    if file.endswith("gz"):
        os.remove(file)

subprocess.run(["cp", "eccDNA_" + args.sample + "_CM.bed", "../temp_Res_" + args.sample])
subprocess.run(["cp", os.path.join(args.output_dir, args.sample, args.sample + ".csv"), "../temp_Res_" + args.sample])

# 运行BedTools处理器脚本
subprocess.run([
    "bash", bedtools_processor_path,
    "-f", str(args.overlap_ratio),
    "-d", "../temp_Res_" + args.sample
], check=True)

# 最终的文件复制操作
final_output_file = "../temp_Res_" + args.sample + "/" + args.sample + ".final_output." + str(args.threshold) + ".csv"
subprocess.run([
    "cp", final_output_file, os.path.join(args.output_dir, args.sample + ".Detected_eccDNA.csv")
], check=True)

# 清理中间文件
if os.path.exists(final_output_file):
    os.remove(final_output_file)

