import os
import glob
import pandas as pd
import subprocess
import argparse
import logging
from pathlib import Path

# 初始化日志记录
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# 设置命令行参数
parser = argparse.ArgumentParser(description='Process eccDNA data.')
parser.add_argument('-i', '--input', required=True, help='Input BED file')
parser.add_argument('-g', '--genome', required=True, help='Genome file')
parser.add_argument('-o', '--output', default=None, help='Output file name (optional)')
parser.add_argument('--extend_outer', type=int, default=100, help="Number of base pairs to extend outward from the ends (default: 100)")
parser.add_argument('--extend_inner', type=int, default=50, help="Number of base pairs to extend inward from the ends (default: 50)")
args = parser.parse_args()

# 验证输入文件是否存在
input_path = Path(args.input)
if not input_path.exists():
    logging.error(f'Input file {args.input} does not exist.')
    exit(1)

# 读取、处理BED文件，并覆盖原文件
bed_data = []
with open(args.input, 'r') as file:
    for line in file:
        chr, start, end, name = line.strip().split()
        start = int(start)  # 转换为整数
        end = int(end)      # 转换为整数
        bed_data.append(f"{chr}\t{start}\t{end}\t{name}\n")

with open(args.input, 'w') as file:
    file.writelines(bed_data)

# 获取输出文件的完整路径
output_file = args.output if args.output else input_path.parent / 'combined_results.csv'

# 构建指向 bin 文件夹的路径
bin_path = Path(__file__).parent / 'bin'

# 执行外部脚本并检查结果
def run_script(script_name, args):
    script_path = bin_path / script_name
    result = subprocess.run(['python', script_path] + args)
    if result.returncode != 0:
        logging.error(f"Script {script_path} failed with return code {result.returncode}")
        exit(1)

# 准备额外参数
extra_args = []
if args.extend_outer:
    extra_args += ['--extend_outer', str(args.extend_outer)]
if args.extend_inner:
    extra_args += ['--extend_inner', str(args.extend_inner)]

# 执行第一个脚本
run_script('eccDNA_Blast_Analyzer.py', [str(args.genome), str(args.input), 'intermediate_direct_results.csv', 'intermediate_inverted_results.csv'] + extra_args)

# 执行其他Python脚本
run_script('process_direct_repeat_results.py', ['intermediate_direct_results.csv', 'processed_direct_results.csv'])
run_script('process_inverted_repeat_results.py', ['intermediate_inverted_results.csv', 'processed_inverted_results.csv'])

# 读取并处理CSV文件
try:
    df_direct = pd.read_csv('processed_direct_results.csv')
    df_inverted = pd.read_csv('processed_inverted_results.csv')
    df_direct['Repeat_Class'] = 'Direct-Repeat'
    df_inverted['Repeat_Class'] = 'Inverted-Repeat'
    df_direct['q_real_chr'] = df_direct['s_real_chr'] = df_direct['ecc_chr']
    df_inverted['q_real_chr'] = df_inverted['s_real_chr'] = df_inverted['ecc_chr']
    df_combined = pd.concat([df_direct, df_inverted])
    columns_to_keep = ['eccDNA_name', 'identity', 'alignment_length', 'q_real_chr', 'q_real_start', 'q_real_end', 's_real_chr', 's_real_start', 's_real_end', 'Repeat_Class']
    df_combined = df_combined[columns_to_keep]
    df_combined.to_csv(output_file, index=False)
    logging.info(f"Combined results saved to {output_file}")
except Exception as e:
    logging.error(f"Error processing CSV files: {e}")
    exit(1)

# 删除过程文件
for file in ['intermediate_direct_results.csv', 'processed_direct_results.csv', 'intermediate_inverted_results.csv', 'processed_inverted_results.csv']:
    if os.path.exists(file):
        os.remove(file)
        logging.info(f"Deleted process file: {file}")

