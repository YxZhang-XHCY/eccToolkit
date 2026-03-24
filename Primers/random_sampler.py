#!/usr/bin/env python3
import pandas as pd
import argparse
import sys
from pathlib import Path

def sample_csv(input_file: str, output_file: str, n_samples: int, seed: int = None) -> None:
    """
    从CSV文件中随机抽取指定数量的行，并筛选特定列
    
    参数:
    input_file (str): 输入CSV文件路径
    output_file (str): 输出CSV文件路径
    n_samples (int): 需要抽取的样本数量
    seed (int): 随机数种子，用于复现结果
    """
    try:
        # 读取CSV文件，使用逗号作为分隔符
        df = pd.read_csv(input_file, delimiter=',')
        
        # 检查必要的列是否存在
        required_columns = {'eName', 'eChr', 'eStart', 'eEnd', 'eLength'}
        missing_columns = required_columns - set(df.columns)
        if missing_columns:
            raise ValueError(f"输入文件缺少必要的列: {missing_columns}")
        
        # 检查样本数量是否合理
        if n_samples > len(df):
            print(f"警告: 请求的样本数({n_samples})大于文件中的总行数({len(df)})。将返回所有行。")
            n_samples = len(df)
        
        # 随机抽样
        sampled_df = df.sample(n=n_samples, random_state=seed)
        
        # 创建新的eName-eLength列
        sampled_df['name_length'] = sampled_df['eName'] + '-' + sampled_df['eLength'].astype(str)
        
        # 选择需要的列并重新排序
        result_df = sampled_df[['eChr', 'eStart', 'eEnd', 'name_length']]
        
        # 保存结果，不包含表头和索引，使用制表符分隔
        result_df.to_csv(output_file, index=False, header=False, sep='\t')
        print(f"已成功将{n_samples}个随机样本保存到 {output_file}")
        
    except Exception as e:
        print(f"错误: {str(e)}", file=sys.stderr)
        sys.exit(1)

def main():
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='从CSV文件中随机抽取指定数量的行并筛选特定列')
    
    # 添加命令行参数
    parser.add_argument('-i', '--input', required=True,
                      help='输入文件路径（逗号分隔的CSV文件）')
    parser.add_argument('-o', '--output', required=True,
                      help='输出文件路径（将生成制表符分隔的文件）')
    parser.add_argument('-n', '--number', type=int, required=True,
                      help='需要抽取的样本数量')
    parser.add_argument('-s', '--seed', type=int, default=None,
                      help='随机数种子（可选）')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not Path(args.input).exists():
        print(f"错误: 输入文件 '{args.input}' 不存在", file=sys.stderr)
        sys.exit(1)
    
    # 执行抽样
    sample_csv(args.input, args.output, args.number, args.seed)

if __name__ == '__main__':
    main()
