#!/usr/bin/env python3
"""
TSV文件处理脚本
功能：
1. 删除第一行（标题行）
2. 筛选出Chr1和Chr2都为chr1的行
3. 保存处理后的结果

使用方法：
python process_tsv.py input_file.tsv output_file.tsv
"""

import sys
import os

def process_tsv_file(input_file, output_file=None):
    """
    处理TSV文件
    
    Args:
        input_file (str): 输入文件路径
        output_file (str): 输出文件路径，如果为None则自动生成
    
    Returns:
        tuple: (处理的行数, 筛选出的行数)
    """
    
    # 检查输入文件是否存在
    if not os.path.exists(input_file):
        print(f"错误：文件 '{input_file}' 不存在")
        return None, None
    
    # 如果没有指定输出文件，自动生成
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_filtered.tsv"
    
    try:
        # 读取文件
        with open(input_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        print(f"读取文件: {input_file}")
        print(f"原始行数: {len(lines)}")
        
        # 删除第一行（标题行）
        if len(lines) <= 1:
            print("警告：文件只有标题行或为空")
            return 0, 0
        
        data_lines = lines[1:]  # 跳过第一行
        print(f"删除标题行后: {len(data_lines)} 行")
        
        # 筛选出Chr1和Chr2都为chr1的行
        filtered_lines = []
        for line_num, line in enumerate(data_lines, start=2):  # 从第2行开始计数
            line = line.strip()
            if not line:  # 跳过空行
                continue
                
            columns = line.split('\t')
            if len(columns) >= 4:  # 确保有足够的列
                chr1 = columns[0].strip()
                chr2 = columns[3].strip()
                
                if chr1 == 'chr1' and chr2 == 'chr1':
                    filtered_lines.append(line)
                    print(f"找到匹配行 {line_num}: {chr1} -> {chr2}")
        
        print(f"\n筛选结果: {len(filtered_lines)} 行符合条件")
        
        # 保存处理后的结果
        if filtered_lines:
            with open(output_file, 'w', encoding='utf-8') as f:
                for line in filtered_lines:
                    f.write(line + '\n')
            print(f"结果已保存到: {output_file}")
        else:
            print("没有符合条件的数据，未创建输出文件")
        
        # 显示前几行结果
        if filtered_lines:
            print("\n前5行结果预览:")
            for i, line in enumerate(filtered_lines[:5], 1):
                print(f"{i}: {line}")
            if len(filtered_lines) > 5:
                print(f"... 还有 {len(filtered_lines) - 5} 行")
        
        return len(data_lines), len(filtered_lines)
        
    except Exception as e:
        print(f"处理文件时出错: {e}")
        return None, None

def main():
    """主函数"""
    if len(sys.argv) < 2:
        print("使用方法:")
        print("  python process_tsv.py <输入文件> [输出文件]")
        print("  python process_tsv.py HeLa_rep1.Mecc.Links.RandomColor.tsv")
        print("  python process_tsv.py HeLa_rep1.Mecc.Links.RandomColor.tsv filtered_output.tsv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    print("=" * 50)
    print("TSV文件处理脚本")
    print("=" * 50)
    
    processed_lines, filtered_lines = process_tsv_file(input_file, output_file)
    
    if processed_lines is not None:
        print("\n" + "=" * 50)
        print("处理完成！")
        print(f"处理的数据行数: {processed_lines}")
        print(f"筛选出的行数: {filtered_lines}")
        if filtered_lines > 0:
            print(f"筛选比例: {filtered_lines/processed_lines*100:.1f}%")
        print("=" * 50)
    else:
        print("处理失败！")
        sys.exit(1)

if __name__ == "__main__":
    main()
