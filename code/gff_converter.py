#!/usr/bin/env python3
"""
GFF文件批量转换为CSV格式的脚本
将多个GFF文件合并为一个CSV文件，并添加样本和类别信息
"""

import pandas as pd
import os
import glob
import re
from pathlib import Path

def parse_gff_line(line):
    """解析GFF格式的单行数据"""
    if line.startswith('#') or line.strip() == '':
        return None
    
    parts = line.strip().split('\t')
    if len(parts) != 9:
        return None
    
    return {
        'seqname': parts[0],
        'source': parts[1],
        'feature': parts[2],
        'start': int(parts[3]),
        'end': int(parts[4]),
        'score': parts[5],
        'strand': parts[6],
        'frame': parts[7],
        'attribute': parts[8]
    }

def extract_sample_and_class(filename):
    """从文件名中提取样本名和类别"""
    basename = os.path.basename(filename)
    # 移除.gff扩展名
    name_without_ext = basename.replace('.gff', '')
    
    # 使用正则表达式匹配模式：样本名.类别
    match = re.match(r'^(.+)\.(Mecc|Uecc)$', name_without_ext)
    if match:
        sample = match.group(1)
        class_type = match.group(2)
        return sample, class_type
    else:
        # 如果匹配失败，尝试其他方式
        parts = name_without_ext.split('.')
        if len(parts) >= 2:
            sample = '.'.join(parts[:-1])
            class_type = parts[-1]
            return sample, class_type
        else:
            return name_without_ext, 'Unknown'

def process_gff_files(input_dir='.', output_csv='merged_gff_data.csv'):
    """处理所有GFF文件并合并为CSV"""
    all_data = []
    
    # 查找所有GFF文件
    gff_files = glob.glob(os.path.join(input_dir, '*.gff'))
    
    if not gff_files:
        print(f"在目录 {input_dir} 中未找到GFF文件")
        return
    
    print(f"找到 {len(gff_files)} 个GFF文件")
    
    for gff_file in gff_files:
        print(f"处理文件: {gff_file}")
        
        # 从文件名提取样本和类别信息
        sample, class_type = extract_sample_and_class(gff_file)
        
        try:
            with open(gff_file, 'r', encoding='utf-8') as f:
                file_data = []
                for line_num, line in enumerate(f, 1):
                    parsed = parse_gff_line(line)
                    if parsed:
                        # 添加样本和类别信息
                        parsed['Sample'] = sample
                        parsed['Class'] = class_type
                        parsed['source_file'] = os.path.basename(gff_file)
                        file_data.append(parsed)
                
                all_data.extend(file_data)
                print(f"  - 从 {gff_file} 读取了 {len(file_data)} 条记录")
                
        except Exception as e:
            print(f"处理文件 {gff_file} 时出错: {e}")
            continue
    
    if not all_data:
        print("没有有效的数据可以处理")
        return
    
    # 转换为DataFrame
    df = pd.DataFrame(all_data)
    
    # 重新排列列的顺序，将Sample和Class放在前面
    columns_order = ['Sample', 'Class', 'seqname', 'source', 'feature', 
                    'start', 'end', 'score', 'strand', 'frame', 'attribute', 'source_file']
    df = df[columns_order]
    
    # 保存为CSV
    df.to_csv(output_csv, index=False, encoding='utf-8')
    print(f"\n合并完成！总共处理了 {len(df)} 条记录")
    print(f"结果已保存到: {output_csv}")
    
    return df

def generate_report(df):
    """生成简单的数据报告"""
    print("\n" + "="*60)
    print("数据处理报告")
    print("="*60)
    
    print(f"总记录数: {len(df):,}")
    print(f"总文件数: {df['source_file'].nunique()}")
    
    print("\n按样本统计:")
    sample_counts = df['Sample'].value_counts()
    for sample, count in sample_counts.items():
        print(f"  {sample}: {count:,} 条记录")
    
    print("\n按类别统计:")
    class_counts = df['Class'].value_counts()
    for class_type, count in class_counts.items():
        print(f"  {class_type}: {count:,} 条记录")
    
    print("\n按样本和类别统计:")
    sample_class_counts = df.groupby(['Sample', 'Class']).size().reset_index(name='count')
    for _, row in sample_class_counts.iterrows():
        print(f"  {row['Sample']}.{row['Class']}: {row['count']:,} 条记录")
    
    print("\n特征类型统计:")
    feature_counts = df['feature'].value_counts()
    for feature, count in feature_counts.head(10).items():
        print(f"  {feature}: {count:,} 条记录")
    if len(feature_counts) > 10:
        print(f"  ... 还有 {len(feature_counts) - 10} 种其他特征类型")
    
    print("\n数据源统计:")
    source_counts = df['source'].value_counts()
    for source, count in source_counts.items():
        print(f"  {source}: {count:,} 条记录")
    
    print("\n链方向统计:")
    strand_counts = df['strand'].value_counts()
    for strand, count in strand_counts.items():
        print(f"  {strand}: {count:,} 条记录")
    
    print("\n数据质量检查:")
    print(f"  缺失值检查:")
    missing_values = df.isnull().sum()
    for col, missing in missing_values.items():
        if missing > 0:
            print(f"    {col}: {missing} 个缺失值")
    if missing_values.sum() == 0:
        print("    无缺失值")
    
    print(f"\n  坐标范围:")
    print(f"    起始位置范围: {df['start'].min():,} - {df['start'].max():,}")
    print(f"    结束位置范围: {df['end'].min():,} - {df['end'].max():,}")
    print(f"    平均序列长度: {(df['end'] - df['start'] + 1).mean():.1f}")
    
    print("\n处理完成！")

def main():
    """主函数"""
    print("GFF文件批量转换工具")
    print("="*40)
    
    # 设置输入目录和输出文件
    input_directory = '.'  # 当前目录，可以修改为你的GFF文件所在目录
    output_file = 'merged_gff_data.csv'
    
    # 处理GFF文件
    df = process_gff_files(input_directory, output_file)
    
    if df is not None and not df.empty:
        # 生成报告
        generate_report(df)
        
        # 显示前几行数据作为预览
        print(f"\n数据预览 (前5行):")
        print(df.head().to_string(index=False))
    else:
        print("没有数据可以处理")

if __name__ == "__main__":
    main()
