import pandas as pd
import numpy as np

def random_sample_data(input_file, output_file, random_seed=42):
    # 读取CSV文件
    df = pd.read_csv(input_file, sep='\t')
    
    # 按Sample分组
    sample_groups = df['Sample'].unique()
    
    final_samples = []
    
    # 对每个Sample组进行处理
    for sample in sample_groups:
        sample_df = df[df['Sample'] == sample]
        
        # 对每个Repeat_Class进行抽样
        for repeat_class in sample_df['Repeat_Class'].unique():
            class_df = sample_df[sample_df['Repeat_Class'] == repeat_class]
            
            # 如果数据量少于100，取全部数据；否则随机抽取100个
            if len(class_df) <= 200:
                sampled = class_df
            else:
                sampled = class_df.sample(n=200, random_state=random_seed)
            
            final_samples.append(sampled)
    
    # 合并所有抽样结果
    result = pd.concat(final_samples, ignore_index=True)
    
    # 保存结果到新的CSV文件
    result.to_csv(output_file, index=False)
    
    # 打印统计信息
    print(f"Total samples selected: {len(result)}")
    print("\nSamples per group:")
    for sample in sample_groups:
        sample_df = result[result['Sample'] == sample]
        print(f"\n{sample}:")
        print(sample_df['Repeat_Class'].value_counts())

# 使用示例
if __name__ == "__main__":
    input_file = "input.csv"  # 替换为你的输入文件名
    output_file = "sampled_output.csv"  # 替换为期望的输出文件名
    random_sample_data(input_file, output_file)
