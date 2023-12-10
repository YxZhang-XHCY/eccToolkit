import pandas as pd
import sys

def process_blast_results(csv_input_file, csv_output_file):
    # 读取CSV文件
    df = pd.read_csv(csv_input_file)

    # 确保相关列是数值类型
    df['identity'] = pd.to_numeric(df['identity'], errors='coerce')
    df['q_end'] = pd.to_numeric(df['q_end'], errors='coerce')
    df['q_start'] = pd.to_numeric(df['q_start'], errors='coerce')
    df['s_end'] = pd.to_numeric(df['s_end'], errors='coerce')
    df['s_start'] = pd.to_numeric(df['s_start'], errors='coerce')

    # 删除不满足条件的行，并创建一个副本以避免SettingWithCopyWarning
    df_filtered = df[(df['identity'] > 90) & (df['q_end'] - df['q_start'] > 0) & (df['s_end'] - df['s_start'] > 0)].copy()

    # 提取起始位置的数值并进行转换
    df_filtered['q_real_start'] = df_filtered['q_start'] + pd.to_numeric(df_filtered['query_id'].str.extract(r':(\d+)-')[0], errors='coerce')
    df_filtered['q_real_end'] = df_filtered['q_end'] + pd.to_numeric(df_filtered['query_id'].str.extract(r':(\d+)-')[0], errors='coerce')
    df_filtered['s_real_start'] = df_filtered['s_start'] + pd.to_numeric(df_filtered['subject_id'].str.extract(r':(\d+)-')[0], errors='coerce')
    df_filtered['s_real_end'] = df_filtered['s_end'] + pd.to_numeric(df_filtered['subject_id'].str.extract(r':(\d+)-')[0], errors='coerce')

    # 分割eccDNA_name为三列
    ecc_cols = df_filtered['eccDNA_name'].str.split('_', expand=True)
    df_filtered['ecc_chr'] = ecc_cols[0]
    df_filtered['ecc_start'] = pd.to_numeric(ecc_cols[1], errors='coerce')
    df_filtered['ecc_end'] = pd.to_numeric(ecc_cols[2], errors='coerce')

    # 筛选接近ecc_start或ecc_end的行
    close_to_start = (df_filtered['q_real_start'] - df_filtered['ecc_start']).abs() < 2
    close_to_end = (df_filtered['s_real_end'] - df_filtered['ecc_end']).abs() < 2
    df_final = df_filtered[close_to_start | close_to_end]

    # 保存处理后的数据
    df_final.to_csv(csv_output_file, index=False)

def main():
    if len(sys.argv) != 3:
        print("Usage: python process_direct_repeat_results.py <input_file.csv> <output_file.csv>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    process_blast_results(input_file, output_file)

if __name__ == "__main__":
    main()
