import argparse
import logging
import os
import glob
import pandas as pd
from typing import List

# --- Configuration / 配置 ---
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(asctime)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# --- Helper Functions / 辅助函数 ---

def extract_sample_name(file_path: str, suffix: str) -> str:
    """
    Extracts a sample name from a file path by removing a known suffix.
    中文：通过移除已知后缀，从文件路径中提取样本名。

    Args:
        file_path (str): The full path to the file. / 文件的完整路径。
        suffix (str): The suffix to remove from the filename. / 需要从文件名中移除的后缀。

    Returns:
        str: The extracted sample name. / 提取出的样本名。
    """
    base_name = os.path.basename(file_path)
    if base_name.endswith(suffix):
        return base_name[:-len(suffix)]
    # Fallback for unexpected naming / 针对未能预期的命名方式，提供备用方案
    return base_name.split('.')[0]

def process_fled_output(file_path: str) -> pd.DataFrame:
    """
    Reads and processes a single FLED output file (OnesegJunction.out).
    中文：读取并处理单个 FLED 输出文件 (OnesegJunction.out)。

    Args:
        file_path (str): Path to the FLED output file. / FLED 输出文件的路径。

    Returns:
        pd.DataFrame: A processed DataFrame, or an empty one if an error occurs. / 处理后的 DataFrame，如果出错则返回空的 DataFrame。
    """
    logging.info(f"Processing FLED file: {file_path}")
    try:
        # Define column names for FLED output / 定义 FLED 输出文件的列名
        columns = [
            'chrom', 'start', 'end', 'strand', 'label', 'Tag', 
            'Nfullpass', 'Nbreakpoint', 'L_Pvalue', 'R_Pvalue',
            'L_covRatio', 'R_covRatio', 'readID'
        ]
        df = pd.read_csv(file_path, sep='\t', header=None, names=columns)

        # Vectorized operations for new columns for better performance / 使用向量化操作创建新列以获得更好性能
        df['name'] = df['chrom'].astype(str) + '-' + df['start'].astype(str) + '-' + df['end'].astype(str)
        df['Length'] = df['end'] - df['start'] + 1
        df['sample'] = extract_sample_name(file_path, '.DiGraph.OnesegJunction.out')

        # Filter out rows with no full pass reads / 筛选掉 Nfullpass 为 0 的记录
        df_filtered = df[df['Nfullpass'] > 0].copy() # Use .copy() to avoid SettingWithCopyWarning / 使用.copy()以避免SettingWithCopyWarning

        # Clean up readID column / 清理 readID 列
        df_filtered['readID'] = (
            df_filtered['readID']
            .str.strip('[]')
            .str.replace(r"['\"]", "", regex=True)
            .str.replace(r",\s*", ",", regex=True)
        )
        
        # Round float columns / 将浮点数相关的列四舍五入
        float_cols = ['L_Pvalue', 'R_Pvalue', 'L_covRatio', 'R_covRatio']
        for col in float_cols:
            df_filtered[col] = df_filtered[col].round(2)
        
        logging.info(f"  Found {len(df_filtered):,} valid records in sample {df_filtered['sample'].iloc[0]}")
        return df_filtered

    except Exception as e:
        logging.error(f"Could not process file {file_path}: {e}")
        return pd.DataFrame()

def process_circlemap_output(file_path: str) -> pd.DataFrame:
    """
    Reads and processes a single Circle-Map output file (.bed).
    中文：读取并处理单个 Circle-Map 输出文件 (.bed)。

    Args:
        file_path (str): Path to the Circle-Map output file. / Circle-Map 输出文件的路径。

    Returns:
        pd.DataFrame: A processed and filtered DataFrame. / 处理并筛选后的 DataFrame。
    """
    logging.info(f"Processing Circle-Map file: {file_path}")
    try:
        # Define column names for Circle-Map output / 定义 Circle-Map 输出文件的列名
        columns = [
            'Chromosome', 'Start', 'End', 'Discordants', 'Split_reads', 
            'Circle_score', 'Mean_coverage', 'Standard_deviation',
            'Coverage_increase_start', 'Coverage_increase_end', 'Coverage_continuity'
        ]
        df = pd.read_csv(file_path, sep='\t', header=None, names=columns)

        # Vectorized operations for new columns / 使用向量化操作创建新列
        df['name'] = df['Chromosome'].astype(str) + '-' + df['Start'].astype(str) + '-' + df['End'].astype(str)
        df['Length'] = df['End'] - df['Start'] + 1
        df['sample'] = extract_sample_name(file_path, '_CM.bed')
        
        # Filter based on established Circle-Map criteria / 根据既定的 Circle-Map 标准进行筛选
        df_filtered = df[
            (df['Circle_score'] > 50) & 
            (df['Split_reads'] > 2) & 
            (df['Discordants'] > 2) & 
            (df['Coverage_increase_start'] > 0.33) & 
            (df['Coverage_increase_end'] > 0.33) & 
            (df['Length'] < 1e7)
        ].copy() # Use .copy() to avoid SettingWithCopyWarning / 使用.copy()以避免SettingWithCopyWarning

        # Round float columns / 将浮点数相关的列四舍五入
        float_cols = [
            'Circle_score', 'Mean_coverage', 'Standard_deviation',
            'Coverage_increase_start', 'Coverage_increase_end', 'Coverage_continuity'
        ]
        for col in float_cols:
            df_filtered[col] = df_filtered[col].round(2)

        logging.info(f"  Found {len(df_filtered):,} valid records in sample {df_filtered['sample'].iloc[0]}")
        return df_filtered

    except Exception as e:
        logging.error(f"Could not process file {file_path}: {e}")
        return pd.DataFrame()

def main() -> None:
    """
    Main function to drive the analysis.
    中文: 驱动分析的主函数。
    """
    parser = argparse.ArgumentParser(
        description="Process and combine eccDNA results from FLED and Circle-Map. / 处理并整合来自 FLED 和 Circle-Map 的 eccDNA 结果。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--fled-dir", required=True, help="Directory containing FLED output files. / 包含 FLED 输出文件的目录。")
    parser.add_argument("--circlemap-dir", required=True, help="Directory containing Circle-Map output files. / 包含 Circle-Map 输出文件的目录。")
    parser.add_argument("--fled-pattern", default="*.DiGraph.OnesegJunction.out", help="Filename pattern for FLED files. / 用于查找 FLED 文件的文件名模式。")
    parser.add_argument("--circlemap-pattern", default="*_CM.bed", help="Filename pattern for Circle-Map files. / 用于查找 Circle-Map 文件的文件名模式。")
    parser.add_argument("--output-dir", required=True, help="Directory to save the final result files. / 用于保存最终结果文件的目录。")
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist / 如果输出目录不存在则创建
    os.makedirs(args.output_dir, exist_ok=True)
    
    # --- Process FLED files / 处理 FLED 文件 ---
    fled_files = glob.glob(os.path.join(args.fled_dir, args.fled_pattern)) # Find all files matching the FLED pattern / 查找所有匹配 FLED 模式的文件
    if not fled_files:
        logging.warning(f"No FLED files found in '{args.fled_dir}' with pattern '{args.fled_pattern}'")
        fled_df = pd.DataFrame()
    else:
        fled_df = pd.concat([process_fled_output(f) for f in fled_files], ignore_index=True)

    # --- Process Circle-Map files / 处理 Circle-Map 文件 ---
    circlemap_files = glob.glob(os.path.join(args.circlemap_dir, args.circlemap_pattern)) # Find all files matching the Circle-Map pattern / 查找所有匹配 Circle-Map 模式的文件
    if not circlemap_files:
        logging.warning(f"No Circle-Map files found in '{args.circlemap_dir}' with pattern '{args.circlemap_pattern}'")
        circlemap_df = pd.DataFrame()
    else:
        circlemap_df = pd.concat([process_circlemap_output(f) for f in circlemap_files], ignore_index=True)

    # Exit if no data could be processed / 如果没有任何数据被处理则退出
    if fled_df.empty and circlemap_df.empty:
        logging.error("No valid data could be processed from any input files. Exiting.")
        return
        
    # --- Create simplified merged table / 创建简化的合并表格 ---
    fled_simple = pd.DataFrame()
    if not fled_df.empty:
        fled_simple = fled_df[['name', 'sample', 'Length', 'chrom', 'start', 'end']].rename(columns={'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'})
        fled_simple['Source'] = 'FLED'
    
    circlemap_simple = pd.DataFrame()
    if not circlemap_df.empty:
        circlemap_simple = circlemap_df[['name', 'sample', 'Length', 'Chromosome', 'Start', 'End']]
        circlemap_simple['Source'] = 'Circle-Map'
    
    combined_df = pd.concat([fled_simple, circlemap_simple], ignore_index=True).sort_values(['Chromosome', 'Start'])
    
    # --- Save results / 保存结果 ---
    # Define primary columns for reordering / 定义重排的列顺序
    primary_cols = ['name', 'sample', 'Length', 'Chromosome', 'Start', 'End', 'Source']
    
    if not fled_df.empty:
        other_fled_cols = [c for c in fled_df.columns if c not in primary_cols[:3]]
        fled_df_final = fled_df[primary_cols[:3] + other_fled_cols]
        fled_df_final.to_csv(os.path.join(args.output_dir, "FLED_filtered_results.csv"), index=False)
        logging.info(f"FLED results saved. Total records: {len(fled_df_final):,}")

    if not circlemap_df.empty:
        other_circlemap_cols = [c for c in circlemap_df.columns if c not in primary_cols[:3]]
        circlemap_df_final = circlemap_df[primary_cols[:3] + other_circlemap_cols]
        circlemap_df_final.to_csv(os.path.join(args.output_dir, "CircleMap_filtered_results.csv"), index=False, float_format='%.2f')
        logging.info(f"Circle-Map results saved. Total records: {len(circlemap_df_final):,}")

    if not combined_df.empty:
        combined_df.to_csv(os.path.join(args.output_dir, "FLED_and_CircleMap_combined.csv"), index=False)
        logging.info(f"Combined results saved. Total records: {len(combined_df):,}")

    logging.info(f"✔ Processing complete! All output files are saved in: {args.output_dir}")

if __name__ == "__main__":
    main()