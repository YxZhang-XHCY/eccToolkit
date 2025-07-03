import os
import logging
import pandas as pd
from typing import List

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(asctime)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def extract_sample_name(file_path: str) -> str:
    """Extract sample name from the given file path."""
    base_name = os.path.basename(file_path)
    sample_name = base_name.split('.DiGraph')[0]
    return sample_name

def round_float_columns(df: pd.DataFrame, columns: List[str], decimals: int = 2) -> pd.DataFrame:
    """Round specified float columns in a DataFrame to the given number of decimals."""
    for col in columns:
        if col in df.columns:
            df[col] = df[col].round(decimals)
    return df

def reorder_columns(df: pd.DataFrame, primary_cols: List[str]) -> pd.DataFrame:
    """Reorder DataFrame columns so that primary_cols appear in the front in the given order."""
    other_cols = [c for c in df.columns if c not in primary_cols]
    df = df[primary_cols + other_cols]
    return df

def read_fled_output(file_path: str) -> pd.DataFrame:
    """
    Read FLED OnesegJunction.out output file, add relevant columns, 
    filter data, and return a processed DataFrame.
    """
    columns = [
        'chrom', 'start', 'end', 'strand', 'label', 'Tag', 
        'Nfullpass', 'Nbreakpoint', 'L_Pvalue', 'R_Pvalue',
        'L_covRatio', 'R_covRatio', 'readID'
    ]

    # Check file existence
    if not os.path.isfile(file_path):
        logging.warning(f"File not found: {file_path}")
        return pd.DataFrame()

    df = pd.read_csv(file_path, sep='\t', header=None, names=columns)

    # Clean up readID and convert to comma-separated string
    df['readID'] = (
        df['readID']
        .str.strip('[]')
        .str.split(',')
        .apply(lambda x: ','.join(i.strip().strip("'\"") for i in x))
    )

    # Add 'sample' column
    df['sample'] = extract_sample_name(file_path)

    # Create name column
    df['name'] = df.apply(lambda row: f"{row['chrom']}-{row['start']}-{row['end']}", axis=1)

    # Calculate length
    df['Length'] = df['end'] - df['start'] + 1

    # Filter out rows where Nfullpass == 0
    df = df[df['Nfullpass'] > 0]

    # Log number of records
    if not df.empty:
        logging.info(f"{df['sample'].iloc[0]}: {len(df):,} records after filtering.")

    # Round float columns
    float_columns = ['L_Pvalue', 'R_Pvalue', 'L_covRatio', 'R_covRatio']
    df = round_float_columns(df, float_columns, 2)

    # Reorder columns so that name, sample, and Length are in front
    primary_cols = ['name', 'sample', 'Length']
    df = reorder_columns(df, primary_cols)

    return df

def process_circlemap_output(file_list: List[str]) -> pd.DataFrame:
    """
    Process Circle-Map output files, apply filters, and return a combined DataFrame.
    """
    columns = [
        'Chromosome', 'Start', 'End', 'Discordants', 'Split_reads', 
        'Circle_score', 'Mean_coverage', 'Standard_deviation',
        'Coverage_increase_start', 'Coverage_increase_end', 'Coverage_continuity'
    ]
    
    float_columns = [
        'Circle_score', 'Mean_coverage', 'Standard_deviation',
        'Coverage_increase_start', 'Coverage_increase_end', 'Coverage_continuity'
    ]
    
    all_data = []

    for file in file_list:
        if not os.path.isfile(file):
            logging.warning(f"File not found: {file}")
            continue

        df = pd.read_csv(file, sep='\t', header=None, names=columns)
        sample_name = os.path.basename(file).split('.')[0]

        # Create name column
        df['name'] = df.apply(lambda row: f"{row['Chromosome']}-{row['Start']}-{row['End']}", axis=1)

        # Round float columns
        df = round_float_columns(df, float_columns, 2)

        # Add sample name
        df['sample'] = sample_name

        # Calculate length
        df['Length'] = df['End'] - df['Start'] + 1

        # Filter based on Circle-Map criteria
        filtered_df = df[
            (df['Circle_score'] > 50) & 
            (df['Split_reads'] > 2) & 
            (df['Discordants'] > 2) & 
            (df['Coverage_increase_start'] > 0.33) & 
            (df['Coverage_increase_end'] > 0.33) & 
            (df['Length'] < 1e7)
        ]

        logging.info(f"{sample_name}: {len(filtered_df):,} records after filtering.")
        all_data.append(filtered_df)

    # Combine all data
    if all_data:
        final_df = pd.concat(all_data, ignore_index=True)
    else:
        logging.warning("No valid data was found among Circle-Map files.")
        return pd.DataFrame()

    # Reorder columns
    primary_cols = ['name', 'sample', 'Length']
    final_df = reorder_columns(final_df, primary_cols)

    return final_df

def create_simplified_merge(fled_df: pd.DataFrame, circlemap_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a simplified merged table from FLED and Circle-Map DataFrames.
    """
    # Extract needed columns from FLED
    fled_simple = fled_df.copy()
    fled_simple['Chromosome'] = fled_simple['chrom']
    fled_simple['Start'] = fled_simple['start']
    fled_simple['End'] = fled_simple['end']
    fled_simple['Source'] = 'FLED'
    fled_simple = fled_simple[['name', 'sample', 'Length', 'Chromosome', 'Start', 'End', 'Source']]

    # Extract needed columns from Circle-Map
    circlemap_simple = circlemap_df.copy()
    circlemap_simple['Source'] = 'Circle-Map'
    circlemap_simple = circlemap_simple[['name', 'sample', 'Length', 'Chromosome', 'Start', 'End', 'Source']]

    # Merge two DataFrames
    merged_df = pd.concat([fled_simple, circlemap_simple], ignore_index=True)
    merged_df.sort_values(['Chromosome', 'Start', 'End'], inplace=True)
    
    return merged_df

def main() -> None:
    # List of FLED files
    fled_files = [
        "3SEP_LR_25_1.DiGraph.OnesegJunction.out",
        "3SEP_LR_25_2.DiGraph.OnesegJunction.out",
        "3SEP_LR_25_3.DiGraph.OnesegJunction.out",
        "Circel-Seq_LR_25_1.DiGraph.OnesegJunction.out",
        "Circel-Seq_LR_25_2.DiGraph.OnesegJunction.out",
        "Circel-Seq_LR_25_3.DiGraph.OnesegJunction.out"
    ]
    
    # List of Circle-Map files
    circlemap_files = [
        "eccDNA_Circel-Seq_NGS_25_1_CM.bed",
        "eccDNA_Circel-Seq_NGS_25_2_CM.bed",
        "eccDNA_Circel-Seq_NGS_25_3_CM.bed"
    ]

    logging.info("Processing FLED outputs:")
    fled_dfs = []
    for file in fled_files:
        df = read_fled_output(file)
        if not df.empty:
            fled_dfs.append(df)

    if fled_dfs:
        merged_fled = pd.concat(fled_dfs, ignore_index=True)
    else:
        logging.warning("No valid data from FLED outputs. Exiting...")
        return

    logging.info("Processing Circle-Map outputs:")
    merged_circlemap = process_circlemap_output(circlemap_files)
    if merged_circlemap.empty:
        logging.warning("No valid data from Circle-Map outputs. Exiting...")
        return

    # Create simplified merged table
    simplified_merge = create_simplified_merge(merged_fled, merged_circlemap)

    # Save all results
    merged_fled.to_csv("FLED_Merged_Results.csv", index=False)
    merged_fled.to_csv("FLED_Merged_Results.tsv", sep='\t', index=False)
    merged_circlemap.to_csv("Circle_Map_Merged_Results.csv", index=False, float_format='%.2f')
    merged_circlemap.to_csv("Circle_Map_Merged_Results.tsv", sep='\t', index=False, float_format='%.2f')
    simplified_merge.to_csv("Combined_Basic_Results.csv", index=False)
    simplified_merge.to_csv("Combined_Basic_Results.tsv", sep='\t', index=False)

    logging.info("Processing complete!")
    logging.info(f"FLED results: {len(merged_fled):,} total records.")
    logging.info("Output files: FLED_Merged_Results.csv and FLED_Merged_Results.tsv.")
    logging.info(f"Circle-Map results: {len(merged_circlemap):,} total records.")
    logging.info("Output files: Circle_Map_Merged_Results.csv and Circle_Map_Merged_Results.tsv.")
    logging.info(f"Combined results: {len(simplified_merge):,} total records.")
    logging.info("Output files: Combined_Basic_Results.csv and Combined_Basic_Results.tsv.")

if __name__ == "__main__":
    main()
