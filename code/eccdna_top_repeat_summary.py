#!/usr/bin/env python3
"""
eccDNA Top Repeat Summary Analysis Script
Analyzes eccDNA data to summarize top repeat motifs
Outputs separate tables for each Species-Class combination
"""

import pandas as pd
import argparse
import sys
import os
from collections import defaultdict

def get_species(sample):
    """Determine species based on sample name"""
    sample_upper = sample.upper()
    
    # Remove trailing numbers and underscores for species identification
    if sample_upper.startswith('HELA'):
        return 'HeLa'
    elif sample_upper.startswith('ZM'):
        return 'ZM'
    elif sample_upper.startswith('U87'):
        return 'U87'
    elif sample_upper.startswith('AT'):
        return 'AT'
    # Add more species mappings as needed
    return sample

def analyze_motif_distribution(df):
    """Analyze motif distribution for a given dataframe"""
    if len(df) == 0:
        return pd.DataFrame()
    
    # Count occurrences of each motif
    motif_counts = df['motif'].value_counts()
    total_count = len(df)
    
    # Calculate percentages
    motif_percentages = (motif_counts / total_count * 100).round(2)
    
    # Sort by percentage (descending)
    motif_percentages = motif_percentages.sort_values(ascending=False)
    
    # Find top motifs that sum up to >= 70%
    cumulative_percentage = 0
    top_motifs = []
    other_count = 0
    other_percentage = 0
    
    for idx, (motif_name, percentage) in enumerate(motif_percentages.items()):
        if idx < 10:  # Maximum 10 named motifs
            cumulative_percentage += percentage
            count = motif_counts[motif_name]
            top_motifs.append({
                'Motif': motif_name,
                'Count': count,
                'Percentage': percentage,
                'Cumulative_Percentage': cumulative_percentage
            })
            
            if cumulative_percentage >= 70:
                # Calculate "other" for remaining motifs
                other_count = motif_counts.iloc[idx+1:].sum()
                other_percentage = motif_percentages.iloc[idx+1:].sum()
                break
        else:
            other_count += motif_counts[motif_name]
            other_percentage += percentage
    
    # If we've gone through top 10 and still < 70%, add remaining as "other"
    if len(top_motifs) == 10 and cumulative_percentage < 70:
        other_count = motif_counts.iloc[10:].sum()
        other_percentage = motif_percentages.iloc[10:].sum()
    
    # Add "other" if there are any
    if other_percentage > 0:
        cumulative_percentage += other_percentage
        top_motifs.append({
            'Motif': 'other',
            'Count': int(other_count),
            'Percentage': round(other_percentage, 2),
            'Cumulative_Percentage': round(cumulative_percentage, 2)
        })
    
    return pd.DataFrame(top_motifs)

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Analyze eccDNA repeat motifs by species and class')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file')
    parser.add_argument('-o', '--output', required=True, help='Output directory for result files')
    
    args = parser.parse_args()
    
    try:
        # Read the input CSV file
        print(f"Reading input file: {args.input}")
        df = pd.read_csv(args.input)
        print(f"Total rows: {len(df)}")
        
        # Create output directory if it doesn't exist
        output_dir = args.output
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
        
        # Add Species column
        df['Species'] = df['Sample'].apply(get_species)
        
        # Get unique species and classes
        species_list = sorted(df['Species'].unique())
        classes_list = sorted(df['Class'].unique())
        
        print(f"Species found: {', '.join(species_list)}")
        print(f"Classes found: {', '.join(classes_list)}")
        
        # Create a summary of all outputs
        summary_info = []
        
        # Process each species-class combination
        for species in species_list:
            species_df = df[df['Species'] == species]
            
            for class_name in classes_list:
                species_class_df = species_df[species_df['Class'] == class_name]
                
                if len(species_class_df) == 0:
                    continue
                
                print(f"\n{'='*60}")
                print(f"Processing: {species} - {class_name} ({len(species_class_df)} rows)")
                print(f"{'='*60}")
                
                # Analyze motif distribution
                result_df = analyze_motif_distribution(species_class_df)
                
                if len(result_df) > 0:
                    # Create output filename
                    output_file = os.path.join(output_dir, f"{species}_{class_name}_motif_summary.csv")
                    
                    # Save to file
                    result_df.to_csv(output_file, index=False)
                    print(f"Saved to: {output_file}")
                    
                    # Print the table
                    print(f"\n{species} - {class_name} Motif Distribution:")
                    print(result_df.to_string(index=False))
                    
                    # Add to summary
                    summary_info.append({
                        'Species': species,
                        'Class': class_name,
                        'Total_Rows': len(species_class_df),
                        'Output_File': os.path.basename(output_file)
                    })
        
        # Create a master summary file
        summary_df = pd.DataFrame(summary_info)
        summary_file = os.path.join(output_dir, "analysis_summary.csv")
        summary_df.to_csv(summary_file, index=False)
        
        print(f"\n{'='*60}")
        print("Analysis complete!")
        print(f"{'='*60}")
        print(f"\nSummary of all outputs:")
        print(summary_df.to_string(index=False))
        print(f"\nAll output files saved to: {output_dir}")
        print(f"Master summary saved to: {summary_file}")
        
    except FileNotFoundError:
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
