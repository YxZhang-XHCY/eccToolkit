#!/usr/bin/env python3
"""
eccDNA Hotspot Detection and Classification Script
Author: Your Name
Date: 2024
Description: Multi-scale detection of eccDNA hotspots with reproducibility classification
"""

import pandas as pd
import numpy as np
from scipy import stats
import argparse
from pathlib import Path
import warnings
import sys
warnings.filterwarnings('ignore')

class eccDNAHotspotDetector:
    def __init__(self, fold_threshold=3, fdr_threshold=0.01, overlap_threshold=0.5):
        """
        Initialize the detector with parameters
        
        Parameters:
        -----------
        fold_threshold : float
            Minimum fold enrichment over background (default: 3)
        fdr_threshold : float
            FDR threshold for significance (default: 0.01)
        overlap_threshold : float
            Minimum overlap fraction for peak merging (default: 0.5)
        """
        self.fold_threshold = fold_threshold
        self.fdr_threshold = fdr_threshold
        self.overlap_threshold = overlap_threshold
        
    def load_data(self, file_paths):
        """Load multi-scale data from CSV files"""
        data = {}
        for scale, path in file_paths.items():
            df = pd.read_csv(path)
            # Get sample columns (excluding coordinate columns)
            sample_cols = [col for col in df.columns if col not in ['Chrom', 'WindowStart', 'WindowEnd']]
            data[scale] = {
                'df': df,
                'samples': sample_cols
            }
        return data
    
    def calculate_background(self, values, method='median'):
        """Calculate background signal"""
        if method == 'median':
            # Use median of non-zero values as background
            non_zero = values[values > 0]
            if len(non_zero) > 0:
                return np.median(non_zero)
            else:
                return 0.1  # Pseudocount
        elif method == 'percentile':
            return np.percentile(values[values > 0], 25) if np.any(values > 0) else 0.1
            
    def detect_peaks_single_sample(self, df, sample_col, scale):
        """Detect peaks in a single sample"""
        # Calculate background
        background = self.calculate_background(df[sample_col].values)
        
        # Calculate fold enrichment
        df['fold_enrichment'] = (df[sample_col] + 1) / (background + 1)
        
        # Identify peaks
        df['is_peak'] = df['fold_enrichment'] >= self.fold_threshold
        
        # Get peak regions
        peaks = []
        in_peak = False
        peak_start = None
        
        for idx, row in df.iterrows():
            if row['is_peak'] and not in_peak:
                # Start of a new peak
                in_peak = True
                peak_start = idx
            elif not row['is_peak'] and in_peak:
                # End of current peak
                in_peak = False
                peak_end = idx - 1
                peaks.append({
                    'chrom': df.loc[peak_start, 'Chrom'],
                    'start': df.loc[peak_start, 'WindowStart'],
                    'end': df.loc[peak_end, 'WindowEnd'],
                    'max_fold': df.loc[peak_start:peak_end, 'fold_enrichment'].max(),
                    'mean_signal': df.loc[peak_start:peak_end, sample_col].mean(),
                    'peak_width': df.loc[peak_end, 'WindowEnd'] - df.loc[peak_start, 'WindowStart'],
                    'scale': scale,
                    'sample': sample_col
                })
        
        # Handle case where peak extends to end of chromosome
        if in_peak:
            peak_end = len(df) - 1
            peaks.append({
                'chrom': df.loc[peak_start, 'Chrom'],
                'start': df.loc[peak_start, 'WindowStart'],
                'end': df.loc[peak_end, 'WindowEnd'],
                'max_fold': df.loc[peak_start:peak_end, 'fold_enrichment'].max(),
                'mean_signal': df.loc[peak_start:peak_end, sample_col].mean(),
                'peak_width': df.loc[peak_end, 'WindowEnd'] - df.loc[peak_start, 'WindowStart'],
                'scale': scale,
                'sample': sample_col
            })
            
        return peaks
    
    def calculate_overlap(self, peak1, peak2):
        """Calculate overlap fraction between two peaks"""
        if peak1['chrom'] != peak2['chrom']:
            return 0
        
        overlap_start = max(peak1['start'], peak2['start'])
        overlap_end = min(peak1['end'], peak2['end'])
        
        if overlap_start >= overlap_end:
            return 0
        
        overlap_length = overlap_end - overlap_start
        peak1_length = peak1['end'] - peak1['start']
        peak2_length = peak2['end'] - peak2['start']
        
        # Return reciprocal overlap fraction
        return min(overlap_length / peak1_length, overlap_length / peak2_length)
    
    def classify_hotspots(self, all_peaks, n_samples=3):
        """Classify hotspots based on reproducibility"""
        # Group peaks by chromosome and position
        hotspots = []
        
        # First, collect all unique peak regions
        unique_regions = []
        for peaks in all_peaks.values():
            for peak in peaks:
                # Check if this peak overlaps with any existing region
                merged = False
                for region in unique_regions:
                    if self.calculate_overlap(peak, region) >= self.overlap_threshold:
                        # Merge regions
                        region['start'] = min(region['start'], peak['start'])
                        region['end'] = max(region['end'], peak['end'])
                        region['peaks'].append(peak)
                        merged = True
                        break
                
                if not merged:
                    # Create new region
                    unique_regions.append({
                        'chrom': peak['chrom'],
                        'start': peak['start'],
                        'end': peak['end'],
                        'peaks': [peak]
                    })
        
        # Classify each region
        for region in unique_regions:
            # Count samples that have peaks in this region
            samples_with_peak = set([p['sample'] for p in region['peaks']])
            n_samples_with_peak = len(samples_with_peak)
            
            # Calculate statistics
            all_signals = [p['mean_signal'] for p in region['peaks']]
            all_folds = [p['max_fold'] for p in region['peaks']]
            
            # Determine classification
            if n_samples_with_peak == n_samples:
                classification = 'core'
            elif n_samples_with_peak >= 2:
                classification = 'recurrent'
            else:
                classification = 'variable'
            
            # Find the best scale (where signal is highest)
            scales_signals = {}
            for peak in region['peaks']:
                scale = peak['scale']
                if scale not in scales_signals:
                    scales_signals[scale] = []
                scales_signals[scale].append(peak['mean_signal'])
            
            best_scale = max(scales_signals.keys(), 
                           key=lambda s: np.mean(scales_signals[s]))
            
            hotspot = {
                'chrom': region['chrom'],
                'start': region['start'],
                'end': region['end'],
                'width': region['end'] - region['start'],
                'classification': classification,
                'n_samples': n_samples_with_peak,
                'mean_signal': np.mean(all_signals),
                'max_signal': np.max(all_signals),
                'signal_cv': np.std(all_signals) / np.mean(all_signals) if np.mean(all_signals) > 0 else 0,
                'mean_fold': np.mean(all_folds),
                'max_fold': np.max(all_folds),
                'best_scale': best_scale,
                'samples': list(samples_with_peak)
            }
            
            hotspots.append(hotspot)
        
        return hotspots
    
    def run_analysis(self, file_paths):
        """Run the complete hotspot detection pipeline"""
        print("Loading data...")
        data = self.load_data(file_paths)
        
        # Detect peaks at each scale
        all_peaks = {}
        for scale, scale_data in data.items():
            print(f"\nAnalyzing {scale} scale...")
            df = scale_data['df']
            samples = scale_data['samples']
            
            scale_peaks = []
            for sample in samples:
                print(f"  - Processing {sample}")
                peaks = self.detect_peaks_single_sample(df.copy(), sample, scale)
                scale_peaks.extend(peaks)
                print(f"    Found {len(peaks)} peaks")
            
            all_peaks[scale] = scale_peaks
        
        # Classify hotspots
        print("\nClassifying hotspots...")
        hotspots = self.classify_hotspots(all_peaks, n_samples=len(data[list(data.keys())[0]]['samples']))
        
        # Sort by classification and signal
        hotspots_df = pd.DataFrame(hotspots)
        if not hotspots_df.empty:
            hotspots_df = hotspots_df.sort_values(
                ['classification', 'mean_signal'], 
                ascending=[True, False]
            )
        
        return hotspots_df
    
    def generate_report(self, hotspots_df, output_prefix):
        """Generate summary report and output files"""
        # Save full results
        hotspots_df.to_csv(f"{output_prefix}_all_hotspots.csv", index=False)
        
        # Generate summary statistics
        summary = []
        for classification in ['core', 'recurrent', 'variable']:
            class_df = hotspots_df[hotspots_df['classification'] == classification]
            summary.append({
                'Classification': classification,
                'Count': len(class_df),
                'Mean_Width': class_df['width'].mean() if len(class_df) > 0 else 0,
                'Mean_Signal': class_df['mean_signal'].mean() if len(class_df) > 0 else 0,
                'Mean_Fold': class_df['mean_fold'].mean() if len(class_df) > 0 else 0
            })
        
        summary_df = pd.DataFrame(summary)
        summary_df.to_csv(f"{output_prefix}_summary.csv", index=False)
        
        # Generate BED files for each classification
        for classification in ['core', 'recurrent', 'variable']:
            class_df = hotspots_df[hotspots_df['classification'] == classification]
            if not class_df.empty:
                bed_df = class_df[['chrom', 'start', 'end', 'classification', 'mean_signal']]
                bed_df.to_csv(f"{output_prefix}_{classification}_hotspots.bed", 
                            sep='\t', header=False, index=False)
        
        # Print summary
        print("\n" + "="*60)
        print("HOTSPOT DETECTION SUMMARY")
        print("="*60)
        print(summary_df.to_string(index=False))
        print("\nTop 10 Core Hotspots:")
        if 'core' in hotspots_df['classification'].values:
            core_df = hotspots_df[hotspots_df['classification'] == 'core'].head(10)
            print(core_df[['chrom', 'start', 'end', 'mean_signal', 'mean_fold']].to_string(index=False))
        else:
            print("No core hotspots found")

def main():
    parser = argparse.ArgumentParser(description='Detect eccDNA hotspots from multi-scale data')
    parser.add_argument('--input-10kb', required=True, help='10kb window matrix CSV file')
    parser.add_argument('--input-50kb', required=True, help='50kb window matrix CSV file')
    parser.add_argument('--input-100kb', required=True, help='100kb window matrix CSV file')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')
    parser.add_argument('--fold-threshold', type=float, default=3, 
                       help='Minimum fold enrichment (default: 3)')
    parser.add_argument('--fdr-threshold', type=float, default=0.01,
                       help='FDR threshold (default: 0.01)')
    
    args = parser.parse_args()
    
    # Create file paths dictionary
    file_paths = {
        '10kb': args.input_10kb,
        '50kb': args.input_50kb,
        '100kb': args.input_100kb
    }
    
    # Initialize detector
    detector = eccDNAHotspotDetector(
        fold_threshold=args.fold_threshold,
        fdr_threshold=args.fdr_threshold
    )
    
    # Run analysis
    hotspots_df = detector.run_analysis(file_paths)
    
    # Generate report
    detector.generate_report(hotspots_df, args.output_prefix)

if __name__ == "__main__":
    # If running directly without command line arguments
    if len(sys.argv) == 1:
        # Example usage
        file_paths = {
            '10kb': 'U87_window_10kb_matrix.csv',
            '50kb': 'U87_window_50kb_matrix.csv',
            '100kb': 'U87_window_100kb_matrix.csv'
        }
        
        detector = eccDNAHotspotDetector(fold_threshold=3)
        hotspots_df = detector.run_analysis(file_paths)
        detector.generate_report(hotspots_df, 'U87_hotspots')
    else:
        main()
