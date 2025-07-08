#!/usr/bin/env python3
"""
eccDNA Hotspot Detection and Classification Script with Permutation Test
Author: Your Name
Date: 2024
Description: Multi-scale detection of eccDNA hotspots with reproducibility classification
             and global permutation significance test
"""

import pandas as pd
import numpy as np
from scipy import stats
import argparse
from pathlib import Path
import warnings
import sys
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')

class eccDNAHotspotDetector:
    def __init__(self, fold_threshold=3, overlap_threshold=0.5, n_permutations=1000, seed=0):
        """
        Initialize the detector with parameters
        """
        self.fold_threshold = fold_threshold
        self.overlap_threshold = overlap_threshold
        self.n_permutations = n_permutations
        self.rng = np.random.default_rng(seed)
        
    def load_data(self, file_paths):
        data = {}
        for scale, path in file_paths.items():
            df = pd.read_csv(path)
            sample_cols = [col for col in df.columns if col not in ['Chrom', 'WindowStart', 'WindowEnd']]
            data[scale] = {
                'df': df,
                'samples': sample_cols
            }
        return data
    
    def calculate_background(self, values):
        non_zero = values[values > 0]
        if len(non_zero) > 0:
            return np.median(non_zero)
        else:
            return 0.1  # Pseudocount

    def detect_peaks_single_sample(self, df, sample_col, scale):
        background = self.calculate_background(df[sample_col].values)
        df['fold_enrichment'] = (df[sample_col] + 1) / (background + 1)
        df['is_peak'] = df['fold_enrichment'] >= self.fold_threshold
        
        peaks = []
        in_peak = False
        peak_start = None
        
        for idx, row in df.iterrows():
            if row['is_peak'] and not in_peak:
                in_peak = True
                peak_start = idx
            elif not row['is_peak'] and in_peak:
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
        if peak1['chrom'] != peak2['chrom']:
            return 0
        overlap_start = max(peak1['start'], peak2['start'])
        overlap_end = min(peak1['end'], peak2['end'])
        if overlap_start >= overlap_end:
            return 0
        overlap_length = overlap_end - overlap_start
        peak1_length = peak1['end'] - peak1['start']
        peak2_length = peak2['end'] - peak2['start']
        return min(overlap_length / peak1_length, overlap_length / peak2_length)
    
    def classify_hotspots(self, all_peaks, n_samples=3):
        hotspots = []
        unique_regions = []
        for peaks in all_peaks.values():
            for peak in peaks:
                merged = False
                for region in unique_regions:
                    if self.calculate_overlap(peak, region) >= self.overlap_threshold:
                        region['start'] = min(region['start'], peak['start'])
                        region['end'] = max(region['end'], peak['end'])
                        region['peaks'].append(peak)
                        merged = True
                        break
                if not merged:
                    unique_regions.append({
                        'chrom': peak['chrom'],
                        'start': peak['start'],
                        'end': peak['end'],
                        'peaks': [peak]
                    })
        for region in unique_regions:
            samples_with_peak = set([p['sample'] for p in region['peaks']])
            n_samples_with_peak = len(samples_with_peak)
            all_signals = [p['mean_signal'] for p in region['peaks']]
            all_folds = [p['max_fold'] for p in region['peaks']]
            if n_samples_with_peak == n_samples:
                classification = 'core'
            elif n_samples_with_peak >= 2:
                classification = 'recurrent'
            else:
                classification = 'variable'
            scales_signals = {}
            for peak in region['peaks']:
                scale = peak['scale']
                if scale not in scales_signals:
                    scales_signals[scale] = []
                scales_signals[scale].append(peak['mean_signal'])
            best_scale = max(scales_signals.keys(), key=lambda s: np.mean(scales_signals[s]))
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
    
    def permutation_pval_for_hotspot(self, hotspot, window_df, sample_cols, n_perm=1000):
        """
        For each hotspot, do permutation test: random sample the same length window N times.
        Return mean_signal, null distribution, pval.
        """
        chrom = hotspot['chrom']
        length = hotspot['width']
        # restrict to the same chrom
        dfc = window_df[window_df['Chrom']==chrom].copy()
        window_size = dfc['WindowEnd'] - dfc['WindowStart']
        # mean_signal in all samples/peaks contributing to this hotspot, in best scale
        # Here, we use the mean signal of all involved windows, all samples.
        dfc['in_hotspot'] = (dfc['WindowStart'] < hotspot['end']) & (dfc['WindowEnd'] > hotspot['start'])
        observed = dfc.loc[dfc['in_hotspot'], sample_cols].values.mean()
        # collect eligible start indices for permutation
        valid_starts = []
        for idx in range(len(dfc)):
            win_start = dfc.iloc[idx]['WindowStart']
            win_end = win_start + length
            if win_end > dfc['WindowEnd'].max():
                continue
            # check if all windows in this interval are in the same chrom
            overlapped = (dfc['WindowStart'] < win_end) & (dfc['WindowEnd'] > win_start)
            if overlapped.sum() > 0:
                valid_starts.append(win_start)
        # random sample
        null_dist = []
        if not valid_starts:
            return observed, [observed], 1.0
        for _ in range(n_perm):
            start = self.rng.choice(valid_starts)
            end = start + length
            overlapped = (dfc['WindowStart'] < end) & (dfc['WindowEnd'] > start)
            null_val = dfc.loc[overlapped, sample_cols].values.mean()
            null_dist.append(null_val)
        null_dist = np.array(null_dist)
        # pval: >= observed (conservative, one-tailed)
        pval = (np.sum(null_dist >= observed) + 1) / (n_perm + 1)
        return observed, null_dist, pval
    
    def run_analysis(self, file_paths, n_permutations=1000):
        print("Loading data...")
        data = self.load_data(file_paths)
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
        print("\nClassifying hotspots...")
        n_samples = len(data[list(data.keys())[0]]['samples'])
        hotspots = self.classify_hotspots(all_peaks, n_samples=n_samples)
        hotspots_df = pd.DataFrame(hotspots)
        if hotspots_df.empty:
            return hotspots_df
        # ==== Permutation test ====
        print("\nRunning permutation tests...")
        # Use best scale window matrix as background
        all_window_dfs = {scale: data[scale]['df'] for scale in data}
        all_sample_cols = {scale: data[scale]['samples'] for scale in data}
        # Add pval and FDR columns
        pvals = []
        for i, row in hotspots_df.iterrows():
            best_scale = row['best_scale']
            win_df = all_window_dfs[best_scale]
            sample_cols = all_sample_cols[best_scale]
            obs, null, pval = self.permutation_pval_for_hotspot(row, win_df, sample_cols, n_perm=n_permutations)
            pvals.append(pval)
        hotspots_df['pval'] = pvals
        # FDR correction
        _, qvals, _, _ = multipletests(hotspots_df['pval'], alpha=0.05, method='fdr_bh')
        hotspots_df['qval'] = qvals
        hotspots_df['significant'] = hotspots_df['qval'] < 0.05
        # Sort
        hotspots_df = hotspots_df.sort_values(['classification', 'mean_signal'], ascending=[True, False])
        return hotspots_df
    
    def generate_report(self, hotspots_df, output_prefix):
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
                'Mean_Fold': class_df['mean_fold'].mean() if len(class_df) > 0 else 0,
                'FDR<0.05': class_df['significant'].sum()
            })
        summary_df = pd.DataFrame(summary)
        summary_df.to_csv(f"{output_prefix}_summary.csv", index=False)
        # Generate BED files for each classification (significant only)
        for classification in ['core', 'recurrent', 'variable']:
            class_df = hotspots_df[(hotspots_df['classification'] == classification) & (hotspots_df['significant'])]
            if not class_df.empty:
                bed_df = class_df[['chrom', 'start', 'end', 'classification', 'mean_signal']]
                bed_df.to_csv(f"{output_prefix}_{classification}_hotspots.bed", 
                            sep='\t', header=False, index=False)
        print("\n" + "="*60)
        print("HOTSPOT DETECTION SUMMARY")
        print("="*60)
        print(summary_df.to_string(index=False))
        print("\nTop 10 Core Hotspots:")
        if 'core' in hotspots_df['classification'].values:
            core_df = hotspots_df[(hotspots_df['classification'] == 'core') & (hotspots_df['significant'])].head(10)
            print(core_df[['chrom', 'start', 'end', 'mean_signal', 'mean_fold', 'pval', 'qval']].to_string(index=False))
        else:
            print("No core hotspots found")
            
def main():
    parser = argparse.ArgumentParser(description='Detect eccDNA hotspots from multi-scale data with permutation FDR test')
    parser.add_argument('--input-10kb', required=True, help='10kb window matrix CSV file')
    parser.add_argument('--input-50kb', required=True, help='50kb window matrix CSV file')
    parser.add_argument('--input-100kb', required=True, help='100kb window matrix CSV file')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')
    parser.add_argument('--fold-threshold', type=float, default=3, 
                       help='Minimum fold enrichment (default: 3)')
    parser.add_argument('--n-permutations', type=int, default=1000,
                       help='Number of permutations for permutation test (default: 1000)')
    args = parser.parse_args()
    file_paths = {
        '10kb': args.input_10kb,
        '50kb': args.input_50kb,
        '100kb': args.input_100kb
    }
    detector = eccDNAHotspotDetector(
        fold_threshold=args.fold_threshold,
        n_permutations=args.n_permutations
    )
    hotspots_df = detector.run_analysis(file_paths, n_permutations=args.n_permutations)
    detector.generate_report(hotspots_df, args.output_prefix)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        file_paths = {
            '10kb': 'U87_window_10kb_matrix.csv',
            '50kb': 'U87_window_50kb_matrix.csv',
            '100kb': 'U87_window_100kb_matrix.csv'
        }
        detector = eccDNAHotspotDetector(fold_threshold=3, n_permutations=500)
        hotspots_df = detector.run_analysis(file_paths, n_permutations=500)
        detector.generate_report(hotspots_df, 'U87_hotspots')
    else:
        main()
