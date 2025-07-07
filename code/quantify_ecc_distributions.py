# -*- coding: utf-8 -*-
"""
quantify_ecc_distributions.py (v5 - Data-Driven Bins)

Description:
This is the final version of the main analysis script.
It now uses a data-driven binning strategy based on observed period clusters
from the user's specific dataset to generate a clean, interpretable
periodicity fingerprint for each sample.
"""
import argparse
import os
import warnings
from itertools import product
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from statsmodels.tsa.stattools import acf
from statsmodels.nonparametric.smoothers_lowess import lowess
from tqdm import tqdm

# Suppress common warnings for cleaner output
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message="The behavior of DataFrame concatenation with empty or all-NA entries is deprecated.")

def analyze_single_sample(df_sample, cutoff, period_range, long_tail_min_len, tolerance=5):
    """
    Performs the core analysis, returning raw period intensities and long-tail stats.
    The internal logic of this function remains unchanged.
    """
    # Step 1: Create a 1-bp resolution signal vector (histogram)
    data = df_sample[df_sample["eLength"] <= cutoff]["eLength"]
    signal = np.bincount(data.astype(int), minlength=cutoff + 1)

    # Step 2: Pre-process the signal (Smooth and Detrend)
    smoothed_signal = gaussian_filter1d(signal.astype(float), sigma=2)
    trend = lowess(smoothed_signal, np.arange(len(smoothed_signal)), frac=0.1, it=0)[:, 1]
    detrended_signal = smoothed_signal - trend

    min_lag, max_lag = period_range

    # Step 3, Method A: Autocorrelation Function (ACF)
    try:
        acf_vals = acf(detrended_signal, nlags=max_lag, fft=True)
        acf_peaks_indices, acf_props = find_peaks(acf_vals[min_lag:], prominence=0.01, height=0.01, distance=5)
        acf_peaks = acf_peaks_indices + min_lag
        acf_heights = acf_vals[min_lag:][acf_peaks_indices]
        acf_prominences = acf_props.get('prominences', np.zeros_like(acf_peaks, dtype=float))
    except Exception:
        acf_peaks, acf_heights, acf_prominences = np.array([]), np.array([]), np.array([])

    # Step 3, Method B: Fast Fourier Transform (FFT)
    try:
        n_points = len(detrended_signal)
        fft_vals = np.fft.rfft(detrended_signal)
        fft_power = np.abs(fft_vals) ** 2
        freqs = np.fft.rfftfreq(n_points, d=1)
        periods = np.divide(1.0, freqs, out=np.full_like(freqs, np.inf), where=(freqs != 0))
        valid_indices = (periods > min_lag) & (periods < max_lag)
        fft_peaks_indices, _ = find_peaks(
            fft_power[valid_indices],
            prominence=np.quantile(fft_power[valid_indices], 0.8) if any(valid_indices) else 1,
            height=np.quantile(fft_power[valid_indices], 0.75) if any(valid_indices) else 1,
            distance=5
        )
        fft_periods = periods[valid_indices][fft_peaks_indices]
    except Exception:
        fft_periods = np.array([])
    
    # Step 4: Cross-validate to find common periods and their intensities
    period_intensities = {}
    for i, p_acf in enumerate(acf_peaks):
        if np.any(np.abs(p_acf - fft_periods) <= tolerance):
            closest_fft_peak = fft_periods[np.argmin(np.abs(p_acf - fft_periods))]
            common_period = round((p_acf + closest_fft_peak) / 2.0, 1)
            period_intensities[common_period] = acf_heights[i]

    # Calculate ~200bp noise strength separately for the other analysis
    noise_strength_200bp = 0.0
    for i, p in enumerate(acf_peaks):
        if 180 <= p <= 220:
            noise_strength_200bp = acf_prominences[i] * acf_heights[i]
            break

    # Calculate long-tail statistics
    long_tail_data = df_sample[df_sample["eLength"] > long_tail_min_len]["eLength"]
    total_eccdnas = len(df_sample)
    if len(long_tail_data) > 1:
        long_tail_stats = {
            "tail_mean": long_tail_data.mean(),
            "tail_median": long_tail_data.median(),
            "tail_abundance_ratio": len(long_tail_data) / total_eccdnas if total_eccdnas > 0 else 0,
        }
    else:
        long_tail_stats = {"tail_mean": 0, "tail_median": 0, "tail_abundance_ratio": 0}
        
    return {
        "period_intensities": period_intensities,
        "noise_strength_200bp": noise_strength_200bp,
        "long_tail_stats": long_tail_stats
    }

def main():
    parser = argparse.ArgumentParser(description="Quantify eccDNA Distributions (v5 - Data-Driven Bins).")
    parser.add_argument("--input", required=True, help="Path to the input CSV file.")
    parser.add_argument("--output", required=True, help="Path to the output directory.")
    parser.add_argument("--cutoffs", default="2000", help="Periodicity analysis length cutoff.")
    parser.add_argument("--long_tail_thresholds", default="5000", help="Thresholds for defining the long tail.")
    parser.add_argument("--min_records", type=int, default=1000, help="Minimum number of records per sample.")
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    df = pd.read_csv(args.input)
    cutoffs = [int(c) for c in args.cutoffs.split(',')]
    long_tail_thresholds = [int(t) for t in args.long_tail_thresholds.split(',')]
    samples = df['sample'].unique()
    
    long_tail_results = []
    period_results_long_format = []

    # =========== UPDATED: Data-Driven Period Bins ===========
    period_bins = {
        '009-012bp (Helix)': (9, 12),
        '018-024bp': (18, 24),
        '028-034bp': (28, 34),
        '040-048bp': (40, 48),
        '050-055bp': (50, 55),
        '060-075bp': (60, 75),
        '080-100bp': (80, 100),
        '110-135bp': (110, 135),
        '165-180bp': (165, 180),
        '195-205bp (Nucleosome)': (195, 205),
        '330-340bp': (330, 340)
    }
    # =======================================================

    analysis_combinations = list(product(samples, cutoffs, long_tail_thresholds))
    progress_bar = tqdm(analysis_combinations, desc="Quantifying & Binning Fingerprints")
    
    for sample_name, cutoff, long_tail_thresh in progress_bar:
        df_sample = df[df['sample'] == sample_name]
        if len(df_sample) < args.min_records: continue
            
        try:
            results = analyze_single_sample(df_sample, cutoff, (5, 500), long_tail_thresh)
            
            summary_row = {"sample": sample_name, "cutoff": cutoff, "long_tail_threshold": long_tail_thresh, "noise_strength_200bp": results['noise_strength_200bp'], **results['long_tail_stats']}
            long_tail_results.append(summary_row)
            
            binned_intensities = {bin_label: 0.0 for bin_label in period_bins}
            for period, intensity in results['period_intensities'].items():
                for bin_label, (start, end) in period_bins.items():
                    if start <= period < end:
                        binned_intensities[bin_label] += intensity
                        break
            
            for bin_label, total_intensity in binned_intensities.items():
                if total_intensity > 0:
                    period_results_long_format.append({"sample": sample_name, "period_bin": bin_label, "intensity": total_intensity})
        except Exception as e:
            tqdm.write(f"ERROR analyzing {sample_name}: {e}")
            
    if long_tail_results:
        pd.DataFrame(long_tail_results).to_csv(os.path.join(args.output, "all_samples_multithresh_summary.csv"), index=False)
        print("\nLong-tail summary saved.")

    if period_results_long_format:
        period_df = pd.DataFrame(period_results_long_format)
        binned_matrix = period_df.pivot_table(index='sample', columns='period_bin', values='intensity', fill_value=0)
        
        # Ensure the columns (bins) are sorted for consistent heatmaps
        sorted_bins = sorted(binned_matrix.columns)
        binned_matrix = binned_matrix[sorted_bins]

        matrix_filename = os.path.join(args.output, "binned_period_intensity_matrix.csv")
        binned_matrix.to_csv(matrix_filename)
        print(f"Binned period intensity matrix saved to: '{matrix_filename}'")

if __name__ == "__main__":
    main()
