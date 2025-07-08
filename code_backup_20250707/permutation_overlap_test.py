#!/usr/bin/env python3
"""
Permutation Test for Genomic Region Overlap Analysis
Uses bedtools shuffle/intersect to test if overlap between two sets of genomic regions
is significantly higher than expected by chance
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import pandas as pd
import numpy as np
import subprocess
import sys
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool, cpu_count
import tempfile
import time
import logging
from typing import Tuple, Dict, List

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PermutationOverlapTest:
    def __init__(self, args):
        self.input_file = args.input
        self.annotation_files = args.annotations
        self.genome_file = args.genome
        self.n_permutations = args.permutations
        self.n_cores = args.cores or max(1, cpu_count() - 1)
        self.output_dir = Path(args.output)
        self.keep_temp = args.keep_temp
        self.report_prefix = args.report_prefix
        self.exclude_sex_chr = args.exclude_sex_chr
        
        # Input specifications
        self.chr_col = args.chr_col
        self.start_col = args.start_col
        self.end_col = args.end_col
        self.name_col = args.name_col
        self.header = not args.no_header
        
        # Create output directories
        self.results_dir = self.output_dir / 'results'
        self.figures_dir = self.output_dir / 'figures'
        self.temp_dir = self.output_dir / 'temp'
        
        for d in [self.results_dir, self.figures_dir, self.temp_dir]:
            d.mkdir(exist_ok=True, parents=True)
        
        logger.info(f"Using {self.n_cores} CPU cores")
        logger.info(f"Output directory: {self.output_dir}")
        
        # Initialize data containers
        self.input_bed = None
        self.input_bed_file = None
        self.genome_sizes = None

    def load_input_regions(self) -> Tuple[pd.DataFrame, Path]:
        """Load input regions from CSV file and convert to BED format"""
        logger.info(f"Loading input regions from {self.input_file}...")
        
        # Load CSV
        if self.header:
            df = pd.read_csv(self.input_file)
        else:
            df = pd.read_csv(self.input_file, header=None)
            # Assign column names based on position
            col_names = {}
            col_names[self.chr_col] = 'chr'
            col_names[self.start_col] = 'start'
            col_names[self.end_col] = 'end'
            if self.name_col is not None:
                col_names[self.name_col] = 'name'
            df.rename(columns=col_names, inplace=True)
        
        # Extract relevant columns
        if self.header:
            bed_df = df[[self.chr_col, self.start_col, self.end_col]].copy()
            bed_df.columns = ['chr', 'start', 'end']
            if self.name_col and self.name_col in df.columns:
                bed_df['name'] = df[self.name_col]
            else:
                bed_df['name'] = [f"region_{i}" for i in range(len(df))]
        else:
            bed_df = df[['chr', 'start', 'end']].copy()
            if 'name' not in bed_df.columns:
                bed_df['name'] = [f"region_{i}" for i in range(len(df))]
        
        # Filter out mitochondrial and unplaced chromosomes
        original_count = len(bed_df)
        bed_df = bed_df[~bed_df['chr'].str.contains('chrM|_', na=False)]
        
        # Optional: filter out sex chromosomes
        if self.exclude_sex_chr:
            bed_df = bed_df[~bed_df['chr'].str.contains('chrX|chrY', na=False)]
        
        filtered_count = len(bed_df)
        
        logger.info(f"  Loaded {original_count} regions")
        if original_count != filtered_count:
            filters = ["chrM and unplaced"]
            if self.exclude_sex_chr:
                filters.append("chrX/Y")
            logger.info(f"  Filtered to {filtered_count} regions (removed {', '.join(filters)})")
        
        # Save as BED file
        bed_file = self.temp_dir / 'input_regions.bed'
        bed_df[['chr', 'start', 'end', 'name']].to_csv(
            bed_file, sep='\t', index=False, header=False
        )
        
        self.input_bed = bed_df
        self.input_bed_file = bed_file
        
        return bed_df, bed_file

    def create_genome_file(self) -> Path:
        """Create genome file for bedtools shuffle (hg38)"""
        logger.info("Creating hg38 genome file...")
        
        self.genome_sizes = {
            'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
            'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
            'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
            'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
            'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
            'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
            'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
            'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
        }
        
        genome_file = self.temp_dir / 'hg38.genome'
        with open(genome_file, 'w') as f:
            for chrom, size in self.genome_sizes.items():
                f.write(f"{chrom}\t{size}\n")
        
        return genome_file

    def run_bedtools_intersect(self, a_file: Path, b_file: Path) -> int:
        """Run bedtools intersect and count unique overlapping regions from A"""
        cmd = f"bedtools intersect -a {a_file} -b {b_file} -wa | sort -u | wc -l"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"Bedtools error: {result.stderr}")
            return 0
        
        return int(result.stdout.strip())

    def single_shuffle_intersect(self, args: Tuple) -> int:
        """Single shuffle and intersect operation for parallel processing"""
        iteration, input_file, annotation_file, genome_file, temp_base = args
        
        # Create temporary file for shuffled regions
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False, dir=temp_base) as tmp:
            shuffled_file = tmp.name
        
        try:
            # Shuffle maintaining chromosome distribution
            shuffle_cmd = (
                f"bedtools shuffle -i {input_file} -g {genome_file} "
                f"-chrom -seed {iteration} > {shuffled_file} 2>/dev/null"
            )
            result = subprocess.run(shuffle_cmd, shell=True)
            
            if result.returncode != 0:
                logger.error(f"Shuffle failed for iteration {iteration}")
                return 0
            
            # Count overlaps
            count = self.run_bedtools_intersect(shuffled_file, annotation_file)
            
        finally:
            # Clean up
            if os.path.exists(shuffled_file):
                os.unlink(shuffled_file)
        
        return count

    def run_permutation_test(self, annotation_file: Path, annotation_name: str) -> Dict:
        """Run permutation test for a single annotation file"""
        logger.info(f"\nAnalyzing overlap with {annotation_name}...")
        
        # Get observed overlap count
        observed_count = self.run_bedtools_intersect(self.input_bed_file, annotation_file)
        logger.info(f"  Observed overlaps: {observed_count}")
        
        # Prepare genome file
        if self.genome_file:
            genome_file = self.genome_file
        else:
            genome_file = self.create_genome_file()
        
        # Run permutations in parallel
        logger.info(f"  Running {self.n_permutations} permutations...")
        start_time = time.time()
        
        with tempfile.TemporaryDirectory(dir=self.temp_dir) as temp_base:
            # Prepare arguments for parallel processing
            args_list = [
                (i, self.input_bed_file, annotation_file, genome_file, temp_base)
                for i in range(self.n_permutations)
            ]
            
            # Run in parallel
            with Pool(processes=self.n_cores) as pool:
                chunksize = max(1, self.n_permutations // (self.n_cores * 4))
                null_counts = pool.map(self.single_shuffle_intersect, args_list, chunksize=chunksize)
        
        null_counts = np.array(null_counts)
        elapsed_time = time.time() - start_time
        logger.info(f"  Completed in {elapsed_time:.1f} seconds")
        
        # Calculate statistics
        null_mean = np.mean(null_counts)
        null_std = np.std(null_counts)
        
        # Calculate z-score
        if null_std == 0:
            z_score = np.inf if observed_count > null_mean else -np.inf
        else:
            z_score = (observed_count - null_mean) / null_std
        
        # Calculate fold enrichment
        fold_enrichment = observed_count / null_mean if null_mean > 0 else np.inf
        
        # Calculate empirical p-value
        p_value = np.sum(null_counts >= observed_count) / len(null_counts)
        if p_value == 0:
            p_value = 1 / (len(null_counts) + 1)
        
        # Calculate percentile
        percentile = np.sum(null_counts < observed_count) / len(null_counts) * 100
        
        results = {
            'annotation': annotation_name,
            'observed_count': observed_count,
            'expected_mean': null_mean,
            'expected_std': null_std,
            'fold_enrichment': fold_enrichment,
            'z_score': z_score,
            'p_value': p_value,
            'percentile': percentile,
            'n_permutations': self.n_permutations,
            'null_distribution': null_counts
        }
        
        return results

    def plot_null_distribution(self, results: Dict, output_name: str):
        """Plot null distribution with observed value"""
        plt.figure(figsize=(10, 6))
        
        null_counts = results['null_distribution']
        observed = results['observed_count']
        
        # Plot histogram
        plt.hist(null_counts, bins=50, density=True, alpha=0.7, color='skyblue', 
                edgecolor='black', label='Null distribution')
        
        # Add observed value line
        plt.axvline(observed, color='red', linestyle='--', linewidth=2, 
                   label=f'Observed ({observed})')
        
        # Add expected mean line
        plt.axvline(results['expected_mean'], color='green', linestyle=':', 
                   linewidth=2, label=f'Expected ({results["expected_mean"]:.1f})')
        
        # Add text box with statistics
        textstr = (
            f'Fold enrichment: {results["fold_enrichment"]:.2f}\n'
            f'Z-score: {results["z_score"]:.2f}\n'
            f'P-value: {results["p_value"]:.4f}\n'
            f'Percentile: {results["percentile"]:.1f}%'
        )
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        plt.text(0.98, 0.97, textstr, transform=plt.gca().transAxes, 
                fontsize=10, verticalalignment='top', horizontalalignment='right',
                bbox=props)
        
        plt.xlabel('Number of overlapping regions')
        plt.ylabel('Density')
        plt.title(f'Permutation Test: {results["annotation"]}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Save figure
        plt.savefig(self.figures_dir / f'{output_name}_null_distribution.png', dpi=300)
        plt.close()

    def plot_summary(self, all_results: List[Dict]):
        """Create summary plots for all annotations"""
        if len(all_results) < 2:
            return
        
        # Prepare data for plotting
        annotations = [r['annotation'] for r in all_results]
        fold_enrichments = [r['fold_enrichment'] for r in all_results]
        p_values = [r['p_value'] for r in all_results]
        z_scores = [r['z_score'] for r in all_results]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        # Bar plot of fold enrichments
        bars = ax1.bar(range(len(annotations)), fold_enrichments)
        
        # Color bars by significance
        colors = ['red' if p < 0.001 else 'orange' if p < 0.01 else 'yellow' if p < 0.05 else 'lightblue' 
                 for p in p_values]
        for bar, color in zip(bars, colors):
            bar.set_color(color)
        
        ax1.set_xticks(range(len(annotations)))
        ax1.set_xticklabels(annotations, rotation=45, ha='right')
        ax1.set_ylabel('Fold Enrichment')
        ax1.set_title('Fold Enrichment by Annotation')
        ax1.axhline(y=1, color='black', linestyle='--', alpha=0.5)
        ax1.grid(True, alpha=0.3)
        
        # Add significance stars
        for i, (fe, p) in enumerate(zip(fold_enrichments, p_values)):
            if p < 0.001:
                ax1.text(i, fe + 0.1, '***', ha='center', va='bottom')
            elif p < 0.01:
                ax1.text(i, fe + 0.1, '**', ha='center', va='bottom')
            elif p < 0.05:
                ax1.text(i, fe + 0.1, '*', ha='center', va='bottom')
        
        # Volcano plot
        neg_log_p = [-np.log10(p) if p > 0 else -np.log10(1/(all_results[0]['n_permutations']+1)) 
                     for p in p_values]
        
        ax2.scatter(z_scores, neg_log_p, s=100, alpha=0.7)
        
        # Add labels for significant points
        for i, (z, nlp, ann) in enumerate(zip(z_scores, neg_log_p, annotations)):
            if p_values[i] < 0.05:
                ax2.annotate(ann, (z, nlp), xytext=(5, 5), textcoords='offset points', 
                           fontsize=8, alpha=0.7)
        
        # Add significance thresholds
        ax2.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, label='p=0.05')
        ax2.axhline(y=-np.log10(0.01), color='gray', linestyle=':', alpha=0.5, label='p=0.01')
        
        ax2.set_xlabel('Z-score')
        ax2.set_ylabel('-log10(p-value)')
        ax2.set_title('Statistical Significance of Overlaps')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'summary_analysis.png', dpi=300)
        plt.close()

    def generate_report(self, all_results: List[Dict]):
        """Generate markdown report"""
        timestamp = pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')
        
        report = f"""# Permutation Test for Genomic Region Overlap Analysis

## Analysis Summary
- **Date**: {timestamp}
- **Input file**: {self.input_file}
- **Input regions**: {len(self.input_bed):,}
- **Permutations**: {self.n_permutations:,}
- **CPU cores used**: {self.n_cores}
- **Output directory**: {self.output_dir}

## Methods
This analysis uses permutation testing to assess whether the overlap between input regions 
and various genomic annotations is greater than expected by chance. The method:

1. Counts observed overlaps between input regions and each annotation
2. Randomly shuffles input regions {self.n_permutations} times (maintaining chromosome distribution)
3. Counts overlaps for each shuffled set to build null distribution
4. Calculates enrichment statistics and empirical p-values

## Results

| Annotation | Observed | Expected (mean±std) | Fold Enrichment | Z-score | P-value | Significance |
|------------|----------|-------------------|-----------------|---------|---------|--------------|
"""
        
        for r in all_results:
            sig = '***' if r['p_value'] < 0.001 else '**' if r['p_value'] < 0.01 else '*' if r['p_value'] < 0.05 else 'ns'
            report += (
                f"| {r['annotation']} | {r['observed_count']} | "
                f"{r['expected_mean']:.1f}±{r['expected_std']:.1f} | "
                f"{r['fold_enrichment']:.2f} | {r['z_score']:.2f} | "
                f"{r['p_value']:.4f} | {sig} |\n"
            )
        
        report += """
### Significance levels:
- *** : p < 0.001 (highly significant)
- **  : p < 0.01 (very significant)  
- *   : p < 0.05 (significant)
- ns  : p ≥ 0.05 (not significant)

## Interpretation

### Fold Enrichment:
- Values > 1 indicate enrichment (more overlap than expected)
- Values < 1 indicate depletion (less overlap than expected)
- Values ≈ 1 indicate no enrichment/depletion

### Z-score:
- Measures how many standard deviations the observed value is from the expected mean
- |Z| > 2 generally indicates significant deviation from random expectation

### P-value:
- Empirical p-value from permutation test
- Represents the fraction of permutations with overlap ≥ observed
- Lower values indicate stronger evidence against the null hypothesis

## Output Files

### Results:
- `results/permutation_test_results.csv`: Detailed statistics for all annotations
- `results/analysis_report.md`: This report

### Figures:
- `figures/*_null_distribution.png`: Null distribution plots for each annotation
- `figures/summary_analysis.png`: Summary plots (if multiple annotations analyzed)

## Software Information
- Python version: {sys.version.split()[0]}
- bedtools: Required for overlap calculations
- Random seed: Set per iteration for reproducibility

---
*Report generated by Permutation Test for Genomic Region Overlap Analysis*
"""
        
        # Save report with custom prefix
        report_file = f'{self.report_prefix}_report.md' if self.report_prefix else 'analysis_report.md'
        with open(self.results_dir / report_file, 'w') as f:
            f.write(report)
        
        logger.info(f"Report saved to: {self.results_dir / report_file}")

    def run(self):
        """Run the complete analysis"""
        logger.info("="*60)
        logger.info("Permutation Test for Genomic Region Overlap Analysis")
        logger.info("="*60)
        
        start_time = time.time()
        
        try:
            # Load input regions
            input_df, input_bed_file = self.load_input_regions()
            
            # Process each annotation file
            all_results = []
            
            for i, annotation_file in enumerate(self.annotation_files):
                annotation_path = Path(annotation_file)
                if not annotation_path.exists():
                    logger.warning(f"Annotation file not found: {annotation_file}")
                    continue
                
                # Use filename as annotation name
                annotation_name = annotation_path.stem
                
                # Run permutation test
                results = self.run_permutation_test(annotation_path, annotation_name)
                
                # Create plot
                safe_name = annotation_name.replace(' ', '_').replace('/', '_')
                self.plot_null_distribution(results, safe_name)
                
                # Remove null distribution from results before appending
                results_summary = {k: v for k, v in results.items() if k != 'null_distribution'}
                all_results.append(results_summary)
            
            if not all_results:
                logger.error("No annotation files could be processed!")
                return
            
            # Save results
            results_df = pd.DataFrame(all_results)
            results_file = f'{self.report_prefix}_results.csv' if self.report_prefix else 'permutation_test_results.csv'
            results_df.to_csv(self.results_dir / results_file, index=False)
            
            # Create summary plots if multiple annotations
            if len(all_results) > 1:
                self.plot_summary(all_results)
            
            # Generate report
            self.generate_report(all_results)
            
            # Clean up temp files unless requested to keep
            if not self.keep_temp:
                import shutil
                shutil.rmtree(self.temp_dir)
                logger.info("Temporary files cleaned up")
            
            elapsed_time = time.time() - start_time
            logger.info("="*60)
            logger.info(f"Analysis completed in {elapsed_time:.1f} seconds")
            logger.info(f"Results saved to: {self.output_dir}")
            logger.info("="*60)
            
        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            raise

def check_dependencies():
    """Check if required dependencies are installed"""
    # Check bedtools
    try:
        result = subprocess.run(['bedtools', '--version'], capture_output=True, text=True)
        logger.info(f"Found {result.stdout.strip()}")
    except:
        logger.error("bedtools not found! Please install bedtools (v2.30.0+)")
        logger.error("Installation: conda install -c bioconda bedtools")
        sys.exit(1)
    
    # Check Python packages
    required_packages = ['pandas', 'numpy', 'matplotlib', 'seaborn']
    missing = []
    
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing.append(package)
    
    if missing:
        logger.error(f"Missing Python packages: {', '.join(missing)}")
        logger.error(f"Installation: pip install {' '.join(missing)}")
        sys.exit(1)

def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Permutation test for genomic region overlap analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with eccDNA CSV and enhancer regions
  %(prog)s -i eccDNA.csv -a enhancers.bed -o results/
  
  # Multiple annotation files
  %(prog)s -i eccDNA.csv -a enhancers.bed promoters.bed ATAC_peaks.bed -o results/
  
  # Custom column names for CSV input
  %(prog)s -i regions.csv --chr-col chromosome --start-col begin --end-col stop -a annotations.bed
  
  # More permutations and cores
  %(prog)s -i eccDNA.csv -a enhancers.bed -p 10000 --cores 32 -o results/
  
  # Use custom genome file
  %(prog)s -i eccDNA.csv -a enhancers.bed -g hg19.genome -o results/
  
  # Exclude sex chromosomes and use custom report prefix
  %(prog)s -i eccDNA.csv -a enhancers.bed --exclude-sex-chr --report-prefix eccDNA_enhancer -o results/
        """
    )
    
    # Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='Input CSV file with genomic regions')
    parser.add_argument('-a', '--annotations', nargs='+', required=True,
                       help='Annotation BED file(s) to test for overlap enrichment')
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory for results')
    
    # Optional arguments
    parser.add_argument('-p', '--permutations', type=int, default=1000,
                       help='Number of permutations (default: 1000)')
    parser.add_argument('-g', '--genome', default=None,
                       help='Genome file for bedtools shuffle (default: built-in hg38)')
    parser.add_argument('--cores', type=int, default=None,
                       help='Number of CPU cores to use (default: auto-detect)')
    
    # CSV format options
    parser.add_argument('--chr-col', default='eChr',
                       help='Column name for chromosome (default: eChr)')
    parser.add_argument('--start-col', default='eStart',
                       help='Column name for start position (default: eStart)')
    parser.add_argument('--end-col', default='eEnd',
                       help='Column name for end position (default: eEnd)')
    parser.add_argument('--name-col', default='eName',
                       help='Column name for region names (default: eName)')
    parser.add_argument('--no-header', action='store_true',
                       help='Input CSV has no header row')
    
    # Other options
    parser.add_argument('--keep-temp', action='store_true',
                       help='Keep temporary files for debugging')
    parser.add_argument('--exclude-sex-chr', action='store_true',
                       help='Exclude chrX and chrY from analysis')
    parser.add_argument('--report-prefix', default='',
                       help='Prefix for output files (default: none)')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    # Check dependencies
    check_dependencies()
    
    # Run analysis
    analyzer = PermutationOverlapTest(args)
    analyzer.run()

if __name__ == "__main__":
    main()
