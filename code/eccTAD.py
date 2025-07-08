#!/usr/bin/env python3
"""
eccTAD - eccDNA TAD Boundary Enrichment Analysis
Version 1.0
Analyzes enrichment of eccDNA at TAD boundaries using permutation tests
"""

# Set environment variables to prevent thread oversubscription
import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import pandas as pd
import numpy as np
import subprocess
import sys
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool, cpu_count
import tempfile
import time
import warnings
import logging

warnings.filterwarnings("ignore")
np.random.seed(42)

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class eccTADAnalysis:
    def __init__(self, args):
        self.eccdna_files = args.input
        self.tad_file = args.tad
        self.genome_file = args.genome
        self.n_cores = args.cores or max(1, cpu_count() - 2)
        self.n_shuffles = args.nshuffle
        self.prefix = args.prefix
        self.work_dir = Path(args.workdir)
        self.keep_sex_chr = (
            args.keep_sex_chr
        )  # New parameter for sex chromosome filtering

        # Create directories with prefix if provided
        if self.prefix:
            self.results_dir = self.work_dir / f"{self.prefix}_TAD_results"
            self.figures_dir = self.work_dir / f"{self.prefix}_TAD_figures"
            self.plot_data_dir = self.work_dir / f"{self.prefix}_TAD_plot_data"
        else:
            self.results_dir = self.work_dir / "TAD_results"
            self.figures_dir = self.work_dir / "TAD_figures"
            self.plot_data_dir = self.work_dir / "TAD_plot_data"

        self.temp_dir = self.work_dir / "temp"

        # Create directories
        for d in [
            self.results_dir,
            self.figures_dir,
            self.temp_dir,
            self.plot_data_dir,
        ]:
            d.mkdir(exist_ok=True, parents=True)

        logger.info(f"eccTAD Analysis v1.0")
        logger.info(f"Using {self.n_cores} CPU cores for parallel processing")
        logger.info(f"Working directory: {self.work_dir}")
        if self.prefix:
            logger.info(f"Output prefix: {self.prefix}")

        self.eccdna_data = None
        self.genome_sizes = None

    def load_eccdna_data(self):
        """Load and merge eccDNA data from input files"""
        logger.info("Loading eccDNA data...")
        all_eccdna = []

        for i, file in enumerate(self.eccdna_files):
            filepath = Path(file)
            if not filepath.exists():
                logger.error(f"File not found: {filepath}")
                continue

            try:
                df = pd.read_csv(filepath)
                df["source"] = filepath.stem
                df["replicate"] = f"rep{i+1}"
                all_eccdna.append(df)
                logger.info(f"  Loaded {len(df)} eccDNAs from {filepath.name}")
            except Exception as e:
                logger.error(f"Error loading {filepath}: {e}")
                continue

        if not all_eccdna:
            raise ValueError("No eccDNA files could be loaded!")

        eccdna_df = pd.concat(all_eccdna, ignore_index=True)

        # Check required columns
        required_cols = ["eChr", "eStart", "eEnd", "eName"]
        missing_cols = [col for col in required_cols if col not in eccdna_df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")

        # Calculate eLength if missing
        if "eLength" not in eccdna_df.columns:
            eccdna_df["eLength"] = eccdna_df["eEnd"] - eccdna_df["eStart"]

        logger.info(f"Total eccDNA before filtering: {len(eccdna_df)}")

        # Filter mitochondrial DNA
        eccdna_df = eccdna_df[~eccdna_df["eChr"].str.contains("chrM")]

        # Filter sex chromosomes if requested
        if not self.keep_sex_chr:
            eccdna_df = eccdna_df[~eccdna_df["eChr"].str.contains("chrX|chrY")]
            logger.info(f"Filtered out sex chromosomes (chrX, chrY)")

        # Filter other non-standard chromosomes
        eccdna_df = eccdna_df[~eccdna_df["eChr"].str.contains("_")]

        logger.info(f"Total eccDNA after filtering: {len(eccdna_df)}")

        # Save BED file
        eccdna_bed = eccdna_df[["eChr", "eStart", "eEnd", "eName"]].copy()
        eccdna_bed_file = self.temp_dir / "eccdna_all.bed"
        eccdna_bed.to_csv(eccdna_bed_file, sep="\t", index=False, header=False)

        self.eccdna_data = eccdna_bed
        return eccdna_df, eccdna_bed_file

    def load_tad_boundaries(self):
        """Load TAD boundary data"""
        logger.info("Loading TAD boundary data...")

        if not Path(self.tad_file).exists():
            raise FileNotFoundError(f"TAD file not found: {self.tad_file}")

        try:
            tad_df = pd.read_csv(
                self.tad_file,
                sep="\t",
                header=None,
                names=["chr", "start", "end", "boundary_type", "score"],
            )
        except Exception as e:
            logger.error(f"Error loading TAD file: {e}")
            raise

        # Filter chromosomes
        tad_df = tad_df[~tad_df["chr"].str.contains("chrM")]

        # Filter sex chromosomes if requested
        if not self.keep_sex_chr:
            tad_df = tad_df[~tad_df["chr"].str.contains("chrX|chrY")]

        # Filter other non-standard chromosomes
        tad_df = tad_df[~tad_df["chr"].str.contains("_")]

        logger.info(f"Total TAD boundaries: {len(tad_df)}")
        logger.info(
            f"Strong boundaries: {len(tad_df[tad_df['boundary_type'] == 'Strong'])}"
        )
        logger.info(
            f"Weak boundaries: {len(tad_df[tad_df['boundary_type'] == 'Weak'])}"
        )

        return tad_df

    def create_genome_file(self):
        """Create genome file for bedtools shuffle (hg38)"""
        self.genome_sizes = {
            "chr1": 248956422,
            "chr2": 242193529,
            "chr3": 198295559,
            "chr4": 190214555,
            "chr5": 181538259,
            "chr6": 170805979,
            "chr7": 159345973,
            "chr8": 145138636,
            "chr9": 138394717,
            "chr10": 133797422,
            "chr11": 135086622,
            "chr12": 133275309,
            "chr13": 114364328,
            "chr14": 107043718,
            "chr15": 101991189,
            "chr16": 90338345,
            "chr17": 83257441,
            "chr18": 80373285,
            "chr19": 58617616,
            "chr20": 64444167,
            "chr21": 46709983,
            "chr22": 50818468,
            "chrX": 156040895,
            "chrY": 57227415,
        }

        genome_file = self.temp_dir / "hg38.genome"
        with open(genome_file, "w") as f:
            for chrom, size in self.genome_sizes.items():
                # Skip sex chromosomes if not keeping them
                if not self.keep_sex_chr and chrom in ["chrX", "chrY"]:
                    continue
                f.write(f"{chrom}\t{size}\n")

        return genome_file

    def create_tad_windows(self, tad_df, buffer_size, boundary_type="All"):
        """Create windows around TAD boundaries"""
        if boundary_type != "All":
            tad_subset = tad_df[tad_df["boundary_type"] == boundary_type].copy()
        else:
            tad_subset = tad_df.copy()

        tad_windows = tad_subset[["chr", "start", "end"]].copy()
        tad_windows["start"] = tad_windows["start"] - buffer_size
        tad_windows["end"] = tad_windows["end"] + buffer_size
        tad_windows["start"] = tad_windows["start"].clip(lower=0)

        window_file = self.temp_dir / f"tad_windows_{boundary_type}_{buffer_size}.bed"
        tad_windows.to_csv(window_file, sep="\t", index=False, header=False)

        return window_file

    def run_bedtools_intersect(self, a_file, b_file):
        """Run bedtools intersect and count overlaps"""
        cmd = f"bedtools intersect -a {a_file} -b {b_file} -wa | sort -u | wc -l"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Bedtools intersect error: {result.stderr}")
            return 0
        return int(result.stdout.strip())

    def single_shuffle_intersect(self, args):
        """Single shuffle and intersect operation for parallel processing"""
        i, eccdna_file, target_file, genome_file, temp_base = args

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".bed", delete=False, dir=temp_base
        ) as tmp:
            shuffled_file = tmp.name

        try:
            shuffle_cmd = f"bedtools shuffle -i {eccdna_file} -g {genome_file} -chrom -seed {i} > {shuffled_file} 2>/dev/null"
            result = subprocess.run(shuffle_cmd, shell=True)
            if result.returncode != 0:
                logger.error(f"Shuffle failed for iteration {i}")
                return 0

            count = self.run_bedtools_intersect(shuffled_file, target_file)
        finally:
            if os.path.exists(shuffled_file):
                os.unlink(shuffled_file)

        return count

    def parallel_shuffle_and_intersect(
        self, eccdna_file, target_file, genome_file, n_shuffles=1000
    ):
        """Perform parallel permutation test using bedtools shuffle"""
        logger.info(
            f"Running {n_shuffles} permutations in parallel using {self.n_cores} cores..."
        )

        with tempfile.TemporaryDirectory(dir=self.temp_dir) as temp_base:
            args_list = [
                (i, eccdna_file, target_file, genome_file, temp_base)
                for i in range(n_shuffles)
            ]

            start_time = time.time()
            with Pool(processes=self.n_cores) as pool:
                chunksize = max(1, n_shuffles // (self.n_cores * 4))
                null_counts = pool.map(
                    self.single_shuffle_intersect, args_list, chunksize=chunksize
                )

            elapsed = time.time() - start_time
            logger.info(
                f"  Completed {n_shuffles} permutations in {elapsed:.1f} seconds"
            )

        return np.array(null_counts)

    def calculate_enrichment_stats(self, observed, null_distribution):
        """Calculate enrichment statistics"""
        null_mean = np.mean(null_distribution)
        null_sd = np.std(null_distribution)

        if null_sd == 0:
            z_score = np.inf if observed > null_mean else -np.inf
        else:
            z_score = (observed - null_mean) / null_sd

        fold_enrichment = observed / null_mean if null_mean > 0 else np.inf

        # Calculate empirical p-value
        p_value = np.sum(null_distribution >= observed) / len(null_distribution)
        if p_value == 0:
            p_value = 1 / (len(null_distribution) + 1)

        # Calculate percentile
        percentile = np.sum(null_distribution < observed) / len(null_distribution) * 100

        return {
            "observed": observed,
            "null_mean": null_mean,
            "null_sd": null_sd,
            "fold_enrichment": fold_enrichment,
            "z_score": z_score,
            "p_value": p_value,
            "percentile": percentile,
        }

    def analyze_tad_enrichment(self, eccdna_file, tad_df, n_shuffles=1000):
        """Analyze TAD boundary enrichment with different parameters"""
        logger.info("\nAnalyzing TAD boundary enrichment...")

        genome_file = (
            self.genome_file if self.genome_file else self.create_genome_file()
        )
        results = []
        all_null_distributions = {}

        boundary_types = ["Strong", "Weak", "All"]
        buffer_sizes = [25000, 50000, 100000]

        for boundary_type in boundary_types:
            for buffer_size in buffer_sizes:
                logger.info(
                    f"\nProcessing {boundary_type} boundaries with ±{buffer_size/1000:.0f}kb buffer..."
                )

                window_file = self.create_tad_windows(
                    tad_df, buffer_size, boundary_type
                )
                observed = self.run_bedtools_intersect(eccdna_file, window_file)
                logger.info(f"  Observed count: {observed}")

                null_counts = self.parallel_shuffle_and_intersect(
                    eccdna_file, window_file, genome_file, n_shuffles
                )

                stats = self.calculate_enrichment_stats(observed, null_counts)
                stats["boundary_type"] = boundary_type
                stats["buffer_size"] = buffer_size
                results.append(stats)

                # Store null distribution for export
                key = f"TAD_{boundary_type}_{buffer_size}"
                all_null_distributions[key] = {
                    "observed": observed,
                    "null_counts": null_counts,
                }

                self.plot_null_distribution(null_counts, observed, key)

        results_df = pd.DataFrame(results)
        _, results_df["q_value"] = multipletests(
            results_df["p_value"], method="fdr_bh"
        )[:2]

        # Save with prefix
        output_file = "TAD_enrichment_results.csv"
        if self.prefix:
            output_file = f"{self.prefix}_{output_file}"
        results_df.to_csv(self.results_dir / output_file, index=False)

        # Export null distributions
        self.export_null_distributions(all_null_distributions)

        self.plot_tad_sensitivity(results_df)
        self.plot_enrichment_heatmap(results_df)

        return results_df

    def export_null_distributions(self, null_distributions):
        """Export null distributions for external analysis"""
        export_data = {}
        for key, data in null_distributions.items():
            export_data[f"{key}_observed"] = [data["observed"]]
            export_data[f"{key}_null"] = data["null_counts"].tolist()

        # Find max length
        max_len = max(len(v) for v in export_data.values())

        # Pad shorter arrays
        for key in export_data:
            if len(export_data[key]) < max_len:
                export_data[key].extend([np.nan] * (max_len - len(export_data[key])))

        df = pd.DataFrame(export_data)

        # Apply prefix
        filename = "TAD_null_distributions.csv"
        if self.prefix:
            filename = f"{self.prefix}_{filename}"
        df.to_csv(self.plot_data_dir / filename, index=False)
        logger.info(f"Exported null distributions to {filename}")

    def plot_null_distribution(self, null_counts, observed, title):
        """Plot null distribution with observed value"""
        plt.figure(figsize=(8, 6))

        plt.hist(
            null_counts,
            bins=30,
            density=True,
            alpha=0.7,
            color="skyblue",
            edgecolor="black",
            label="Null distribution",
        )
        plt.axvline(
            observed,
            color="red",
            linestyle="--",
            linewidth=2,
            label=f"Observed ({observed})",
        )

        null_mean = np.mean(null_counts)
        fold_enrich = observed / null_mean if null_mean > 0 else np.inf
        p_value = np.sum(null_counts >= observed) / len(null_counts)

        plt.text(
            0.95,
            0.95,
            f"Fold enrichment: {fold_enrich:.2f}\np-value: {p_value:.4f}",
            transform=plt.gca().transAxes,
            ha="right",
            va="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        )

        plt.xlabel("Number of intersections")
        plt.ylabel("Density")
        plt.title(f"Permutation Test: {title}")
        plt.legend()
        plt.tight_layout()

        safe_title = title.replace(" ", "_").replace("/", "_")
        # Apply prefix to figure name
        if self.prefix:
            safe_title = f"{self.prefix}_{safe_title}"
        plt.savefig(self.figures_dir / f"null_dist_{safe_title}.png", dpi=300)
        plt.close()

    def plot_tad_sensitivity(self, results_df):
        """Plot TAD boundary sensitivity analysis"""
        plt.figure(figsize=(10, 8))

        buffer_sizes = [25, 50, 100]  # in kb
        colors = {"Strong": "red", "Weak": "blue", "All": "green"}

        # Export data for external plotting with prefix
        plot_data = results_df[
            [
                "boundary_type",
                "buffer_size",
                "fold_enrichment",
                "null_mean",
                "null_sd",
                "p_value",
                "q_value",
            ]
        ].copy()
        plot_data["buffer_kb"] = plot_data["buffer_size"] / 1000

        output_file = "TAD_sensitivity_data.csv"
        if self.prefix:
            output_file = f"{self.prefix}_{output_file}"
        plot_data.to_csv(self.plot_data_dir / output_file, index=False)

        for boundary_type in ["Strong", "Weak", "All"]:
            subset = results_df[results_df["boundary_type"] == boundary_type]
            fold_enrichments = subset["fold_enrichment"].values

            ci_lower = []
            ci_upper = []

            for i, fe in enumerate(fold_enrichments):
                se = subset.iloc[i]["null_sd"] / subset.iloc[i]["null_mean"]
                ci_lower.append(max(0, fe - 1.96 * se * fe))
                ci_upper.append(fe + 1.96 * se * fe)

            plt.errorbar(
                buffer_sizes,
                fold_enrichments,
                yerr=[fold_enrichments - ci_lower, ci_upper - fold_enrichments],
                marker="o",
                markersize=8,
                linewidth=2,
                capsize=5,
                label=boundary_type,
                color=colors[boundary_type],
            )

        plt.xlabel("Buffer size (kb)")
        plt.ylabel("Fold enrichment")
        plt.title("TAD Boundary Enrichment Sensitivity Analysis")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        # Apply prefix to figure name
        output_file = "TAD_sensitivity_analysis.png"
        if self.prefix:
            output_file = f"{self.prefix}_{output_file}"
        plt.savefig(self.figures_dir / output_file, dpi=300)
        plt.close()

    def plot_enrichment_heatmap(self, results_df):
        """Create enrichment heatmap"""
        plt.figure(figsize=(10, 8))

        # Prepare data for heatmap
        heatmap_data = []
        labels = []
        p_values = []

        for _, row in results_df.iterrows():
            heatmap_data.append(row["fold_enrichment"])
            labels.append(f"{row['boundary_type']} - {int(row['buffer_size']/1000)}kb")
            p_values.append(row["q_value"])

        # Export heatmap data with prefix
        heatmap_df = pd.DataFrame(
            {"category": labels, "fold_enrichment": heatmap_data, "q_value": p_values}
        )
        output_file = "TAD_enrichment_summary_data.csv"
        if self.prefix:
            output_file = f"{self.prefix}_{output_file}"
        heatmap_df.to_csv(self.plot_data_dir / output_file, index=False)

        # Create heatmap
        data_matrix = np.array(heatmap_data).reshape(-1, 1)

        # Create significance annotations
        annot_data = []
        for i, (fe, pval) in enumerate(zip(heatmap_data, p_values)):
            if pval < 0.001:
                annot_data.append(f"{fe:.2f}***")
            elif pval < 0.01:
                annot_data.append(f"{fe:.2f}**")
            elif pval < 0.05:
                annot_data.append(f"{fe:.2f}*")
            else:
                annot_data.append(f"{fe:.2f}")

        annot_matrix = np.array(annot_data).reshape(-1, 1)

        sns.heatmap(
            data_matrix,
            annot=annot_matrix,
            fmt="",
            cmap="RdYlBu_r",
            yticklabels=labels,
            xticklabels=["Fold Enrichment"],
            cbar_kws={"label": "Fold Enrichment"},
        )

        plt.title("TAD Boundary Enrichment Summary\n(* q<0.05, ** q<0.01, *** q<0.001)")
        plt.tight_layout()

        # Apply prefix to figure name
        output_file = "TAD_enrichment_heatmap.png"
        if self.prefix:
            output_file = f"{self.prefix}_{output_file}"
        plt.savefig(self.figures_dir / output_file, dpi=300)
        plt.close()

    def plot_quality_control(self, eccdna_df):
        """Create quality control plots"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # eccDNA chromosome distribution
        chr_counts = eccdna_df["eChr"].value_counts().sort_index()
        chr_order = [f"chr{i}" for i in range(1, 23)]
        if self.keep_sex_chr:
            chr_order += ["chrX", "chrY"]
        chr_order = [c for c in chr_order if c in chr_counts.index]

        # Export chromosome distribution data with prefix
        chr_dist_data = pd.DataFrame(
            {
                "chromosome": chr_order,
                "count": [chr_counts.get(c, 0) for c in chr_order],
            }
        )
        output_file = "chromosome_distribution.csv"
        if self.prefix:
            output_file = f"{self.prefix}_{output_file}"
        chr_dist_data.to_csv(self.plot_data_dir / output_file, index=False)

        ax1.bar(range(len(chr_order)), [chr_counts.get(c, 0) for c in chr_order])
        ax1.set_xticks(range(len(chr_order)))
        ax1.set_xticklabels(chr_order, rotation=45, ha="right")
        ax1.set_xlabel("Chromosome")
        ax1.set_ylabel("eccDNA count")
        ax1.set_title("eccDNA Distribution Across Chromosomes")

        # eccDNA size distribution
        if "eLength" in eccdna_df.columns:
            # Export size distribution data with prefix
            size_data = pd.DataFrame({"eLength": eccdna_df["eLength"]})
            output_file = "size_distribution.csv"
            if self.prefix:
                output_file = f"{self.prefix}_{output_file}"
            size_data.to_csv(self.plot_data_dir / output_file, index=False)

            ax2.hist(eccdna_df["eLength"], bins=50, edgecolor="black", alpha=0.7)
            ax2.set_xlabel("eccDNA length (bp)")
            ax2.set_ylabel("Count")
            ax2.set_title("eccDNA Size Distribution")
            ax2.set_yscale("log")

            median_size = eccdna_df["eLength"].median()
            ax2.axvline(
                median_size,
                color="red",
                linestyle="--",
                label=f"Median: {median_size:.0f} bp",
            )
            ax2.legend()
        else:
            ax2.text(
                0.5,
                0.5,
                "eLength column not found",
                transform=ax2.transAxes,
                ha="center",
                va="center",
            )
            ax2.set_title("eccDNA Size Distribution (N/A)")

        plt.tight_layout()

        # Apply prefix to figure name
        output_file = "quality_control.png"
        if self.prefix:
            output_file = f"{self.prefix}_{output_file}"
        plt.savefig(self.figures_dir / output_file, dpi=300)
        plt.close()

    def generate_report(self, results_df):
        """Generate markdown report"""
        timestamp = pd.Timestamp.now().strftime("%Y-%m-%d %H:%M")
        total_eccdna = len(
            pd.read_csv(self.temp_dir / "eccdna_all.bed", sep="\t", header=None)
        )

        report = f"""# eccTAD Analysis Report

## Analysis Overview
- **Date**: {timestamp}
- **Tool**: eccTAD v1.0
- **Working Directory**: {self.work_dir}"""

        if self.prefix:
            report += f"\n- **Output Prefix**: {self.prefix}"

        report += f"""
- **CPU cores used**: {self.n_cores}
- **Permutations**: {self.n_shuffles}

## Input Files
- **eccDNA files**: {', '.join([Path(f).name for f in self.eccdna_files])}
- **TAD boundary file**: {Path(self.tad_file).name}

## Data Summary
- **Total eccDNA regions analyzed**: {total_eccdna:,}
- **TAD boundary types**: Strong, Weak, and All combined
- **Buffer sizes tested**: 25kb, 50kb, 100kb

## Key Findings

### TAD Boundary Enrichment Results
"""
        sig_count = sum(results_df["q_value"] < 0.05)
        report += f"\n**Significant enrichments (q < 0.05)**: {sig_count} out of {len(results_df)}\n\n"

        for _, row in results_df.iterrows():
            significance = (
                "***"
                if row["q_value"] < 0.001
                else (
                    "**"
                    if row["q_value"] < 0.01
                    else "*" if row["q_value"] < 0.05 else ""
                )
            )
            report += f"- **{row['boundary_type']} boundaries (±{row['buffer_size']/1000:.0f}kb)**: "
            report += f"{row['fold_enrichment']:.2f}-fold enrichment "
            report += (
                f"(p={row['p_value']:.4f}, q={row['q_value']:.4f}){significance}\n"
            )

        # Summary statistics
        max_enrichment = results_df.loc[results_df["fold_enrichment"].idxmax()]
        report += f"\n### Summary Statistics\n"
        report += (
            f"- **Highest enrichment**: {max_enrichment['fold_enrichment']:.2f}-fold "
        )
        report += f"({max_enrichment['boundary_type']} boundaries, ±{max_enrichment['buffer_size']/1000:.0f}kb)\n"

        report += f"""
## Statistical Methods
- **Permutation test**: {self.n_shuffles} iterations using bedtools shuffle
- **Shuffle parameters**: -chrom (maintains chromosome distribution)
- **Multiple testing correction**: Benjamini-Hochberg FDR
- **Overlap criterion**: ≥1bp
- **Parallel processing**: Enabled ({self.n_cores} cores)
- **Random seed**: 42 (for reproducibility)

## Output Files
"""

        if self.prefix:
            report += f"\nAll output files are prefixed with: **{self.prefix}_**\n"

        report += """
### Results
- `*TAD_enrichment_results.csv`: Enrichment statistics for all conditions
- `*TAD_analysis_report.md`: This report

### Figures
- `*quality_control.png`: eccDNA distribution and size plots
- `*TAD_sensitivity_analysis.png`: Enrichment vs buffer size
- `*TAD_enrichment_heatmap.png`: Summary heatmap with significance
- `*null_dist_TAD_*.png`: Permutation test distributions

### Plot Data (for external visualization)
- `*chromosome_distribution.csv`: eccDNA counts per chromosome
- `*size_distribution.csv`: eccDNA size data
- `*TAD_sensitivity_data.csv`: Enrichment vs buffer size data
- `*TAD_enrichment_summary_data.csv`: All enrichment values
- `*TAD_null_distributions.csv`: Null distribution data

## Interpretation Guidelines

### P-value significance levels:
- *** : q < 0.001 (highly significant)
- **  : q < 0.01 (very significant)
- *   : q < 0.05 (significant)
- no asterisk : q ≥ 0.05 (not significant)

### Fold enrichment interpretation:
- > 2.0 : Strong enrichment
- 1.5-2.0 : Moderate enrichment
- 1.2-1.5 : Weak enrichment
- 0.8-1.2 : No enrichment/depletion
- < 0.8 : Depletion

## Software Requirements
- Python 3.7+
- bedtools 2.30.0+
- Required Python packages: pandas, numpy, matplotlib, seaborn, scipy, statsmodels

## Citation
If you use this analysis in your research, please cite:
- Bedtools: Quinlan AR and Hall IM, 2010. Bioinformatics.
- eccTAD: TAD Boundary Enrichment Analysis Pipeline v1.0

---
*Report generated automatically by eccTAD*
"""

        # Apply prefix to report filename
        output_file = "TAD_analysis_report.md"
        if self.prefix:
            output_file = f"{self.prefix}_{output_file}"

        with open(self.results_dir / output_file, "w") as f:
            f.write(report)

        report_path = self.results_dir / output_file
        logger.info(f"\nReport saved to: {report_path}")

    def run_analysis(self):
        """Run complete analysis pipeline"""
        print("=" * 60)
        print("eccTAD Analysis Pipeline v1.0")
        print("=" * 60)

        start_total = time.time()

        try:
            # Load data
            eccdna_df, eccdna_file = self.load_eccdna_data()
            tad_df = self.load_tad_boundaries()

            # Quality control plots
            self.plot_quality_control(eccdna_df)

            # Run analysis
            results_df = self.analyze_tad_enrichment(
                eccdna_file, tad_df, self.n_shuffles
            )

            # Generate report
            self.generate_report(results_df)

            total_time = time.time() - start_total
            print("\n" + "=" * 60)
            print(f"Analysis complete in {total_time:.1f} seconds!")
            print("Check directories for outputs:")
            print(f"  - Results: {self.results_dir}")
            print(f"  - Figures: {self.figures_dir}")
            print(f"  - Plot data: {self.plot_data_dir}")
            print("=" * 60)

        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            raise


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="eccTAD - eccDNA TAD Boundary Enrichment Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis
  %(prog)s -i sample1.csv sample2.csv -t TAD_boundaries.bed
  
  # Keep sex chromosomes in analysis
  %(prog)s -i eccDNA.csv -t TAD.bed --keep-sex-chr
  
  # With custom output prefix
  %(prog)s -i eccDNA.csv -t TAD.bed --prefix HeLa_rep1
  
  # Custom genome file and more permutations
  %(prog)s -i eccDNA.csv -t TAD.bed -g hg19.genome -n 10000
  
  # Specify output directory, CPU cores, and prefix
  %(prog)s -i eccDNA.csv -t TAD.bed -w ./output --cores 32 --prefix experiment_001
        """,
    )

    # Required arguments
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Input eccDNA CSV file(s). Multiple files will be merged.",
    )
    parser.add_argument(
        "-t",
        "--tad",
        required=True,
        help="TAD boundary BED file (format: chr start end boundary_type score)",
    )

    # Optional arguments
    parser.add_argument(
        "-g",
        "--genome",
        default=None,
        help="Genome file for bedtools shuffle (default: built-in hg38)",
    )
    parser.add_argument(
        "-w",
        "--workdir",
        default=".",
        help="Working directory for output (default: current directory)",
    )
    parser.add_argument(
        "-n",
        "--nshuffle",
        default=1000,
        type=int,
        help="Number of permutations for enrichment testing (default: 1000)",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=None,
        help="Number of CPU cores to use (default: auto-detect)",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default=None,
        help="Prefix for all output files and directories (e.g., HeLa_rep1)",
    )
    parser.add_argument(
        "--keep-sex-chr",
        action="store_true",
        help="Keep sex chromosomes (chrX, chrY) in analysis (default: filter out)",
    )
    parser.add_argument("--version", action="version", version="%(prog)s 1.0")

    return parser.parse_args()


def check_dependencies():
    """Check if required dependencies are installed"""
    # Check bedtools
    try:
        result = subprocess.run(
            ["bedtools", "--version"], capture_output=True, text=True
        )
        bedtools_version = result.stdout.strip()
        logger.info(f"Found {bedtools_version}")
    except Exception as e:
        logger.error(
            "bedtools not found. Please install bedtools (version 2.30.0 or higher)"
        )
        logger.error("Installation: conda install -c bioconda bedtools")
        sys.exit(1)

    # Check Python packages
    required_packages = [
        "pandas",
        "numpy",
        "matplotlib",
        "seaborn",
        "scipy",
        "statsmodels",
    ]
    missing_packages = []

    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing_packages.append(package)

    if missing_packages:
        logger.error(f"Missing Python packages: {', '.join(missing_packages)}")
        logger.error(f"Installation: pip install {' '.join(missing_packages)}")
        sys.exit(1)


def main():
    """Main function"""
    args = parse_args()

    # Check dependencies
    check_dependencies()

    # Create analyzer and run
    analyzer = eccTADAnalysis(args)
    analyzer.run_analysis()


if __name__ == "__main__":
    main()
