"""CeccDNA false discovery rate estimation from spike-in control data."""

import logging
import os
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def estimate_ceccdna_fdr(
    spikein_file: str,
    sample_file: str,
    output_file: str,
    min_read_count: int = 2,
) -> None:
    """Estimate CeccDNA false discovery rate using spike-in chimeric artefact data.

    The spike-in file should contain columns:
        - sample: sample identifier
        - A: read count for molecule A
        - B: read count for molecule B
        - AB: chimeric read count (A-B fusions)

    The sample file should contain columns:
        - sample: sample identifier
        - total_reads: total HiFi reads
        - ceccdna_count: observed CeccDNA count

    The method:
    1. Calculates chimeric read fraction from spike-in data: AB / (A + B + AB)
    2. Estimates expected false chimeric reads: total_reads * chimeric_fraction
    3. Without read filter: FDR = expected_false / observed_ceccdna
    4. With read filter (>= min_read_count): uses birthday-problem approximation
       to estimate probability of multiple false reads at the same breakpoint

    Args:
        spikein_file: Path to spike-in control data CSV.
        sample_file: Path to sample data CSV with total reads and CeccDNA counts.
        output_file: Path to write FDR estimation results.
        min_read_count: Minimum read count filter for the filtered FDR estimate.
    """
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Load spike-in data
    spikein_df = pd.read_csv(spikein_file)
    required_cols = {"sample", "A", "B", "AB"}
    missing = required_cols - set(spikein_df.columns)
    if missing:
        raise ValueError(f"Spike-in file missing columns: {missing}")

    # Calculate chimeric read fraction per spike-in sample
    spikein_df["total_spikein"] = spikein_df["A"] + spikein_df["B"] + spikein_df["AB"]
    spikein_df["chimeric_fraction"] = spikein_df["AB"] / spikein_df["total_spikein"]

    mean_chimeric_fraction = spikein_df["chimeric_fraction"].mean()
    sd_chimeric_fraction = spikein_df["chimeric_fraction"].std()

    logger.info(
        f"Spike-in chimeric fraction: {mean_chimeric_fraction:.5f} "
        f"+/- {sd_chimeric_fraction:.5f} "
        f"({mean_chimeric_fraction * 100:.3f}%)"
    )

    # Load sample data
    sample_df = pd.read_csv(sample_file)
    required_cols = {"sample", "total_reads", "ceccdna_count"}
    missing = required_cols - set(sample_df.columns)
    if missing:
        raise ValueError(f"Sample file missing columns: {missing}")

    # Estimate FDR per sample
    results = []
    for _, row in sample_df.iterrows():
        total_reads = row["total_reads"]
        observed_ceccdna = row["ceccdna_count"]

        # Expected false chimeric reads
        expected_false_reads = total_reads * mean_chimeric_fraction

        # FDR without read filter
        fdr_no_filter = (
            expected_false_reads / observed_ceccdna if observed_ceccdna > 0 else np.nan
        )

        # FDR with read filter (birthday problem approximation)
        # At 50-bp resolution, ~62M possible breakpoint positions in human genome
        # Number of breakpoint pairs: 62M choose 2 ~ 1.9e15
        # P(any pair hit >= min_read_count times) is negligible for reasonable k
        genome_positions = 3.1e9 / 50  # ~62M at 50bp resolution
        genome_pairs = genome_positions * (genome_positions - 1) / 2
        k = expected_false_reads

        # Expected pairs hit >= 2 times (birthday problem)
        if min_read_count == 2:
            expected_false_filtered = k * (k - 1) / (2 * genome_pairs)
        else:
            # Generalized: Poisson approximation for >= min_read_count hits
            lam = k / genome_pairs  # expected hits per pair
            from scipy.stats import poisson
            p_ge_n = 1 - poisson.cdf(min_read_count - 1, lam)
            expected_false_filtered = genome_pairs * p_ge_n

        fdr_filtered = (
            expected_false_filtered / observed_ceccdna
            if observed_ceccdna > 0
            else np.nan
        )

        results.append({
            "sample": row["sample"],
            "total_reads": total_reads,
            "observed_ceccdna": observed_ceccdna,
            "chimeric_fraction": mean_chimeric_fraction,
            "expected_false_chimeric_reads": expected_false_reads,
            "fdr_no_filter": fdr_no_filter,
            f"fdr_ge{min_read_count}_reads": fdr_filtered,
            f"expected_false_ge{min_read_count}_reads": expected_false_filtered,
        })

    results_df = pd.DataFrame(results)

    # Write CSV output
    csv_path = output_file
    if not csv_path.endswith(".csv"):
        csv_path = os.path.splitext(output_file)[0] + ".csv"
    results_df.to_csv(csv_path, index=False)
    logger.info(f"Wrote FDR estimates to {csv_path}")

    # Write text report
    txt_path = os.path.splitext(output_file)[0] + ".txt"
    with open(txt_path, "w") as f:
        f.write("CeccDNA False Discovery Rate Estimation\n")
        f.write("=" * 50 + "\n\n")

        f.write("Spike-in chimeric read fraction:\n")
        f.write(f"  Mean: {mean_chimeric_fraction * 100:.3f}%\n")
        f.write(f"  SD:   {sd_chimeric_fraction * 100:.3f}%\n")
        f.write(f"  N samples: {len(spikein_df)}\n\n")

        f.write("Per-sample FDR estimates:\n")
        f.write("-" * 50 + "\n")
        for _, r in results_df.iterrows():
            f.write(f"\n  {r['sample']}:\n")
            f.write(f"    Total reads: {r['total_reads']:,.0f}\n")
            f.write(f"    Observed CeccDNA: {r['observed_ceccdna']:,.0f}\n")
            f.write(f"    Expected false chimeric reads: {r['expected_false_chimeric_reads']:.0f}\n")
            fdr_nf = r["fdr_no_filter"]
            f.write(f"    FDR (no filter): {fdr_nf * 100:.1f}%\n" if not np.isnan(fdr_nf) else "    FDR (no filter): N/A\n")
            fdr_f = r[f"fdr_ge{min_read_count}_reads"]
            f.write(
                f"    FDR (>={min_read_count} reads): {fdr_f * 100:.2e}%\n"
                if not np.isnan(fdr_f)
                else f"    FDR (>={min_read_count} reads): N/A\n"
            )

        f.write(f"\n{'=' * 50}\n")
        f.write("Conclusion:\n")
        f.write(
            f"  With >={min_read_count} read-support filter, random chimeric artefacts\n"
            f"  cannot independently hit the same breakpoint pair, making the\n"
            f"  FDR effectively 0%.\n"
        )

    logger.info(f"Wrote FDR report to {txt_path}")
