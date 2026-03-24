"""Saturation curve analysis for eccDNA detection."""

import logging
import os
import subprocess
import time
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from ecctoolkit.utils.subprocess_utils import require_tools, run_command
from ecctoolkit.utils.validation import validate_file_exists

logger = logging.getLogger(__name__)


def _exp_saturation(x: np.ndarray, n: float, c: float) -> np.ndarray:
    """Exponential saturation model: y = N * (1 - exp(-c * x))."""
    return n * (1 - np.exp(-c * x))


def _downsample(
    input_file: str,
    output_file: str,
    fraction: float,
    seed: int,
) -> None:
    """Downsample a FASTA/FASTQ file using seqkit sample."""
    cmd = [
        "seqkit", "sample",
        "-p", str(fraction),
        "-s", str(seed),
        input_file,
    ]
    logger.info(f"Downsampling to {fraction*100:.0f}%: {output_file}")
    with open(output_file, "w") as fout:
        subprocess.run(cmd, check=True, stdout=fout, stderr=subprocess.PIPE)


def _count_eccdna(output_dir: str, prefix: str) -> int:
    """Count eccDNA from detector output files."""
    # Try common output file patterns
    patterns = [
        os.path.join(output_dir, f"{prefix}*summary*.csv"),
        os.path.join(output_dir, f"{prefix}*.bed"),
        os.path.join(output_dir, f"{prefix}*eccDNA*.txt"),
    ]

    import glob
    for pattern in patterns:
        matches = glob.glob(pattern)
        if matches:
            result_file = matches[0]
            with open(result_file) as f:
                # Count non-header lines
                lines = [line for line in f if line.strip() and not line.startswith("#")]
            # Subtract header if CSV
            if result_file.endswith(".csv"):
                return max(0, len(lines) - 1)
            return len(lines)

    logger.warning(f"No eccDNA result files found for prefix {prefix}")
    return 0


def _run_detector(
    input_file: str,
    reference: str,
    output_dir: str,
    prefix: str,
    threads: int,
    detector: str,
) -> None:
    """Run eccDNA detection tool on input file."""
    if detector == "circleseeker":
        cmd = [
            "CircleSeeker",
            "-i", input_file,
            "-p", prefix,
            "-r", reference,
            "-o", output_dir,
            "-t", str(threads),
            "--enable_X",
        ]
    elif detector == "circlemap":
        raise ValueError(
            "CircleMap requires paired-end reads. "
            "Use circleseeker for saturation analysis on single-end data."
        )
    else:
        raise ValueError(f"Unsupported detector: {detector}")

    logger.info(f"Running {detector}: {' '.join(cmd)}")
    run_command(cmd, shell=False, check=True, capture_output=False)


def _fit_saturation_model(
    fractions: np.ndarray,
    counts: np.ndarray,
) -> Optional[dict]:
    """Fit exponential saturation model and return estimated parameters."""
    if len(fractions) < 3 or np.all(counts == 0):
        return None

    try:
        p0 = [counts[-1] * 1.5, 2.0]
        popt, pcov = curve_fit(
            _exp_saturation, fractions, counts,
            p0=p0, maxfev=10000,
        )
        n_est, c_est = popt
        perr = np.sqrt(np.diag(pcov))

        residuals = counts - _exp_saturation(fractions, *popt)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((counts - np.mean(counts)) ** 2)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

        return {
            "N_estimated": n_est,
            "N_stderr": perr[0],
            "c_estimated": c_est,
            "c_stderr": perr[1],
            "r_squared": r_squared,
            "saturation_pct": counts[-1] / n_est * 100 if n_est > 0 else 0.0,
        }
    except (RuntimeError, ValueError) as e:
        logger.warning(f"Failed to fit saturation model: {e}")
        return None


def run_saturation(
    input_file: str,
    reference: str,
    output_dir: str,
    prefix: str,
    threads: int = 8,
    fractions: Optional[List[float]] = None,
    seed: int = 42,
    detector: str = "circleseeker",
) -> None:
    """
    Generate saturation curve data for eccDNA detection.

    Downsamples input reads at multiple fractions, runs eccDNA detection
    on each, and produces a saturation curve CSV. Optionally fits an
    exponential saturation model to estimate the total eccDNA pool size.

    Args:
        input_file: Input FASTA/FASTQ file
        reference: Reference genome FASTA
        output_dir: Output directory
        prefix: Sample prefix
        threads: Number of threads
        fractions: Sampling fractions (default: 0.1 to 1.0 by 0.1)
        seed: Random seed for reproducibility
        detector: Detection tool to use (default: circleseeker)
    """
    validate_file_exists(input_file, "Input file")
    validate_file_exists(reference, "Reference genome")

    tools = ["seqkit"]
    if detector == "circleseeker":
        tools.append("CircleSeeker")
    require_tools(tools)

    os.makedirs(output_dir, exist_ok=True)

    if fractions is None:
        fractions = [p / 100.0 for p in range(10, 110, 10)]

    logger.info("=== Saturation Curve Analysis Started ===")
    logger.info(f"Input: {input_file}")
    logger.info(f"Reference: {reference}")
    logger.info(f"Prefix: {prefix}")
    logger.info(f"Fractions: {fractions}")
    logger.info(f"Detector: {detector}")

    results = []

    for frac in fractions:
        label = int(frac * 100)
        start_time = time.time()

        if frac >= 1.0:
            sampled_file = input_file
        else:
            sampled_file = os.path.join(
                output_dir, f"{prefix}_{label}pct.fasta"
            )
            _downsample(input_file, sampled_file, frac, seed)

        # Count input reads
        count_cmd = ["seqkit", "stats", "-T", sampled_file]
        count_result = subprocess.run(
            count_cmd, capture_output=True, text=True
        )
        n_reads = 0
        if count_result.returncode == 0:
            lines = count_result.stdout.strip().split("\n")
            if len(lines) >= 2:
                fields = lines[1].split("\t")
                if len(fields) >= 4:
                    n_reads = int(fields[3].replace(",", ""))

        # Run detector
        sub_prefix = f"{prefix}_{label}pct"
        sub_outdir = os.path.join(output_dir, sub_prefix)
        os.makedirs(sub_outdir, exist_ok=True)

        _run_detector(
            sampled_file, reference, sub_outdir,
            sub_prefix, threads, detector,
        )

        n_eccdna = _count_eccdna(sub_outdir, sub_prefix)
        elapsed = time.time() - start_time

        results.append({
            "fraction": frac,
            "n_reads": n_reads,
            "n_eccdna": n_eccdna,
            "time_seconds": round(elapsed, 2),
        })

        logger.info(
            f"  {label}%: {n_reads:,} reads -> {n_eccdna:,} eccDNA "
            f"({elapsed:.1f}s)"
        )

    # Save saturation curve data
    df = pd.DataFrame(results)
    csv_path = os.path.join(output_dir, f"{prefix}_saturation_curve.csv")
    df.to_csv(csv_path, index=False)
    logger.info(f"Saturation curve saved: {csv_path}")

    # Fit saturation model
    frac_arr = np.array(df["fraction"])
    count_arr = np.array(df["n_eccdna"], dtype=float)

    fit_result = _fit_saturation_model(frac_arr, count_arr)
    if fit_result is not None:
        logger.info("Exponential saturation model fit:")
        logger.info(
            f"  N_total = {fit_result['N_estimated']:,.0f} "
            f"+/- {fit_result['N_stderr']:,.0f}"
        )
        logger.info(f"  Capture rate c = {fit_result['c_estimated']:.4f}")
        logger.info(f"  R^2 = {fit_result['r_squared']:.6f}")
        logger.info(f"  Saturation = {fit_result['saturation_pct']:.1f}%")

        # Save model parameters
        model_path = os.path.join(
            output_dir, f"{prefix}_saturation_model.csv"
        )
        pd.DataFrame([fit_result]).to_csv(model_path, index=False)
        logger.info(f"Model parameters saved: {model_path}")

    logger.info("=== Saturation Curve Analysis Finished ===")
