"""CReSIL eccDNA detection pipeline for long-read data."""

import logging
import os
import subprocess
from pathlib import Path
from typing import Optional

from ecctoolkit.utils.subprocess_utils import require_tools, run_command
from ecctoolkit.utils.validation import validate_file_exists

logger = logging.getLogger(__name__)


def _run_pipeline_command(cmd: list, log_file: str) -> None:
    """Execute a pipeline command with logging to file."""
    logger.info(f"Executing: {' '.join(cmd)}")
    with open(log_file, "w") as log_handle:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        stdout, stderr = process.communicate()

        log_handle.write(stdout)
        if stderr:
            log_handle.write(stderr)

        if process.returncode != 0:
            logger.error(f"Command failed (exit {process.returncode}): {' '.join(cmd)}")
            if stderr:
                logger.error(f"stderr: {stderr[:500]}")
            raise subprocess.CalledProcessError(process.returncode, cmd)


def run_cresil_pipeline(
    input_fastq: str,
    reference_mmi: str,
    reference_fasta: str,
    output_dir: str,
    rmsk_bed: Optional[str] = None,
    cpg_bed: Optional[str] = None,
    gene_bed: Optional[str] = None,
    threads: int = 8,
) -> None:
    """
    Run CReSIL pipeline for chimeric eccDNA detection from long reads.

    Executes the CReSIL workflow: trim -> identify -> annotate (optional).

    Args:
        input_fastq: Path to input FASTQ file (long reads)
        reference_mmi: Path to minimap2 index (.mmi)
        reference_fasta: Path to reference FASTA file
        output_dir: Output directory
        rmsk_bed: Path to RepeatMasker BED file (for annotation)
        cpg_bed: Path to CpG islands BED file (for annotation)
        gene_bed: Path to gene annotation BED file (for annotation)
        threads: Number of threads
    """
    # Validate required inputs
    validate_file_exists(input_fastq, "Input FASTQ")
    validate_file_exists(reference_mmi, "Reference MMI")
    validate_file_exists(reference_fasta, "Reference FASTA")

    # Check for reference .fai index
    reference_fai = reference_fasta + ".fai"
    validate_file_exists(reference_fai, "Reference FAI index")

    # Validate optional annotation files if provided
    if rmsk_bed:
        validate_file_exists(rmsk_bed, "RepeatMasker BED")
    if cpg_bed:
        validate_file_exists(cpg_bed, "CpG BED")
    if gene_bed:
        validate_file_exists(gene_bed, "Gene BED")

    # Check tools
    require_tools(["cresil"])

    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)

    logger.info("=== CReSIL Pipeline Started ===")
    logger.info(f"Input: {input_fastq}")
    logger.info(f"Reference MMI: {reference_mmi}")
    logger.info(f"Reference FASTA: {reference_fasta}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Threads: {threads}")

    # Step 1: Trim
    logger.info("Step 1: Running CReSIL trim...")
    trim_output = os.path.join(output_dir, "trim_result")
    os.makedirs(trim_output, exist_ok=True)

    trim_cmd = [
        "cresil", "trim",
        "-t", str(threads),
        "-fq", input_fastq,
        "-r", reference_mmi,
        "-o", trim_output,
    ]
    _run_pipeline_command(
        trim_cmd, os.path.join(log_dir, "cresil_trim.log"),
    )

    trim_file = os.path.join(trim_output, "trim.txt")
    if not os.path.exists(trim_file):
        raise FileNotFoundError(
            f"Expected trim output not found: {trim_file}. "
            "Check CReSIL trim logs."
        )
    logger.info(f"Trim completed: {trim_file}")

    # Step 2: Identify
    logger.info("Step 2: Running CReSIL identify...")
    identify_cmd = [
        "cresil", "identify",
        "-t", str(threads),
        "-fa", reference_fasta,
        "-fai", reference_fai,
        "-fq", input_fastq,
        "-trim", trim_file,
    ]
    _run_pipeline_command(
        identify_cmd, os.path.join(log_dir, "cresil_identify.log"),
    )

    identify_file = os.path.join(output_dir, "eccDNA_final.txt")
    if not os.path.exists(identify_file):
        raise FileNotFoundError(
            f"Expected identify output not found: {identify_file}. "
            "Check CReSIL identify logs."
        )
    logger.info(f"Identify completed: {identify_file}")

    # Step 3: Annotate (optional, requires all annotation files)
    if rmsk_bed and cpg_bed and gene_bed:
        logger.info("Step 3: Running CReSIL annotate...")
        annotate_cmd = [
            "cresil", "annotate",
            "-t", str(threads),
            "-rp", rmsk_bed,
            "-cg", cpg_bed,
            "-gb", gene_bed,
            "-identify", identify_file,
        ]
        _run_pipeline_command(
            annotate_cmd, os.path.join(log_dir, "cresil_annotate.log"),
        )
        logger.info("Annotation completed.")
    else:
        logger.info(
            "Step 3: Skipping annotation (rmsk_bed, cpg_bed, and "
            "gene_bed are all required)."
        )

    logger.info("=== CReSIL Pipeline Finished ===")
    logger.info(f"Results: {output_dir}")
