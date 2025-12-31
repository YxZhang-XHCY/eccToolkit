"""Circle-Map eccDNA identification pipeline."""

import logging
from pathlib import Path

from ecctoolkit.utils.subprocess_utils import require_tools, run_command
from ecctoolkit.utils.io import create_output_dirs

logger = logging.getLogger(__name__)


def run_circlemap_pipeline(
    fastq1: str,
    fastq2: str,
    reference: str,
    output_dir: str,
    sample_name: str,
    threads: int = 8,
    verbose: bool = False,
) -> None:
    """
    Run the Circle-Map eccDNA identification pipeline.

    Args:
        fastq1: Path to FASTQ file (Read 1)
        fastq2: Path to FASTQ file (Read 2)
        reference: Path to reference genome FASTA
        output_dir: Output directory
        sample_name: Sample name prefix
        threads: Number of threads
        verbose: Enable verbose output
    """
    # Check dependencies
    require_tools(["fastp", "bwa", "samtools", "Circle-Map"])

    # Create output directories
    output_path = Path(output_dir)
    dirs = create_output_dirs(output_path, ["qc", "alignment", "circlemap", "logs"])

    logger.info(f"Starting Circle-Map pipeline for {sample_name}")
    logger.info(f"Input: {fastq1}, {fastq2}")
    logger.info(f"Reference: {reference}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Threads: {threads}")

    # Step 1: Quality control with fastp
    logger.info("Step 1: Running fastp QC...")
    clean_r1 = dirs["qc"] / f"{sample_name}_clean_R1.fq.gz"
    clean_r2 = dirs["qc"] / f"{sample_name}_clean_R2.fq.gz"
    fastp_cmd = (
        f"fastp -i {fastq1} -I {fastq2} "
        f"-o {clean_r1} -O {clean_r2} "
        f"-h {dirs['qc']}/{sample_name}_fastp.html "
        f"-j {dirs['qc']}/{sample_name}_fastp.json "
        f"-w {threads}"
    )
    run_command(fastp_cmd, check=True)

    # Step 2: BWA alignment
    logger.info("Step 2: Running BWA alignment...")
    sam_file = dirs["alignment"] / f"{sample_name}.sam"
    bwa_cmd = f"bwa mem -t {threads} {reference} {clean_r1} {clean_r2} > {sam_file}"
    run_command(bwa_cmd, check=True)

    # Step 3: Sort and convert to BAM
    logger.info("Step 3: Sorting and converting to BAM...")
    sorted_bam = dirs["alignment"] / f"{sample_name}_sorted.bam"
    sort_cmd = f"samtools sort -@ {threads} -o {sorted_bam} {sam_file}"
    run_command(sort_cmd, check=True)

    # Index BAM
    run_command(f"samtools index {sorted_bam}", check=True)

    # Step 4: Extract discordant reads
    logger.info("Step 4: Extracting discordant reads...")
    disc_bam = dirs["alignment"] / f"{sample_name}_discordant.bam"
    disc_cmd = f"samtools view -bF 2 {sorted_bam} > {disc_bam}"
    run_command(disc_cmd, check=True)
    run_command(f"samtools index {disc_bam}", check=True)

    # Step 5: Run Circle-Map ReadExtractor
    logger.info("Step 5: Running Circle-Map ReadExtractor...")
    extracted_bam = dirs["circlemap"] / f"{sample_name}_extracted.bam"
    extract_cmd = (
        f"Circle-Map ReadExtractor -i {sorted_bam} "
        f"-o {extracted_bam}"
    )
    run_command(extract_cmd, check=True)
    run_command(f"samtools sort -@ {threads} -o {extracted_bam}.sorted {extracted_bam}", check=True)
    run_command(f"mv {extracted_bam}.sorted {extracted_bam}", check=True)
    run_command(f"samtools index {extracted_bam}", check=True)

    # Step 6: Run Circle-Map Realign
    logger.info("Step 6: Running Circle-Map Realign...")
    circlemap_out = dirs["circlemap"] / f"{sample_name}_circle.bed"
    realign_cmd = (
        f"Circle-Map Realign -i {extracted_bam} "
        f"-qbam {sorted_bam} "
        f"-sbam {disc_bam} "
        f"-fasta {reference} "
        f"-o {circlemap_out}"
    )
    run_command(realign_cmd, check=True)

    # Cleanup intermediate files
    logger.info("Cleaning up intermediate files...")
    sam_file.unlink(missing_ok=True)

    logger.info(f"Pipeline completed. Results: {circlemap_out}")
