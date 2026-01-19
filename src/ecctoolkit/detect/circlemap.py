"""Circle-Map eccDNA identification pipeline."""

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from ecctoolkit.utils.subprocess_utils import require_tools, run_command
from ecctoolkit.utils.io import create_output_dirs
from ecctoolkit.utils.validation import validate_file_exists

logger = logging.getLogger(__name__)


def _resolve_circlemap_executable() -> Optional[str]:
    """Resolve Circle-Map executable name across common installations."""
    candidates = ("Circle-Map", "Circle-Map.py", "circle-map", "circle-map.py")
    for name in candidates:
        if shutil.which(name):
            return name
    return None


def _has_bwa_index(reference: Path) -> bool:
    """Return True if all BWA index sidecar files exist for reference."""
    required_exts = (".amb", ".ann", ".bwt", ".pac", ".sa")
    return all(Path(str(reference) + ext).exists() for ext in required_exts)


def _ensure_reference_indices(
    reference: Path,
    threads: int,
    log_path: Path,
    auto_index: bool,
) -> None:
    """Ensure BWA index and FASTA index exist (optionally auto-generate)."""
    if reference.suffix == ".gz":
        raise ValueError(
            "Circle-Map pipeline requires an uncompressed reference FASTA (not .gz). "
            f"Got: {reference}"
        )

    bwa_index_ok = _has_bwa_index(reference)
    fai_path = Path(str(reference) + ".fai")
    fai_ok = fai_path.exists()

    if bwa_index_ok and fai_ok:
        return

    if not auto_index:
        missing = []
        if not bwa_index_ok:
            missing.append("BWA index (.amb/.ann/.bwt/.pac/.sa)")
        if not fai_ok:
            missing.append("FASTA index (.fai)")
        raise RuntimeError(
            "Missing reference index files: "
            + ", ".join(missing)
            + f". Run `bwa index {reference}` and/or `samtools faidx {reference}`, "
            "or rerun with --auto-index."
        )

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as log_handle:
        if not bwa_index_ok:
            logger.info("Building BWA index...")
            run_command(
                ["bwa", "index", str(reference)],
                check=True,
                shell=False,
                stdout=log_handle,
                stderr=log_handle,
            )
        if not fai_ok:
            logger.info("Building FASTA index (.fai)...")
            run_command(
                ["samtools", "faidx", str(reference)],
                check=True,
                shell=False,
                stdout=log_handle,
                stderr=log_handle,
            )


def run_circlemap_pipeline(
    fastq1: str,
    fastq2: str,
    reference: str,
    output_dir: str,
    sample_name: str,
    threads: int = 8,
    skip_fastp: bool = False,
    min_mapq: int = 0,
    auto_index: bool = False,
    keep_intermediate: bool = True,
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
        skip_fastp: Skip fastp QC/trim and use input FASTQ as-is
        min_mapq: Minimum MAPQ for discordant BAM extraction
        auto_index: Auto-generate BWA and FASTA indices if missing
        keep_intermediate: Keep intermediate BAMs (discordant, extracted, etc.)
        verbose: Enable verbose output
    """
    # Validate inputs early for better UX.
    validate_file_exists(fastq1, "FASTQ R1")
    validate_file_exists(fastq2, "FASTQ R2")
    validate_file_exists(reference, "Reference FASTA")

    threads = max(1, int(threads))
    min_mapq = max(0, int(min_mapq))

    root_logger = logging.getLogger()
    if not root_logger.handlers:
        logging.basicConfig(
            level=logging.DEBUG if verbose else logging.INFO,
            format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
        )

    circlemap_exe = _resolve_circlemap_executable()
    if circlemap_exe is None:
        raise RuntimeError(
            "Missing required tool: Circle-Map. "
            "Tried: Circle-Map, Circle-Map.py, circle-map, circle-map.py. "
            "Install Circle-Map and ensure it is in your PATH."
        )

    # Check dependencies
    tools = ["bwa", "samtools", circlemap_exe]
    if not skip_fastp:
        tools.append("fastp")
    require_tools(tools)

    # Create output directories
    output_path = Path(output_dir)
    dirs = create_output_dirs(output_path, ["qc", "alignment", "circlemap", "logs"])

    logger.info(f"Starting Circle-Map pipeline for {sample_name}")
    logger.info(f"Input: {fastq1}, {fastq2}")
    logger.info(f"Reference: {reference}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Threads: {threads}")

    reference_path = Path(reference)
    _ensure_reference_indices(
        reference_path,
        threads=threads,
        log_path=dirs["logs"] / f"{sample_name}_index.log",
        auto_index=auto_index,
    )

    # Step 1: Quality control with fastp
    if skip_fastp:
        logger.info("Step 1: Skipping fastp QC (--skip-fastp).")
        clean_r1 = Path(fastq1)
        clean_r2 = Path(fastq2)
    else:
        logger.info("Step 1: Running fastp QC...")
        clean_r1 = dirs["qc"] / f"{sample_name}_clean_R1.fq.gz"
        clean_r2 = dirs["qc"] / f"{sample_name}_clean_R2.fq.gz"
        fastp_cmd = [
            "fastp",
            "-i",
            str(fastq1),
            "-I",
            str(fastq2),
            "-o",
            str(clean_r1),
            "-O",
            str(clean_r2),
            "-h",
            str(dirs["qc"] / f"{sample_name}_fastp.html"),
            "-j",
            str(dirs["qc"] / f"{sample_name}_fastp.json"),
            "-w",
            str(threads),
        ]
        fastp_log = dirs["logs"] / f"{sample_name}_fastp.log"
        with open(fastp_log, "w") as log_handle:
            run_command(
                fastp_cmd,
                check=True,
                shell=False,
                stdout=log_handle,
                stderr=log_handle,
            )

    # Step 2: BWA alignment
    sorted_bam = dirs["alignment"] / f"{sample_name}_sorted.bam"
    logger.info("Step 2: Running BWA alignment + samtools sort (streaming)...")
    bwa_cmd = [
        "bwa",
        "mem",
        "-t",
        str(threads),
        str(reference_path),
        str(clean_r1),
        str(clean_r2),
    ]
    sort_cmd = [
        "samtools",
        "sort",
        "-@",
        str(threads),
        "-o",
        str(sorted_bam),
        "-",
    ]
    bwa_log = dirs["logs"] / f"{sample_name}_bwa_sort.log"
    with open(bwa_log, "w") as log_handle:
        bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=log_handle)
        try:
            sort_proc = subprocess.Popen(
                sort_cmd,
                stdin=bwa_proc.stdout,
                stdout=subprocess.DEVNULL,
                stderr=log_handle,
            )
        finally:
            if bwa_proc.stdout is not None:
                bwa_proc.stdout.close()

        sort_rc = sort_proc.wait()
        bwa_rc = bwa_proc.wait()
        if bwa_rc != 0:
            raise subprocess.CalledProcessError(bwa_rc, bwa_cmd)
        if sort_rc != 0:
            raise subprocess.CalledProcessError(sort_rc, sort_cmd)

    # Step 3: Index BAM
    logger.info("Step 3: Indexing BAM...")
    index_log = dirs["logs"] / f"{sample_name}_samtools_index.log"
    with open(index_log, "w") as log_handle:
        run_command(
            ["samtools", "index", "-@", str(threads), str(sorted_bam)],
            check=True,
            shell=False,
            stdout=log_handle,
            stderr=log_handle,
        )

    # Step 4: Extract discordant reads
    logger.info("Step 4: Extracting discordant reads...")
    disc_bam = dirs["alignment"] / f"{sample_name}_discordant.bam"
    disc_log = dirs["logs"] / f"{sample_name}_discordant.log"
    # Typical Circle-Map discordant extraction:
    # -f 1: paired
    # -F 1294: exclude proper pair, unmapped read, unmapped mate, secondary, duplicate
    disc_cmd = [
        "samtools",
        "view",
        "-@",
        str(threads),
        "-b",
        "-f",
        "1",
        "-F",
        "1294",
    ]
    if min_mapq > 0:
        disc_cmd.extend(["-q", str(min_mapq)])
    disc_cmd.append(str(sorted_bam))
    with open(disc_bam, "wb") as disc_out, open(disc_log, "w") as log_handle:
        run_command(disc_cmd, check=True, shell=False, stdout=disc_out, stderr=log_handle)

    with open(index_log, "a") as log_handle:
        run_command(
            ["samtools", "index", "-@", str(threads), str(disc_bam)],
            check=True,
            shell=False,
            stdout=log_handle,
            stderr=log_handle,
        )

    # Step 5: Run Circle-Map ReadExtractor
    logger.info("Step 5: Running Circle-Map ReadExtractor...")
    extracted_bam = dirs["circlemap"] / f"{sample_name}_extracted.bam"
    extractor_log = dirs["logs"] / f"{sample_name}_circlemap_readextractor.log"
    extract_cmd = [
        circlemap_exe,
        "ReadExtractor",
        "-i",
        str(sorted_bam),
        "-o",
        str(extracted_bam),
    ]
    with open(extractor_log, "w") as log_handle:
        run_command(extract_cmd, check=True, shell=False, stdout=log_handle, stderr=log_handle)

    extracted_sorted_bam = dirs["circlemap"] / f"{sample_name}_extracted.sorted.bam"
    with open(index_log, "a") as log_handle:
        run_command(
            [
                "samtools",
                "sort",
                "-@",
                str(threads),
                "-o",
                str(extracted_sorted_bam),
                str(extracted_bam),
            ],
            check=True,
            shell=False,
            stdout=log_handle,
            stderr=log_handle,
        )
    extracted_sorted_bam.replace(extracted_bam)
    with open(index_log, "a") as log_handle:
        run_command(
            ["samtools", "index", "-@", str(threads), str(extracted_bam)],
            check=True,
            shell=False,
            stdout=log_handle,
            stderr=log_handle,
        )

    # Step 6: Run Circle-Map Realign
    logger.info("Step 6: Running Circle-Map Realign...")
    circlemap_out = dirs["circlemap"] / f"{sample_name}_circle.bed"
    realign_log = dirs["logs"] / f"{sample_name}_circlemap_realign.log"
    realign_cmd = [
        circlemap_exe,
        "Realign",
        "-i",
        str(extracted_bam),
        "-qbam",
        str(sorted_bam),
        "-sbam",
        str(disc_bam),
        "-fasta",
        str(reference_path),
        "-o",
        str(circlemap_out),
    ]
    with open(realign_log, "w") as log_handle:
        run_command(realign_cmd, check=True, shell=False, stdout=log_handle, stderr=log_handle)

    if not keep_intermediate:
        logger.info("Cleaning up intermediate files (keep_intermediate=False)...")
        disc_bam.unlink(missing_ok=True)
        Path(str(disc_bam) + ".bai").unlink(missing_ok=True)
        extracted_bam.unlink(missing_ok=True)
        Path(str(extracted_bam) + ".bai").unlink(missing_ok=True)

    logger.info(f"Pipeline completed. Results: {circlemap_out}")
