"""Bedtools wrapper utilities for eccToolkit."""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Union

logger = logging.getLogger(__name__)


def run_command(
    cmd: str,
    capture_output: bool = True,
    check: bool = False,
    shell: bool = True,
) -> subprocess.CompletedProcess:
    """
    Run a shell command.

    Args:
        cmd: Command string to execute
        capture_output: Capture stdout/stderr
        check: Raise exception on non-zero return
        shell: Run through shell

    Returns:
        CompletedProcess instance
    """
    logger.debug(f"Running: {cmd}")
    result = subprocess.run(
        cmd,
        shell=shell,
        capture_output=capture_output,
        text=True,
    )
    if result.returncode != 0 and result.stderr:
        logger.warning(f"Command stderr: {result.stderr.strip()}")
    return result


def intersect(
    a_file: Union[str, Path],
    b_file: Union[str, Path],
    unique: bool = True,
    count_only: bool = True,
) -> Union[int, str]:
    """
    Run bedtools intersect.

    Args:
        a_file: First BED file (-a)
        b_file: Second BED file (-b)
        unique: Return unique entries only
        count_only: Return count instead of output

    Returns:
        Count of intersections or intersection output
    """
    cmd = f"bedtools intersect -a {a_file} -b {b_file} -wa"

    if unique:
        cmd += " | sort -u"

    if count_only:
        cmd += " | wc -l"
        result = run_command(cmd)
        if result.returncode != 0:
            logger.error(f"Bedtools intersect failed: {result.stderr}")
            return 0
        return int(result.stdout.strip())
    else:
        result = run_command(cmd)
        return result.stdout


def shuffle(
    input_file: Union[str, Path],
    genome_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    seed: Optional[int] = None,
    chrom_only: bool = True,
) -> Path:
    """
    Run bedtools shuffle.

    Args:
        input_file: Input BED file
        genome_file: Genome sizes file
        output_file: Output file (default: temp file)
        seed: Random seed
        chrom_only: Shuffle within chromosomes only

    Returns:
        Path to output file
    """
    if output_file is None:
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".bed", delete=False
        )
        output_file = tmp.name
        tmp.close()

    cmd = f"bedtools shuffle -i {input_file} -g {genome_file}"

    if chrom_only:
        cmd += " -chrom"

    if seed is not None:
        cmd += f" -seed {seed}"

    cmd += f" > {output_file} 2>/dev/null"

    result = run_command(cmd)
    if result.returncode != 0:
        logger.error(f"Bedtools shuffle failed: {result.stderr}")

    return Path(output_file)


def create_genome_file(
    output_path: Union[str, Path],
    assembly: str = "hg38",
    exclude_sex: bool = False,
    exclude_mito: bool = True,
) -> Path:
    """
    Create genome sizes file for bedtools.

    Args:
        output_path: Output file path
        assembly: Genome assembly (hg38 or hg19)
        exclude_sex: Exclude sex chromosomes
        exclude_mito: Exclude mitochondrial chromosome

    Returns:
        Path to genome file
    """
    from ecctoolkit.utils.config import get_genome_sizes

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    genome_sizes = get_genome_sizes(assembly)

    with open(output_path, "w") as f:
        for chrom, size in genome_sizes.items():
            if exclude_sex and chrom in ["chrX", "chrY"]:
                continue
            if exclude_mito and chrom in ["chrM", "chrMT"]:
                continue
            f.write(f"{chrom}\t{size}\n")

    logger.info(f"Created genome file: {output_path}")
    return output_path


def make_windows(
    genome_file: Union[str, Path],
    window_size: int,
    output_file: Optional[Union[str, Path]] = None,
) -> Path:
    """
    Create genomic windows using bedtools makewindows.

    Args:
        genome_file: Genome sizes file
        window_size: Window size in bp
        output_file: Output file path

    Returns:
        Path to output file
    """
    if output_file is None:
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".bed", delete=False
        )
        output_file = tmp.name
        tmp.close()

    cmd = f"bedtools makewindows -g {genome_file} -w {window_size} > {output_file}"
    result = run_command(cmd)

    if result.returncode != 0:
        logger.error(f"Bedtools makewindows failed: {result.stderr}")

    return Path(output_file)


def coverage(
    a_file: Union[str, Path],
    b_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
) -> Path:
    """
    Run bedtools coverage.

    Args:
        a_file: Reference BED file
        b_file: Query BED file
        output_file: Output file path

    Returns:
        Path to output file
    """
    if output_file is None:
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".bed", delete=False
        )
        output_file = tmp.name
        tmp.close()

    cmd = f"bedtools coverage -a {a_file} -b {b_file} > {output_file}"
    result = run_command(cmd)

    if result.returncode != 0:
        logger.error(f"Bedtools coverage failed: {result.stderr}")

    return Path(output_file)
