"""eccDNA validation from amplicon sequencing."""

import csv
import logging
import os
import re
import subprocess
from pathlib import Path
from typing import List, Optional

import pandas as pd
import pysam

from ecctoolkit.utils.subprocess_utils import require_tools, run_command
from ecctoolkit.utils.validation import validate_file_exists

logger = logging.getLogger(__name__)


def _build_junction_reference(
    input_file: str,
    reference: str,
    output_fasta: str,
    inside_bp: int,
    flank_bp: int,
) -> None:
    """Build junction-centered reference for each eccDNA candidate.

    For each candidate, extracts inside (inside_bp) and outside (flank_bp)
    flanking sequences at both start and end, concatenating them to create
    a junction-spanning reference:
      [tail_inside | tail_outside | head_outside | head_inside]
    """
    ref = pysam.FastaFile(reference)
    with open(input_file) as fh, open(output_fasta, "w") as out:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            chrom = row["Chr"]
            start = int(row["Start"])
            end = int(row["End"])
            name = row.get("Name") or f"ecc_{row.get('No', 'unknown')}"
            chrlen = ref.get_reference_length(chrom)

            # Dynamic trim: shrink inside if eccDNA is too short
            inside = min(inside_bp, max(1, (end - start + 1) // 2))

            # Tail: inside/outside coordinates (1-based inclusive -> 0-based)
            tail_in_s = max(1, end - inside + 1)
            tail_in_e = end
            tail_out_s = end + 1
            tail_out_e = min(chrlen, end + flank_bp)

            # Head: outside/inside coordinates
            head_out_s = max(1, start - flank_bp)
            head_out_e = start - 1
            head_in_s = start
            head_in_e = min(chrlen, start + inside - 1)

            seq = (
                ref.fetch(chrom, tail_in_s - 1, tail_in_e)
                + ref.fetch(chrom, tail_out_s - 1, tail_out_e)
                + ref.fetch(chrom, head_out_s - 1, head_out_e)
                + ref.fetch(chrom, head_in_s - 1, head_in_e)
            )

            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i : i + 60] + "\n")

    ref.close()
    logger.info(f"Built junction reference: {output_fasta}")


def _bwa_align(
    ref_fasta: str,
    fastq1: str,
    fastq2: Optional[str],
    threads: int,
    output_prefix: str,
) -> str:
    """Index reference with BWA, align reads, sort and index BAM."""
    # BWA index
    if not os.path.exists(ref_fasta + ".bwt"):
        run_command(
            ["bwa", "index", ref_fasta],
            shell=False, check=True, capture_output=False,
        )

    bam = output_prefix + ".sorted.bam"

    bwa_cmd = ["bwa", "mem", "-t", str(threads), ref_fasta, fastq1]
    if fastq2:
        bwa_cmd.append(fastq2)

    sort_cmd = ["samtools", "sort", f"-@{threads}", "-o", bam]

    logger.info("Running BWA alignment...")
    bwa_proc = subprocess.Popen(
        bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    subprocess.check_call(sort_cmd, stdin=bwa_proc.stdout)
    bwa_proc.stdout.close()
    if bwa_proc.wait() != 0:
        raise subprocess.CalledProcessError(bwa_proc.returncode, bwa_cmd)

    run_command(
        ["samtools", "index", bam],
        shell=False, check=True, capture_output=False,
    )

    return bam


def _crosses_center(aln: pysam.AlignedSegment, center: int, window: int) -> bool:
    """Check if alignment directly spans the junction center +/- window."""
    return aln.reference_start < center - window and aln.reference_end > center + window


def _sa_cross(aln: pysam.AlignedSegment, center: int, window: int) -> bool:
    """Check if SA tag indicates a split read crossing the junction."""
    if not aln.has_tag("SA"):
        return False
    for block in aln.get_tag("SA").split(";"):
        if not block:
            continue
        fields = block.split(",")
        if len(fields) < 2:
            continue
        sa_pos = int(fields[1]) - 1  # 0-based
        if (aln.reference_end < center - window and sa_pos > center + window) or (
            aln.reference_start > center + window and sa_pos < center - window
        ):
            return True
    return False


def _analyze_bam(
    bamfile: str,
    window_bp: int,
    min_crossing: int,
    min_unique: int,
    min_split: int,
) -> pd.DataFrame:
    """Analyze BAM to classify eccDNA candidates."""
    bam = pysam.AlignmentFile(bamfile)
    results = []

    for ref in bam.references:
        unique_reads = set()
        split_reads = 0

        ref_len = bam.get_reference_length(ref)
        center = ref_len // 2
        crossing_reads = set()

        for aln in bam.fetch(ref):
            if aln.is_unmapped or aln.is_duplicate:
                continue

            unique_reads.add(aln.query_name)

            if aln.has_tag("SA"):
                split_reads += 1
            else:
                cig = aln.cigarstring or ""
                if re.match(r"^\d+S", cig) or cig.endswith("S"):
                    split_reads += 1

            if _crosses_center(aln, center, window_bp) or _sa_cross(
                aln, center, window_bp
            ):
                crossing_reads.add(aln.query_name)

        n_unique = len(unique_reads)
        n_crossing = len(crossing_reads)

        # Two-tier detection strategy
        basic_pass = n_unique >= min_unique and split_reads >= min_split
        crossing_pass = n_crossing >= min_crossing

        if basic_pass and crossing_pass:
            status = "Detected"
            confidence = "High"
        elif basic_pass or crossing_pass:
            status = "Detected"
            confidence = "Medium"
        else:
            status = "Not_detected"
            confidence = "Low"

        results.append(
            {
                "ecc_id": ref,
                "unique_reads": n_unique,
                "split_reads": split_reads,
                "crossing_reads": n_crossing,
                "ref_length": ref_len,
                "status": status,
                "confidence": confidence,
            }
        )

        logger.debug(
            f"{ref}: {n_unique} unique, {split_reads} split, "
            f"{n_crossing} crossing -> {status} ({confidence})"
        )

    bam.close()
    return pd.DataFrame(results)


def run_validation(
    input_file: str,
    fastq1: str,
    fastq2: str,
    reference: str,
    output_dir: str,
    threads: int = 8,
    inside_bp: int = 300,
    flank_bp: int = 200,
    window_bp: int = 20,
    min_crossing: int = 3,
    min_unique: int = 20,
    min_split: int = 3,
) -> None:
    """
    Validate eccDNA candidates from amplicon sequencing.

    Builds junction-spanning reference sequences for each candidate,
    aligns amplicon reads, and analyzes junction-crossing evidence.

    Args:
        input_file: TSV file with eccDNA candidates (Chr, Start, End columns)
        fastq1: Path to FASTQ file (Read 1)
        fastq2: Path to FASTQ file (Read 2, or None for single-end)
        reference: Path to indexed reference genome FASTA
        output_dir: Output directory
        threads: Number of threads
        inside_bp: Base pairs to keep inside eccDNA at each side
        flank_bp: Base pairs to keep outside eccDNA at each side
        window_bp: Half-window around center for junction crossing
        min_crossing: Minimum crossing reads for junction support
        min_unique: Minimum unique reads for initial detection
        min_split: Minimum split reads for initial detection
    """
    validate_file_exists(input_file, "Candidates TSV")
    validate_file_exists(fastq1, "FASTQ R1")
    if fastq2:
        validate_file_exists(fastq2, "FASTQ R2")
    validate_file_exists(reference, "Reference FASTA")

    require_tools(["bwa", "samtools"])

    os.makedirs(output_dir, exist_ok=True)
    prefix = os.path.join(output_dir, Path(input_file).stem)

    logger.info("=== eccDNA Validation Pipeline ===")
    logger.info(f"Candidates: {input_file}")
    logger.info(f"Reads: {fastq1}, {fastq2}")
    logger.info(
        f"Parameters: inside={inside_bp}bp, flank={flank_bp}bp, "
        f"window={window_bp}bp"
    )
    logger.info(
        f"Thresholds: unique>={min_unique}, split>={min_split}, "
        f"crossing>={min_crossing}"
    )

    # Step 1: Build junction reference
    logger.info("Step 1: Building junction-spanning reference...")
    amplicon_fa = prefix + ".amplicon.fa"
    _build_junction_reference(
        input_file, reference, amplicon_fa, inside_bp, flank_bp,
    )

    # Step 2: BWA alignment
    logger.info("Step 2: Aligning amplicon reads...")
    bam_file = _bwa_align(amplicon_fa, fastq1, fastq2, threads, prefix)

    # Step 3: Analyze BAM
    logger.info("Step 3: Analyzing junction evidence...")
    df = _analyze_bam(bam_file, window_bp, min_crossing, min_unique, min_split)

    # Save results
    detail_file = prefix + ".validation_detail.tsv"
    df.to_csv(detail_file, sep="\t", index=False)

    summary_file = prefix + ".validation_summary.tsv"
    summary = df[
        ["ecc_id", "status", "confidence", "unique_reads", "crossing_reads"]
    ].copy()
    summary.to_csv(summary_file, sep="\t", index=False)

    # Report
    n_total = len(df)
    n_detected = len(df[df["status"] == "Detected"])
    n_high = len(df[df["confidence"] == "High"])
    n_medium = len(df[df["confidence"] == "Medium"])
    pct = (n_detected / n_total * 100) if n_total > 0 else 0.0

    logger.info("=== Validation Complete ===")
    logger.info(f"  Total candidates: {n_total}")
    logger.info(f"  Detected: {n_detected} ({pct:.1f}%)")
    logger.info(f"  High confidence: {n_high}")
    logger.info(f"  Medium confidence: {n_medium}")
    logger.info(f"  Detail: {detail_file}")
    logger.info(f"  Summary: {summary_file}")
