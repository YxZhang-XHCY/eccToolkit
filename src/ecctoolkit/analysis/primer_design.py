"""Outward-facing primer design for eccDNA junction validation.

Designs primers that span the circularization junction of eccDNA,
with BLAST-based specificity validation against the reference genome.

Algorithm:
  1. Extract eccDNA sequence from reference genome
  2. Concatenate tail + head to create a linear template spanning the junction
  3. Use Primer3 to design primer pairs flanking the junction
  4. Validate specificity with BLAST (only keep primers with exactly 1 perfect match)
  5. Retry with relaxed parameters if needed
"""

import logging
import os
import random
import shutil
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Default Primer3 parameters
DEFAULT_PRIMER3_ARGS = {
    "PRIMER_NUM_RETURN": 5,
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_TM": 58.0,
    "PRIMER_MAX_TM": 62.0,
    "PRIMER_MIN_GC": 40.0,
    "PRIMER_MAX_GC": 60.0,
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_DNA_CONC": 50.0,
    "PRIMER_MAX_NS_ACCEPTED": 0,
    "PRIMER_MAX_SELF_ANY": 8,
    "PRIMER_MAX_SELF_END": 3,
    "PRIMER_PAIR_MAX_COMPL_ANY": 8,
    "PRIMER_PAIR_MAX_COMPL_END": 3,
    "PRIMER_GC_CLAMP": 2,
}


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def _parse_input(input_file: str) -> pd.DataFrame:
    """Parse eccDNA input into DataFrame with name, chr, start, end."""
    import re

    is_bed = input_file.endswith((".bed",))
    sep = "\t" if input_file.endswith((".bed", ".tsv")) else ","

    if is_bed:
        # BED files never have headers
        df = pd.read_csv(input_file, sep="\t", comment="#", header=None)
        if df.shape[1] >= 4:
            df.columns = ["chr", "start", "end", "name"] + [
                f"col{i}" for i in range(4, df.shape[1])
            ]
        elif df.shape[1] == 3:
            df.columns = ["chr", "start", "end"]
        else:
            raise ValueError(f"BED file must have at least 3 columns, got {df.shape[1]}")
    else:
        try:
            df = pd.read_csv(input_file, sep=sep, comment="#")
        except Exception:
            df = pd.read_csv(input_file, sep="\t", comment="#", header=None)

        col_map = {}
        for c in df.columns:
            cl = str(c).lower().strip()
            if cl in ("chr", "chrom", "chromosome", "#chrom", "#chr"):
                col_map[c] = "chr"
            elif cl == "start":
                col_map[c] = "start"
            elif cl == "end":
                col_map[c] = "end"
            elif cl in ("name", "eccdna_name", "eccdna_id", "id"):
                col_map[c] = "name"
            elif cl == "seqname":
                col_map[c] = "seqname"
        df = df.rename(columns=col_map)

        if "chr" not in df.columns and df.shape[1] >= 3:
            df.columns = ["chr", "start", "end"] + [
                f"col{i}" for i in range(3, df.shape[1])
            ]

    if "chr" not in df.columns and "seqname" in df.columns:
        pattern = re.compile(r"^(.+?):(\d+)-(\d+)")
        parsed = df["seqname"].str.extract(pattern)
        df["chr"] = parsed[0]
        df["start"] = parsed[1].astype(int)
        df["end"] = parsed[2].astype(int)

    if "chr" not in df.columns:
        raise ValueError("Cannot find chr/start/end columns.")

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["length"] = df["end"] - df["start"]

    if "name" not in df.columns:
        df["name"] = df.apply(
            lambda r: f"{r['chr']}:{r['start']}-{r['end']}", axis=1
        )

    return df[["name", "chr", "start", "end", "length"]].copy()


# ---------------------------------------------------------------------------
# Sequence preparation
# ---------------------------------------------------------------------------

def _prepare_junction_template(sequence: str, junction_flank: int = 300) -> Tuple[str, int]:
    """Create a linear template spanning the eccDNA junction.

    For outward-facing primers, we concatenate the tail and head of the
    eccDNA sequence to create a template that spans the circularization point.

    Args:
        sequence: The eccDNA sequence.
        junction_flank: Size of flanking region around junction.

    Returns:
        template: The junction-spanning template sequence.
        junction_pos: Position of the junction point in the template.
    """
    length = len(sequence)
    if length <= junction_flank * 2:
        # Short eccDNA: duplicate entire sequence
        template = sequence + sequence
        junction_pos = length
    else:
        # Take tail + head flanking the junction
        tail = sequence[-junction_flank:]
        head = sequence[:junction_flank]
        template = tail + head
        junction_pos = junction_flank

    return template, junction_pos


# ---------------------------------------------------------------------------
# Primer3 design
# ---------------------------------------------------------------------------

def _design_primers_primer3(
    template: str,
    junction_pos: int,
    seq_id: str,
    primer3_args: Dict,
    product_size_range: List[int],
) -> Optional[Dict]:
    """Design primers using Primer3, targeting a product spanning the junction.

    Args:
        template: Junction-spanning template sequence.
        junction_pos: Position of the junction in the template.
        seq_id: Sequence identifier.
        primer3_args: Primer3 global arguments.
        product_size_range: [min_size, max_size] for product.

    Returns:
        Primer3 result dict, or None if no primers found.
    """
    import primer3

    args = primer3_args.copy()
    args["PRIMER_PRODUCT_SIZE_RANGE"] = [product_size_range]

    seq_args = {
        "SEQUENCE_ID": seq_id,
        "SEQUENCE_TEMPLATE": template,
        # Force the product to span the junction point
        "SEQUENCE_TARGET": [junction_pos, 1],
    }

    try:
        result = primer3.bindings.design_primers(seq_args, args)
        if result is None or result.get("PRIMER_PAIR_NUM_RETURNED", 0) == 0:
            return None
        return result
    except Exception as e:
        logger.warning(f"Primer3 error for {seq_id}: {e}")
        return None


def _relax_parameters(primer3_args: Dict, attempt: int) -> Dict:
    """Progressively relax Primer3 parameters for retry."""
    args = primer3_args.copy()

    if attempt >= 1:
        args["PRIMER_MAX_SELF_ANY"] = args.get("PRIMER_MAX_SELF_ANY", 8) + 1
        args["PRIMER_MAX_SELF_END"] = args.get("PRIMER_MAX_SELF_END", 3) + 1
        args["PRIMER_PAIR_MAX_COMPL_ANY"] = args.get("PRIMER_PAIR_MAX_COMPL_ANY", 8) + 1
        args["PRIMER_PAIR_MAX_COMPL_END"] = args.get("PRIMER_PAIR_MAX_COMPL_END", 3) + 1

    if attempt >= 2:
        args["PRIMER_MIN_TM"] = args.get("PRIMER_MIN_TM", 58.0) - 1.0
        args["PRIMER_MAX_TM"] = args.get("PRIMER_MAX_TM", 62.0) + 1.0
        args["PRIMER_MIN_GC"] = args.get("PRIMER_MIN_GC", 40.0) - 5.0
        args["PRIMER_MAX_GC"] = args.get("PRIMER_MAX_GC", 60.0) + 5.0

    return args


# ---------------------------------------------------------------------------
# BLAST specificity check
# ---------------------------------------------------------------------------

def _check_blast_db(reference: str) -> str:
    """Find or create BLAST database for the reference genome."""
    db_base = os.path.splitext(reference)[0]

    # Check if BLAST DB already exists
    for ext in [".ndb", ".nin", ".nsq"]:
        if os.path.exists(db_base + ext):
            return db_base

    # Create BLAST DB
    logger.info("Creating BLAST database (first run only)...")
    cmd = ["makeblastdb", "-in", reference, "-dbtype", "nucl", "-out", db_base]
    subprocess.run(cmd, capture_output=True, check=True)
    logger.info(f"BLAST database created: {db_base}")
    return db_base


def _blast_check_specificity(
    primer_seq: str,
    db_name: str,
    threads: int = 1,
) -> bool:
    """Check if a primer has exactly 1 perfect match in the genome.

    Returns True if specific (exactly 1 perfect hit), False otherwise.
    """
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fa", delete=False
    ) as tmp:
        tmp.write(f">query\n{primer_seq}\n")
        tmp_path = tmp.name

    try:
        cmd = [
            "blastn", "-task", "blastn-short",
            "-db", db_name,
            "-query", tmp_path,
            "-outfmt", "6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore",
            "-num_threads", str(threads),
            "-evalue", "10",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        hits = result.stdout.strip().split("\n")

        perfect_matches = 0
        for hit in hits:
            if not hit.strip():
                continue
            fields = hit.split("\t")
            if len(fields) < 4:
                continue
            identity = float(fields[2])
            aln_len = int(fields[3])
            if identity == 100.0 and aln_len == len(primer_seq):
                perfect_matches += 1

        return perfect_matches == 1

    except (subprocess.TimeoutExpired, Exception) as e:
        logger.warning(f"BLAST check failed for primer: {e}")
        return False
    finally:
        os.unlink(tmp_path)


# ---------------------------------------------------------------------------
# Single eccDNA processing
# ---------------------------------------------------------------------------

def _process_single_eccdna(
    name: str,
    chrom: str,
    start: int,
    end: int,
    genome_fh,
    db_name: Optional[str],
    primer3_args: Dict,
    product_size_range: List[int],
    junction_flank: int,
    max_attempts: int,
    blast_threads: int,
    skip_blast: bool,
) -> List[Dict]:
    """Design primers for a single eccDNA.

    Returns list of primer result dicts.
    """
    # Extract sequence
    try:
        sequence = genome_fh.fetch(chrom, start, end).upper()
    except Exception as e:
        logger.warning(f"Cannot fetch sequence for {name} ({chrom}:{start}-{end}): {e}")
        return [_empty_result(name, chrom, start, end, "fetch_failed")]

    if len(sequence) < 50:
        logger.warning(f"Sequence too short for {name}: {len(sequence)}bp")
        return [_empty_result(name, chrom, start, end, "too_short")]

    # Prepare junction template
    template, junction_pos = _prepare_junction_template(sequence, junction_flank)

    seq_id = f"{chrom}:{start}-{end}"
    all_results = []

    for attempt in range(max_attempts):
        current_args = _relax_parameters(primer3_args, attempt)

        p3_result = _design_primers_primer3(
            template, junction_pos, seq_id, current_args, product_size_range
        )

        if p3_result is None:
            continue

        n_returned = p3_result.get("PRIMER_PAIR_NUM_RETURNED", 0)

        for i in range(n_returned):
            fwd = p3_result.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
            rev = p3_result.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
            fwd_pos = p3_result.get(f"PRIMER_LEFT_{i}", (0, 0))
            rev_pos = p3_result.get(f"PRIMER_RIGHT_{i}", (0, 0))
            fwd_tm = p3_result.get(f"PRIMER_LEFT_{i}_TM", 0)
            rev_tm = p3_result.get(f"PRIMER_RIGHT_{i}_TM", 0)
            fwd_gc = p3_result.get(f"PRIMER_LEFT_{i}_GC_PERCENT", 0)
            rev_gc = p3_result.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", 0)
            product_size = p3_result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0)

            # BLAST check
            fwd_specific = True
            rev_specific = True
            if not skip_blast and db_name:
                fwd_specific = _blast_check_specificity(fwd, db_name, blast_threads)
                rev_specific = _blast_check_specificity(rev, db_name, blast_threads)

            result = {
                "eccDNA_name": name,
                "chr": chrom,
                "start": start,
                "end": end,
                "eccDNA_length": end - start,
                "primer_pair_id": f"{name}_pair{i + 1}",
                "forward_primer": fwd,
                "reverse_primer": rev,
                "forward_tm": round(fwd_tm, 1),
                "reverse_tm": round(rev_tm, 1),
                "forward_gc": round(fwd_gc, 1),
                "reverse_gc": round(rev_gc, 1),
                "product_size": product_size,
                "forward_specific": fwd_specific,
                "reverse_specific": rev_specific,
                "both_specific": fwd_specific and rev_specific,
                "attempt": attempt + 1,
                "status": "success",
            }
            all_results.append(result)

        # Stop if we found at least one fully specific pair
        if any(r["both_specific"] for r in all_results):
            break

    if not all_results:
        return [_empty_result(name, chrom, start, end, "no_primers")]

    return all_results


def _empty_result(name: str, chrom: str, start: int, end: int, reason: str) -> Dict:
    """Create an empty result for failed primer design."""
    return {
        "eccDNA_name": name,
        "chr": chrom,
        "start": start,
        "end": end,
        "eccDNA_length": end - start,
        "primer_pair_id": "",
        "forward_primer": "",
        "reverse_primer": "",
        "forward_tm": 0,
        "reverse_tm": 0,
        "forward_gc": 0,
        "reverse_gc": 0,
        "product_size": 0,
        "forward_specific": False,
        "reverse_specific": False,
        "both_specific": False,
        "attempt": 0,
        "status": reason,
    }


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_primer_design(
    input_file: str,
    genome_file: str,
    output_file: str,
    product_size: int = 200,
    junction_flank: int = 300,
    max_attempts: int = 3,
    threads: int = 1,
    sample_n: Optional[int] = None,
    seed: int = 42,
    skip_blast: bool = False,
    best_only: bool = False,
) -> None:
    """Design outward-facing primers for eccDNA junction validation.

    For each eccDNA, designs primers that flank the circularization junction.
    These outward-facing primers will only produce a PCR product if the
    DNA is circular (or contains the junction sequence).

    Args:
        input_file: eccDNA regions file (CSV/TSV/BED with chr, start, end).
        genome_file: Reference genome FASTA (must have .fai index).
        output_file: Output CSV/TSV file for primer results.
        product_size: Target PCR product size in bp (default: 200).
        junction_flank: Flanking size around junction for template (default: 300).
        max_attempts: Maximum design attempts with relaxed params (default: 3).
        threads: Threads for BLAST specificity check (default: 1).
        sample_n: Randomly sample N eccDNA for primer design (default: all).
        seed: Random seed for sampling (default: 42).
        skip_blast: Skip BLAST specificity check (default: False).
        best_only: Only output the best (first specific) pair per eccDNA (default: False).
    """
    import pysam

    # Check tools
    if not skip_blast:
        if not shutil.which("blastn"):
            raise RuntimeError(
                "blastn not found. Install BLAST+ or use --skip-blast."
            )
        if not shutil.which("makeblastdb"):
            raise RuntimeError("makeblastdb not found. Install BLAST+.")

    # Load input
    logger.info(f"Loading eccDNA regions from {input_file}")
    df = _parse_input(input_file)
    logger.info(f"Loaded {len(df)} eccDNA regions")

    # Random sampling
    if sample_n is not None and sample_n < len(df):
        logger.info(f"Randomly sampling {sample_n} eccDNA (seed={seed})")
        df = df.sample(n=sample_n, random_state=seed).reset_index(drop=True)
        logger.info(f"Selected {len(df)} eccDNA for primer design")

    # Open genome
    genome_fh = pysam.FastaFile(genome_file)

    # BLAST database
    db_name = None
    if not skip_blast:
        db_name = _check_blast_db(genome_file)

    # Primer3 args
    primer3_args = DEFAULT_PRIMER3_ARGS.copy()
    product_min = max(80, product_size - 100)
    product_max = product_size + 100
    product_size_range = [product_min, product_max]

    # Process each eccDNA
    all_results = []
    n_success = 0
    n_total = len(df)

    for idx, row in df.iterrows():
        if (idx + 1) % 5 == 0 or idx == 0:
            logger.info(f"Processing {idx + 1}/{n_total}: {row['name']}")

        results = _process_single_eccdna(
            name=row["name"],
            chrom=row["chr"],
            start=int(row["start"]),
            end=int(row["end"]),
            genome_fh=genome_fh,
            db_name=db_name,
            primer3_args=primer3_args,
            product_size_range=product_size_range,
            junction_flank=junction_flank,
            max_attempts=max_attempts,
            blast_threads=threads,
            skip_blast=skip_blast,
        )

        if best_only:
            # Keep only the best specific pair, or first result if none specific
            specific = [r for r in results if r.get("both_specific", False)]
            if specific:
                all_results.append(specific[0])
                n_success += 1
            else:
                all_results.append(results[0])
        else:
            all_results.extend(results)
            if any(r.get("both_specific", False) for r in results):
                n_success += 1

    genome_fh.close()

    # Output
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    result_df = pd.DataFrame(all_results)
    sep = "\t" if output_file.endswith(".tsv") else ","
    result_df.to_csv(output_file, index=False, sep=sep)

    logger.info(f"Results written to {output_file}")
    logger.info(f"=== Primer Design Summary ===")
    logger.info(f"Total eccDNA: {n_total}")
    logger.info(f"Successful (specific primers): {n_success}")
    logger.info(f"Success rate: {n_success / n_total * 100:.1f}%")

    if not skip_blast and "both_specific" in result_df.columns:
        n_specific = result_df["both_specific"].sum()
        logger.info(f"Total specific primer pairs: {n_specific}")

    logger.info("Primer design completed!")
