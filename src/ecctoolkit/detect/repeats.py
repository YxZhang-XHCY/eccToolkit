"""eccDNA terminal repeat identification using BLAST."""

import logging
import os
import shutil
import subprocess
import tempfile
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

# BLAST output column names (outfmt 6 + qlen/slen)
BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen",
]

OUTPUT_COLUMNS = [
    "eccDNA_name", "chr", "start", "end",
    "repeat_type", "repeat_length", "pct_identity", "e_value",
    "direct_length", "inverted_length",
]


def _load_chrom_lengths(genome_file: str) -> Dict[str, int]:
    """Load chromosome lengths from .fai index."""
    fai_path = genome_file + ".fai"
    chrom_lengths = {}
    if not os.path.exists(fai_path):
        logger.warning(f"FAI index not found: {fai_path}. Coordinate clamping will be limited.")
        return chrom_lengths
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                chrom_lengths[parts[0]] = int(parts[1])
    logger.info(f"Loaded {len(chrom_lengths)} chromosome lengths from {fai_path}")
    return chrom_lengths


def _get_chrom_length(chrom: str, chrom_lengths: Dict[str, int]) -> Optional[int]:
    """Get chromosome length, trying both chr-prefixed and non-prefixed variants."""
    if chrom in chrom_lengths:
        return chrom_lengths[chrom]
    alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    if alt in chrom_lengths:
        return chrom_lengths[alt]
    return None


def _clamp_coords(
    chrom: str, start: int, end: int, chrom_lengths: Dict[str, int]
) -> Tuple[int, int, bool]:
    """Clamp coordinates to chromosome boundaries. Returns (start, end, is_valid)."""
    cstart = max(0, start)
    cend = end
    chrom_len = _get_chrom_length(chrom, chrom_lengths)
    if chrom_len is not None:
        cend = min(chrom_len, end)
    return cstart, cend, cstart < cend


def _parse_input(input_file: str) -> pd.DataFrame:
    """Parse input file with eccDNA coordinates.

    Supports formats:
    - TSV/CSV with columns: name, chr, start, end (or chr, start, end)
    - TSV/CSV with a 'seqname' column containing chr:start-end format
    """
    sep = "\t"
    if input_file.endswith(".csv"):
        sep = ","
    df = pd.read_csv(input_file, sep=sep, comment="#")

    # Normalize column names
    col_map = {}
    for c in df.columns:
        cl = c.lower().strip()
        if cl in ("chrom", "chromosome", "chr"):
            col_map[c] = "chr"
        elif cl in ("start", "chromstart", "begin"):
            col_map[c] = "start"
        elif cl in ("end", "chromend", "stop"):
            col_map[c] = "end"
        elif cl in ("name", "ename", "eccDNA_name", "eccdna_name", "id"):
            col_map[c] = "name"
        elif cl in ("seqname", "seq_name"):
            col_map[c] = "seqname"
    df.rename(columns=col_map, inplace=True)

    # Handle seqname format: chr:start-end
    if "seqname" in df.columns and "chr" not in df.columns:
        parsed = df["seqname"].str.extract(r"([^:]+):(\d+)-(\d+)")
        df["chr"] = parsed[0]
        df["start"] = pd.to_numeric(parsed[1])
        df["end"] = pd.to_numeric(parsed[2])

    if "name" not in df.columns:
        df["name"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)

    required = ["chr", "start", "end", "name"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Required column '{col}' not found. Available: {list(df.columns)}")

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df.dropna(subset=["chr", "start", "end"], inplace=True)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    # Filter invalid
    valid = df["end"] > df["start"]
    if not valid.all():
        n_bad = (~valid).sum()
        logger.warning(f"Removed {n_bad} records with end <= start")
        df = df[valid].copy()

    logger.info(f"Loaded {len(df)} eccDNA records from {input_file}")
    return df[required].reset_index(drop=True)


def _reverse_complement(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(complement)[::-1]


def _revcomp_fasta(inpath: str, outpath: str) -> bool:
    """Write reverse complement of a FASTA file."""
    try:
        with open(inpath) as fin, open(outpath, "w") as fout:
            header = ""
            seq_parts = []
            for line in fin:
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        fout.write(header + "\n")
                        fout.write(_reverse_complement("".join(seq_parts)) + "\n")
                    header = line
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if header:
                fout.write(header + "\n")
                fout.write(_reverse_complement("".join(seq_parts)) + "\n")
        return True
    except Exception as e:
        logger.error(f"Error creating reverse complement FASTA: {e}")
        return False


def _extract_flanks_batch(
    records: List[Tuple[str, str, int, int]],
    genome_file: str,
    flank_size: int,
    chrom_lengths: Dict[str, int],
    temp_dir: str,
) -> Tuple[Optional[str], Optional[str]]:
    """Extract start and end flank FASTA for a batch of eccDNA records.

    Returns paths to (start_flanks.fa, end_flanks.fa).
    Each record: (name, chr, start, end) in BED-like 0-based half-open coords.
    """
    uid = uuid.uuid4().hex[:8]
    start_bed = os.path.join(temp_dir, f"start_flanks_{uid}.bed")
    end_bed = os.path.join(temp_dir, f"end_flanks_{uid}.bed")
    start_fa = os.path.join(temp_dir, f"start_flanks_{uid}.fa")
    end_fa = os.path.join(temp_dir, f"end_flanks_{uid}.fa")

    valid_names = []
    with open(start_bed, "w") as sf, open(end_bed, "w") as ef:
        for name, chrom, start, end in records:
            # Start flank: region around the start breakpoint
            # (start - flank_size) to (start + flank_size)
            s_start, s_end, s_valid = _clamp_coords(
                chrom, start - flank_size, start + flank_size, chrom_lengths
            )
            # End flank: region around the end breakpoint
            # (end - flank_size) to (end + flank_size)
            e_start, e_end, e_valid = _clamp_coords(
                chrom, end - flank_size, end + flank_size, chrom_lengths
            )
            if s_valid and e_valid:
                # Use name as the BED name field for matching later
                sf.write(f"{chrom}\t{s_start}\t{s_end}\t{name}\n")
                ef.write(f"{chrom}\t{e_start}\t{e_end}\t{name}\n")
                valid_names.append(name)

    if not valid_names:
        return None, None

    # Run bedtools getfasta with -name flag to keep record names
    for bed_path, fa_path in [(start_bed, start_fa), (end_bed, end_fa)]:
        cmd = [
            "bedtools", "getfasta",
            "-fi", genome_file,
            "-bed", bed_path,
            "-fo", fa_path,
            "-name",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"bedtools getfasta failed: {result.stderr.strip()}")
            return None, None

    # Clean up BED files
    for p in [start_bed, end_bed]:
        try:
            os.remove(p)
        except OSError:
            pass

    return start_fa, end_fa


def _split_fasta_by_name(fasta_path: str) -> Dict[str, str]:
    """Parse a multi-record FASTA into {name: sequence} dict.

    bedtools getfasta with -name produces headers like >name::chr:start-end
    """
    sequences = {}
    current_name = None
    seq_parts = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(seq_parts)
                # Header format: >name::chr:start-end
                header = line[1:]
                current_name = header.split("::")[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if current_name is not None:
        sequences[current_name] = "".join(seq_parts)

    return sequences


def _run_blast_pair(
    query_fa: str, subject_fa: str, temp_dir: str
) -> Optional[str]:
    """Run blastn-short on query vs subject, return output path or None."""
    uid = uuid.uuid4().hex[:8]
    out_path = os.path.join(temp_dir, f"blast_out_{uid}.txt")

    cmd = [
        "blastn",
        "-query", query_fa,
        "-subject", subject_fa,
        "-out", out_path,
        "-outfmt", "6 std qlen slen",
        "-task", "blastn-short",
        "-word_size", "7",
        "-evalue", "10",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        # Check if it's just "no hits"
        combined = (result.stdout + result.stderr).lower()
        if "no hits found" in combined or "not found" in combined:
            return None
        logger.warning(f"blastn failed: {result.stderr.strip()}")
        return None

    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return out_path
    return None


def _parse_blast_output(
    blast_path: str,
    min_identity: float,
    min_length: int,
) -> pd.DataFrame:
    """Parse BLAST tabular output and filter by identity/length."""
    try:
        df = pd.read_csv(blast_path, sep="\t", header=None, names=BLAST_COLUMNS)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

    df = df[(df["pident"] >= min_identity) & (df["length"] >= min_length)]
    return df


def _process_single_eccdna(args_tuple):
    """Process a single eccDNA record for repeat detection.

    This is the worker function for parallel processing.
    args_tuple: (name, chrom, start, end, genome_file, flank_size,
                 blast_identity, min_repeat_length, chrom_lengths)
    """
    (name, chrom, start, end, genome_file, flank_size,
     blast_identity, min_repeat_length, chrom_lengths) = args_tuple

    temp_dir = tempfile.mkdtemp(prefix="ecc_repeat_")
    direct_best = None
    inverted_best = None

    try:
        # Create single-record BED and extract flanks
        records = [(name, chrom, start, end)]
        start_fa, end_fa = _extract_flanks_batch(
            records, genome_file, flank_size, chrom_lengths, temp_dir
        )
        if not start_fa or not end_fa:
            return _make_result(name, chrom, start, end, None, None)

        # Direct repeat: start_flank vs end_flank
        direct_out = _run_blast_pair(start_fa, end_fa, temp_dir)
        if direct_out:
            df_direct = _parse_blast_output(direct_out, blast_identity, min_repeat_length)
            if not df_direct.empty:
                # Take best hit by bitscore
                best_idx = df_direct["bitscore"].idxmax()
                direct_best = df_direct.loc[best_idx]

        # Inverted repeat: start_flank vs revcomp(end_flank)
        end_rc_fa = os.path.join(temp_dir, "end_revcomp.fa")
        if _revcomp_fasta(end_fa, end_rc_fa):
            inv_out = _run_blast_pair(start_fa, end_rc_fa, temp_dir)
            if inv_out:
                df_inv = _parse_blast_output(inv_out, blast_identity, min_repeat_length)
                if not df_inv.empty:
                    best_idx = df_inv["bitscore"].idxmax()
                    inverted_best = df_inv.loc[best_idx]

    except Exception as e:
        logger.error(f"Error processing {name}: {e}")
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    return _make_result(name, chrom, start, end, direct_best, inverted_best)


def _make_result(
    name: str, chrom: str, start: int, end: int,
    direct_best: Optional[pd.Series], inverted_best: Optional[pd.Series],
) -> dict:
    """Create a result dict for one eccDNA."""
    has_direct = direct_best is not None
    has_inverted = inverted_best is not None

    if has_direct and has_inverted:
        repeat_type = "Both"
    elif has_direct:
        repeat_type = "Direct"
    elif has_inverted:
        repeat_type = "Inverted"
    else:
        repeat_type = "No_Repeat"

    # Pick the best hit for summary fields
    best = None
    if has_direct and has_inverted:
        best = direct_best if direct_best["bitscore"] >= inverted_best["bitscore"] else inverted_best
    elif has_direct:
        best = direct_best
    elif has_inverted:
        best = inverted_best

    return {
        "eccDNA_name": name,
        "chr": chrom,
        "start": start,
        "end": end,
        "repeat_type": repeat_type,
        "repeat_length": int(best["length"]) if best is not None else 0,
        "pct_identity": float(best["pident"]) if best is not None else 0.0,
        "e_value": float(best["evalue"]) if best is not None else 0.0,
        "direct_length": int(direct_best["length"]) if has_direct else 0,
        "inverted_length": int(inverted_best["length"]) if has_inverted else 0,
    }


def find_terminal_repeats(
    input_file: str,
    genome_file: str,
    output_file: str,
    threads: int = 8,
    breakpoint_tolerance: int = 5,
    flank_size: int = 100,
    blast_identity: float = 80.0,
    min_repeat_length: int = 10,
    batch_size: int = 100,
) -> None:
    """Find terminal repeats at eccDNA breakpoints using BLAST.

    For each eccDNA, extracts flanking sequences around the start and end
    breakpoints, then runs bidirectional BLAST to detect direct and inverted
    repeats at the junction.

    Args:
        input_file: TSV/CSV file with eccDNA coordinates
        genome_file: Reference genome FASTA (with .fai index)
        output_file: Output CSV file
        threads: Number of parallel workers
        breakpoint_tolerance: bp tolerance for breakpoint proximity
        flank_size: Size of flanking region to extract (bp)
        blast_identity: Minimum BLAST percent identity
        min_repeat_length: Minimum repeat alignment length (bp)
        batch_size: Batch size for parallel processing
    """
    # Check dependencies
    for tool in ["bedtools", "blastn"]:
        if shutil.which(tool) is None:
            raise RuntimeError(f"Required tool '{tool}' not found in PATH.")

    # Load inputs
    chrom_lengths = _load_chrom_lengths(genome_file)
    df_input = _parse_input(input_file)

    if df_input.empty:
        logger.warning("No valid eccDNA records to process.")
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(output_file, index=False)
        return

    # Build argument tuples
    work_items = [
        (row["name"], row["chr"], row["start"], row["end"],
         genome_file, flank_size, blast_identity, min_repeat_length, chrom_lengths)
        for _, row in df_input.iterrows()
    ]

    results = []
    n_total = len(work_items)
    logger.info(f"Processing {n_total} eccDNA records with {threads} worker(s)")

    if threads <= 1:
        for i, item in enumerate(work_items):
            results.append(_process_single_eccdna(item))
            if (i + 1) % batch_size == 0 or (i + 1) == n_total:
                logger.info(f"Progress: {i + 1}/{n_total}")
    else:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {
                executor.submit(_process_single_eccdna, item): i
                for i, item in enumerate(work_items)
            }
            completed = 0
            for future in as_completed(futures):
                try:
                    results.append(future.result())
                except Exception as e:
                    idx = futures[future]
                    logger.error(f"Worker error for record {idx}: {e}")
                completed += 1
                if completed % batch_size == 0 or completed == n_total:
                    logger.info(f"Progress: {completed}/{n_total}")

    # Build output
    df_out = pd.DataFrame(results, columns=OUTPUT_COLUMNS)
    df_out.to_csv(output_file, index=False)

    # Summary
    counts = df_out["repeat_type"].value_counts()
    logger.info(f"Results saved to {output_file}")
    logger.info(f"Summary: {dict(counts)}")
