"""RepeatMasker .out format conversion utilities."""

import logging
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)


def _load_chrom_alias(
    alias_path: str,
    style: str = "ucsc",
) -> Dict[str, str]:
    """Load chromosome name mapping from UCSC chromAlias file.

    Args:
        alias_path: Path to chromAlias.txt
        style: Name style to use ('ucsc', 'assembly', 'ncbi')

    Returns:
        Dictionary mapping RefSeq names to target names
    """
    style_index = {"ucsc": 4, "assembly": 1, "genbank": 2, "ncbi": 3}
    if style not in style_index:
        raise ValueError(f"Unknown style '{style}'. Use: {list(style_index.keys())}")

    col_idx = style_index[style]
    mapping = {}

    with open(alias_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) <= col_idx:
                continue
            refseq = parts[0]
            target = parts[col_idx]
            if target:
                mapping[refseq] = target

    logger.info(f"Loaded {len(mapping)} chromosome mappings ({style} style)")
    return mapping


def _parse_rmsk_out_to_ucsc(
    input_path: str,
    output_path: str,
    alias_map: Optional[Dict[str, str]],
) -> int:
    """Parse RepeatMasker .out and convert to UCSC rmsk.txt format.

    Returns:
        Number of converted entries
    """
    n_converted = 0
    n_skipped = 0

    with open(input_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("SW") or line.startswith("score"):
                continue

            # Remove trailing * (overlapping annotation marker)
            if line.endswith("*"):
                line = line[:-1].strip()

            fields = line.split()
            if len(fields) < 15:
                continue

            try:
                sw_score = int(fields[0])
                perc_div = float(fields[1])
                perc_del = float(fields[2])
                perc_ins = float(fields[3])
                chrom = fields[4]
                begin = int(fields[5])
                end = int(fields[6])
                left_query = fields[7]
                strand = fields[8]
                rep_name = fields[9]
                class_family = fields[10]

                # Convert strand: C -> -
                strand = "-" if strand == "C" else "+"

                # Convert chromosome name
                if alias_map is not None:
                    if chrom in alias_map:
                        chrom = alias_map[chrom]
                    else:
                        n_skipped += 1
                        continue

                # Split class/family
                if "/" in class_family:
                    rep_class, rep_family = class_family.split("/", 1)
                else:
                    rep_class = class_family
                    rep_family = class_family

                # Convert to milliDiv/Del/Ins (x10)
                milli_div = int(round(perc_div * 10))
                milli_del = int(round(perc_del * 10))
                milli_ins = int(round(perc_ins * 10))

                # 0-based start (RepeatMasker .out uses 1-based)
                geno_start = begin - 1

                # UCSC rmsk.txt format:
                # bin swScore milliDiv milliDel milliIns genoName genoStart
                # genoEnd genoLeft strand repName repClass repFamily
                # repStart repEnd repLeft id
                out_fields = [
                    "0", str(sw_score),
                    str(milli_div), str(milli_del), str(milli_ins),
                    chrom, str(geno_start), str(end), left_query,
                    strand, rep_name, rep_class, rep_family,
                    "0", "0", "(0)", str(n_converted),
                ]
                fout.write("\t".join(out_fields) + "\n")
                n_converted += 1

            except (ValueError, IndexError):
                n_skipped += 1
                continue

    logger.info(f"Converted {n_converted:,} entries, skipped {n_skipped:,}")
    return n_converted


def _parse_rmsk_out_to_bed(
    input_path: str,
    output_path: str,
    alias_map: Optional[Dict[str, str]],
) -> int:
    """Parse RepeatMasker .out and convert to BED format.

    Output columns: chrom, start, end, name, score, strand, class, family

    Returns:
        Number of converted entries
    """
    n_converted = 0
    n_skipped = 0

    with open(input_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("SW") or line.startswith("score"):
                continue

            if line.endswith("*"):
                line = line[:-1].strip()

            fields = line.split()
            if len(fields) < 15:
                continue

            try:
                sw_score = int(fields[0])
                chrom = fields[4]
                begin = int(fields[5])
                end = int(fields[6])
                strand = fields[8]
                rep_name = fields[9]
                class_family = fields[10]

                strand = "-" if strand == "C" else "+"

                if alias_map is not None:
                    if chrom in alias_map:
                        chrom = alias_map[chrom]
                    else:
                        n_skipped += 1
                        continue

                if "/" in class_family:
                    rep_class, rep_family = class_family.split("/", 1)
                else:
                    rep_class = class_family
                    rep_family = class_family

                # 0-based start
                geno_start = begin - 1

                out_fields = [
                    chrom, str(geno_start), str(end),
                    rep_name, str(sw_score), strand,
                    rep_class, rep_family,
                ]
                fout.write("\t".join(out_fields) + "\n")
                n_converted += 1

            except (ValueError, IndexError):
                n_skipped += 1
                continue

    logger.info(f"Converted {n_converted:,} entries, skipped {n_skipped:,}")
    return n_converted


def convert_repeatmasker(
    input_file: str,
    output_file: str,
    chrom_alias: Optional[str] = None,
    output_format: str = "ucsc",
) -> None:
    """
    Convert RepeatMasker .out format to UCSC rmsk.txt or BED format.

    Handles coordinate conversion (1-based to 0-based), strand conversion
    (C to -), milliDiv calculation, and class/family splitting.

    Args:
        input_file: RepeatMasker .out file
        output_file: Output file path
        chrom_alias: UCSC chromAlias.txt for chromosome name conversion.
            If None, chromosome names are kept as-is.
        output_format: Output format ('ucsc' or 'bed')

    Raises:
        FileNotFoundError: If input file or alias file not found
        ValueError: If output_format is unsupported
    """
    if output_format not in ("ucsc", "bed"):
        raise ValueError(
            f"Unsupported output_format '{output_format}'. Use 'ucsc' or 'bed'."
        )

    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    # Load chromosome alias mapping if provided
    alias_map = None
    if chrom_alias is not None:
        alias_path = Path(chrom_alias)
        if not alias_path.exists():
            raise FileNotFoundError(f"Alias file not found: {chrom_alias}")
        alias_map = _load_chrom_alias(chrom_alias)

    # Ensure output directory exists
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Converting RepeatMasker: {input_file} -> {output_file}")
    logger.info(f"Output format: {output_format}")

    if output_format == "ucsc":
        n = _parse_rmsk_out_to_ucsc(input_file, output_file, alias_map)
    else:
        n = _parse_rmsk_out_to_bed(input_file, output_file, alias_map)

    logger.info(f"Output saved: {output_file} ({n:,} entries)")
