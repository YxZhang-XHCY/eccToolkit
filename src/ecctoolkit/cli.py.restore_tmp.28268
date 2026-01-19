"""
eccToolkit CLI - Command Line Interface for eccDNA analysis.

Usage:
    ecc <command> [options]

Each command is an independent analysis tool.
"""

import click

from ecctoolkit import __version__


@click.group()
@click.version_option(version=__version__, prog_name="eccToolkit")
def main():
    """eccToolkit - A comprehensive toolkit for eccDNA analysis.

    Each command is an independent analysis tool. Use 'ecc <command> --help'
    for detailed usage of each command.
    """
    pass


# ============================================================================
# Detection & Validation Commands
# ============================================================================

@main.command()
@click.option("-1", "--fastq1", required=True, help="Input FASTQ file (Read 1)")
@click.option("-2", "--fastq2", required=True, help="Input FASTQ file (Read 2)")
@click.option("-r", "--reference", required=True, help="Reference genome FASTA")
@click.option("-o", "--output", required=True, help="Output directory")
@click.option("-s", "--sample", required=True, help="Sample name")
@click.option("-t", "--threads", default=8, help="Number of threads")
def circlemap(fastq1, fastq2, reference, output, sample, threads):
    """Run Circle-Map eccDNA detection pipeline.

    Complete pipeline: fastp QC -> BWA alignment -> Circle-Map detection.
    """
    from ecctoolkit.detect.circlemap import run_circlemap_pipeline
    run_circlemap_pipeline(fastq1, fastq2, reference, output, sample, threads)


@main.command()
@click.option("-i", "--input", "input_file", required=True, help="eccDNA candidates TSV")
@click.option("-1", "--fastq1", required=True, help="Amplicon reads R1")
@click.option("-2", "--fastq2", required=True, help="Amplicon reads R2")
@click.option("-r", "--reference", required=True, help="Reference genome FASTA")
@click.option("-o", "--output", required=True, help="Output directory")
@click.option("-t", "--threads", default=8, help="Number of threads")
def validate(input_file, fastq1, fastq2, reference, output, threads):
    """Validate eccDNA candidates from amplicon sequencing.

    Builds junction-centered reference and aligns reads to verify eccDNA.
    """
    from ecctoolkit.detect.validate import run_validation
    run_validation(input_file, fastq1, fastq2, reference, output, threads)


@main.command("find-repeats")
@click.option("-i", "--input", "input_file", required=True, help="eccDNA candidates TSV")
@click.option("-g", "--genome", required=True, help="Reference genome FASTA")
@click.option("-o", "--output", required=True, help="Output CSV file")
@click.option("-t", "--threads", default=8, help="Number of threads")
@click.option("--tolerance", default=5, help="Breakpoint tolerance (bp)")
def find_repeats(input_file, genome, output, threads, tolerance):
    """Find terminal repeats at eccDNA breakpoints using BLAST."""
    from ecctoolkit.detect.repeats import find_terminal_repeats
    find_terminal_repeats(input_file, genome, output, threads, tolerance)


@main.command()
@click.option("-i", "--input", "input_file", required=True, help="Assembled reads FASTA")
@click.option("-r", "--reference", required=True, help="Reference genome FASTA")
@click.option("-o", "--output", required=True, help="Output directory")
@click.option("-p", "--prefix", required=True, help="Sample prefix")
@click.option("-t", "--threads", default=8, help="Number of threads")
def saturation(input_file, reference, output, prefix, threads):
    """Generate saturation curve by subsampling and running CircleSeeker."""
    from ecctoolkit.detect.saturation import run_saturation
    run_saturation(input_file, reference, output, prefix, threads)


# ============================================================================
# Enrichment Analysis Commands (Permutation Tests)
# ============================================================================

@main.command("enrich-cnv")
@click.option("-i", "--input", "input_files", required=True, multiple=True, help="eccDNA CSV files")
@click.option("--cnv", required=True, help="CNV regions BED file")
@click.option("-o", "--output", required=True, help="Output directory")
@click.option("-g", "--genome", help="Genome sizes file (auto-generated if not provided)")
@click.option("--prefix", help="Output prefix")
@click.option("-n", "--nshuffle", default=1000, help="Number of permutations")
@click.option("-t", "--threads", default=8, help="Number of threads")
@click.option("--keep-sex", is_flag=True, help="Keep sex chromosomes")
def enrich_cnv(input_files, cnv, output, genome, prefix, nshuffle, threads, keep_sex):
    """Analyze eccDNA enrichment in CNV gain/loss/neutral regions.

    Uses permutation test to assess significance.
    """
    from ecctoolkit.enrich.cnv import run_cnv_enrichment
    run_cnv_enrichment(list(input_files), cnv, output, genome, prefix, nshuffle, threads, keep_sex)


@main.command("enrich-tad")
@click.option("-i", "--input", "input_files", required=True, multiple=True, help="eccDNA CSV files")
@click.option("--tad", required=True, help="TAD boundaries BED file")
@click.option("-o", "--output", required=True, help="Output directory")
@click.option("-g", "--genome", help="Genome sizes file")
@click.option("--prefix", help="Output prefix")
@click.option("-n", "--nshuffle", default=1000, help="Number of permutations")
@click.option("-t", "--threads", default=8, help="Number of threads")
@click.option("--keep-sex", is_flag=True, help="Keep sex chromosomes")
def enrich_tad(input_files, tad, output, genome, prefix, nshuffle, threads, keep_sex):
    """Analyze eccDNA enrichment at TAD boundaries.

    Uses permutation test to assess significance.
    """
    from ecctoolkit.enrich.tad import run_tad_enrichment
    run_tad_enrichment(list(input_files), tad, output, genome, prefix, nshuffle, threads, keep_sex)


@main.command("enrich-overlap")
@click.option("-i", "--input", "input_file", required=True, help="Query regions BED/CSV")
@click.option("-a", "--annotation", required=True, multiple=True, help="Annotation BED files")
@click.option("-o", "--output", required=True, help="Output directory")
@click.option("-g", "--genome", required=True, help="Genome sizes file")
@click.option("-n", "--nshuffle", default=1000, help="Number of permutations")
@click.option("-t", "--threads", default=8, help="Number of threads")
def enrich_overlap(input_file, annotation, output, genome, nshuffle, threads):
    """General permutation test for genomic region overlap."""
    from ecctoolkit.enrich.overlap import run_overlap_test
    run_overlap_test(input_file, list(annotation), output, genome, nshuffle, threads)


# ============================================================================
# Hotspot Analysis Commands
# ============================================================================

@main.command("hotspot-matrix")
@click.option("-i", "--input", "input_file", required=True, help="eccDNA CSV file")
@click.option("-f", "--fai", required=True, help="Reference genome FAI index")
@click.option("-o", "--output", required=True, help="Output directory")
@click.option("--windows", default="10000,50000,100000", help="Window sizes (comma-separated)")
def hotspot_matrix(input_file, fai, output, windows):
    """Generate multi-scale window count matrices using bedtools.

    Creates eccDNA density matrices at different resolutions (default: 10kb, 50kb, 100kb).
    """
    from ecctoolkit.hotspot.matrix import generate_multiscale_matrix
    window_list = [int(w) for w in windows.split(",")]
    generate_multiscale_matrix(input_file, fai, output, window_list)


@main.command("hotspot-detect")
@click.option("-i", "--input", "input_dir", required=True, help="Multi-scale matrix directory")
@click.option("-e", "--eccdna", required=True, help="Original eccDNA CSV")
@click.option("-o", "--output", required=True, help="Output directory")
def hotspot_detect(input_dir, eccdna, output):
    """Detect hotspots with precise boundaries (SHARP method).

    Identifies candidate regions from multi-scale data and refines boundaries.
    """
    from ecctoolkit.hotspot.detect import run_sharp_detection
    run_sharp_detection(input_dir, eccdna, output)


@main.command("hotspot-test")
@click.option("-i", "--input", "input_dir", required=True, help="Multi-scale matrix directory")
@click.option("-o", "--output", required=True, help="Output directory")
@click.option("-n", "--nshuffle", default=1000, help="Number of permutations")
@click.option("-t", "--threads", default=8, help="Number of threads")
def hotspot_test(input_dir, output, nshuffle, threads):
    """Test hotspot significance using permutation (parallel version)."""
    from ecctoolkit.hotspot.permtest import run_hotspot_permtest
    run_hotspot_permtest(input_dir, output, nshuffle, threads)


@main.command("hotspot-refine")
@click.option("-i", "--input", "input_file", required=True, help="Candidate hotspots file")
@click.option("-e", "--eccdna", required=True, help="Original eccDNA CSV")
@click.option("-o", "--output", required=True, help="Output file")
def hotspot_refine(input_file, eccdna, output):
    """Refine hotspot boundaries at high resolution."""
    from ecctoolkit.hotspot.refine import refine_hotspot_boundaries
    refine_hotspot_boundaries(input_file, eccdna, output)


# ============================================================================
# Transposable Element (TE) Analysis Commands
# ============================================================================

@main.command("te-analyze")
@click.option("-i", "--input", "input_file", required=True, help="GFF annotation file")
@click.option("-o", "--output", required=True, help="Output directory")
@click.option("-t", "--threads", default=8, help="Number of threads")
def te_analyze(input_file, output, threads):
    """Analyze TE composition from GFF, compare Mecc vs Uecc."""
    from ecctoolkit.te.analyze import run_te_analysis
    run_te_analysis(input_file, output, threads)


@main.command("te-classify")
@click.option("-i", "--input", "input_file", required=True, help="TE annotation CSV")
@click.option("-o", "--output", required=True, help="Output directory")
def te_classify(input_file, output):
    """Classify eccDNA as single or composite TE based on motif count."""
    from ecctoolkit.te.classify import classify_te_composition
    classify_te_composition(input_file, output)


@main.command("te-distribution")
@click.option("-i", "--input", "input_file", required=True, help="TE annotation CSV")
@click.option("-o", "--output", required=True, help="Output directory")
def te_distribution(input_file, output):
    """Analyze TE percentage distribution in 5 bins (0-20%, 20-40%, etc.)."""
    from ecctoolkit.te.distribution import analyze_te_distribution
    analyze_te_distribution(input_file, output)


@main.command("te-composition")
@click.option("-i", "--input", "input_file", required=True, help="TE annotation CSV")
@click.option("-o", "--output", required=True, help="Output directory")
def te_composition(input_file, output):
    """Analyze composite TE composition and motif combinations."""
    from ecctoolkit.te.composition import analyze_te_composition
    analyze_te_composition(input_file, output)


@main.command("te-process")
@click.option("-i", "--input", "input_file", required=True, help="TE annotation CSV")
@click.option("-r", "--reference", help="Reference table for seq_length")
@click.option("-o", "--output", required=True, help="Output CSV file")
def te_process(input_file, reference, output):
    """Process TE data: fill missing values, recalculate percentages."""
    from ecctoolkit.te.process import process_te_data
    process_te_data(input_file, reference, output)


# ============================================================================
# Expression Correlation Commands
# ============================================================================

@main.command("expr-correlate")
@click.option("-i", "--input", "input_file", required=True, help="eccDNA enrichment CSV")
@click.option("-d", "--deg", required=True, help="DEG results TSV")
@click.option("-o", "--output", required=True, help="Output directory")
@click.option("--mode", type=click.Choice(["single", "gradient"]), default="gradient",
              help="Analysis mode: single threshold or gradient FC")
def expr_correlate(input_file, deg, output, mode):
    """Analyze correlation between eccDNA enrichment and DEGs (Fisher test)."""
    from ecctoolkit.expression.correlate import run_expression_correlation
    run_expression_correlation(input_file, deg, output, mode)


@main.command("expr-enrich")
@click.option("-i", "--input", "input_file", required=True, help="Gene list CSV")
@click.option("-d", "--deg", required=True, help="DEG results TSV")
@click.option("-o", "--output", required=True, help="Output CSV")
def expr_enrich(input_file, deg, output):
    """Analyze DEG enrichment at multiple FC thresholds."""
    from ecctoolkit.expression.enrich import run_deg_enrichment
    run_deg_enrichment(input_file, deg, output)


# ============================================================================
# Visualization Commands
# ============================================================================

@main.command()
@click.option("-i", "--input", "input_file", required=True, help="eccDNA CSV with eName")
@click.option("-o", "--output", required=True, help="Output TSV file")
@click.option("--chromosomes", help="Chromosomes to include (comma-separated)")
@click.option("--color-mode", type=click.Choice(["random", "highlight-max"]), default="random")
def genolink(input_file, output, chromosomes, color_mode):
    """Generate genomic link data for Circos-style visualization.

    Groups by eName and creates pairwise links with RGB colors.
    """
    from ecctoolkit.visualize.genolink import generate_genolink
    chrom_list = chromosomes.split(",") if chromosomes else None
    generate_genolink(input_file, output, chrom_list, color_mode)


# ============================================================================
# Data Processing Commands
# ============================================================================

@main.command()
@click.option("--fled", required=True, help="FLED output directory")
@click.option("--circlemap", required=True, help="CircleMap output directory")
@click.option("-o", "--output", required=True, help="Output directory")
def combine(fled, circlemap, output):
    """Combine and merge FLED and CircleMap detection results."""
    from ecctoolkit.process.combine import combine_fled_circlemap
    combine_fled_circlemap(fled, circlemap, output)


@main.command()
@click.option("-i", "--input", "input_dir", required=True, help="Input directory")
@click.option("-o", "--output", required=True, help="Output CSV file")
@click.option("--pattern", default="*.csv", help="File pattern to match")
def merge(input_dir, output, pattern):
    """Merge multiple CSV files and add sample column."""
    from ecctoolkit.process.merge import merge_files
    merge_files(input_dir, output, pattern)


@main.command()
@click.option("-i", "--input", "input_file", required=True, help="Input eccDNA CSV")
@click.option("-o", "--output", required=True, help="Output CSV file")
def parse(input_file, output):
    """Parse eccDNA CSV: extract chr/start/end from seqname, calculate motif_percent."""
    from ecctoolkit.process.parse import parse_eccdna
    parse_eccdna(input_file, output)


@main.command()
@click.option("-i", "--input", "input_file", required=True, help="Input CSV file")
@click.option("-o", "--output", required=True, help="Output CSV file")
@click.option("--min-percent", default=80.0, help="Minimum anno_Percent threshold")
def filter(input_file, output, min_percent):
    """Filter eccDNA by annotation percentage (default: >= 80%)."""
    from ecctoolkit.process.filter import filter_by_annotation
    filter_by_annotation(input_file, output, min_percent)


@main.command("convert-gff")
@click.option("-i", "--input", "input_dir", required=True, help="Directory with GFF files")
@click.option("-o", "--output", required=True, help="Output CSV file")
def convert_gff(input_dir, output):
    """Batch convert GFF files to merged CSV."""
    from ecctoolkit.process.convert import convert_gff_to_csv
    convert_gff_to_csv(input_dir, output)


@main.command("convert-cnv")
@click.option("-i", "--input", "input_file", required=True, help="CNV TSV file (1-based)")
@click.option("-o", "--output", required=True, help="Output directory")
def convert_cnv(input_file, output):
    """Convert CNV TSV to per-cell-line BED files (0-based)."""
    from ecctoolkit.process.convert import convert_cnv_to_bed
    convert_cnv_to_bed(input_file, output)


@main.command()
@click.option("-i", "--input", "input_file", required=True, help="Combined CSV file")
@click.option("-o", "--output", required=True, help="Output prefix")
def report(input_file, output):
    """Generate summary reports (sample counts, length distribution)."""
    from ecctoolkit.process.report import generate_reports
    generate_reports(input_file, output)


@main.command()
@click.option("-i", "--input", "input_file", required=True, help="eccDNA CSV with eLength")
@click.option("-o", "--output", required=True, help="Output directory")
def quantify(input_file, output):
    """Quantify eccDNA length distribution (periodicity, ACF/FFT analysis)."""
    from ecctoolkit.process.quantify import run_quantification
    run_quantification(input_file, output)


# ============================================================================
# Simulation Commands
# ============================================================================

@main.command("sim-region")
@click.option("-r", "--reference", required=True, help="Reference genome FASTA")
@click.option("-o", "--output", required=True, help="Output prefix")
@click.option("-u", "--num-unique", default=1000, help="Number of UeccDNA (Unique)")
@click.option("-m", "--num-multi", default=0, help="Number of MeccDNA (Multi-mapped)")
@click.option("-c", "--num-chimeric", default=0, help="Number of CeccDNA (Chimeric)")
@click.option("-t", "--threads", default=8, help="Number of threads")
@click.option("--seed", default=42, help="Random seed")
@click.option("--mode", default=400.0, help="Lognormal mode/peak in bp")
@click.option("--sigma", default=0.8, help="Lognormal sigma")
@click.option("--tail-weight", default=0.05, help="Tail distribution weight")
@click.option("--tail-min", default=5000, help="Minimum length for tail distribution")
@click.option("--tail-max", default=500000, help="Maximum length for tail distribution")
@click.option("--min-length", default=100, help="Minimum eccDNA length")
@click.option("--max-length", default=500000, help="Maximum eccDNA length")
@click.option("--identity", default=99.0, help="Identity threshold (%)")
@click.option("--min-coverage", default=90.0, help="Minimum coverage threshold (%)")
@click.option("--max-secondary", default=50, help="Max secondary alignments")
@click.option("--no-hit-policy", type=click.Choice(["skip", "unique"]), default="skip",
              help="Handle no-hit regions: skip or count as unique")
@click.option("--split-by-length/--no-split-by-length", default=True,
              help="Split candidates for short/long minimap2 presets")
@click.option("--split-length", default=5000, help="Length threshold for preset split")
@click.option("--minimap-preset-short", default="sr", help="Minimap2 preset for short")
@click.option("--minimap-preset-long", default="asm5", help="Minimap2 preset for long")
@click.option("--multiplier-u", default=2.0, help="Candidate multiplier for Unique")
@click.option("--multiplier-m", default=10.0, help="Candidate multiplier for Multi")
@click.option("--keep-tmp", is_flag=True, help="Keep temporary minimap2 files")
@click.option("-v", "--verbose", is_flag=True, help="Verbose output")
def sim_region(reference, output, num_unique, num_multi, num_chimeric,
               threads, seed, mode, sigma, tail_weight, tail_min, tail_max,
               min_length, max_length, identity, min_coverage, max_secondary,
               no_hit_policy, split_by_length, split_length, minimap_preset_short,
               minimap_preset_long, multiplier_u, multiplier_m, keep_tmp, verbose):
    """Simulate eccDNA regions from reference genome.

    Generates three types:
    - UeccDNA: Unique mapping (one genomic location)
    - MeccDNA: Multi-mapping (multiple locations)
    - CeccDNA: Chimeric (2-5 fragments joined)

    Uses minimap2 for classification.
    """
    from ecctoolkit.simulate.region import run_region_simulation
    run_region_simulation(
        reference=reference,
        output_prefix=output,
        num_unique=num_unique,
        num_multi=num_multi,
        num_chimeric=num_chimeric,
        threads=threads,
        seed=seed,
        mode=mode,
        sigma=sigma,
        tail_weight=tail_weight,
        tail_min=tail_min,
        tail_max=tail_max,
        min_length=min_length,
        max_length=max_length,
        identity=identity,
        min_coverage=min_coverage,
        max_secondary=max_secondary,
        no_hit_policy=no_hit_policy,
        split_by_length=split_by_length,
        split_length=split_length,
        minimap_preset_short=minimap_preset_short,
        minimap_preset_long=minimap_preset_long,
        candidate_multiplier_u=multiplier_u,
        candidate_multiplier_m=multiplier_m,
        keep_tmp=keep_tmp,
        verbose=verbose,
    )


@main.command("sim-reads")
@click.option("-i", "--input", "input_file", required=True, help="eccDNA FASTA file")
@click.option("-o", "--output", required=True, help="Output directory")
@click.option("--cov-ngs", default=10000.0, help="NGS output (reads or coverage)")
@click.option("--cov-hifi", default=1000.0, help="HiFi output (reads or coverage)")
@click.option("--cov-ont", default=1000.0, help="ONT output (reads or coverage)")
@click.option("--output-mode", type=click.Choice(["read_count", "coverage"]),
              default="read_count", help="Output scale mode")
@click.option("--force-ecc-coverage/--no-force-ecc-coverage", default=True,
              help="Force each eccDNA to appear in reads when output is large")
@click.option("--force-ecc-min-reads-factor", default=2.0, type=float,
              help="Threshold factor for forcing eccDNA coverage")
@click.option("--rca-mode", type=click.Choice(["empirical", "kinetic"]),
              help="RCA repeat model")
@click.option("--reaction-hours", type=float,
              help="Reaction time in hours (kinetic mode)")
@click.option("--branch-length-mode", type=click.Choice(["ratio_trunk", "ratio_ecc"]),
              default="ratio_trunk", help="Branch length mode")
@click.option("--chimera-breakpoint-mode", type=click.Choice(["uniform", "beta"]),
              default="uniform", help="Chimera breakpoint sampling mode")
@click.option("--chimera-breakpoint-alpha", default=0.5, type=float,
              help="Chimera breakpoint alpha (beta mode)")
@click.option("--chimera-breakpoint-beta", default=0.5, type=float,
              help="Chimera breakpoint beta (beta mode)")
@click.option("-r", "--ref", "--reference", "reference",
              help="Reference genome FASTA for background linear DNA")
@click.option("--linear-ratio", default=0.0,
              help="Background linear DNA ratio (0-1, default 0=disabled)")
@click.option("--linear-min-len", default=200, help="Min background fragment length")
@click.option("--linear-max-len", default=10000, help="Max background fragment length")
@click.option("--chimera-background-ratio", default=0.0, type=float,
              help="Chimera donor fraction from background pool (0-1)")
@click.option("--chimera-background-pool", default=200, type=int,
              help="Background donor pool size for chimera")
@click.option("-t", "--threads", default=1, help="Number of threads (0=auto)")
@click.option("--seed", type=int, help="Random seed")
@click.option("--config", "config_file", help="Config file (YAML/JSON)")
@click.option("--skip-ngs", is_flag=True, help="Skip NGS generation")
@click.option("--skip-hifi", is_flag=True, help="Skip HiFi generation")
@click.option("--skip-ont", is_flag=True, help="Skip ONT generation")
@click.option("--compress", is_flag=True, help="Compress output (gzip)")
@click.option("-sample", "--sample", default="", help="Sample prefix for output files")
@click.option("-v", "--verbose", is_flag=True, help="Verbose output")
def sim_reads(input_file, output, cov_ngs, cov_hifi, cov_ont, output_mode,
              force_ecc_coverage, force_ecc_min_reads_factor,
              rca_mode, reaction_hours, branch_length_mode,
              chimera_breakpoint_mode, chimera_breakpoint_alpha, chimera_breakpoint_beta,
              reference, linear_ratio, linear_min_len, linear_max_len,
              chimera_background_ratio, chimera_background_pool,
              threads, seed, config_file, skip_ngs, skip_hifi,
              skip_ont, compress, sample, verbose):
    """Simulate sequencing reads from eccDNA using RCA model.

    Generates NGS (paired-end), HiFi, and ONT reads with realistic
    RCA (Rolling Circle Amplification) artifacts.

    Output modes:
    - read_count: Specify exact number of reads
    - coverage: Specify target coverage (auto-calculates reads)

    Background DNA:
    Use --ref and --linear-ratio to add background linear DNA reads
    that simulate genomic contamination (not from eccDNA).
    """
    from ecctoolkit.simulate.reads import run_read_simulation
    run_read_simulation(
        input_file=input_file,
        output_dir=output,
        cov_ngs=cov_ngs,
        cov_hifi=cov_hifi,
        cov_ont=cov_ont,
        output_mode=output_mode,
        force_ecc_coverage=force_ecc_coverage,
        force_ecc_min_reads_factor=force_ecc_min_reads_factor,
        rca_mode=rca_mode,
        reaction_hours=reaction_hours,
        branch_length_mode=branch_length_mode,
        chimera_breakpoint_mode=chimera_breakpoint_mode,
        chimera_breakpoint_alpha=chimera_breakpoint_alpha,
        chimera_breakpoint_beta=chimera_breakpoint_beta,
        reference=reference,
        linear_ratio=linear_ratio,
        linear_min_len=linear_min_len,
        linear_max_len=linear_max_len,
        chimera_background_ratio=chimera_background_ratio,
        chimera_background_pool=chimera_background_pool,
        threads=threads,
        seed=seed,
        config_file=config_file,
        skip_ngs=skip_ngs,
        skip_hifi=skip_hifi,
        skip_ont=skip_ont,
        compress=compress,
        sample=sample,
        verbose=verbose,
    )


if __name__ == "__main__":
    main()
