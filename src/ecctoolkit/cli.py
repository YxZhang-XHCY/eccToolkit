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
@click.option("--skip-fastp", is_flag=True, help="Skip fastp QC and use input FASTQ as-is")
@click.option("--min-mapq", default=0, type=int, help="Min MAPQ for discordant BAM (default: 0)")
@click.option("--auto-index", is_flag=True, help="Auto-generate BWA/samtools indices if missing")
@click.option(
    "--keep-intermediate/--cleanup-intermediate",
    default=True,
    help="Keep intermediate BAMs (default: keep)",
)
@click.option("-v", "--verbose", is_flag=True, help="Verbose output (debug logging)")
def circlemap(
    fastq1,
    fastq2,
    reference,
    output,
    sample,
    threads,
    skip_fastp,
    min_mapq,
    auto_index,
    keep_intermediate,
    verbose,
):
    """Run Circle-Map eccDNA detection pipeline.

    Complete pipeline: fastp QC -> BWA alignment -> Circle-Map detection.
    """
    from ecctoolkit.detect.circlemap import run_circlemap_pipeline
    run_circlemap_pipeline(
        fastq1=fastq1,
        fastq2=fastq2,
        reference=reference,
        output_dir=output,
        sample_name=sample,
        threads=threads,
        skip_fastp=skip_fastp,
        min_mapq=min_mapq,
        auto_index=auto_index,
        keep_intermediate=keep_intermediate,
        verbose=verbose,
    )


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
@click.option("-o", "--output", required=True, help="Output prefix (filename prefix, e.g., 'sample1')")
@click.option("-d", "--output-dir", default=None, help="Output directory (if not set, creates dir from -o)")
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
@click.option("--max-length-multi", default=10000, help="Maximum MeccDNA length")
@click.option("--identity", default=99.0, help="Identity threshold (%)")
@click.option("--min-coverage", default=90.0, help="Minimum coverage threshold (%)")
@click.option("--max-secondary", default=50, help="Max secondary alignments")
@click.option("--no-hit-policy", type=click.Choice(["skip", "unique"]), default="skip",
              help="Handle no-hit regions: skip or count as unique")
@click.option("--split-by-length/--no-split-by-length", default=False,
              help="Split candidates for short/long minimap2 presets (default: off)")
@click.option("--split-length", default=10000, help="Length threshold for preset split")
@click.option("--minimap-preset-short", default="map-hifi", help="Minimap2 preset for short (default: map-hifi)")
@click.option("--minimap-preset-long", default="map-hifi", help="Minimap2 preset for long (default: map-hifi)")
@click.option("--multiplier-u", default=2.0, help="Candidate multiplier for Unique")
@click.option("--multiplier-m", default=10.0, help="Candidate multiplier for Multi")
@click.option("--keep-tmp", is_flag=True, help="Keep temporary minimap2 files")
@click.option("-v", "--verbose", is_flag=True, help="Verbose output")
@click.option("--generate-ml-data", is_flag=True, help="Generate ML training data (doubled sequences + alignment TSV)")
@click.option("--ml-identity", default=99.0, help="Identity threshold for ML training data (%)")
def sim_region(reference, output, output_dir, num_unique, num_multi, num_chimeric,
               threads, seed, mode, sigma, tail_weight, tail_min, tail_max,
               min_length, max_length, max_length_multi, identity, min_coverage, max_secondary,
               no_hit_policy, split_by_length, split_length, minimap_preset_short,
               minimap_preset_long, multiplier_u, multiplier_m, keep_tmp, verbose,
               generate_ml_data, ml_identity):
    """Simulate eccDNA regions from reference genome.

    Generates three types:
    - UeccDNA: Unique mapping (one genomic location)
    - MeccDNA: Multi-mapping (multiple locations)
    - CeccDNA: Chimeric (2-5 fragments joined)

    Uses minimap2 for classification.

    \b
    Examples:
      # Output to specific directory (recommended for scripts)
      ecc sim-region -r hg38.fa -o sample1 -d /path/to/output/

      # Auto-create directory from prefix
      ecc sim-region -r hg38.fa -o /path/to/output/sample1
    """
    from ecctoolkit.simulate.region import run_region_simulation
    run_region_simulation(
        reference=reference,
        output_prefix=output,
        output_dir=output_dir,
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
        max_length_multi=max_length_multi,
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
        generate_ml_data=generate_ml_data,
        ml_identity_threshold=ml_identity,
    )


# ============================================================================
# Read Simulation Command (readsim)
# ============================================================================

@main.command("readsim")
# === Required parameters ===
@click.option("--sample", required=True, help="Sample name (output prefix)")
@click.option("-o", "--output", "path", required=True, help="Output directory")
@click.option("-t", "--threads", required=True, type=int, help="Number of threads")
# === Common optional parameters ===
@click.option("--meancov", default=25.0, help="Mean coverage (default=25)")
@click.option("--seed", type=int, help="Random seed")
# === Platform selection ===
@click.option("--skip-sr", is_flag=True, help="Skip short-read (NGS) simulation")
@click.option("--skip-hifi", is_flag=True, help="Skip HiFi simulation")
@click.option("--skip-ont", is_flag=True, help="Skip ONT simulation")
# === Input options ===
@click.option("-r", "--reference", help="Reference genome FASTA")
@click.option("-i", "--eccdna-fasta", help="eccDNA FASTA from sim-region (e.g., sample.all.fa)")
@click.option("--pool-csv", help="Pre-generated molecule pool CSV (skip seqsim/libsim)")
# === Sequence generation options ===
@click.option("--circular-number", default=5000, help="Circular DNA count (default=5000)")
@click.option("--linear-number", default=5000, help="Linear DNA count (default=5000)")
@click.option("--amp", default=5000, help="RCA amplification length bp (default=5000)")
@click.option("--simple-ratio", type=float, help="Simple eccDNA ratio")
@click.option("--simple-template", help="Simple eccDNA template BED")
@click.option("--chimeric-template", help="Chimeric eccDNA template BED")
# === NGS platform options ===
@click.option("--sr-platform", default="HS25", help="ART platform (default=HS25)")
@click.option("--sr-mean", default=400.0, help="Illumina insert mean (default=400)")
@click.option("--sr-std", default=125.0, help="Illumina insert std (default=125)")
@click.option("--sr-readlen", default=150.0, help="Illumina read length (default=150)")
# === ONT platform options ===
@click.option("--ont-model", default="R94", help="PBSIM2 ONT model (default=R94)")
@click.option("--ont-mean", default=3000.0, help="ONT mean read length (default=3000)")
@click.option("--ont-std", default=2500.0, help="ONT read length std (default=2500)")
# === HiFi platform options ===
@click.option("--hifi-sample-fastq", help="HiFi sample FASTQ for PBSIM2")
@click.option("--hifi-mode", default="auto", type=click.Choice(["auto", "sampling", "simple"]),
              help="HiFi simulation mode (default=auto)")
@click.option("--hifi-profile-id", help="HiFi profile ID for PBSIM2")
@click.option("--hifi-profile-root", help="HiFi profile directory")
@click.option("--hifi-len-min", default=5000, help="HiFi min length (default=5000)")
@click.option("--hifi-len-peak-min", default=10000, help="HiFi peak min (default=10000)")
@click.option("--hifi-len-peak-max", default=25000, help="HiFi peak max (default=25000)")
@click.option("--hifi-len-max", default=60000, help="HiFi max length (default=60000)")
@click.option("--hifi-qmin", default=20, help="HiFi min quality (default=20)")
@click.option("--hifi-qmean", default=30, help="HiFi mean quality (default=30)")
@click.option("--hifi-qsd", default=0.0, help="HiFi quality std (default=0.0)")
@click.option("--hifi-total-reads", type=int, help="Total HiFi reads (simple-mode)")
# === Output options ===
@click.option("-v", "--verbose", is_flag=True, help="Verbose output (debug logging)")
def readsim(sample, path, threads, meancov, seed, skip_sr, skip_hifi, skip_ont,
            reference, eccdna_fasta, pool_csv, circular_number, linear_number, amp,
            simple_ratio, simple_template, chimeric_template,
            sr_platform, sr_mean, sr_std, sr_readlen,
            ont_model, ont_mean, ont_std,
            hifi_sample_fastq, hifi_mode, hifi_profile_id, hifi_profile_root,
            hifi_len_min, hifi_len_peak_min, hifi_len_peak_max, hifi_len_max,
            hifi_qmin, hifi_qmean, hifi_qsd, hifi_total_reads,
            verbose):
    """Simulate sequencing reads from eccDNA.

    Uses external tools (ART, PBSIM2) for read simulation.

    \b
    Pipeline: seqsim → libsim → fqsim
    External tools required: art_illumina, pbsim, seqkit

    \b
    Examples:
      # Use eccDNA from sim-region (recommended)
      ecc readsim -r hg38.fa -i sample.all.fa -o output --sample test -t 4

      # Random generation with reference genome
      ecc readsim -r hg38.fa -o output --sample test -t 4

      # Use pre-generated molecule pool
      ecc readsim --pool-csv pool.csv -o output --sample test -t 4

      # Skip certain platforms
      ecc readsim -r hg38.fa -i sample.all.fa -o output --sample test -t 4 --skip-ont
    """
    import logging
    import os

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    logger = logging.getLogger(__name__)

    if not pool_csv and not reference:
        raise click.ClickException("-r/--reference is required unless --pool-csv is provided")

    from ecctoolkit.simulate.readsim import seqsim, libsim, fqsim

    # 所有输出都放在 path/sample/ 目录下
    sample_dir = os.path.join(path, sample)
    os.makedirs(sample_dir, exist_ok=True)

    if pool_csv:
        csv_path = pool_csv
    elif eccdna_fasta:
        # 使用 sim-region 生成的 eccDNA FASTA
        logger.info(f"Using eccDNA FASTA from sim-region: {eccdna_fasta}")
        from pathlib import Path
        from ecctoolkit.simulate.pipeline import SimulatePipeline

        # 复用 pipeline 的转换逻辑
        eccdna_path = Path(eccdna_fasta)
        if not eccdna_path.exists():
            raise click.ClickException(f"eccDNA FASTA not found: {eccdna_fasta}")

        # 创建一个临时的 pipeline 来执行转换
        class _TempConfig:
            def __init__(self):
                self.sample = sample
                self.seed = seed
                self.threads = threads

        temp_config = _TempConfig()
        temp_pipeline = SimulatePipeline.__new__(SimulatePipeline)
        temp_pipeline.config = temp_config
        temp_pipeline.reference = reference
        temp_pipeline.sequencing_dir = Path(sample_dir)

        # 转换格式
        temp_pipeline._convert_region_format(eccdna_path)

        # 运行 libsim (RCA 扩增)
        lib = libsim(
            sample=sample,
            reference=reference,
            path=sample_dir,
            seed=seed,
            meancov=meancov,
            amp=amp,
            threads=threads,
        )
        csv_path = os.path.join(sample_dir, sample + '.lib.csv')
    else:
        # Run seqsim (随机生成 eccDNA 和线性 DNA)
        logger.info("Generating random eccDNA from reference genome...")
        seq = seqsim(
            sample=sample,
            reference=reference,
            path=path,
            linear_number=linear_number,
            circular_number=circular_number,
            seed=seed,
            simple_ratio=simple_ratio,
            simple_template=simple_template,
            chimeric_template=chimeric_template,
        )

        # Run libsim (RCA 扩增，合并分子库)
        lib = libsim(
            sample=sample,
            reference=reference,
            path=sample_dir,
            seed=seed,
            meancov=meancov,
            amp=amp,
            threads=threads,
        )
        csv_path = os.path.join(sample_dir, sample + '.lib.csv')

    # Run fqsim (生成 FASTQ reads)
    fqsim(
        sample=sample,
        csv=csv_path,
        path=sample_dir,
        seed=seed,
        thread=threads,
        skip_sr=skip_sr,
        skip_hifi=skip_hifi,
        skip_ont=skip_ont,
        ont_model=ont_model,
        ont_mean=ont_mean,
        ont_std=ont_std,
        hifi_sample_fastq=hifi_sample_fastq,
        hifi_mode=hifi_mode,
        hifi_profile_id=hifi_profile_id,
        hifi_profile_root=hifi_profile_root,
        hifi_len_min=hifi_len_min,
        hifi_len_peak_min=hifi_len_peak_min,
        hifi_len_peak_max=hifi_len_peak_max,
        hifi_len_max=hifi_len_max,
        hifi_qmin=hifi_qmin,
        hifi_qmean=hifi_qmean,
        hifi_qsd=hifi_qsd,
        hifi_total_reads=hifi_total_reads,
        sr_platform=sr_platform,
        sr_mean=sr_mean,
        sr_std=sr_std,
        sr_readlen=sr_readlen,
    )

    logger.info("Read simulation completed!")


# ============================================================================
# Unified Simulation Command (simulate)
# ============================================================================

@main.command("simulate")
# === Required parameters ===
@click.option("-r", "--reference", required=True,
              type=click.Path(exists=True),
              help="Reference genome FASTA file")
@click.option("-o", "--output", required=True,
              type=click.Path(),
              help="Output directory")
# === eccDNA counts ===
@click.option("-u", "--num-unique", default=1000, type=int,
              help="Number of UeccDNA (Unique) to generate (default: 1000)")
@click.option("-m", "--num-multi", default=0, type=int,
              help="Number of MeccDNA (Multi-mapped) to generate (default: 0)")
@click.option("-c", "--num-chimeric", default=0, type=int,
              help="Number of CeccDNA (Chimeric) to generate (default: 0)")
# === Read simulation control ===
@click.option("--skip-readsim", is_flag=True,
              help="Skip read simulation (only generate eccDNA regions)")
# === Coverage ===
@click.option("--cov", "--coverage", "coverage", default=(30.0,), type=float,
              multiple=True,
              help="Target coverage(s) for read simulation (default: 30). "
                   "Use multiple times for multi-coverage: --cov 10 --cov 25 --cov 50")
# === Execution options ===
@click.option("-t", "--threads", default=8, type=int,
              help="Number of threads (default: 8)")
@click.option("--seed", default=42, type=int,
              help="Random seed (default: 42)")
@click.option("--sample", default="sim_ecc",
              help="Sample name/output prefix (default: sim_ecc)")
# === Platform selection ===
@click.option("--skip-sr", is_flag=True,
              help="Skip short-read (NGS) simulation")
@click.option("--skip-hifi", is_flag=True,
              help="Skip HiFi simulation")
@click.option("--skip-ont", is_flag=True,
              help="Skip ONT simulation")
# === Advanced options ===
@click.option("--config", "config_file", type=click.Path(exists=True),
              help="Optional YAML config file (CLI options override config)")
@click.option("--skip-region", is_flag=True,
              help="Skip region generation, use --input-eccdna instead")
@click.option("--input-eccdna", type=click.Path(exists=True),
              help="Use existing eccDNA FASTA (implies --skip-region)")
@click.option("--compress", is_flag=True,
              help="Compress output files (gzip)")
@click.option("--dry-run", is_flag=True,
              help="Show execution plan without running")
@click.option("-v", "--verbose", is_flag=True,
              help="Verbose output")
def simulate(reference, output, num_unique, num_multi, num_chimeric,
             skip_readsim, coverage, threads, seed, sample,
             skip_sr, skip_hifi, skip_ont,
             config_file, skip_region, input_eccdna, compress, dry_run, verbose):
    """Run complete eccDNA simulation pipeline.

    Combines region generation (sim-region) with read simulation (readsim).

    \b
    Examples:
      # Basic usage: 1000 UeccDNA with read simulation
      ecc simulate -r hg38.fa -o output/

      # Generate 10000 UeccDNA, 1000 MeccDNA, 1000 CeccDNA
      ecc simulate -r hg38.fa -o output/ -u 10000 -m 1000 -c 1000

      # Custom coverage
      ecc simulate -r hg38.fa -o output/ -u 5000 --cov 50

      # Only generate eccDNA regions (no reads)
      ecc simulate -r hg38.fa -o output/ -u 5000 --skip-readsim

      # Use existing eccDNA FASTA
      ecc simulate -r hg38.fa -o output/ --input-eccdna my_eccdna.fa

      # Preview execution plan
      ecc simulate -r hg38.fa -o output/ -u 10000 -m 1000 -c 1000 --dry-run
    """
    import logging

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    logger = logging.getLogger(__name__)

    from ecctoolkit.simulate.unified_config import UnifiedSimulateConfig
    from ecctoolkit.simulate.pipeline import SimulatePipeline, show_simulation_plan

    # Load config from file or use defaults
    if config_file:
        logger.info(f"Loading config from: {config_file}")
        config = UnifiedSimulateConfig.from_yaml(config_file)
    else:
        config = UnifiedSimulateConfig()

    # Apply CLI parameters (override config)
    config.sample = sample
    config.threads = threads
    config.seed = seed

    # Region parameters
    config.region.num_unique = num_unique
    config.region.num_multi = num_multi
    config.region.num_chimeric = num_chimeric

    # Readsim mode
    config.readsim.mode = "none" if skip_readsim else "on"

    # Coverage setting - convert tuple to list, use first for single mode
    coverage_list = list(coverage)
    config.readsim.params.meancov = coverage_list[0]  # First coverage for config

    # Handle input_eccdna
    if input_eccdna:
        skip_region = True

    # Validate config
    warnings = config.validate()
    if warnings:
        for w in warnings:
            logger.warning(f"Config warning: {w}")

    # Dry run - show plan and exit
    if dry_run:
        show_simulation_plan(config, output, reference, skip_region, input_eccdna, coverage_list)
        return

    # Run pipeline
    pipeline = SimulatePipeline(
        config=config,
        output_dir=output,
        reference=reference,
        skip_region=skip_region,
        input_eccdna=input_eccdna,
        skip_sr=skip_sr,
        skip_hifi=skip_hifi,
        skip_ont=skip_ont,
        compress=compress,
        verbose=verbose,
        coverage_list=coverage_list,
    )
    pipeline.run()

    logger.info("Simulation pipeline completed successfully!")


# ============================================================================
# Benchmark Command
# ============================================================================

@main.command("benchmark")
@click.option("-r", "--result", required=True,
              type=click.Path(exists=True),
              help="Detection result file (CircleSeeker merged_output.csv)")
@click.option("--truth-dir",
              type=click.Path(exists=True),
              help="Directory containing truth BED files (*.unique.bed, *.multi.bed, *.chimeric.bed)")
@click.option("--truth-unique",
              type=click.Path(exists=True),
              help="UeccDNA ground truth BED file")
@click.option("--truth-multi",
              type=click.Path(exists=True),
              help="MeccDNA ground truth BED file")
@click.option("--truth-chimeric",
              type=click.Path(exists=True),
              help="CeccDNA ground truth BED file")
@click.option("-o", "--output",
              type=click.Path(),
              help="Output directory for reports (default: current directory)")
@click.option("--prefix", default="benchmark",
              help="Output file prefix (default: benchmark)")
@click.option("--overlap", default=0.8, type=float,
              help="Minimum overlap ratio for matching (default: 0.8)")
@click.option("--format", "output_format",
              type=click.Choice(["text", "json", "csv", "all"]),
              default="all",
              help="Output format (default: all)")
@click.option("-v", "--verbose", is_flag=True,
              help="Verbose output")
def benchmark(result, truth_dir, truth_unique, truth_multi, truth_chimeric,
              output, prefix, overlap, output_format, verbose):
    """Benchmark eccDNA detection results against ground truth.

    Compare CircleSeeker (or other tool) outputs against simulated ground truth
    to calculate precision, recall, F1 scores, and other metrics.

    \b
    Examples:
      # Using truth directory (auto-detect BED files)
      ecc benchmark -r merged_output.csv --truth-dir simulation_output/

      # Using individual BED files
      ecc benchmark -r merged_output.csv \\
          --truth-unique sim.unique.bed \\
          --truth-multi sim.multi.bed \\
          --truth-chimeric sim.chimeric.bed

      # Custom output
      ecc benchmark -r merged_output.csv --truth-dir sim/ -o reports/ --prefix my_test
    """
    import logging
    import os

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    logger = logging.getLogger(__name__)

    # Validate inputs
    if not truth_dir and not (truth_unique or truth_multi or truth_chimeric):
        raise click.ClickException(
            "Either --truth-dir or at least one of --truth-unique/--truth-multi/--truth-chimeric is required"
        )

    from ecctoolkit.benchmark.evaluator import BenchmarkEvaluator
    from ecctoolkit.benchmark.report import BenchmarkReporter

    # Create evaluator
    evaluator = BenchmarkEvaluator(overlap_threshold=overlap)

    # Load truth data
    evaluator.load_truth(
        unique_bed=truth_unique,
        multi_bed=truth_multi,
        chimeric_bed=truth_chimeric,
        truth_dir=truth_dir,
    )

    # Load detection results
    evaluator.load_detected(result)

    # Run evaluation
    logger.info("Running benchmark evaluation...")
    metrics = evaluator.evaluate()

    # Generate reports
    reporter = BenchmarkReporter(metrics)

    # Print summary to console
    reporter.print_summary()

    # Save reports
    if output:
        os.makedirs(output, exist_ok=True)
        output_dir = output
    else:
        output_dir = "."

    if output_format in ["text", "all"]:
        txt_path = os.path.join(output_dir, f"{prefix}_report.txt")
        reporter.save_text_report(txt_path)

    if output_format in ["json", "all"]:
        json_path = os.path.join(output_dir, f"{prefix}_report.json")
        reporter.save_json_report(json_path)

    if output_format in ["csv", "all"]:
        csv_path = os.path.join(output_dir, f"{prefix}_summary.csv")
        reporter.save_csv_summary(csv_path)

    logger.info("Benchmark completed!")


if __name__ == "__main__":
    main()
