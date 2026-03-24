#!/bin/bash
# ============================================================================
# eccToolkit Simulation Benchmark Experiment
# For iMeta paper: CircleSeeker + eccToolkit
#
# Usage: bash scripts/benchmark_experiment.sh <REFERENCE> <OUTPUT_DIR> [THREADS]
#
# Requires: ecc (eccToolkit), minimap2, samtools
# Optional: CircleSeeker, CReSIL, Circle-Map, ecc_finder
#           (tools that are not installed will be skipped)
# ============================================================================

set -euo pipefail

# ============================================================================
# Arguments
# ============================================================================

if [[ $# -lt 2 ]]; then
    echo "Usage: bash $0 <REFERENCE> <OUTPUT_DIR> [THREADS]"
    echo ""
    echo "  REFERENCE   Reference genome FASTA file"
    echo "  OUTPUT_DIR  Base output directory for all experiments"
    echo "  THREADS     Number of threads (default: 8)"
    exit 1
fi

REFERENCE=$(realpath "$1")
OUTPUT_DIR=$(realpath "$2")
THREADS=${3:-8}
SEED=42

# ============================================================================
# Tool detection
# ============================================================================

declare -A TOOL_AVAILABLE

check_tool() {
    local name="$1"
    local cmd="$2"
    if command -v "$cmd" &>/dev/null; then
        TOOL_AVAILABLE[$name]=1
        echo "[OK] $name ($cmd)"
    else
        TOOL_AVAILABLE[$name]=0
        echo "[SKIP] $name ($cmd not found)"
    fi
}

echo "============================================================"
echo "Checking available tools..."
echo "============================================================"
check_tool "ecc" "ecc"
check_tool "CircleSeeker" "CircleSeeker"
check_tool "CReSIL" "CReSIL"
check_tool "Circle-Map" "Circle-Map"
check_tool "ecc_finder" "ecc_finder"
check_tool "minimap2" "minimap2"
check_tool "samtools" "samtools"
echo ""

# ecc (eccToolkit) is required
if [[ "${TOOL_AVAILABLE[ecc]}" -eq 0 ]]; then
    echo "ERROR: ecc (eccToolkit) is required but not found."
    exit 1
fi

# ============================================================================
# Helper functions
# ============================================================================

log_section() {
    echo ""
    echo "============================================================"
    echo "$1"
    echo "============================================================"
}

log_step() {
    echo ">>> $1"
}

# Format coverage value for directory naming (remove trailing zeros)
fmt_cov() {
    local cov="$1"
    # Remove trailing .0 for integer coverages
    echo "$cov" | sed 's/\.0$//'
}

# Run benchmark for a single tool against ground truth
run_benchmark() {
    local tool_name="$1"
    local result_file="$2"
    local truth_dir="$3"
    local output_dir="$4"
    local prefix="$5"

    if [[ ! -f "$result_file" ]]; then
        echo "  [WARN] Result file not found: $result_file (skipping $tool_name)"
        return
    fi

    log_step "Benchmarking $tool_name: $prefix"
    mkdir -p "$output_dir"
    ecc benchmark \
        -r "$result_file" \
        --truth-dir "$truth_dir" \
        -o "$output_dir" \
        --prefix "${prefix}" \
        --format all \
        -v
}

# Align reads to reference with minimap2
align_reads() {
    local platform="$1"
    local reads="$2"
    local output_bam="$3"
    local preset=""

    case "$platform" in
        hifi) preset="map-hifi" ;;
        ont)  preset="map-ont"  ;;
        ngs)  preset="sr"      ;;
        *)    echo "Unknown platform: $platform"; return 1 ;;
    esac

    log_step "Aligning $platform reads -> $output_bam"
    minimap2 -a -x "$preset" -t "$THREADS" "$REFERENCE" "$reads" \
        | samtools sort -@ "$THREADS" -o "$output_bam"
    samtools index -@ "$THREADS" "$output_bam"
}

# ============================================================================
# Detection tool wrappers
# Each function runs a specific tool and produces output in a standard location.
# Modify these functions to match your tool installations and parameters.
# ============================================================================

run_circleseeker() {
    local bam="$1"
    local outdir="$2"
    local platform="$3"

    if [[ "${TOOL_AVAILABLE[CircleSeeker]}" -eq 0 ]]; then
        echo "  [SKIP] CircleSeeker not available"
        return
    fi

    log_step "Running CircleSeeker ($platform)"
    mkdir -p "$outdir"
    # CircleSeeker command - adjust parameters as needed
    CircleSeeker -i "$bam" -o "$outdir" -r "$REFERENCE" -t "$THREADS"
}

run_cresil() {
    local bam="$1"
    local outdir="$2"

    if [[ "${TOOL_AVAILABLE[CReSIL]}" -eq 0 ]]; then
        echo "  [SKIP] CReSIL not available"
        return
    fi

    log_step "Running CReSIL"
    mkdir -p "$outdir"
    # CReSIL command - adjust parameters as needed
    CReSIL -i "$bam" -o "$outdir" -r "$REFERENCE" -t "$THREADS"
}

run_circle_map() {
    local bam="$1"
    local outdir="$2"

    if [[ "${TOOL_AVAILABLE[Circle-Map]}" -eq 0 ]]; then
        echo "  [SKIP] Circle-Map not available"
        return
    fi

    log_step "Running Circle-Map"
    mkdir -p "$outdir"
    # Circle-Map command - adjust parameters as needed
    Circle-Map ReadExtractor -i "$bam" -o "$outdir/candidates.bam"
    Circle-Map Realign \
        -i "$outdir/candidates.bam" \
        -qbam "$bam" \
        -sbam "$bam" \
        -fasta "$REFERENCE" \
        -o "$outdir/circle_map_output.bed" \
        -t "$THREADS"
}

run_ecc_finder() {
    local bam="$1"
    local outdir="$2"
    local platform="$3"

    if [[ "${TOOL_AVAILABLE[ecc_finder]}" -eq 0 ]]; then
        echo "  [SKIP] ecc_finder not available"
        return
    fi

    log_step "Running ecc_finder ($platform)"
    mkdir -p "$outdir"
    # ecc_finder command - adjust parameters as needed
    if [[ "$platform" == "ngs" ]]; then
        ecc_finder -i "$bam" -o "$outdir" -r "$REFERENCE" -t "$THREADS" --short-read
    else
        ecc_finder -i "$bam" -o "$outdir" -r "$REFERENCE" -t "$THREADS" --long-read
    fi
}


# ############################################################################
#
# EXPERIMENT 1: Standard Benchmark (Main Figure)
#
# Ground truth: 10000 U + 1000 M + 500 C
# Platforms: HiFi (30x), NGS (30x), ONT (30x)
# Tools: CircleSeeker, CReSIL, Circle-Map, ecc_finder
#
# ############################################################################

run_experiment_standard() {
    local EXP_DIR="${OUTPUT_DIR}/exp1_standard"
    local GT_DIR="${EXP_DIR}/ground_truth"
    local READS_DIR="${EXP_DIR}/reads"
    local ALIGN_DIR="${EXP_DIR}/alignments"
    local DETECT_DIR="${EXP_DIR}/detection"
    local BENCH_DIR="${EXP_DIR}/benchmark"

    log_section "Experiment 1: Standard Benchmark"

    # ------------------------------------------------------------------
    # Step 1: Generate ground truth eccDNA regions
    # ------------------------------------------------------------------
    log_step "Step 1: Generate ground truth (U=10000, M=1000, C=500)"

    if [[ -d "$GT_DIR" && -f "$GT_DIR/standard.all.fa" ]]; then
        echo "  Ground truth already exists, skipping."
    else
        ecc sim-region \
            -r "$REFERENCE" \
            -o standard \
            -d "$GT_DIR" \
            -u 10000 -m 1000 -c 500 \
            -t "$THREADS" \
            --seed "$SEED"
    fi

    # ------------------------------------------------------------------
    # Step 2: Simulate reads for each platform
    # ------------------------------------------------------------------
    log_step "Step 2: Simulate reads (30x coverage)"

    # Use the unified simulate command with --skip-region to just do read sim
    # But since we already have the eccDNA FASTA, use sim-reads directly

    local ECCDNA_FASTA="${GT_DIR}/standard.all.fa"

    # HiFi reads
    if [[ ! -d "${READS_DIR}/hifi" ]]; then
        log_step "  Simulating HiFi reads (30x)"
        ecc sim-reads \
            -r "$REFERENCE" \
            -i "$ECCDNA_FASTA" \
            -o "${READS_DIR}/hifi" \
            --sample hifi_30x \
            -t "$THREADS" \
            --meancov 30 \
            --skip-sr --skip-ont \
            --seed "$SEED"
    else
        echo "  HiFi reads already exist, skipping."
    fi

    # NGS reads
    if [[ ! -d "${READS_DIR}/ngs" ]]; then
        log_step "  Simulating NGS reads (30x)"
        ecc sim-reads \
            -r "$REFERENCE" \
            -i "$ECCDNA_FASTA" \
            -o "${READS_DIR}/ngs" \
            --sample ngs_30x \
            -t "$THREADS" \
            --meancov 30 \
            --skip-hifi --skip-ont \
            --seed "$SEED"
    else
        echo "  NGS reads already exist, skipping."
    fi

    # ONT reads
    if [[ ! -d "${READS_DIR}/ont" ]]; then
        log_step "  Simulating ONT reads (30x)"
        ecc sim-reads \
            -r "$REFERENCE" \
            -i "$ECCDNA_FASTA" \
            -o "${READS_DIR}/ont" \
            --sample ont_30x \
            -t "$THREADS" \
            --meancov 30 \
            --skip-sr --skip-hifi \
            --seed "$SEED"
    else
        echo "  ONT reads already exist, skipping."
    fi

    # ------------------------------------------------------------------
    # Step 3: Align reads to reference
    # ------------------------------------------------------------------
    log_step "Step 3: Align reads to reference"
    mkdir -p "$ALIGN_DIR"

    # Find the simulated FASTQ files and align them
    # Adjust glob patterns based on actual sim-reads output structure
    for platform in hifi ngs ont; do
        local bam="${ALIGN_DIR}/${platform}_30x.sorted.bam"
        if [[ -f "$bam" ]]; then
            echo "  $platform alignment already exists, skipping."
            continue
        fi

        # Find FASTQ files from sim-reads output
        local reads_found=0
        local fq_files=""

        # sim-reads output structure: <output>/<sample>/*.fastq or *.fq.gz
        local sample_dir="${READS_DIR}/${platform}/${platform}_30x"
        if [[ -d "$sample_dir" ]]; then
            fq_files=$(find "$sample_dir" -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" 2>/dev/null | head -1)
        fi

        if [[ -n "$fq_files" ]]; then
            align_reads "$platform" "$fq_files" "$bam"
            reads_found=1
        fi

        if [[ "$reads_found" -eq 0 ]]; then
            echo "  [WARN] No FASTQ found for $platform in ${sample_dir}"
        fi
    done

    # ------------------------------------------------------------------
    # Step 4: Run detection tools
    # ------------------------------------------------------------------
    log_step "Step 4: Run detection tools"
    mkdir -p "$DETECT_DIR"

    # HiFi-based tools
    local hifi_bam="${ALIGN_DIR}/hifi_30x.sorted.bam"
    if [[ -f "$hifi_bam" ]]; then
        run_circleseeker "$hifi_bam" "${DETECT_DIR}/circleseeker_hifi" "hifi"
        run_cresil "$hifi_bam" "${DETECT_DIR}/cresil_hifi"
        run_ecc_finder "$hifi_bam" "${DETECT_DIR}/eccfinder_hifi" "hifi"
    fi

    # NGS-based tools
    local ngs_bam="${ALIGN_DIR}/ngs_30x.sorted.bam"
    if [[ -f "$ngs_bam" ]]; then
        run_circleseeker "$ngs_bam" "${DETECT_DIR}/circleseeker_ngs" "ngs"
        run_circle_map "$ngs_bam" "${DETECT_DIR}/circlemap_ngs"
        run_ecc_finder "$ngs_bam" "${DETECT_DIR}/eccfinder_ngs" "ngs"
    fi

    # ONT-based tools
    local ont_bam="${ALIGN_DIR}/ont_30x.sorted.bam"
    if [[ -f "$ont_bam" ]]; then
        run_circleseeker "$ont_bam" "${DETECT_DIR}/circleseeker_ont" "ont"
    fi

    # ------------------------------------------------------------------
    # Step 5: Benchmark evaluation
    # ------------------------------------------------------------------
    log_step "Step 5: Benchmark evaluation"
    mkdir -p "$BENCH_DIR"

    # Benchmark each tool's results against ground truth
    # CircleSeeker outputs: merged_output.csv
    # CReSIL outputs: adjust based on actual output
    # Circle-Map outputs: circle_map_output.bed
    # ecc_finder outputs: adjust based on actual output

    # HiFi benchmarks
    run_benchmark "CircleSeeker_HiFi" \
        "${DETECT_DIR}/circleseeker_hifi/merged_output.csv" \
        "$GT_DIR" "$BENCH_DIR" "cs_hifi"

    run_benchmark "CReSIL_HiFi" \
        "${DETECT_DIR}/cresil_hifi/result.txt" \
        "$GT_DIR" "$BENCH_DIR" "cresil_hifi"

    run_benchmark "ecc_finder_HiFi" \
        "${DETECT_DIR}/eccfinder_hifi/result.txt" \
        "$GT_DIR" "$BENCH_DIR" "eccfinder_hifi"

    # NGS benchmarks
    run_benchmark "CircleSeeker_NGS" \
        "${DETECT_DIR}/circleseeker_ngs/merged_output.csv" \
        "$GT_DIR" "$BENCH_DIR" "cs_ngs"

    run_benchmark "Circle-Map_NGS" \
        "${DETECT_DIR}/circlemap_ngs/circle_map_output.bed" \
        "$GT_DIR" "$BENCH_DIR" "circlemap_ngs"

    run_benchmark "ecc_finder_NGS" \
        "${DETECT_DIR}/eccfinder_ngs/result.txt" \
        "$GT_DIR" "$BENCH_DIR" "eccfinder_ngs"

    # ONT benchmarks
    run_benchmark "CircleSeeker_ONT" \
        "${DETECT_DIR}/circleseeker_ont/merged_output.csv" \
        "$GT_DIR" "$BENCH_DIR" "cs_ont"

    echo ""
    echo "Experiment 1 completed. Results in: $BENCH_DIR"
}


# ############################################################################
#
# EXPERIMENT 2: Coverage Gradient (Supplementary)
#
# Same ground truth as Exp 1, vary coverage: 5x, 10x, 30x, 50x, 100x
# Platform: HiFi
# Tools: CircleSeeker, CReSIL
#
# ############################################################################

run_experiment_coverage() {
    local EXP_DIR="${OUTPUT_DIR}/exp2_coverage"
    local GT_DIR="${EXP_DIR}/ground_truth"
    local BENCH_DIR="${EXP_DIR}/benchmark"

    log_section "Experiment 2: Coverage Gradient"

    # ------------------------------------------------------------------
    # Step 1: Generate ground truth (same params as standard)
    # ------------------------------------------------------------------
    log_step "Step 1: Generate ground truth"

    if [[ -d "$GT_DIR" && -f "$GT_DIR/covgrad.all.fa" ]]; then
        echo "  Ground truth already exists, skipping."
    else
        ecc sim-region \
            -r "$REFERENCE" \
            -o covgrad \
            -d "$GT_DIR" \
            -u 10000 -m 1000 -c 500 \
            -t "$THREADS" \
            --seed "$SEED"
    fi

    local ECCDNA_FASTA="${GT_DIR}/covgrad.all.fa"

    # ------------------------------------------------------------------
    # Step 2-5: For each coverage level
    # ------------------------------------------------------------------
    for COV in 5 10 30 50 100; do
        local COV_STR=$(fmt_cov "$COV")
        local COV_DIR="${EXP_DIR}/cov_${COV_STR}x"
        local READS_DIR="${COV_DIR}/reads"
        local ALIGN_DIR="${COV_DIR}/alignments"
        local DETECT_DIR="${COV_DIR}/detection"

        log_step "Coverage ${COV_STR}x"

        # Simulate HiFi reads
        if [[ ! -d "${READS_DIR}/hifi" ]]; then
            ecc sim-reads \
                -r "$REFERENCE" \
                -i "$ECCDNA_FASTA" \
                -o "${READS_DIR}/hifi" \
                --sample "hifi_${COV_STR}x" \
                -t "$THREADS" \
                --meancov "$COV" \
                --skip-sr --skip-ont \
                --seed "$SEED"
        fi

        # Align reads
        local bam="${ALIGN_DIR}/hifi_${COV_STR}x.sorted.bam"
        if [[ ! -f "$bam" ]]; then
            mkdir -p "$ALIGN_DIR"
            local sample_dir="${READS_DIR}/hifi/hifi_${COV_STR}x"
            local fq_files=$(find "$sample_dir" -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" 2>/dev/null | head -1)
            if [[ -n "$fq_files" ]]; then
                align_reads "hifi" "$fq_files" "$bam"
            else
                echo "  [WARN] No FASTQ found for HiFi ${COV_STR}x"
                continue
            fi
        fi

        # Run tools
        run_circleseeker "$bam" "${DETECT_DIR}/circleseeker" "hifi"
        run_cresil "$bam" "${DETECT_DIR}/cresil"

        # Benchmark
        mkdir -p "$BENCH_DIR"
        run_benchmark "CircleSeeker_${COV_STR}x" \
            "${DETECT_DIR}/circleseeker/merged_output.csv" \
            "$GT_DIR" "$BENCH_DIR" "cs_${COV_STR}x"

        run_benchmark "CReSIL_${COV_STR}x" \
            "${DETECT_DIR}/cresil/result.txt" \
            "$GT_DIR" "$BENCH_DIR" "cresil_${COV_STR}x"
    done

    echo ""
    echo "Experiment 2 completed. Results in: $BENCH_DIR"
}


# ############################################################################
#
# EXPERIMENT 3: Complexity Gradient (Supplementary)
#
# Vary eccDNA count: 1K, 5K, 10K, 50K (keeping U:M:C ratio ~ 80:15:5)
# Coverage: 30x, Platform: HiFi
# Tool: CircleSeeker
#
# ############################################################################

run_experiment_complexity() {
    local EXP_DIR="${OUTPUT_DIR}/exp3_complexity"
    local BENCH_DIR="${EXP_DIR}/benchmark"

    log_section "Experiment 3: Complexity Gradient"

    # Define complexity levels: total (U M C)
    local -a COMPLEXITIES=(
        "1000:800:150:50"
        "5000:4000:750:250"
        "10000:8000:1500:500"
        "50000:40000:7500:2500"
    )

    for entry in "${COMPLEXITIES[@]}"; do
        IFS=':' read -r TOTAL U M C <<< "$entry"
        local LABEL="${TOTAL}"
        local COMP_DIR="${EXP_DIR}/n_${LABEL}"
        local GT_DIR="${COMP_DIR}/ground_truth"
        local READS_DIR="${COMP_DIR}/reads"
        local ALIGN_DIR="${COMP_DIR}/alignments"
        local DETECT_DIR="${COMP_DIR}/detection"

        log_step "Complexity N=${TOTAL} (U=${U}, M=${M}, C=${C})"

        # Generate ground truth
        if [[ ! -d "$GT_DIR" || ! -f "$GT_DIR/comp_${LABEL}.all.fa" ]]; then
            ecc sim-region \
                -r "$REFERENCE" \
                -o "comp_${LABEL}" \
                -d "$GT_DIR" \
                -u "$U" -m "$M" -c "$C" \
                -t "$THREADS" \
                --seed "$SEED"
        fi

        local ECCDNA_FASTA="${GT_DIR}/comp_${LABEL}.all.fa"

        # Simulate HiFi reads at 30x
        if [[ ! -d "${READS_DIR}/hifi" ]]; then
            ecc sim-reads \
                -r "$REFERENCE" \
                -i "$ECCDNA_FASTA" \
                -o "${READS_DIR}/hifi" \
                --sample "hifi_30x" \
                -t "$THREADS" \
                --meancov 30 \
                --skip-sr --skip-ont \
                --seed "$SEED"
        fi

        # Align reads
        local bam="${ALIGN_DIR}/hifi_30x.sorted.bam"
        if [[ ! -f "$bam" ]]; then
            mkdir -p "$ALIGN_DIR"
            local sample_dir="${READS_DIR}/hifi/hifi_30x"
            local fq_files=$(find "$sample_dir" -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" 2>/dev/null | head -1)
            if [[ -n "$fq_files" ]]; then
                align_reads "hifi" "$fq_files" "$bam"
            else
                echo "  [WARN] No FASTQ found for complexity N=${TOTAL}"
                continue
            fi
        fi

        # Run CircleSeeker
        run_circleseeker "$bam" "${DETECT_DIR}/circleseeker" "hifi"

        # Benchmark
        mkdir -p "$BENCH_DIR"
        run_benchmark "CircleSeeker_N${TOTAL}" \
            "${DETECT_DIR}/circleseeker/merged_output.csv" \
            "$GT_DIR" "$BENCH_DIR" "cs_n${LABEL}"
    done

    echo ""
    echo "Experiment 3 completed. Results in: $BENCH_DIR"
}


# ############################################################################
#
# Main execution
#
# ############################################################################

log_section "eccToolkit Simulation Benchmark Experiment"
echo "Reference: $REFERENCE"
echo "Output:    $OUTPUT_DIR"
echo "Threads:   $THREADS"
echo "Seed:      $SEED"

mkdir -p "$OUTPUT_DIR"

# Record experiment metadata
cat > "${OUTPUT_DIR}/experiment_info.txt" <<EOF
eccToolkit Simulation Benchmark
================================
Date:      $(date -Iseconds)
Reference: $REFERENCE
Output:    $OUTPUT_DIR
Threads:   $THREADS
Seed:      $SEED

Tool availability:
$(for tool in "${!TOOL_AVAILABLE[@]}"; do
    if [[ "${TOOL_AVAILABLE[$tool]}" -eq 1 ]]; then
        echo "  [OK]   $tool"
    else
        echo "  [SKIP] $tool"
    fi
done)
EOF

# Run all experiments
run_experiment_standard
run_experiment_coverage
run_experiment_complexity

log_section "All experiments completed!"
echo "Results are in: $OUTPUT_DIR"
echo ""
echo "To generate figures, run:"
echo "  python scripts/benchmark_plot.py $OUTPUT_DIR"
