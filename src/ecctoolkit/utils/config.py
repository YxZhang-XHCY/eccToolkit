"""Configuration constants and environment setup for eccToolkit."""

import os

# Default parameters
DEFAULT_CORES = 8
DEFAULT_PERMUTATIONS = 1000
DEFAULT_RANDOM_SEED = 42

# Human genome sizes (hg38/GRCh38)
HG38_GENOME_SIZES = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
}

# Human genome sizes (hg19/GRCh37)
HG19_GENOME_SIZES = {
    "chr1": 249250621,
    "chr2": 243199373,
    "chr3": 198022430,
    "chr4": 191154276,
    "chr5": 180915260,
    "chr6": 171115067,
    "chr7": 159138663,
    "chr8": 146364022,
    "chr9": 141213431,
    "chr10": 135534747,
    "chr11": 135006516,
    "chr12": 133851895,
    "chr13": 115169878,
    "chr14": 107349540,
    "chr15": 102531392,
    "chr16": 90354753,
    "chr17": 81195210,
    "chr18": 78077248,
    "chr19": 59128983,
    "chr20": 63025520,
    "chr21": 48129895,
    "chr22": 51304566,
    "chrX": 155270560,
    "chrY": 59373566,
}

# Mapping for assembly names
GENOME_SIZES = {
    "hg38": HG38_GENOME_SIZES,
    "GRCh38": HG38_GENOME_SIZES,
    "hg19": HG19_GENOME_SIZES,
    "GRCh37": HG19_GENOME_SIZES,
}

# Standard chromosome lists
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
SEX_CHROMOSOMES = ["chrX", "chrY"]
MITOCHONDRIAL = ["chrM", "chrMT"]


def setup_thread_limits(n_threads: int = 1) -> None:
    """
    Set environment variables to prevent thread oversubscription.

    This should be called before importing numpy/scipy to take effect.

    Args:
        n_threads: Number of threads to allow (default: 1)
    """
    thread_vars = [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]
    for var in thread_vars:
        os.environ[var] = str(n_threads)


def get_genome_sizes(assembly: str = "hg38") -> dict:
    """
    Get chromosome sizes for a given assembly.

    Args:
        assembly: Genome assembly name (hg38, hg19, GRCh38, GRCh37)

    Returns:
        Dictionary mapping chromosome names to sizes

    Raises:
        ValueError: If assembly is not supported
    """
    if assembly not in GENOME_SIZES:
        raise ValueError(
            f"Unsupported assembly: {assembly}. "
            f"Supported assemblies: {list(GENOME_SIZES.keys())}"
        )
    return GENOME_SIZES[assembly].copy()


def get_chromosomes(
    include_sex: bool = False,
    include_mito: bool = False,
) -> list:
    """
    Get list of chromosomes to analyze.

    Args:
        include_sex: Include sex chromosomes (chrX, chrY)
        include_mito: Include mitochondrial chromosome

    Returns:
        List of chromosome names
    """
    chroms = AUTOSOMES.copy()
    if include_sex:
        chroms.extend(SEX_CHROMOSOMES)
    if include_mito:
        chroms.extend(MITOCHONDRIAL)
    return chroms
