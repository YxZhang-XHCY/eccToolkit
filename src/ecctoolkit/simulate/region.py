"""
Simulate eccDNA regions from reference genome.

Generates three types of eccDNA:
- UeccDNA (Unique): Only one genomic location matches at >=99% identity
- MeccDNA (Multi-mapped): Multiple genomic locations match at >=99% identity
- CeccDNA (Chimeric): Composed of 2-5 fragments from different locations

Uses minimap2 for fast alignment and classification.
"""

import logging
from typing import Optional

logger = logging.getLogger(__name__)


def run_region_simulation(
    reference: str,
    output_prefix: str,
    num_unique: int = 1000,
    num_multi: int = 0,
    num_chimeric: int = 0,
    threads: int = 8,
    seed: int = 42,
    # Length distribution parameters
    mode: float = 400.0,
    sigma: float = 0.8,
    tail_weight: float = 0.05,
    tail_min: int = 5000,
    tail_max: int = 500000,
    min_length: int = 100,
    max_length: int = 500000,
    # Alignment parameters
    identity: float = 99.0,
    min_coverage: float = 90.0,
    max_secondary: int = 50,
    no_hit_policy: str = "skip",
    split_by_length: bool = True,
    split_length: int = 5000,
    minimap_preset_short: str = "sr",
    minimap_preset_long: str = "asm5",
    candidate_multiplier_u: float = 2.0,
    candidate_multiplier_m: float = 10.0,
    # Output options
    output_dir: Optional[str] = None,
    keep_tmp: bool = False,
    verbose: bool = False,
) -> None:
    """
    Simulate eccDNA regions from reference genome.

    Args:
        reference: Reference genome FASTA file (.fa/.fasta/.fna or .gz)
        output_prefix: Output prefix (also used as directory name)
        num_unique: Number of UeccDNA (Unique) to generate
        num_multi: Number of MeccDNA (Multi-mapped) to generate
        num_chimeric: Number of CeccDNA (Chimeric) to generate
        threads: Number of threads for minimap2
        seed: Random seed for reproducibility
        mode: Lognormal distribution mode/peak in bp
        sigma: Lognormal distribution sigma
        tail_weight: Weight of tail distribution (long eccDNA)
        tail_min: Minimum length for tail distribution
        tail_max: Maximum length for tail distribution
        min_length: Global minimum eccDNA length
        max_length: Global maximum eccDNA length
        identity: Minimum identity threshold for classification (%)
        min_coverage: Minimum coverage threshold (%)
        max_secondary: Maximum secondary alignments to consider
        no_hit_policy: How to handle no-hit regions ("skip" or "unique")
        split_by_length: Split candidates for short/long minimap2 presets
        split_length: Length threshold for short/long presets
        minimap_preset_short: Minimap2 preset for short sequences
        minimap_preset_long: Minimap2 preset for long sequences
        candidate_multiplier_u: Candidate multiplier for Unique regions
        candidate_multiplier_m: Candidate multiplier for Multi regions
        output_dir: Custom output directory (default: create from prefix)
        keep_tmp: Keep temporary minimap2 files
        verbose: Enable verbose logging

    Outputs:
        - {prefix}.all.bed: All eccDNA regions in BED format
        - {prefix}.all.fa: All eccDNA sequences in FASTA format
        - {prefix}.unique.bed/fa: UeccDNA only
        - {prefix}.multi.bed/fa: MeccDNA only
        - {prefix}.chimeric.bed/fa: CeccDNA only
        - {prefix}.qc.log: Quality control report
        - {prefix}.lengths.tsv: Length distribution data
        - {prefix}.multi_hits.tsv: Multi-mapping locations for MeccDNA
    """
    logger.info("eccDNA region simulation")
    logger.info(f"Reference: {reference}")
    logger.info(f"Output prefix: {output_prefix}")
    logger.info(f"Target: {num_unique} Unique, {num_multi} Multi, {num_chimeric} Chimeric")
    logger.info(f"Threads: {threads}, Seed: {seed}")

    # Import and run the simulator
    try:
        from .eccDNA_simulator import SimulationConfig, EccDNASimulator

        # Create configuration
        config = SimulationConfig(
            reference=reference,
            output_prefix=output_prefix,
            output_dir=output_dir,
            num_unique=num_unique,
            num_multi=num_multi,
            num_chimeric=num_chimeric,
            seed=seed,
            mode=mode,
            sigma=sigma,
            tail_weight=tail_weight,
            tail_min=tail_min,
            tail_max=tail_max,
            min_length=min_length,
            max_length=max_length,
            threads=threads,
            identity_threshold=identity,
            min_coverage=min_coverage,
            max_secondary=max_secondary,
            no_hit_policy=no_hit_policy,
            split_by_length=split_by_length,
            split_length=split_length,
            minimap_preset_short=minimap_preset_short,
            minimap_preset_long=minimap_preset_long,
            candidate_multiplier_u=candidate_multiplier_u,
            candidate_multiplier_m=candidate_multiplier_m,
            keep_tmp=keep_tmp,
        )

        # Run simulation
        simulator = EccDNASimulator(config)
        simulator.run()

    except ImportError as e:
        logger.error(f"Failed to import simulator: {e}")
        logger.error("Please ensure eccToolkit is installed and modules are available")
        raise NotImplementedError(
            "Region simulation requires ecctoolkit.simulate.eccDNA_simulator. "
            "Please reinstall eccToolkit and retry."
        )
