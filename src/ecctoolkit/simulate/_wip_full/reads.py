"""
Simulate sequencing reads from eccDNA sequences.

Implements RCA (Rolling Circle Amplification) simulation to generate
realistic NGS, HiFi, and ONT reads from eccDNA templates.

Pipeline:
1. Compartment generation - group eccDNA instances
2. RCA molecule generation - create branched/hyperbranched structures
3. Chimera injection - simulate inter-molecule chimeras
4. Debranching/linearization - convert to linear molecules
5. Read generation - use external tools (ART/PBSIM) via fqsim
"""

import logging
from typing import Optional, Tuple

logger = logging.getLogger(__name__)


def run_read_simulation(
    input_file: str,
    output_dir: str,
    # Output scale (coverage mode only)
    cov_ngs: float = 30.0,
    cov_hifi: float = 30.0,
    cov_ont: float = 30.0,
    force_ecc_coverage: bool = True,
    force_ecc_min_reads_factor: float = 2.0,
    # Adaptive coverage sampling
    adaptive_sampling: bool = True,
    min_coverage: float = 3.0,
    # RCA parameters
    mu_r: float = 3.0,
    sigma_r: float = 0.5,
    k_len: float = 0.5,
    rca_mode: Optional[str] = None,
    reaction_hours: Optional[float] = None,
    # Branch parameters
    b_rate: float = 0.5,
    branch_length_mode: str = "ratio_trunk",
    # Chimera parameters
    chimera_breakpoint_mode: str = "uniform",
    chimera_breakpoint_alpha: float = 0.5,
    chimera_breakpoint_beta: float = 0.5,
    chimera_background_ratio: float = 0.0,
    chimera_background_pool: int = 200,
    # Debranch efficiency
    d_eff: float = 0.9,
    # NGS parameters
    i_mu: float = 400,
    # Background DNA parameters
    reference: Optional[str] = None,
    linear_ratio: float = 0.0,
    linear_mode: str = "additive",
    linear_min_len: int = 200,
    linear_max_len: int = 10000,
    # Parallel processing
    threads: int = 1,
    # Options
    seed: Optional[int] = None,
    config_file: Optional[str] = None,
    skip_ngs: bool = False,
    skip_hifi: bool = False,
    skip_ont: bool = False,
    compress: bool = False,
    truth_format: str = "tsv",
    sample: str = "",
    verbose: bool = False,
    # External tools parameters (for fqsim)
    ont_model: str = "R94",
    ont_mean: float = 3000,
    ont_std: float = 2500,
    hifi_profile_id: Optional[str] = None,
    sr_platform: str = "HS25",
    sr_readlen: int = 150,
    sr_insert_mean: float = 400,
    sr_insert_std: float = 125,
) -> None:
    """
    Simulate sequencing reads from eccDNA FASTA using RCA model.

    Uses external tools (ART/PBSIM) for read generation via fqsim.

    Args:
        input_file: Input eccDNA FASTA file
        output_dir: Output directory
        cov_ngs: NGS target coverage depth
        cov_hifi: HiFi target coverage depth
        cov_ont: ONT target coverage depth
        force_ecc_coverage: Force each eccDNA to appear in reads when output is large
        force_ecc_min_reads_factor: Threshold factor for forcing eccDNA coverage
        mu_r: ln(R) mean at reference length
        sigma_r: ln(R) standard deviation
        k_len: Length correction exponent
        rca_mode: "empirical" or "kinetic"
        reaction_hours: Reaction time in hours (kinetic mode)
        b_rate: Branch rate (/kb)
        branch_length_mode: Branch length mode ("ratio_trunk" or "ratio_ecc")
        chimera_breakpoint_mode: Chimera breakpoint sampling mode ("uniform" or "beta")
        chimera_breakpoint_alpha: Beta distribution alpha for chimera breakpoints
        chimera_breakpoint_beta: Beta distribution beta for chimera breakpoints
        chimera_background_ratio: Fraction of chimera donors from background pool
        chimera_background_pool: Background donor pool size for chimera
        d_eff: Debranching efficiency (0-1)
        i_mu: NGS insert size mean
        reference: Reference genome FASTA for background linear DNA
        linear_ratio: Ratio of background linear DNA reads (0-1)
        linear_mode: How to apply linear_ratio ("fraction" or "additive")
        linear_min_len: Minimum length of background DNA fragments
        linear_max_len: Maximum length of background DNA fragments
        threads: Number of parallel threads (0=auto, 1=single)
        seed: Random seed
        config_file: YAML/JSON config file (overrides other params)
        skip_ngs: Skip NGS read generation
        skip_hifi: Skip HiFi read generation
        skip_ont: Skip ONT read generation
        compress: Compress output files (gzip)
        truth_format: Truth output format ("tsv" or "jsonl")
        sample: Prefix for output files
        verbose: Enable verbose logging
        ont_model: ONT error model for PBSIM (R94, R103, etc.)
        ont_mean: ONT read length mean
        ont_std: ONT read length std
        hifi_profile_id: HiFi profile ID for PBSIM3
        sr_platform: NGS platform for ART (HS25, etc.)
        sr_readlen: NGS read length
        sr_insert_mean: NGS insert size mean
        sr_insert_std: NGS insert size std

    Outputs:
        - {sample}.NGS.R1.fastq, {sample}.NGS.R2.fastq: Paired-end NGS reads
        - {sample}.HiFi.fastq: PacBio HiFi reads
        - {sample}.ONT.fastq: Oxford Nanopore reads
        - {sample}_pool.lib.csv: Molecule pool with eccDNA source tracking
    """
    logger.info("RCA-based read simulation (using external tools)")
    logger.info(f"Input: {input_file}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Target coverage - NGS: {cov_ngs}x, HiFi: {cov_hifi}x, ONT: {cov_ont}x")

    try:
        from pathlib import Path
        import numpy as np
        import os
        import tempfile
        import shutil
        import atexit

        from .rca_readsim.config import SimConfig, get_default_config
        from .rca_readsim.models import (
            EccDNA,
            SegmentType,
            LinearMolecule,
            register_ecc_length,
        )
        from .rca_readsim.io_utils import parse_fasta, build_ecc_db, summarize_eccdnas
        from .rca_readsim.compartment import CompartmentGenerator, estimate_total_instances
        from .rca_readsim.rca_engine import RCAEngine
        from .rca_readsim.chimera import ChimeraInjector
        from .rca_readsim.debranch import Debrancher, LinearMoleculePool
        from .rca_readsim.parallel import (
            get_optimal_workers,
            parallel_generate_rca,
            parallel_inject_chimera,
            parallel_debranch,
        )

        # Load config
        if config_file:
            config_path = Path(config_file)
            if config_path.suffix in ['.yaml', '.yml']:
                config = SimConfig.from_yaml(config_file)
            else:
                config = SimConfig.from_json(config_file)
        else:
            config = get_default_config()

        # Override with CLI parameters
        if seed is not None:
            config.seed = seed
        config.rca.mu_R = mu_r
        config.rca.sigma_R = sigma_r
        config.rca.k_len = k_len
        if rca_mode:
            config.rca.mode = rca_mode
        if reaction_hours is not None:
            config.kinetics.reaction_time_min = reaction_hours * 60.0
        config.branch.B_rate = b_rate
        config.branch.branch_length_mode = branch_length_mode
        config.debranch.D_eff = d_eff
        config.ngs_library.I_mu = i_mu
        config.output_scale.Cov_NGS = cov_ngs
        config.output_scale.Cov_HiFi = cov_hifi
        config.output_scale.Cov_ONT = cov_ont
        config.output_scale.adaptive_sampling = adaptive_sampling
        config.output_scale.min_coverage = min_coverage
        config.constants.force_ecc_coverage = force_ecc_coverage
        config.constants.force_ecc_min_reads_factor = force_ecc_min_reads_factor
        config.constants.chimera_breakpoint_mode = chimera_breakpoint_mode
        config.constants.chimera_breakpoint_alpha = chimera_breakpoint_alpha
        config.constants.chimera_breakpoint_beta = chimera_breakpoint_beta
        config.constants.chimera_background_ratio = chimera_background_ratio
        config.constants.chimera_background_pool = chimera_background_pool

        if config.rca.mode == "kinetic":
            logger.info(
                "RCA mode: kinetic "
                f"(time={config.kinetics.reaction_time_min} min, "
                f"v={config.kinetics.v_nt_per_min} nt/min, "
                f"len_model={config.kinetics.len_scale_mode}, "
                f"penalty={config.kinetics.length_penalty_mode}, "
                f"resource={config.kinetics.resource_mode})"
            )
        else:
            logger.info(
                "RCA mode: empirical "
                f"(mu_R={config.rca.mu_R}, sigma_R={config.rca.sigma_R}, "
                f"k_len={config.rca.k_len})"
            )

        # Background DNA parameters
        if reference is not None:
            config.background.fasta_path = reference
        if reference is not None and linear_ratio > 0:
            config.background.enabled = True
            config.background.ratio = linear_ratio
            config.background.min_length = linear_min_len
            config.background.max_length = linear_max_len
            logger.info(f"Background DNA enabled: ratio={linear_ratio}, ref={reference}")
        elif linear_ratio > 0 and reference is None:
            logger.warning("linear_ratio set but no reference provided; background DNA disabled")

        # Setup
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        prefix = sample or ""
        if prefix and not prefix.endswith(("_", "-", ".")):
            prefix = f"{prefix}_"

        # 设置临时目录在输出目录下（避免挤爆系统临时目录）
        temp_dir = output_path / ".tmp"
        temp_dir.mkdir(parents=True, exist_ok=True)
        os.environ['TMPDIR'] = str(temp_dir)
        os.environ['TEMP'] = str(temp_dir)
        os.environ['TMP'] = str(temp_dir)
        tempfile.tempdir = str(temp_dir)
        logger.info(f"Temp directory: {temp_dir}")

        def cleanup_temp():
            if temp_dir.exists():
                try:
                    shutil.rmtree(temp_dir)
                    logger.info(f"Cleaned up temp directory: {temp_dir}")
                except Exception as e:
                    logger.warning(f"Failed to clean temp directory: {e}")

        atexit.register(cleanup_temp)

        rng = np.random.default_rng(config.seed) if config.seed else np.random.default_rng()

        # Load eccDNA
        eccdnas = parse_fasta(Path(input_file))
        if not eccdnas:
            raise ValueError("No eccDNAs found in input file")

        ecc_db = build_ecc_db(eccdnas)
        for ecc in eccdnas:
            register_ecc_length(ecc.id, ecc.length)

        logger.info(summarize_eccdnas(eccdnas))
        total_ecc_size = sum(ecc.length for ecc in eccdnas)

        # Compute target reads
        def compute_target(platform: str) -> int:
            """计算目标 reads 数量 (仅支持 coverage 模式)"""
            scale = config.output_scale
            if platform == 'NGS':
                value = scale.Cov_NGS
                read_len = config.constants.read_length
            elif platform == 'HiFi':
                value = scale.Cov_HiFi
                read_len = config.constants.hifi_target_size_mean
            else:
                value = scale.Cov_ONT
                read_len = 10000

            # 仅支持 coverage 模式: value 表示覆盖度倍数
            return max(100, int(value * total_ecc_size / read_len))

        requested_ngs = 0 if skip_ngs else compute_target("NGS")
        requested_hifi = 0 if skip_hifi else compute_target("HiFi")
        requested_ont = 0 if skip_ont else compute_target("ONT")

        def _evenize(n: int) -> int:
            return (n // 2) * 2

        def _split_target_reads(platform: str, requested_reads: int) -> Tuple[int, int, int]:
            if requested_reads <= 0:
                return 0, 0, 0

            ratio = config.background.ratio if config.background.enabled else 0.0
            mode = (linear_mode or "additive").lower()
            if mode not in {"fraction", "additive"}:
                raise ValueError(f"Invalid linear_mode: {linear_mode!r}")

            if ratio <= 0:
                if platform == "NGS":
                    requested_reads = _evenize(requested_reads)
                return requested_reads, 0, requested_reads

            if not (0 < ratio < 1):
                raise ValueError(f"linear_ratio must be in (0,1); got {ratio}")

            if platform == "NGS":
                requested_pairs = _evenize(requested_reads) // 2
                if requested_pairs <= 0:
                    return 0, 0, 0

                if mode == "fraction":
                    total_pairs = requested_pairs
                    bg_pairs = int(round(total_pairs * ratio))
                    bg_pairs = max(0, min(bg_pairs, total_pairs))
                    ecc_pairs = total_pairs - bg_pairs
                else:
                    ecc_pairs = requested_pairs
                    bg_pairs = int(round(ecc_pairs * ratio / (1 - ratio)))
                    bg_pairs = max(0, bg_pairs)

                total_pairs = ecc_pairs + bg_pairs
                return ecc_pairs * 2, bg_pairs * 2, total_pairs * 2

            if mode == "fraction":
                total_reads = int(requested_reads)
                bg_reads = int(round(total_reads * ratio))
                bg_reads = max(0, min(bg_reads, total_reads))
                ecc_reads = total_reads - bg_reads
                return ecc_reads, bg_reads, total_reads

            ecc_reads = int(requested_reads)
            bg_reads = int(round(ecc_reads * ratio / (1 - ratio)))
            bg_reads = max(0, bg_reads)
            return ecc_reads, bg_reads, ecc_reads + bg_reads

        target_ngs_ecc, target_ngs_bg, target_ngs_total = _split_target_reads("NGS", requested_ngs)
        target_hifi_ecc, target_hifi_bg, target_hifi_total = _split_target_reads("HiFi", requested_hifi)
        target_ont_ecc, target_ont_bg, target_ont_total = _split_target_reads("ONT", requested_ont)

        logger.info(
            "Target reads (total/ecc/background) - "
            f"NGS: {target_ngs_total}/{target_ngs_ecc}/{target_ngs_bg}, "
            f"HiFi: {target_hifi_total}/{target_hifi_ecc}/{target_hifi_bg}, "
            f"ONT: {target_ont_total}/{target_ont_ecc}/{target_ont_bg}"
        )

        num_workers = get_optimal_workers(threads)
        use_parallel = num_workers > 1

        def should_parallel(item_count: int, min_per_worker: int) -> bool:
            return use_parallel and item_count >= num_workers * min_per_worker

        # Estimate instances
        max_ecc_target = max(target_ngs_ecc, target_hifi_ecc, target_ont_ecc)
        estimate_platform = "NGS"
        if max_ecc_target == target_hifi_ecc:
            estimate_platform = "HiFi"
        elif max_ecc_target == target_ont_ecc:
            estimate_platform = "ONT"
        target_instances = estimate_total_instances(
            eccdnas, config, max_ecc_target, platform=estimate_platform
        )

        # Optional background donor pool for chimera
        background_generator = None
        chimera_background_donors = []
        if config.constants.chimera_background_ratio > 0 and config.constants.chimera_background_pool > 0:
            if not config.background.fasta_path:
                logger.warning("chimera_background_ratio set but no reference provided; skipping background donors")
            else:
                from .rca_readsim.background import BackgroundDNAGenerator

                background_generator = BackgroundDNAGenerator(
                    fasta_path=config.background.fasta_path,
                    params=config.background,
                    seed=config.seed
                )
                for i in range(config.constants.chimera_background_pool):
                    molecule, sequence, bg_id = background_generator.generate_fragment(i)
                    if bg_id not in ecc_db:
                        ecc_db[bg_id] = EccDNA(id=bg_id, seq=sequence, weight=0.0)
                        register_ecc_length(bg_id, len(sequence))
                    chimera_background_donors.append((bg_id, len(sequence), len(sequence)))
                logger.info(f"Chimera background donors: {len(chimera_background_donors)}")

        # Pipeline
        logger.info("Step 1: Generating compartments...")
        comp_gen = CompartmentGenerator(config, rng)
        compartments, _ = comp_gen.generate(eccdnas, target_instances)

        logger.info("Step 2: Generating RCA molecules...")
        if should_parallel(len(compartments), 4):
            logger.info(f"Parallel RCA generation: {num_workers} workers")
            molecules, _ = parallel_generate_rca(
                compartments,
                ecc_db,
                config,
                num_workers=num_workers,
                seed=config.seed,
            )
        else:
            rca_engine = RCAEngine(config, rng)
            molecules, _ = rca_engine.generate_all(compartments, ecc_db)

        logger.info("Step 3: Injecting chimeras...")
        if should_parallel(len(compartments), 4):
            logger.info(f"Parallel chimera injection: {num_workers} workers")
            molecules, _ = parallel_inject_chimera(
                molecules,
                config,
                num_workers=num_workers,
                seed=config.seed + 500 if config.seed else None,
                background_donors=chimera_background_donors,
            )
        else:
            chimera_inj = ChimeraInjector(config, rng, background_donors=chimera_background_donors)
            molecules, _ = chimera_inj.process_all(molecules)

        # ==========================================================================
        # Step 4: 为三个平台生成独立的分子池
        # - NGS: 打断 (~400bp)，不需要 debranch
        # - HiFi: 打断 (15-25kb)，不需要 debranch
        # - ONT: 完全 debranch，不打断
        # ==========================================================================

        from .rca_readsim.fragmentation import Fragmenter
        from ..readsim import fqsim

        # 先将 RCA 分子转换为 LinearMolecule（无 debranch）
        # 这里直接从 RCA 分子提取主干作为线性分子
        raw_linear_molecules = []
        for mol in molecules:
            # 主干分子
            trunk_seg = mol.get_trunk_segment()
            trunk_mol = LinearMolecule(
                molecule_id=f"{mol.instance_id}_trunk",
                segments=[trunk_seg],
                source_graph_id=mol.instance_id,
                is_from_branch=False,
                repeat_count=mol.repeat_count,
                source_ecc_length=mol.source_ecc.length,
                has_chimera=len(mol.chimera_junctions) > 0,
            )
            raw_linear_molecules.append(trunk_mol)

        logger.info(f"  Raw linear molecules from RCA: {len(raw_linear_molecules)}")

        # Background DNA generation (optional)
        background_molecules = []
        background_db = {}
        if config.background.enabled:
            logger.info("Step 4.5: Generating background linear DNA...")
            from .rca_readsim.background import BackgroundDNAGenerator

            max_bg_target = max(target_ngs_bg, target_hifi_bg, target_ont_bg)
            bg_count = max_bg_target
            logger.info(f"  Target background molecules: {bg_count}")
            if bg_count > 0:
                if background_generator is None:
                    background_generator = BackgroundDNAGenerator(
                        fasta_path=config.background.fasta_path,
                        params=config.background,
                        seed=config.seed
                    )
                background_molecules, background_db = background_generator.generate_molecules(
                    bg_count, background_db
                )
                logger.info(f"  Generated {len(background_molecules)} background molecules")

        # 合并 ecc_db 和 background_db
        merged_db = dict(ecc_db)
        merged_db.update(background_db)

        # 创建打断器
        fragmenter = Fragmenter(config, rng)

        # ==========================================================================
        # NGS 分子池：打断
        # ==========================================================================
        ngs_pool_csv = None
        if not skip_ngs:
            logger.info("Step 4a: Preparing NGS molecule pool (fragmentation)...")
            ngs_molecules, ngs_frag_stats = fragmenter.fragment_for_ngs(
                raw_linear_molecules + background_molecules,
                merged_db,
                insert_mean=sr_insert_mean,
                insert_std=sr_insert_std,
            )
            logger.info(f"  {ngs_frag_stats.summary()}")

            if ngs_molecules:
                ngs_pool = LinearMoleculePool(ngs_molecules, merged_db, weight_by_length=True)
                ngs_pool_csv = output_path / f"{prefix}ngs_pool.lib.csv"
                ngs_pool.to_csv(
                    str(ngs_pool_csv),
                    include_coverage=True,
                    include_truth=True,
                    target_coverage=cov_ngs,
                    read_length=config.constants.read_length,
                    original_ecc_total_length=total_ecc_size
                )
                logger.info(f"  NGS pool exported: {ngs_pool_csv} ({len(ngs_molecules)} fragments)")

        # ==========================================================================
        # HiFi 分子池：打断
        # ==========================================================================
        hifi_pool_csv = None
        if not skip_hifi:
            logger.info("Step 4b: Preparing HiFi molecule pool (fragmentation)...")
            hifi_molecules, hifi_frag_stats = fragmenter.fragment_for_hifi(
                raw_linear_molecules + background_molecules,
                merged_db,
                target_mean=config.constants.hifi_target_size_mean,
                target_std=config.constants.hifi_target_size_std,
                min_length=config.constants.hifi_min_length,
                max_length=config.constants.hifi_max_length,
            )
            logger.info(f"  {hifi_frag_stats.summary()}")

            if hifi_molecules:
                hifi_pool = LinearMoleculePool(hifi_molecules, merged_db, weight_by_length=True)
                hifi_pool_csv = output_path / f"{prefix}hifi_pool.lib.csv"
                hifi_pool.to_csv(
                    str(hifi_pool_csv),
                    include_coverage=True,
                    include_truth=True,
                    target_coverage=cov_hifi,
                    read_length=config.constants.hifi_target_size_mean,
                    original_ecc_total_length=total_ecc_size
                )
                logger.info(f"  HiFi pool exported: {hifi_pool_csv} ({len(hifi_molecules)} fragments)")

        # ==========================================================================
        # ONT 分子池：完全 debranch
        # ==========================================================================
        ont_pool_csv = None
        if not skip_ont:
            logger.info("Step 4c: Preparing ONT molecule pool (debranching)...")
            if should_parallel(len(molecules), 50):
                logger.info(f"  Parallel debranching: {num_workers} workers")
                ont_linear_molecules, _ = parallel_debranch(
                    molecules,
                    ecc_db,
                    config,
                    num_workers=num_workers,
                    seed=config.seed + 1000 if config.seed else None,
                )
            else:
                debrancher = Debrancher(config, rng)
                ont_linear_molecules, _ = debrancher.process_all(molecules, ecc_db)

            # 添加背景分子
            all_ont_molecules = list(ont_linear_molecules) + background_molecules
            logger.info(f"  ONT molecules: {len(ont_linear_molecules)} eccDNA + {len(background_molecules)} background")

            if all_ont_molecules:
                ont_pool = LinearMoleculePool(all_ont_molecules, merged_db, weight_by_length=True)
                ont_pool_csv = output_path / f"{prefix}ont_pool.lib.csv"
                ont_pool.to_csv(
                    str(ont_pool_csv),
                    include_coverage=True,
                    include_truth=True,
                    target_coverage=cov_ont,
                    read_length=ont_mean,  # ONT 使用平均 read 长度
                    original_ecc_total_length=total_ecc_size
                )
                logger.info(f"  ONT pool exported: {ont_pool_csv} ({len(all_ont_molecules)} molecules)")

        # ==========================================================================
        # Step 5: 为每个平台调用 fqsim
        # ==========================================================================
        logger.info("Step 5: Generating reads using external tools (ART/PBSIM2)...")

        # NGS reads (ART)
        if not skip_ngs and ngs_pool_csv:
            logger.info("  Generating NGS reads (ART)...")
            fqsim(
                sample=sample or "sample",
                csv=str(ngs_pool_csv),
                path=str(output_path),
                seed=seed,
                thread=threads,
                skip_sr=False,
                skip_hifi=True,
                skip_ont=True,
                sr_mean=sr_insert_mean,
                sr_std=sr_insert_std,
                sr_readlen=sr_readlen,
                sr_platform=sr_platform,
                generate_truth=True,  # 生成 NGS truth 文件
            )
            logger.info("  NGS reads completed.")

        # HiFi reads (simple mode - no external PBSIM dependency)
        if not skip_hifi and hifi_pool_csv:
            logger.info("  Generating HiFi reads (simple mode)...")
            fqsim(
                sample=sample or "sample",
                csv=str(hifi_pool_csv),
                path=str(output_path),
                seed=seed + 100 if seed else None,
                thread=threads,
                skip_sr=True,
                skip_hifi=False,
                skip_ont=True,
                # HiFi 专用参数
                hifi_mode='simple',  # 使用内置 simple 模式，无需外部工具
                hifi_len_min=config.constants.hifi_min_length,
                hifi_len_peak_min=config.constants.hifi_target_size_mean - config.constants.hifi_target_size_std,
                hifi_len_peak_max=config.constants.hifi_target_size_mean + config.constants.hifi_target_size_std,
                hifi_len_max=config.constants.hifi_max_length,
                hifi_qmean=30,  # HiFi 高质量
                hifi_qsd=2.0,
                generate_truth=True,
            )
            logger.info("  HiFi reads completed.")

        # ONT reads (PBSIM2)
        if not skip_ont and ont_pool_csv:
            logger.info("  Generating ONT reads (PBSIM2)...")
            fqsim(
                sample=sample or "sample",
                csv=str(ont_pool_csv),
                path=str(output_path),
                seed=seed + 200 if seed else None,
                thread=threads,
                skip_sr=True,
                skip_hifi=True,
                skip_ont=False,
                ont_mean=ont_mean,
                ont_std=ont_std,
                ont_model=ont_model,
                generate_truth=True,  # 生成 ONT truth 文件
            )
            logger.info("  ONT reads completed.")

        logger.info("  All platform reads generation completed.")

        # Summary
        logger.info("=" * 50)
        logger.info("Full mode simulation completed!")
        logger.info(f"Output directory: {output_path}")
        if ngs_pool_csv:
            logger.info(f"NGS pool: {ngs_pool_csv}")
        if hifi_pool_csv:
            logger.info(f"HiFi pool: {hifi_pool_csv}")
        if ont_pool_csv:
            logger.info(f"ONT pool: {ont_pool_csv}")
        logger.info("=" * 50)

    except Exception as e:
        logger.error(f"Simulation failed: {e}")
        raise
