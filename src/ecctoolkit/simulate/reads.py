"""
Simulate sequencing reads from eccDNA sequences.

Implements RCA (Rolling Circle Amplification) simulation to generate
realistic NGS, HiFi, and ONT reads from eccDNA templates.

Pipeline:
1. Compartment generation - group eccDNA instances
2. RCA molecule generation - create branched/hyperbranched structures
3. Chimera injection - simulate inter-molecule chimeras
4. Debranching/linearization - convert to linear molecules
5. Library preparation - generate platform-specific reads
"""

import logging
from typing import Optional

logger = logging.getLogger(__name__)


def run_read_simulation(
    input_file: str,
    output_dir: str,
    # Output scale
    cov_ngs: float = 10000,
    cov_hifi: float = 1000,
    cov_ont: float = 1000,
    output_mode: str = "read_count",
    force_ecc_coverage: bool = True,
    force_ecc_min_reads_factor: float = 2.0,
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
) -> None:
    """
    Simulate sequencing reads from eccDNA FASTA using RCA model.

    Args:
        input_file: Input eccDNA FASTA file
        output_dir: Output directory
        cov_ngs: NGS output scale (read count or coverage)
        cov_hifi: HiFi output scale (read count or coverage)
        cov_ont: ONT output scale (read count or coverage)
        output_mode: "read_count" or "coverage"
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

    Outputs:
        - reads_ngs_R1.fastq, reads_ngs_R2.fastq: Paired-end NGS reads
        - reads_hifi.fastq: PacBio HiFi reads
        - reads_ont.fastq: Oxford Nanopore reads
        - truth.tsv/jsonl: Ground truth for each read (includes source_type: eccDNA/background)
        - config_used.yaml: Configuration used for simulation
    """
    logger.info("RCA-based read simulation")
    logger.info(f"Input: {input_file}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Mode: {output_mode}")
    if output_mode == "read_count":
        logger.info(f"Target reads - NGS: {cov_ngs}, HiFi: {cov_hifi}, ONT: {cov_ont}")
    else:
        logger.info(f"Target coverage - NGS: {cov_ngs}x, HiFi: {cov_hifi}x, ONT: {cov_ont}x")

    # Import and run RCA read simulator
    try:
        from pathlib import Path
        import numpy as np
        import os
        import tempfile
        import shutil
        import atexit
        import gzip

        from .rca_readsim.config import SimConfig, get_default_config
        from .rca_readsim.models import (
            EccDNA,
            SegmentType,
            register_ecc_length,
        )
        from .rca_readsim.io_utils import parse_fasta, build_ecc_db, summarize_eccdnas
        from .rca_readsim.compartment import CompartmentGenerator, estimate_total_instances
        from .rca_readsim.rca_engine import RCAEngine
        from .rca_readsim.chimera import ChimeraInjector
        from .rca_readsim.debranch import Debrancher, LinearMoleculePool
        from .rca_readsim.library import NGSLibraryGenerator, HiFiLibraryGenerator, ONTLibraryGenerator
        from .rca_readsim.truth import TruthWriter, BedWriter
        from .rca_readsim.parallel import (
            get_optimal_workers,
            serialize_pool,
            ParallelReadGenerator,
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
        config.output_scale.mode = output_mode
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

        # Compute target reads
        def compute_target(platform: str) -> int:
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

            if scale.mode == "coverage":
                total_size = sum(ecc.length for ecc in eccdnas)
                return max(100, int(value * total_size / read_len))
            return int(value)

        target_ngs = 0 if skip_ngs else compute_target('NGS')
        target_hifi = 0 if skip_hifi else compute_target('HiFi')
        target_ont = 0 if skip_ont else compute_target('ONT')

        num_workers = get_optimal_workers(threads)
        use_parallel = num_workers > 1

        def should_parallel(item_count: int, min_per_worker: int) -> bool:
            return use_parallel and item_count >= num_workers * min_per_worker

        def _should_force_reads(platform: str, target_reads: int) -> bool:
            if not config.constants.force_ecc_coverage:
                return False
            if not eccdnas:
                return False
            per_ecc = 2 if platform == "NGS" else 1
            required = len(eccdnas) * per_ecc
            threshold = max(required, required * config.constants.force_ecc_min_reads_factor)
            return target_reads >= threshold

        force_platforms = {
            "NGS": (not skip_ngs) and _should_force_reads("NGS", target_ngs),
            "HiFi": (not skip_hifi) and _should_force_reads("HiFi", target_hifi),
            "ONT": (not skip_ont) and _should_force_reads("ONT", target_ont),
        }

        def _force_min_length(platforms: dict) -> int:
            lengths = []
            if platforms.get("NGS"):
                lengths.append(config.constants.read_length * 2)
            if platforms.get("HiFi"):
                lengths.append(config.constants.hifi_min_length)
            if platforms.get("ONT"):
                lengths.append(config.constants.ont_min_length)
            return max(lengths) if lengths else 0

        force_min_length = _force_min_length(force_platforms)

        # Estimate instances
        max_target = max(target_ngs, target_hifi, target_ont)
        target_instances = estimate_total_instances(eccdnas, config, max_target)

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

        logger.info("Step 4: Debranching...")
        if should_parallel(len(molecules), 50):
            logger.info(f"Parallel debranching: {num_workers} workers")
            linear_molecules, _ = parallel_debranch(
                molecules,
                ecc_db,
                config,
                num_workers=num_workers,
                seed=config.seed + 1000 if config.seed else None,
            )
        else:
            debrancher = Debrancher(config, rng)
            linear_molecules, _ = debrancher.process_all(molecules, ecc_db)

        # Step 4.5: Background DNA generation (optional)
        background_molecules = []
        if config.background.enabled:
            logger.info("Step 4.5: Generating background linear DNA...")
            from .rca_readsim.background import BackgroundDNAGenerator, calculate_background_count

            total_target = max(target_ngs, target_hifi, target_ont)
            bg_count = calculate_background_count(total_target, config.background.ratio)
            logger.info(f"  Target background molecules: {bg_count}")

            if background_generator is None:
                background_generator = BackgroundDNAGenerator(
                    fasta_path=config.background.fasta_path,
                    params=config.background,
                    seed=config.seed
                )
            background_molecules, ecc_db = background_generator.generate_molecules(bg_count, ecc_db)
            logger.info(f"  Generated {len(background_molecules)} background molecules")

        def _has_ecc_candidate(molecules, ecc_id: str, min_len: int) -> bool:
            for mol in molecules:
                if mol.is_background:
                    continue
                if mol.total_length < min_len:
                    continue
                for seg in mol.segments:
                    if seg.segment_type != SegmentType.BACKGROUND and seg.ecc_id == ecc_id:
                        return True
            return False

        def _augment_forced_molecules(molecules, min_len: int):
            if min_len <= 0:
                return []
            missing = [ecc.id for ecc in eccdnas if not _has_ecc_candidate(molecules, ecc.id, min_len)]
            if not missing:
                return []

            import math

            logger.warning(
                f"Force coverage: {len(missing)} eccDNA lack length >= {min_len}, "
                "generating extra molecules"
            )

            extra = []
            rca_engine = RCAEngine(config, rng)
            chimera_inj = ChimeraInjector(config, rng, background_donors=chimera_background_donors)
            debrancher = Debrancher(config, rng)
            force_compartment = -1
            force_counter = 0

            for ecc_id in missing:
                ecc = ecc_db.get(ecc_id)
                if ecc is None or ecc.length <= 0:
                    logger.warning(f"Force coverage: eccDNA {ecc_id} has invalid length")
                    continue
                min_repeat = max(1, int(math.ceil(min_len / ecc.length)))
                if min_repeat > config.R_max:
                    logger.warning(
                        f"Force coverage: eccDNA {ecc_id} requires repeat {min_repeat} "
                        f"> R_max={config.R_max}"
                    )
                    min_repeat = config.R_max

                found = False
                for attempt in range(5):
                    base_repeat = rca_engine.sample_repeat_count(ecc.length)
                    repeat = max(min_repeat, base_repeat)
                    if attempt > 0:
                        repeat = min(config.R_max, repeat * (attempt + 1))
                    mol = rca_engine.generate_molecule(ecc, force_compartment, repeat)
                    force_counter += 1
                    mol.instance_id = f"force_{ecc_id}_{force_counter}"
                    mol.compartment_id = force_compartment
                    force_compartment -= 1

                    mols, _ = chimera_inj.process_all([mol])
                    linear_extra, _ = debrancher.process_all(mols, ecc_db)

                    if _has_ecc_candidate(linear_extra, ecc_id, min_len):
                        extra.extend(linear_extra)
                        found = True
                        break

                if not found:
                    logger.warning(
                        f"Force coverage: failed to build molecule for eccDNA {ecc_id}"
                    )

            return extra

        if force_min_length > 0:
            forced_extra = _augment_forced_molecules(linear_molecules, force_min_length)
            if forced_extra:
                linear_molecules.extend(forced_extra)
                logger.info(
                    f"Force coverage: added {len(forced_extra)} extra eccDNA molecules"
                )

        # Create molecule pool (eccDNA + background)
        all_linear_molecules = linear_molecules + background_molecules
        pool = LinearMoleculePool(all_linear_molecules, ecc_db, weight_by_length=True)
        logger.info(f"  Molecule pool: {len(linear_molecules)} eccDNA + {len(background_molecules)} background")

        # Setup parallel processing for read generation
        pool_data = None
        read_pool = None
        if use_parallel:
            logger.info(f"Parallel mode enabled: {num_workers} workers")
            pool_data = serialize_pool(pool)
            if max(target_ngs, target_hifi, target_ont) > 0:
                read_pool = ParallelReadGenerator(pool_data, num_workers)

        def _open_fastq(path: Path):
            opener = gzip.open if compress else open
            mode = 'wt' if compress else 'w'
            return opener(path, mode)

        def _write_fastq_reads(handle, reads):
            for read in reads:
                handle.write(read.to_fastq())

        def _write_paired_reads(handle_r1, handle_r2, reads):
            for read in reads:
                if read.read_number == 1:
                    handle_r1.write(read.to_fastq())
                elif read.read_number == 2:
                    handle_r2.write(read.to_fastq())

        ecc_ids = [ecc.id for ecc in eccdnas]

        # Build forced-coverage candidates from post-RCA molecules (no synthetic reads).
        def _build_force_candidates(molecules, target_ids):
            target_set = set(target_ids)
            candidates = {ecc_id: [] for ecc_id in target_ids}
            for mol in molecules:
                if mol.is_background:
                    continue
                pos = 0
                spans_by_ecc = {}
                for seg in mol.segments:
                    seg_end = pos + seg.length
                    if seg.segment_type != SegmentType.BACKGROUND and seg.ecc_id in target_set:
                        spans_by_ecc.setdefault(seg.ecc_id, []).append((pos, seg_end))
                    pos = seg_end
                for ecc_id, spans in spans_by_ecc.items():
                    candidates[ecc_id].append((mol, spans))
            return candidates

        ecc_force_candidates = _build_force_candidates(linear_molecules, ecc_ids)

        def _pick_overlap_start(seg_start: int, seg_end: int, span_length: int, mol_len: int):
            if span_length <= 0 or span_length > mol_len:
                return None
            min_start = max(0, seg_start - span_length + 1)
            max_start = min(mol_len - span_length, seg_end - 1)
            if min_start > max_start:
                return None
            return int(rng.integers(min_start, max_start + 1))

        def _pick_ngs_insert_start(
            seg_start: int,
            seg_end: int,
            mol_len: int,
            insert_size: int,
            read_length: int,
        ):
            if insert_size > mol_len:
                return None
            # Try to overlap R1
            min_start = max(0, seg_start - read_length + 1)
            max_start = min(mol_len - insert_size, seg_end - 1)
            if min_start <= max_start:
                return int(rng.integers(min_start, max_start + 1))
            # Try to overlap R2
            min_start = max(0, seg_start - insert_size + 1)
            max_start = min(mol_len - insert_size, seg_end - insert_size + read_length - 1)
            if min_start <= max_start:
                return int(rng.integers(min_start, max_start + 1))
            return None

        def _force_ngs_for_ecc(ecc_id: str, generator) -> list:
            candidates = ecc_force_candidates.get(ecc_id, [])
            if not candidates:
                return []
            read_length = config.constants.read_length
            min_insert = 2 * read_length
            candidate_order = rng.permutation(len(candidates))
            for idx in candidate_order:
                mol, spans = candidates[idx]
                mol_seq = pool.get_sequence(mol)
                mol_len = len(mol_seq)
                if mol_len < min_insert:
                    continue
                insert_sizes = [generator.sample_insert_size() for _ in range(5)]
                insert_sizes.append(mol_len)
                for insert_size in insert_sizes:
                    if insert_size < min_insert or insert_size > mol_len:
                        continue
                    span_order = rng.permutation(len(spans))
                    for sidx in span_order:
                        seg_start, seg_end = spans[sidx]
                        insert_start = _pick_ngs_insert_start(
                            seg_start, seg_end, mol_len, insert_size, read_length
                        )
                        if insert_start is None:
                            continue
                        branch_result = generator.check_branch_chimera(
                            mol, insert_start, insert_size
                        )
                        branch_chimera = None
                        if branch_result is not None:
                            branch_info, breakpoint = branch_result
                            branch_chimera = generator.generate_branch_chimera_insert(
                                mol, mol_seq, insert_start, insert_size,
                                branch_info, breakpoint
                            )
                        r1, r2 = generator.generate_paired_reads(
                            mol, mol_seq, insert_start, insert_size, branch_chimera
                        )
                        if ecc_id not in r1.source_ecc_ids and ecc_id not in r2.source_ecc_ids:
                            continue
                        return [r1, r2]
            return []

        def _force_hifi_for_ecc(ecc_id: str, generator) -> list:
            candidates = ecc_force_candidates.get(ecc_id, [])
            if not candidates:
                return []
            min_len = config.constants.hifi_min_length
            max_len = config.constants.hifi_max_length
            target_mean = config.constants.hifi_target_size_mean
            target_std = config.constants.hifi_target_size_std
            candidate_order = rng.permutation(len(candidates))
            for idx in candidate_order:
                mol, spans = candidates[idx]
                mol_seq = pool.get_sequence(mol)
                mol_len = len(mol_seq)
                if mol_len < min_len:
                    continue
                if mol_len <= max_len:
                    frag_length = mol_len
                else:
                    frag_length = int(rng.normal(target_mean, target_std))
                    frag_length = max(min_len, min(max_len, frag_length, mol_len))
                span_order = rng.permutation(len(spans))
                for sidx in span_order:
                    seg_start, seg_end = spans[sidx]
                    frag_start = _pick_overlap_start(seg_start, seg_end, frag_length, mol_len)
                    if frag_start is None:
                        continue
                    frag_mol = generator._create_fragment_molecule(mol, frag_start, frag_length)
                    if frag_mol is None:
                        continue
                    frag_seq = mol_seq[frag_start:frag_start + frag_length]
                    read = generator.generate_read(frag_mol, frag_seq, 0, frag_length)
                    if read is None:
                        continue
                    if ecc_id not in read.source_ecc_ids:
                        continue
                    return [read]
            return []

        def _force_ont_for_ecc(ecc_id: str, generator) -> list:
            candidates = ecc_force_candidates.get(ecc_id, [])
            if not candidates:
                return []
            min_len = config.constants.ont_min_length
            candidate_order = rng.permutation(len(candidates))
            for idx in candidate_order:
                mol, spans = candidates[idx]
                mol_seq = pool.get_sequence(mol)
                mol_len = len(mol_seq)
                if mol_len < min_len:
                    continue
                span_order = rng.permutation(len(spans))
                for sidx in span_order:
                    seg_start, _ = spans[sidx]
                    read_start = min(seg_start, mol_len - min_len)
                    if read_start < 0:
                        continue
                    read_length = mol_len - read_start
                    read = generator.generate_read(
                        mol, mol_seq, read_start=read_start, read_length=read_length
                    )
                    if read is None:
                        continue
                    if ecc_id not in read.source_ecc_ids:
                        continue
                    return [read]
            return []

        def _generate_forced_reads(platform: str, generator, target_reads: int) -> list:
            forced_reads = []
            for ecc_id in ecc_ids:
                if platform == "NGS":
                    if len(forced_reads) + 2 > target_reads:
                        break
                    reads = _force_ngs_for_ecc(ecc_id, generator)
                elif platform == "HiFi":
                    if len(forced_reads) + 1 > target_reads:
                        break
                    reads = _force_hifi_for_ecc(ecc_id, generator)
                else:
                    if len(forced_reads) + 1 > target_reads:
                        break
                    reads = _force_ont_for_ecc(ecc_id, generator)
                if not reads:
                    logger.warning(f"{platform}: failed to force eccDNA {ecc_id} into reads")
                    continue
                forced_reads.extend(reads)
            return forced_reads

        def _generate_reads(
            platform: str,
            generator_class,
            target_reads: int,
            seed_base: Optional[int],
            write_fn,
        ):
            if target_reads <= 0:
                return
            if platform == "NGS":
                target_reads = (target_reads // 2) * 2
                if target_reads == 0:
                    return
            chunk_size = min(50000, max(1000, target_reads // 10))
            generated = 0
            generator = None
            chunk_idx = 0

            if force_platforms.get(platform, False):
                if generator is None:
                    generator = generator_class(config, rng)
                    if platform == "NGS":
                        generator._ecc_db = pool.ecc_db

                forced_reads = _generate_forced_reads(platform, generator, target_reads)

                if forced_reads:
                    write_fn(forced_reads)
                    truth_writer.write_reads(forced_reads)
                    bed_writer.write_reads(forced_reads)
                    generated += len(forced_reads)
                    if generated >= target_reads:
                        return

            while generated < target_reads:
                remaining = target_reads - generated
                chunk_target = min(chunk_size, remaining)
                if platform == "NGS":
                    chunk_target = (chunk_target // 2) * 2
                    if chunk_target == 0:
                        break

                chunk_seed = None
                if seed_base is not None:
                    chunk_seed = seed_base + chunk_idx * 10000

                if use_parallel and chunk_target >= num_workers * 100 and read_pool is not None:
                    reads, _ = read_pool.generate(
                        generator_class,
                        config,
                        chunk_target,
                        seed=chunk_seed,
                        job_id=chunk_idx,
                    )
                else:
                    if generator is None:
                        generator = generator_class(config, rng)
                        if platform == "NGS":
                            generator._ecc_db = pool.ecc_db
                    reads, _ = generator.generate(pool, chunk_target)

                if not reads:
                    logger.warning(f"{platform}: no reads generated for chunk {chunk_idx}")
                    break

                write_fn(reads)
                truth_writer.write_reads(reads)
                bed_writer.write_reads(reads)
                generated += len(reads)
                chunk_idx += 1

        logger.info("Step 5: Generating reads (streaming outputs)...")
        truth_ext = ".jsonl" if truth_format == "jsonl" else ".tsv"
        truth_path = output_path / f"{prefix}truth{truth_ext}"
        try:
            with TruthWriter(truth_path, format=truth_format) as truth_writer, \
                    BedWriter(output_path, prefix=prefix) as bed_writer:
                if not skip_ngs:
                    logger.info("Step 5a: Generating NGS reads...")
                    r1_name = f"{prefix}reads_ngs_R1.fastq"
                    r2_name = f"{prefix}reads_ngs_R2.fastq"
                    if compress:
                        r1_name += ".gz"
                        r2_name += ".gz"
                    r1_path = output_path / r1_name
                    r2_path = output_path / r2_name
                    with _open_fastq(r1_path) as r1_handle, _open_fastq(r2_path) as r2_handle:
                        _generate_reads(
                            "NGS",
                            NGSLibraryGenerator,
                            target_ngs,
                            config.seed,
                            lambda reads: _write_paired_reads(r1_handle, r2_handle, reads),
                        )

                if not skip_hifi:
                    logger.info("Step 5b: Generating HiFi reads...")
                    hifi_name = f"{prefix}reads_hifi.fastq"
                    if compress:
                        hifi_name += ".gz"
                    hifi_path = output_path / hifi_name
                    with _open_fastq(hifi_path) as hifi_handle:
                        _generate_reads(
                            "HiFi",
                            HiFiLibraryGenerator,
                            target_hifi,
                            config.seed + 1000 if config.seed else None,
                            lambda reads: _write_fastq_reads(hifi_handle, reads),
                        )

                if not skip_ont:
                    logger.info("Step 5c: Generating ONT reads...")
                    ont_name = f"{prefix}reads_ont.fastq"
                    if compress:
                        ont_name += ".gz"
                    ont_path = output_path / ont_name
                    with _open_fastq(ont_path) as ont_handle:
                        _generate_reads(
                            "ONT",
                            ONTLibraryGenerator,
                            target_ont,
                            config.seed + 2000 if config.seed else None,
                            lambda reads: _write_fastq_reads(ont_handle, reads),
                        )
        finally:
            if read_pool is not None:
                read_pool.close()

        config.to_yaml(output_path / f"{prefix}config_used.yaml")
        logger.info("Simulation completed successfully!")

    except ImportError as e:
        logger.error(f"Failed to import RCA read simulator: {e}")
        raise NotImplementedError(f"Read simulation import error: {e}")
