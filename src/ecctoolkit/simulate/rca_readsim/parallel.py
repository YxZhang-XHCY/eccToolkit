"""
并行处理模块

提供多进程支持以加速读段生成。
"""

import multiprocessing as mp
from typing import List, Tuple, Dict, Any, Optional
from dataclasses import dataclass
import numpy as np
import logging

logger = logging.getLogger(__name__)

_READ_POOL = None


@dataclass
class ParallelConfig:
    """并行配置"""
    num_workers: int = 1
    chunk_size: int = 1000  # 每个worker处理的read数


def get_optimal_workers(requested: int = 0) -> int:
    """
    获取最优worker数量

    Args:
        requested: 请求的worker数，0表示自动

    Returns:
        实际使用的worker数
    """
    cpu_count = mp.cpu_count()

    if requested <= 0:
        # 自动：使用 CPU 核数 - 1，至少1个
        return max(1, cpu_count - 1)
    else:
        # 限制在 [1, CPU数] 范围内
        return min(max(1, requested), cpu_count)


def parallel_generate_reads(
    generator_class,
    config,
    pool_data: Dict[str, Any],
    target_reads: int,
    num_workers: int = 1,
    seed: Optional[int] = None
) -> Tuple[List, Any]:
    """
    并行生成reads

    Args:
        generator_class: 生成器类 (NGSLibraryGenerator 等)
        config: SimConfig
        pool_data: 序列化的分子池数据
        target_reads: 目标read数
        num_workers: worker数量
        seed: 随机种子基数

    Returns:
        (reads列表, 统计信息)
    """
    if num_workers <= 1:
        # 单进程模式
        return _generate_reads_single(generator_class, config, pool_data, target_reads, seed)

    # 多进程模式
    logger.info(f"Using {num_workers} workers for parallel read generation")

    # 分配任务
    reads_per_worker = target_reads // num_workers
    remainder = target_reads % num_workers

    tasks = []
    for i in range(num_workers):
        worker_target = reads_per_worker + (1 if i < remainder else 0)
        worker_seed = seed + i if seed is not None else None
        tasks.append((generator_class, config, pool_data, worker_target, worker_seed, i))

    # 并行执行
    with mp.Pool(num_workers) as pool:
        results = pool.starmap(_worker_generate_reads, tasks)

    # 合并结果
    all_reads = []
    for reads, stats in results:
        all_reads.extend(reads)

    # 合并统计信息 (简化处理，只返回第一个worker的stats类型)
    combined_stats = results[0][1] if results else None

    return all_reads, combined_stats


class ParallelReadGenerator:
    """Persistent worker pool for read generation."""

    def __init__(self, pool_data: Dict[str, Any], num_workers: int):
        self.num_workers = max(1, num_workers)
        self._pool = mp.Pool(
            self.num_workers,
            initializer=_init_read_worker,
            initargs=(pool_data,)
        )

    def generate(
        self,
        generator_class,
        config,
        target_reads: int,
        seed: Optional[int] = None,
        job_id: Optional[int] = None,
    ):
        reads_per_worker = target_reads // self.num_workers
        remainder = target_reads % self.num_workers

        tasks = []
        for i in range(self.num_workers):
            worker_target = reads_per_worker + (1 if i < remainder else 0)
            worker_seed = seed + i if seed is not None else None
            tasks.append((generator_class, config, worker_target, worker_seed, i, job_id))

        results = self._pool.starmap(_worker_generate_reads_global, tasks)

        all_reads = []
        for reads, _ in results:
            all_reads.extend(reads)

        combined_stats = results[0][1] if results else None
        return all_reads, combined_stats

    def close(self):
        if self._pool is not None:
            self._pool.close()
            self._pool.join()
            self._pool = None

    def terminate(self):
        if self._pool is not None:
            self._pool.terminate()
            self._pool.join()
            self._pool = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def _init_read_worker(pool_data: Dict[str, Any]):
    global _READ_POOL
    _READ_POOL = _rebuild_pool(pool_data)


def _worker_generate_reads_global(generator_class, config, target_reads, seed, worker_id, job_id):
    import numpy as np

    if _READ_POOL is None:
        raise RuntimeError("Read pool not initialized in worker process")

    rng = np.random.default_rng(seed)
    generator = generator_class(config, rng)
    generator._worker_id = worker_id

    reads, stats = generator.generate(_READ_POOL, target_reads)

    prefix = f"w{worker_id}_"
    if job_id is not None:
        prefix = f"j{job_id}_{prefix}"

    for read in reads:
        read.read_id = f"{prefix}{read.read_id}"
        if hasattr(read, 'mate_id') and read.mate_id:
            read.mate_id = f"{prefix}{read.mate_id}"

    return reads, stats


def _split_round_robin(items: List[Any], num_workers: int) -> List[List[Any]]:
    buckets = [[] for _ in range(num_workers)]
    for idx, item in enumerate(items):
        buckets[idx % num_workers].append(item)
    return buckets


def _assign_groups_by_load(
    groups: List[Tuple[int, List[Any]]],
    num_workers: int
) -> List[List[Tuple[int, List[Any]]]]:
    buckets: List[List[Tuple[int, List[Any]]]] = [[] for _ in range(num_workers)]
    loads = [0 for _ in range(num_workers)]
    for group in groups:
        group_size = len(group[1])
        idx = loads.index(min(loads))
        buckets[idx].append(group)
        loads[idx] += group_size
    return buckets


def _merge_rca_stats(stats_list: List[Any]):
    from .rca_engine import RCAStats

    merged = RCAStats()
    for stats in stats_list:
        if stats is None:
            continue
        merged.total_molecules += stats.total_molecules
        merged.total_branches += stats.total_branches
        for key, value in stats.repeat_count_distribution.items():
            merged.repeat_count_distribution[key] = (
                merged.repeat_count_distribution.get(key, 0) + value
            )
        for key, value in stats.branch_count_distribution.items():
            merged.branch_count_distribution[key] = (
                merged.branch_count_distribution.get(key, 0) + value
            )
    return merged


def parallel_generate_rca(
    compartments: List[Any],
    ecc_db: Dict[str, Any],
    config,
    num_workers: int = 1,
    seed: Optional[int] = None
) -> Tuple[List[Any], Any]:
    """
    并行生成RCA分子图
    """
    from .rca_engine import RCAEngine, RCAStats

    if num_workers <= 1:
        rng = np.random.default_rng(seed)
        engine = RCAEngine(config, rng)
        return engine.generate_all(compartments, ecc_db)

    if not compartments:
        return [], RCAStats()

    chunks = _split_round_robin(compartments, num_workers)
    tasks = []
    for idx, chunk in enumerate(chunks):
        if chunk:
            tasks.append((chunk, ecc_db, config, seed, idx))

    if not tasks:
        return [], RCAStats()

    worker_count = min(num_workers, len(tasks))
    with mp.Pool(worker_count) as pool:
        results = pool.starmap(_worker_generate_rca, tasks)

    all_molecules: List[Any] = []
    stats_list = []
    for molecules, stats in results:
        all_molecules.extend(molecules)
        stats_list.append(stats)

    return all_molecules, _merge_rca_stats(stats_list)


def _worker_generate_rca(compartments, ecc_db, config, seed, worker_id):
    from .rca_engine import RCAEngine

    worker_seed = seed + worker_id if seed is not None else None
    rng = np.random.default_rng(worker_seed)
    engine = RCAEngine(config, rng)
    return engine.generate_all(compartments, ecc_db)


def _merge_chimera_stats(stats_list: List[Any]):
    from .chimera import ChimeraStats

    merged = ChimeraStats()
    for stats in stats_list:
        if stats is None:
            continue
        merged.total_chimera_events += stats.total_chimera_events
        merged.molecules_with_chimera += stats.molecules_with_chimera
        for key, value in stats.chimera_by_compartment_size.items():
            merged.chimera_by_compartment_size[key] = (
                merged.chimera_by_compartment_size.get(key, 0) + value
            )
    return merged


def parallel_inject_chimera(
    molecules: List[Any],
    config,
    num_workers: int = 1,
    seed: Optional[int] = None,
    background_donors: Optional[List[Any]] = None
) -> Tuple[List[Any], Any]:
    """
    并行注入chimera事件（按compartment划分）
    """
    from .chimera import ChimeraInjector, ChimeraStats

    if num_workers <= 1:
        rng = np.random.default_rng(seed)
        injector = ChimeraInjector(config, rng, background_donors=background_donors)
        return injector.process_all(molecules)

    if not molecules:
        return [], ChimeraStats()

    by_compartment: Dict[int, List[Any]] = {}
    comp_order: List[int] = []
    for mol in molecules:
        comp_id = mol.compartment_id
        if comp_id not in by_compartment:
            by_compartment[comp_id] = []
            comp_order.append(comp_id)
        by_compartment[comp_id].append(mol)

    groups = [(comp_id, by_compartment[comp_id]) for comp_id in comp_order]
    buckets = _assign_groups_by_load(groups, num_workers)

    tasks = []
    for idx, bucket in enumerate(buckets):
        if bucket:
            tasks.append((bucket, config, seed, idx, background_donors))

    worker_count = min(num_workers, len(tasks))
    with mp.Pool(worker_count) as pool:
        results = pool.starmap(_worker_inject_chimera, tasks)

    processed: Dict[int, List[Any]] = {}
    stats_list = []
    for group_list, stats in results:
        stats_list.append(stats)
        for comp_id, comp_mols in group_list:
            processed[comp_id] = comp_mols

    merged_stats = _merge_chimera_stats(stats_list)
    ordered = []
    for comp_id in comp_order:
        ordered.extend(processed.get(comp_id, []))
    return ordered, merged_stats


def _worker_inject_chimera(compartment_groups, config, seed, worker_id, background_donors):
    from .chimera import ChimeraInjector, ChimeraStats

    worker_seed = seed + worker_id if seed is not None else None
    rng = np.random.default_rng(worker_seed)
    injector = ChimeraInjector(config, rng, background_donors=background_donors)
    stats = ChimeraStats()

    for _, comp_mols in compartment_groups:
        comp_size = len(comp_mols)
        events = injector.process_compartment_molecules(comp_mols, comp_size)
        stats.total_chimera_events += events
        if events > 0:
            stats.chimera_by_compartment_size[comp_size] = (
                stats.chimera_by_compartment_size.get(comp_size, 0) + events
            )

    stats.molecules_with_chimera = sum(
        1 for _, comp_mols in compartment_groups for mol in comp_mols if mol.chimera_junctions
    )
    return compartment_groups, stats


def _merge_debranch_stats(stats_list: List[Any]):
    from .debranch import DebranchStats

    merged = DebranchStats()
    for stats in stats_list:
        if stats is None:
            continue
        merged.total_branches += stats.total_branches
        merged.debranched_count += stats.debranched_count
        merged.retained_count += stats.retained_count
        merged.nick_count += stats.nick_count
        merged.break_count += stats.break_count
        merged.linear_molecules_from_trunk += stats.linear_molecules_from_trunk
        merged.linear_molecules_from_branch += stats.linear_molecules_from_branch
        merged.trunk_fragments_from_breaks += stats.trunk_fragments_from_breaks
    return merged


def parallel_debranch(
    molecules: List[Any],
    ecc_db: Dict[str, Any],
    config,
    num_workers: int = 1,
    seed: Optional[int] = None
) -> Tuple[List[Any], Any]:
    """
    并行去支化与线性化
    """
    from .debranch import Debrancher, DebranchStats

    if num_workers <= 1:
        rng = np.random.default_rng(seed)
        debrancher = Debrancher(config, rng)
        return debrancher.process_all(molecules, ecc_db)

    if not molecules:
        return [], DebranchStats()

    chunks = _split_round_robin(molecules, num_workers)
    tasks = []
    for idx, chunk in enumerate(chunks):
        if chunk:
            tasks.append((chunk, ecc_db, config, seed, idx))

    worker_count = min(num_workers, len(tasks))
    with mp.Pool(worker_count) as pool:
        results = pool.starmap(_worker_debranch, tasks)

    all_linear: List[Any] = []
    stats_list = []
    for linear_molecules, stats in results:
        all_linear.extend(linear_molecules)
        stats_list.append(stats)

    return all_linear, _merge_debranch_stats(stats_list)


def _worker_debranch(molecules, ecc_db, config, seed, worker_id):
    from .debranch import Debrancher

    worker_seed = seed + worker_id if seed is not None else None
    rng = np.random.default_rng(worker_seed)
    debrancher = Debrancher(config, rng)
    return debrancher.process_all(molecules, ecc_db)


def _generate_reads_single(generator_class, config, pool_data, target_reads, seed):
    """单进程生成reads"""
    from .debranch import LinearMoleculePool

    # 重建分子池
    pool = _rebuild_pool(pool_data)

    # 创建生成器
    rng = np.random.default_rng(seed)
    generator = generator_class(config, rng)

    # 生成reads
    return generator.generate(pool, target_reads)


def _worker_generate_reads(generator_class, config, pool_data, target_reads, seed, worker_id):
    """
    Worker进程生成reads

    注意：这个函数在子进程中运行，需要重新导入模块
    """
    import numpy as np

    # 重建分子池
    pool = _rebuild_pool(pool_data)

    # 创建生成器
    rng = np.random.default_rng(seed)
    generator = generator_class(config, rng)

    # 修改read ID前缀以避免冲突
    generator._worker_id = worker_id

    # 生成reads
    reads, stats = generator.generate(pool, target_reads)

    # 更新read ID以包含worker ID
    for read in reads:
        read.read_id = f"w{worker_id}_{read.read_id}"
        if hasattr(read, 'mate_id') and read.mate_id:
            read.mate_id = f"w{worker_id}_{read.mate_id}"

    return reads, stats


def _rebuild_pool(pool_data: Dict[str, Any]):
    """从序列化数据重建分子池"""
    from .debranch import LinearMoleculePool
    from .models import (
        LinearMolecule, Segment, SegmentType, Strand,
        EccDNA, ActiveBranchInfo, register_ecc_length
    )

    # 重建 EccDNA 数据库
    ecc_db = {}
    for ecc_id, ecc_data in pool_data['ecc_db'].items():
        ecc = EccDNA(
            id=ecc_data['id'],
            seq=ecc_data['seq'],
            weight=ecc_data['weight']
        )
        ecc_db[ecc_id] = ecc
        register_ecc_length(ecc_id, len(ecc_data['seq']))

    # 重建分子列表
    molecules = []
    for mol_data in pool_data['molecules']:
        segments = []
        for seg_data in mol_data['segments']:
            seg = Segment(
                ecc_id=seg_data['ecc_id'],
                ecc_offset=seg_data['ecc_offset'],
                length=seg_data['length'],
                strand=Strand(seg_data['strand']),
                segment_type=SegmentType(seg_data['segment_type'])
            )
            segments.append(seg)

        mol = LinearMolecule(
            molecule_id=mol_data['molecule_id'],
            segments=segments,
            source_graph_id=mol_data['source_graph_id'],
            is_from_branch=mol_data.get('is_from_branch', False),
            is_background=mol_data.get('is_background', False),
            background_chrom=mol_data.get('background_chrom'),
            background_start=mol_data.get('background_start'),
            background_end=mol_data.get('background_end'),
            has_chimera=mol_data.get('has_chimera', False),
            chimera_positions=mol_data.get('chimera_positions', []),
            active_branch_anchors=mol_data.get('active_branch_anchors', []),
            debranched_anchors=mol_data.get('debranched_anchors', []),
            junction_covered_possible=mol_data.get('junction_covered_possible', False)
        )
        molecules.append(mol)

    # 创建分子池
    return LinearMoleculePool(molecules, ecc_db, weight_by_length=True)


def serialize_pool(pool) -> Dict[str, Any]:
    """
    序列化分子池以便跨进程传递

    Args:
        pool: LinearMoleculePool

    Returns:
        可序列化的字典
    """
    # 序列化 EccDNA 数据库
    ecc_db_data = {}
    for ecc_id, ecc in pool.ecc_db.items():
        ecc_db_data[ecc_id] = {
            'id': ecc.id,
            'seq': ecc.seq,
            'weight': ecc.weight
        }

    # 序列化分子列表
    molecules_data = []
    for mol in pool.molecules:
        segments_data = []
        for seg in mol.segments:
            segments_data.append({
                'ecc_id': seg.ecc_id,
                'ecc_offset': seg.ecc_offset,
                'length': seg.length,
                'strand': seg.strand.value,
                'segment_type': seg.segment_type.value
            })

        mol_data = {
            'molecule_id': mol.molecule_id,
            'segments': segments_data,
            'source_graph_id': mol.source_graph_id,
            'is_from_branch': mol.is_from_branch,
            'is_background': mol.is_background,
            'background_chrom': mol.background_chrom,
            'background_start': mol.background_start,
            'background_end': mol.background_end,
            'has_chimera': mol.has_chimera,
            'chimera_positions': list(mol.chimera_positions),
            'active_branch_anchors': list(mol.active_branch_anchors),
            'debranched_anchors': list(mol.debranched_anchors),
            'junction_covered_possible': mol.junction_covered_possible
        }
        molecules_data.append(mol_data)

    return {
        'ecc_db': ecc_db_data,
        'molecules': molecules_data
    }


class ProgressTracker:
    """进度跟踪器（用于显示进度）"""

    def __init__(self, total: int, desc: str = "Processing"):
        self.total = total
        self.desc = desc
        self.current = 0
        self._last_percent = -1

    def update(self, n: int = 1):
        """更新进度"""
        self.current += n
        percent = int(100 * self.current / self.total)
        if percent != self._last_percent and percent % 10 == 0:
            logger.info(f"{self.desc}: {percent}% ({self.current}/{self.total})")
            self._last_percent = percent

    def close(self):
        """完成"""
        if self.current < self.total:
            self.current = self.total
        logger.info(f"{self.desc}: 100% ({self.total}/{self.total})")
