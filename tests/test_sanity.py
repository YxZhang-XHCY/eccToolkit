"""Sanity checks for RCA simulation regression."""

import numpy as np
import warnings

from ecctoolkit.simulate.rca_readsim.config import SimConfig
from ecctoolkit.simulate.rca_readsim.models import (
    EccDNA,
    LinearMolecule,
    Segment,
    SegmentType,
    Strand,
    register_ecc_length,
    RCAMoleculeGraph,
)
from ecctoolkit.simulate.rca_readsim.library.hifi import HiFiLibraryGenerator
from ecctoolkit.simulate.rca_readsim.debranch import LinearMoleculePool, Debrancher
from ecctoolkit.simulate.rca_readsim.rca_engine import RCAEngine
from ecctoolkit.simulate.rca_readsim.chimera import ChimeraInjector


def _make_eccdna(ecc_id: str, length: int, weight: float = 1.0) -> EccDNA:
    seq = ("ACGT" * ((length // 4) + 1))[:length]
    return EccDNA(id=ecc_id, seq=seq, weight=weight)


def _make_linear_molecule(
    ecc: EccDNA,
    length: int,
    mol_id: str,
    segment_type: SegmentType = SegmentType.TRUNK,
    is_background: bool = False,
) -> LinearMolecule:
    segment = Segment(
        ecc_id=ecc.id,
        ecc_offset=0,
        length=length,
        strand=Strand.FORWARD,
        segment_type=segment_type,
    )
    repeat_count = max(1, int(length / max(1, ecc.length)))
    return LinearMolecule(
        molecule_id=mol_id,
        segments=[segment],
        source_graph_id=mol_id,
        is_from_branch=False,
        repeat_count=repeat_count,
        source_ecc_length=ecc.length,
        is_background=is_background,
        background_chrom="chr1" if is_background else None,
        background_start=0 if is_background else None,
        background_end=length if is_background else None,
    )


def test_hifi_length_distribution():
    config = SimConfig()
    rng = np.random.default_rng(42)

    ecc = _make_eccdna("ecc_long", 1000)
    register_ecc_length(ecc.id, ecc.length)
    ecc_db = {ecc.id: ecc}

    mol = _make_linear_molecule(ecc, 200_000, "mol_long")
    pool = LinearMoleculePool([mol], ecc_db, weight_by_length=True)
    hifi = HiFiLibraryGenerator(config, rng)

    reads, _ = hifi.generate(pool, target_reads=150)
    lengths = np.array([len(r.sequence) for r in reads], dtype=float)

    mean_len = lengths.mean()
    std_len = lengths.std()
    in_range = ((lengths >= 15000) & (lengths <= 25000)).mean()

    assert abs(mean_len - 20000) <= 500
    assert abs(std_len - 1200) <= 300
    assert in_range >= 0.95


def test_ctcr_ratio_by_ring_length():
    config = SimConfig()
    rng = np.random.default_rng(123)

    ecc_small = _make_eccdna("ecc_3kb", 3000)
    ecc_large = _make_eccdna("ecc_20kb", 20000)
    register_ecc_length(ecc_small.id, ecc_small.length)
    register_ecc_length(ecc_large.id, ecc_large.length)
    ecc_db = {ecc_small.id: ecc_small, ecc_large.id: ecc_large}

    mol_small = _make_linear_molecule(ecc_small, 200_000, "mol_small")
    mol_large = _make_linear_molecule(ecc_large, 200_000, "mol_large")
    pool = LinearMoleculePool([mol_small, mol_large], ecc_db, weight_by_length=False)
    hifi = HiFiLibraryGenerator(config, rng)

    reads, _ = hifi.generate(pool, target_reads=240)

    small_flags = []
    large_flags = []
    for read in reads:
        if len(read.source_ecc_ids) != 1:
            continue
        ecc_id = read.source_ecc_ids[0]
        is_ctcr = read.repeat_count_truth >= 2
        if ecc_id == ecc_small.id:
            small_flags.append(is_ctcr)
        elif ecc_id == ecc_large.id:
            large_flags.append(is_ctcr)

    assert len(small_flags) >= 80
    assert len(large_flags) >= 80

    p_small = np.mean(small_flags)
    p_large = np.mean(large_flags)
    assert p_small > 0.6
    assert p_large < 0.1


def test_background_repeat_truth_zero():
    config = SimConfig()
    rng = np.random.default_rng(7)

    ecc_bg = _make_eccdna("BG_chr1_0_20000", 20000)
    register_ecc_length(ecc_bg.id, ecc_bg.length)
    ecc_db = {ecc_bg.id: ecc_bg}

    mol_bg = _make_linear_molecule(
        ecc_bg,
        20000,
        "bg_mol",
        segment_type=SegmentType.BACKGROUND,
        is_background=True,
    )
    pool = LinearMoleculePool([mol_bg], ecc_db, weight_by_length=True)
    hifi = HiFiLibraryGenerator(config, rng)

    reads, _ = hifi.generate(pool, target_reads=2)
    assert reads
    for read in reads:
        assert read.repeat_count_truth == 0
        assert read.repeat_count_by_source == {}


def test_chimera_rate_and_break_distribution():
    config = SimConfig()
    rng = np.random.default_rng(202)

    ecc = _make_eccdna("ecc_ref", 1000)
    register_ecc_length(ecc.id, ecc.length)
    ecc_db = {ecc.id: ecc}

    engine = RCAEngine(config, rng)
    molecules = []
    num_compartments = 15
    comp_size = 5
    repeat_count = 30

    for comp_id in range(num_compartments):
        for _ in range(comp_size):
            mol = engine.generate_molecule(ecc, comp_id, repeat_count)
            molecules.append(mol)

    chimera = ChimeraInjector(config, rng)
    molecules, chimera_stats = chimera.process_all(molecules)
    chimera_rate = chimera_stats.molecules_with_chimera / max(1, len(molecules))

    assert 0.05 <= chimera_rate <= 0.35

    debrancher = Debrancher(config, rng)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        linear_mols, _ = debrancher.process_all(molecules, ecc_db)

    frag_lengths = [
        mol.total_length
        for mol in linear_mols
        if (not mol.is_from_branch) and "_frag" in mol.molecule_id
    ]
    assert len(frag_lengths) >= 30

    lengths = np.array(frag_lengths, dtype=float)
    cv = lengths.std() / max(1.0, lengths.mean())
    p10 = np.percentile(lengths, 10)
    p90 = np.percentile(lengths, 90)
    ratio = p90 / max(1.0, p10)

    assert 0.4 <= cv <= 1.6
    assert 2.0 <= ratio <= 20.0


def test_reproducibility_with_fixed_seed():
    def _run(seed: int) -> tuple[int, int, int]:
        config = SimConfig()
        rng = np.random.default_rng(seed)

        ecc = _make_eccdna("ecc_small", 3000)
        register_ecc_length(ecc.id, ecc.length)
        ecc_db = {ecc.id: ecc}

        engine = RCAEngine(config, rng)
        molecules = []
        for comp_id in range(10):
            for _ in range(3):
                mol = engine.generate_molecule(ecc, comp_id, 20)
                molecules.append(mol)

        chimera = ChimeraInjector(config, rng)
        molecules, _ = chimera.process_all(molecules)

        debrancher = Debrancher(config, rng)
        linear_mols, _ = debrancher.process_all(molecules, ecc_db)

        pool = LinearMoleculePool(linear_mols, ecc_db, weight_by_length=True)
        hifi = HiFiLibraryGenerator(config, rng)
        reads, _ = hifi.generate(pool, target_reads=40)

        ctcr = sum(1 for r in reads if r.repeat_count_truth >= 2)
        chimera_reads = sum(1 for r in reads if r.has_inter_chimera)
        return len(reads), ctcr, chimera_reads

    stats_a = _run(999)
    stats_b = _run(999)
    assert stats_a == stats_b
