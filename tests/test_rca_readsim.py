"""
RCA Read Simulation 模块单元测试

测试核心模块：
- models: 数据结构
- rca_engine: RCA 分子生成
- chimera: 嵌合体注入
- debranch: 去支化/线性化
- fragmentation: 分子打断
"""

import pytest
import numpy as np
from dataclasses import asdict


# =============================================================================
# Models 测试
# =============================================================================

class TestModels:
    """测试核心数据结构"""

    def test_eccdna_creation(self):
        """测试 EccDNA 对象创建"""
        from ecctoolkit.simulate._wip_full.rca_readsim.models import EccDNA

        ecc = EccDNA(id="ecc_1", seq="ATCGATCG", weight=1.0)
        assert ecc.id == "ecc_1"
        assert ecc.length == 8
        assert ecc.seq == "ATCGATCG"

    def test_eccdna_circular_substr(self):
        """测试环状序列提取"""
        from ecctoolkit.simulate._wip_full.rca_readsim.models import EccDNA

        ecc = EccDNA(id="ecc_1", seq="ATCGATCG", weight=1.0)

        # 基本提取
        assert ecc.get_circular_substr(0, 4) == "ATCG"

        # 跨越边界的提取
        assert ecc.get_circular_substr(6, 4) == "CGAT"  # CG + AT

        # 多次循环
        assert ecc.get_circular_substr(0, 16) == "ATCGATCGATCGATCG"

    def test_segment_creation(self):
        """测试 Segment 对象创建"""
        from ecctoolkit.simulate._wip_full.rca_readsim.models import (
            Segment, Strand, SegmentType
        )

        seg = Segment(
            ecc_id="ecc_1",
            ecc_offset=0,
            length=100,
            strand=Strand.FORWARD,
            segment_type=SegmentType.TRUNK
        )

        assert seg.ecc_id == "ecc_1"
        assert seg.length == 100
        assert seg.strand == Strand.FORWARD

    def test_linear_molecule_total_length(self):
        """测试 LinearMolecule 总长度计算"""
        from ecctoolkit.simulate._wip_full.rca_readsim.models import (
            LinearMolecule, Segment, Strand, SegmentType
        )

        seg1 = Segment("ecc_1", 0, 100, Strand.FORWARD, SegmentType.TRUNK)
        seg2 = Segment("ecc_1", 100, 50, Strand.FORWARD, SegmentType.TRUNK)

        mol = LinearMolecule(
            molecule_id="mol_1",
            segments=[seg1, seg2],
            source_graph_id="graph_1",
            is_from_branch=False
        )

        assert mol.total_length == 150

    def test_ecc_length_cache(self):
        """测试 eccDNA 长度缓存"""
        from ecctoolkit.simulate._wip_full.rca_readsim.models import (
            register_ecc_length, get_ecc_length, clear_ecc_length_cache
        )

        clear_ecc_length_cache()

        register_ecc_length("test_ecc", 500)
        assert get_ecc_length("test_ecc") == 500

        # 未注册的返回默认值
        assert get_ecc_length("unknown_ecc", warn_on_missing=False) == 1000

        clear_ecc_length_cache()


# =============================================================================
# RCA Engine 测试
# =============================================================================

class TestRCAEngine:
    """测试 RCA 分子生成引擎"""

    @pytest.fixture
    def config(self):
        """创建测试配置"""
        from ecctoolkit.simulate._wip_full.rca_readsim.config import SimConfig
        config = SimConfig()
        config.rca.mode = "empirical"
        config.rca.mu_R = 3.0
        config.rca.sigma_R = 0.5
        config.rca.k_len = 0.5
        return config

    @pytest.fixture
    def rng(self):
        """创建固定种子的随机数生成器"""
        return np.random.default_rng(42)

    def test_sample_base_repeat(self, config, rng):
        """测试基础重复次数采样"""
        from ecctoolkit.simulate._wip_full.rca_readsim.rca_engine import RCAEngine

        engine = RCAEngine(config, rng)

        # 测试不同长度的 eccDNA
        repeats_500 = [engine.sample_base_repeat(500) for _ in range(100)]
        repeats_2000 = [engine.sample_base_repeat(2000) for _ in range(100)]

        # 短环应该有更高的平均重复数
        assert np.mean(repeats_500) > np.mean(repeats_2000)

        # 所有重复数应在有效范围内
        assert all(config.rca.R_min <= r <= config.rca.R_max for r in repeats_500)

    def test_sample_base_repeat_invalid_length(self, config, rng):
        """测试无效长度的处理"""
        from ecctoolkit.simulate._wip_full.rca_readsim.rca_engine import RCAEngine

        engine = RCAEngine(config, rng)

        with pytest.raises(ValueError, match="must be positive"):
            engine.sample_base_repeat(0)

        with pytest.raises(ValueError, match="must be positive"):
            engine.sample_base_repeat(-100)

    def test_sample_branch_count(self, config, rng):
        """测试分支数量采样"""
        from ecctoolkit.simulate._wip_full.rca_readsim.rca_engine import RCAEngine

        engine = RCAEngine(config, rng)

        # 测试分支数量
        trunk_length = 10000
        ecc_length = 1000
        repeat_count = 10

        branch_counts = [
            engine.sample_branch_count(trunk_length, ecc_length, repeat_count)
            for _ in range(100)
        ]

        # 分支数应为非负
        assert all(b >= 0 for b in branch_counts)

    def test_generate_molecule(self, config, rng):
        """测试分子生成"""
        from ecctoolkit.simulate._wip_full.rca_readsim.rca_engine import RCAEngine
        from ecctoolkit.simulate._wip_full.rca_readsim.models import (
            EccDNA, register_ecc_length, clear_ecc_length_cache
        )

        clear_ecc_length_cache()

        engine = RCAEngine(config, rng)

        ecc = EccDNA(id="test_ecc", seq="A" * 1000, weight=1.0)
        register_ecc_length(ecc.id, ecc.length)

        mol = engine.generate_molecule(ecc, compartment_id=1, effective_repeat=20)

        assert mol.instance_id is not None
        assert mol.source_ecc.id == "test_ecc"
        assert mol.repeat_count == 20
        assert mol.trunk_length == 20 * 1000

        clear_ecc_length_cache()


# =============================================================================
# Chimera 测试
# =============================================================================

class TestChimera:
    """测试嵌合体注入"""

    @pytest.fixture
    def config(self):
        """创建测试配置"""
        from ecctoolkit.simulate._wip_full.rca_readsim.config import SimConfig
        config = SimConfig()
        config.constants.beta = 0.2
        config.constants.chimera_mode = "per_molecule"
        return config

    @pytest.fixture
    def rng(self):
        return np.random.default_rng(42)

    def test_compute_inter_chimera_prob(self, config, rng):
        """测试 inter-chimera 概率计算"""
        from ecctoolkit.simulate._wip_full.rca_readsim.chimera import ChimeraInjector

        injector = ChimeraInjector(config, rng)

        # N=1 时概率为 0
        assert injector.compute_inter_chimera_prob(1) == 0.0

        # N>1 时概率 > 0
        prob_2 = injector.compute_inter_chimera_prob(2)
        prob_5 = injector.compute_inter_chimera_prob(5)

        assert 0 < prob_2 < 1
        assert 0 < prob_5 < 1
        assert prob_5 > prob_2  # 更多分子，更高概率

    def test_compute_chimera_length(self, config, rng):
        """测试 chimera 长度计算"""
        from ecctoolkit.simulate._wip_full.rca_readsim.chimera import ChimeraInjector

        injector = ChimeraInjector(config, rng)

        # 正常情况
        chim_len = injector.compute_chimera_length(10000)
        assert 100 <= chim_len <= config.constants.chimera_max_len

        # donor 长度为 0
        assert injector.compute_chimera_length(0) == 0

        # 短 donor
        chim_len_short = injector.compute_chimera_length(50)
        assert chim_len_short <= 50  # 不超过 donor 长度


# =============================================================================
# Debranch 测试
# =============================================================================

class TestDebranch:
    """测试去支化/线性化"""

    @pytest.fixture
    def config(self):
        """创建测试配置"""
        from ecctoolkit.simulate._wip_full.rca_readsim.config import SimConfig
        config = SimConfig()
        config.debranch.D_eff = 0.9
        return config

    @pytest.fixture
    def rng(self):
        return np.random.default_rng(42)

    def test_validate_non_overlapping_junctions(self):
        """测试重叠 junction 验证"""
        from ecctoolkit.simulate._wip_full.rca_readsim.debranch import (
            validate_non_overlapping_junctions
        )
        from ecctoolkit.simulate._wip_full.rca_readsim.models import (
            ChimeraJunction, Segment, Strand, SegmentType
        )

        # 创建重叠的 junctions
        seg1 = Segment("ecc_1", 0, 100, Strand.FORWARD, SegmentType.CHIMERA, parent_offset=0)
        seg2 = Segment("ecc_1", 0, 100, Strand.FORWARD, SegmentType.CHIMERA, parent_offset=50)

        j1 = ChimeraJunction(position=0, donor_ecc_id="ecc_1", donor_segment=seg1)
        j2 = ChimeraJunction(position=50, donor_ecc_id="ecc_1", donor_segment=seg2)

        # 验证应移除重叠的 junction
        result = validate_non_overlapping_junctions([j1, j2])
        assert len(result) == 1
        assert result[0].position == 0

    def test_debranch_molecule(self, config, rng):
        """测试分子去支化"""
        from ecctoolkit.simulate._wip_full.rca_readsim.debranch import Debrancher
        from ecctoolkit.simulate._wip_full.rca_readsim.models import (
            EccDNA, RCAMoleculeGraph, BranchNode, Segment,
            Strand, SegmentType, register_ecc_length, clear_ecc_length_cache
        )

        clear_ecc_length_cache()

        debrancher = Debrancher(config, rng)

        ecc = EccDNA(id="test_ecc", seq="A" * 1000, weight=1.0)
        register_ecc_length(ecc.id, ecc.length)

        # 创建带分支的分子
        branch_seg = Segment("test_ecc", 0, 200, Strand.FORWARD, SegmentType.BRANCH, parent_offset=500)
        branch = BranchNode(branch_id=0, anchor_pos=500, segment=branch_seg, debranched=False)

        mol = RCAMoleculeGraph(
            instance_id="mol_1",
            source_ecc=ecc,
            compartment_id=1,
            repeat_count=10,
            trunk_length=10000,
            branches=[branch]
        )

        # 执行去支化（高效率应该去除大部分分支）
        debrancher.debranch_molecule(mol)

        # 验证分支状态已更新
        # 注意：由于随机性，这里只检查结构正确性
        assert len(mol.branches) == 1
        assert isinstance(mol.branches[0].debranched, bool)

        clear_ecc_length_cache()


# =============================================================================
# Fragmentation 测试
# =============================================================================

class TestFragmentation:
    """测试分子打断"""

    @pytest.fixture
    def config(self):
        """创建测试配置"""
        from ecctoolkit.simulate._wip_full.rca_readsim.config import SimConfig
        return SimConfig()

    @pytest.fixture
    def rng(self):
        return np.random.default_rng(42)

    def test_fragment_for_ngs(self, config, rng):
        """测试 NGS 打断"""
        from ecctoolkit.simulate._wip_full.rca_readsim.fragmentation import Fragmenter
        from ecctoolkit.simulate._wip_full.rca_readsim.models import (
            LinearMolecule, Segment, Strand, SegmentType, EccDNA,
            register_ecc_length, clear_ecc_length_cache
        )

        clear_ecc_length_cache()

        fragmenter = Fragmenter(config, rng)

        ecc = EccDNA(id="test_ecc", seq="A" * 10000, weight=1.0)
        register_ecc_length(ecc.id, ecc.length)

        seg = Segment("test_ecc", 0, 10000, Strand.FORWARD, SegmentType.TRUNK)
        mol = LinearMolecule(
            molecule_id="mol_1",
            segments=[seg],
            source_graph_id="graph_1",
            is_from_branch=False
        )

        ecc_db = {"test_ecc": ecc}

        fragments, stats = fragmenter.fragment_for_ngs(
            [mol], ecc_db,
            insert_mean=400, insert_std=100,
            min_insert=150, max_insert=800
        )

        # 应该产生多个片段
        assert len(fragments) > 1
        assert stats.input_molecules == 1
        assert stats.output_fragments == len(fragments)

        # 所有片段长度应在范围内
        for frag in fragments:
            assert 150 <= frag.total_length <= 800

        clear_ecc_length_cache()

    def test_fragment_for_hifi(self, config, rng):
        """测试 HiFi 打断"""
        from ecctoolkit.simulate._wip_full.rca_readsim.fragmentation import Fragmenter
        from ecctoolkit.simulate._wip_full.rca_readsim.models import (
            LinearMolecule, Segment, Strand, SegmentType, EccDNA,
            register_ecc_length, clear_ecc_length_cache
        )

        clear_ecc_length_cache()

        fragmenter = Fragmenter(config, rng)

        # 使用足够长的分子以产生多个 HiFi 片段
        ecc = EccDNA(id="test_ecc", seq="A" * 80000, weight=1.0)
        register_ecc_length(ecc.id, ecc.length)

        seg = Segment("test_ecc", 0, 80000, Strand.FORWARD, SegmentType.TRUNK)
        mol = LinearMolecule(
            molecule_id="mol_1",
            segments=[seg],
            source_graph_id="graph_1",
            is_from_branch=False
        )

        ecc_db = {"test_ecc": ecc}

        # 使用较低的 min_length 以确保能产生片段
        fragments, stats = fragmenter.fragment_for_hifi(
            [mol], ecc_db,
            target_mean=20000, target_std=2000,
            min_length=5000, max_length=40000
        )

        # 应该产生片段
        assert len(fragments) >= 1

        # 所有片段长度应在合理范围内
        for frag in fragments:
            # min_length // 2 = 2500 是短分子保留的阈值
            assert frag.total_length >= 2500

        clear_ecc_length_cache()


# =============================================================================
# Config 测试
# =============================================================================

class TestConfig:
    """测试配置系统"""

    def test_default_config(self):
        """测试默认配置"""
        from ecctoolkit.simulate._wip_full.rca_readsim.config import (
            SimConfig, get_default_config
        )

        config = get_default_config()

        assert config.rca.mu_R == 3.0
        assert config.rca.sigma_R == 0.5
        assert config.debranch.D_eff == 0.8

    def test_config_validation(self):
        """测试配置验证"""
        from ecctoolkit.simulate._wip_full.rca_readsim.config import SimConfig

        config = SimConfig()

        # 默认配置应该无警告
        warnings = config.validate()
        assert len(warnings) == 0

        # 无效配置应该产生警告
        config.debranch.D_eff = 1.5  # 超出 [0,1] 范围
        warnings = config.validate()
        assert len(warnings) > 0
        assert any("D_eff" in w for w in warnings)

    def test_config_serialization(self, tmp_path):
        """测试配置序列化"""
        from ecctoolkit.simulate._wip_full.rca_readsim.config import SimConfig

        config = SimConfig()
        config.rca.mu_R = 4.0
        config.seed = 42

        # YAML 序列化
        yaml_path = tmp_path / "config.yaml"
        config.to_yaml(str(yaml_path))

        loaded = SimConfig.from_yaml(str(yaml_path))
        assert loaded.rca.mu_R == 4.0
        assert loaded.seed == 42

        # JSON 序列化
        json_path = tmp_path / "config.json"
        config.to_json(str(json_path))

        loaded = SimConfig.from_json(str(json_path))
        assert loaded.rca.mu_R == 4.0


# =============================================================================
# 集成测试
# =============================================================================

class TestIntegration:
    """端到端集成测试"""

    def test_full_pipeline_small(self, tmp_path):
        """测试完整流水线（小规模）"""
        from ecctoolkit.simulate._wip_full.rca_readsim.config import SimConfig
        from ecctoolkit.simulate._wip_full.rca_readsim.models import (
            EccDNA, register_ecc_length, clear_ecc_length_cache
        )
        from ecctoolkit.simulate._wip_full.rca_readsim.io_utils import build_ecc_db
        from ecctoolkit.simulate._wip_full.rca_readsim.compartment import CompartmentGenerator
        from ecctoolkit.simulate._wip_full.rca_readsim.rca_engine import RCAEngine
        from ecctoolkit.simulate._wip_full.rca_readsim.chimera import ChimeraInjector
        from ecctoolkit.simulate._wip_full.rca_readsim.debranch import Debrancher

        clear_ecc_length_cache()

        # 配置
        config = SimConfig()
        config.seed = 42
        config.rca.mode = "empirical"
        config.spatial.lambda_comp = 2.0  # 小 compartment

        rng = np.random.default_rng(config.seed)

        # 创建测试 eccDNA
        eccdnas = [
            EccDNA(id=f"ecc_{i}", seq="ATCG" * 250, weight=1.0)  # 1000bp each
            for i in range(5)
        ]

        ecc_db = build_ecc_db(eccdnas)
        for ecc in eccdnas:
            register_ecc_length(ecc.id, ecc.length)

        # 生成 compartments
        comp_gen = CompartmentGenerator(config, rng)
        compartments, _ = comp_gen.generate(eccdnas, target_instances=20)

        assert len(compartments) > 0

        # 生成 RCA 分子
        rca_engine = RCAEngine(config, rng)
        molecules, rca_stats = rca_engine.generate_all(compartments, ecc_db)

        assert len(molecules) > 0
        assert rca_stats.total_molecules == len(molecules)

        # 注入 chimera
        chimera_inj = ChimeraInjector(config, rng)
        molecules, chimera_stats = chimera_inj.process_all(molecules)

        # Chimera 可能为 0（取决于 compartment 大小）
        assert chimera_stats.total_chimera_events >= 0

        # 去支化
        debrancher = Debrancher(config, rng)
        linear_mols, debranch_stats = debrancher.process_all(molecules, ecc_db)

        assert len(linear_mols) > 0
        assert debranch_stats.linear_molecules_from_trunk > 0

        clear_ecc_length_cache()


# =============================================================================
# Background DNA 测试
# =============================================================================

class TestBackgroundDNA:
    """测试背景 DNA 模拟"""

    @pytest.fixture
    def temp_reference(self, tmp_path):
        """创建临时参考基因组"""
        ref_path = tmp_path / "test_ref.fa"

        # 创建简单的参考基因组
        with open(ref_path, 'w') as f:
            f.write(">chr1\n")
            f.write("ATCGATCGATCG" * 1000 + "\n")  # 12000bp
            f.write(">chr2\n")
            f.write("GCTAGCTAGCTA" * 500 + "\n")  # 6000bp

        return str(ref_path)

    def test_background_generator_creation(self, temp_reference):
        """测试 BackgroundDNAGenerator 创建"""
        from ecctoolkit.simulate._wip_full.rca_readsim.background import (
            BackgroundDNAGenerator
        )
        from ecctoolkit.simulate._wip_full.rca_readsim.config import BackgroundDNAParams

        params = BackgroundDNAParams(
            enabled=True,
            min_length=100,
            max_length=500,
            length_distribution="uniform"
        )

        generator = BackgroundDNAGenerator(temp_reference, params, seed=42)

        # 验证染色体加载
        assert len(generator.chromosomes) == 2
        assert generator.total_genome_length == 18000  # 12000 + 6000

    def test_background_fragment_generation(self, temp_reference):
        """测试背景 DNA 片段生成"""
        from ecctoolkit.simulate._wip_full.rca_readsim.background import (
            BackgroundDNAGenerator
        )
        from ecctoolkit.simulate._wip_full.rca_readsim.config import BackgroundDNAParams

        params = BackgroundDNAParams(
            enabled=True,
            min_length=100,
            max_length=500,
            length_distribution="uniform"
        )

        generator = BackgroundDNAGenerator(temp_reference, params, seed=42)

        # 生成一个片段
        molecule, sequence, bg_id = generator.generate_fragment(0)

        # 验证
        assert molecule.is_background
        assert molecule.background_chrom in ["chr1", "chr2"]
        assert len(sequence) >= params.min_length
        assert len(sequence) <= params.max_length
        assert bg_id.startswith("BG_")

    def test_background_molecules_batch(self, temp_reference):
        """测试批量生成背景 DNA 分子"""
        from ecctoolkit.simulate._wip_full.rca_readsim.background import (
            BackgroundDNAGenerator
        )
        from ecctoolkit.simulate._wip_full.rca_readsim.config import BackgroundDNAParams
        from ecctoolkit.simulate._wip_full.rca_readsim.models import (
            clear_ecc_length_cache
        )

        clear_ecc_length_cache()

        params = BackgroundDNAParams(
            enabled=True,
            min_length=100,
            max_length=500,
            length_distribution="lognormal",
            length_mean=250,
            length_std=50
        )

        generator = BackgroundDNAGenerator(temp_reference, params, seed=42)

        ecc_db = {}
        molecules, updated_db = generator.generate_molecules(10, ecc_db)

        # 验证
        assert len(molecules) == 10
        assert len(updated_db) >= 10  # 每个分子对应一个 EccDNA 条目
        assert all(m.is_background for m in molecules)

        clear_ecc_length_cache()

    def test_calculate_background_count(self):
        """测试背景分子数量计算"""
        from ecctoolkit.simulate._wip_full.rca_readsim.background import (
            calculate_background_count
        )

        # 10000 reads, 10% 背景
        count = calculate_background_count(10000, 0.1)
        assert count == 1000

        # 5000 reads, 20% 背景
        count = calculate_background_count(5000, 0.2)
        assert count == 1000


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
