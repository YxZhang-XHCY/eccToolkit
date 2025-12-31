"""Tests for parallel processing functionality."""

import pytest


class TestParallelImports:
    """Test that parallel module can be imported."""

    def test_import_parallel(self):
        """Test parallel module import."""
        from ecctoolkit.simulate.rca_readsim.parallel import (
            ParallelConfig,
            get_optimal_workers,
            serialize_pool,
            parallel_generate_reads,
            ProgressTracker
        )

        assert callable(get_optimal_workers)
        assert callable(serialize_pool)
        assert callable(parallel_generate_reads)


class TestParallelConfig:
    """Test ParallelConfig class."""

    def test_default_config(self):
        """Test default parallel config."""
        from ecctoolkit.simulate.rca_readsim.parallel import ParallelConfig

        config = ParallelConfig()
        assert config.num_workers == 1
        assert config.chunk_size == 1000


class TestGetOptimalWorkers:
    """Test worker count calculation."""

    def test_auto_workers(self):
        """Test automatic worker detection."""
        import multiprocessing as mp
        from ecctoolkit.simulate.rca_readsim.parallel import get_optimal_workers

        # Auto mode (0)
        workers = get_optimal_workers(0)
        cpu_count = mp.cpu_count()
        assert 1 <= workers <= cpu_count

    def test_explicit_workers(self):
        """Test explicit worker count."""
        import multiprocessing as mp
        from ecctoolkit.simulate.rca_readsim.parallel import get_optimal_workers

        # Explicit count
        workers = get_optimal_workers(4)
        assert workers <= mp.cpu_count()
        assert workers >= 1

    def test_single_worker(self):
        """Test single worker mode."""
        from ecctoolkit.simulate.rca_readsim.parallel import get_optimal_workers

        workers = get_optimal_workers(1)
        assert workers == 1


class TestPoolSerialization:
    """Test pool serialization for multiprocessing."""

    def test_serialize_empty_pool(self):
        """Test serializing an empty pool."""
        from ecctoolkit.simulate.rca_readsim.parallel import serialize_pool
        from ecctoolkit.simulate.rca_readsim.debranch import LinearMoleculePool

        # Create empty pool
        pool = LinearMoleculePool([], {}, weight_by_length=False)
        data = serialize_pool(pool)

        assert 'ecc_db' in data
        assert 'molecules' in data
        assert len(data['molecules']) == 0

    def test_serialize_with_molecules(self):
        """Test serializing a pool with molecules."""
        from ecctoolkit.simulate.rca_readsim.parallel import serialize_pool
        from ecctoolkit.simulate.rca_readsim.debranch import LinearMoleculePool
        from ecctoolkit.simulate.rca_readsim.models import (
            LinearMolecule, Segment, SegmentType, Strand,
            EccDNA, register_ecc_length
        )

        # Create a test molecule
        ecc = EccDNA(id="ecc_1", seq="ACGT" * 100, weight=1.0)
        ecc_db = {"ecc_1": ecc}
        register_ecc_length("ecc_1", 400)

        mol = LinearMolecule(
            molecule_id="mol_1",
            segments=[Segment(
                ecc_id="ecc_1",
                ecc_offset=0,
                length=400,
                strand=Strand.FORWARD,
                segment_type=SegmentType.TRUNK
            )],
            source_graph_id="graph_1"
        )

        pool = LinearMoleculePool([mol], ecc_db, weight_by_length=True)
        data = serialize_pool(pool)

        assert len(data['ecc_db']) == 1
        assert len(data['molecules']) == 1
        assert data['molecules'][0]['molecule_id'] == "mol_1"


class TestProgressTracker:
    """Test progress tracking."""

    def test_progress_tracker(self):
        """Test progress tracker."""
        from ecctoolkit.simulate.rca_readsim.parallel import ProgressTracker

        tracker = ProgressTracker(100, "Test")
        assert tracker.total == 100
        assert tracker.current == 0

        tracker.update(50)
        assert tracker.current == 50

        tracker.close()
        assert tracker.current == 100


class TestCLIParallelOptions:
    """Test CLI with parallel options."""

    def test_sim_reads_help_shows_threads(self):
        """Test sim-reads --help shows threads option."""
        import subprocess

        result = subprocess.run(
            ["ecc", "sim-reads", "--help"],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert "--threads" in result.stdout or "-t" in result.stdout


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
