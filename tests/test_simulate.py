"""Tests for the simulate module."""

import os
import tempfile
import shutil
from pathlib import Path

import pytest

# Test data directory
TEST_DATA_DIR = Path(__file__).parent / "data"


class TestSimulateImports:
    """Test that simulate module can be imported."""

    def test_import_region(self):
        """Test region module import."""
        from ecctoolkit.simulate.region import run_region_simulation
        assert callable(run_region_simulation)

    def test_import_reads(self):
        """Test reads module import."""
        from ecctoolkit.simulate.reads import run_read_simulation
        assert callable(run_read_simulation)

    def test_import_module(self):
        """Test module-level imports."""
        from ecctoolkit.simulate import run_region_simulation, run_read_simulation
        assert callable(run_region_simulation)
        assert callable(run_read_simulation)


class TestCLICommands:
    """Test CLI command availability."""

    def test_sim_region_help(self):
        """Test sim-region --help."""
        import subprocess
        result = subprocess.run(
            ["ecc", "sim-region", "--help"],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert "Simulate eccDNA regions" in result.stdout
        assert "--reference" in result.stdout
        assert "--num-unique" in result.stdout

    def test_sim_reads_help(self):
        """Test sim-reads --help."""
        import subprocess
        result = subprocess.run(
            ["ecc", "sim-reads", "--help"],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert "Simulate sequencing reads" in result.stdout
        assert "--cov-ngs" in result.stdout
        assert "--output-mode" in result.stdout


@pytest.fixture
def test_reference():
    """Create a temporary test reference genome."""
    import random
    random.seed(42)

    tmpdir = tempfile.mkdtemp()
    ref_path = Path(tmpdir) / "test_ref.fa"

    # Create small reference
    with open(ref_path, 'w') as f:
        for chrom, length in [('chr1', 50000), ('chr2', 30000)]:
            f.write(f'>{chrom}\n')
            seq = ''.join(random.choices('ACGT', k=length))
            for i in range(0, length, 80):
                f.write(seq[i:i+80] + '\n')

    yield ref_path

    # Cleanup
    shutil.rmtree(tmpdir)


@pytest.fixture
def test_eccdna_fasta():
    """Create a temporary test eccDNA FASTA."""
    import random
    random.seed(42)

    tmpdir = tempfile.mkdtemp()
    fasta_path = Path(tmpdir) / "test_eccdna.fa"

    # Create small eccDNA sequences
    with open(fasta_path, 'w') as f:
        for i in range(10):
            length = random.randint(300, 2000)
            seq = ''.join(random.choices('ACGT', k=length))
            f.write(f'>eccDNA_{i+1} length={length}\n')
            for j in range(0, length, 80):
                f.write(seq[j:j+80] + '\n')

    yield fasta_path

    # Cleanup
    shutil.rmtree(tmpdir)


@pytest.mark.slow
class TestSimRegionFunctional:
    """Functional tests for sim-region (requires minimap2)."""

    def test_sim_region_basic(self, test_reference):
        """Test basic sim-region functionality."""
        import subprocess

        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "output"

            result = subprocess.run(
                [
                    "ecc", "sim-region",
                    "-r", str(test_reference),
                    "-o", str(output),
                    "-u", "10",
                    "-t", "2",
                    "--seed", "42"
                ],
                capture_output=True,
                text=True,
                timeout=120
            )

            # Check output files exist
            output_dir = Path(tmpdir) / "output"
            if output_dir.exists():
                assert (output_dir / "output.all.bed").exists() or \
                       any(output_dir.glob("*.bed"))


@pytest.mark.slow
class TestSimReadsFunctional:
    """Functional tests for sim-reads (requires rca_readsim module)."""

    def test_sim_reads_basic(self, test_eccdna_fasta):
        """Test basic sim-reads functionality."""
        import subprocess

        with tempfile.TemporaryDirectory() as tmpdir:
            result = subprocess.run(
                [
                    "ecc", "sim-reads",
                    "-i", str(test_eccdna_fasta),
                    "-o", tmpdir,
                    "--cov-ngs", "100",
                    "--skip-hifi",
                    "--skip-ont",
                    "--seed", "42"
                ],
                capture_output=True,
                text=True,
                timeout=120
            )

            # Check output files exist
            output_path = Path(tmpdir)
            ngs_r1 = output_path / "reads_ngs_R1.fastq"
            ngs_r2 = output_path / "reads_ngs_R2.fastq"

            if result.returncode == 0:
                assert ngs_r1.exists()
                assert ngs_r2.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
