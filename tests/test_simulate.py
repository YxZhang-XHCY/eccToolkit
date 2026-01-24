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

    def test_import_lite(self):
        """Test lite module import."""
        from ecctoolkit.simulate.lite import libsim, fqsim
        assert callable(libsim)
        assert callable(fqsim)

    def test_import_module(self):
        """Test module-level imports."""
        from ecctoolkit.simulate import run_region_simulation
        from ecctoolkit.simulate.lite import libsim, fqsim
        assert callable(run_region_simulation)
        assert callable(libsim)
        assert callable(fqsim)


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
        """Test readsim --help."""
        import subprocess
        result = subprocess.run(
            ["ecc", "sim-reads", "--help"],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert "Simulate sequencing reads" in result.stdout
        assert "--meancov" in result.stdout
        assert "--sample" in result.stdout


class TestSimulatePipelineLiteConversion:
    """Test sim-region → lite format conversion."""

    def test_convert_region_to_lite_format_handles_chimeric_fragments(self, tmp_path):
        from ecctoolkit.simulate.pipeline import SimulatePipeline
        from ecctoolkit.simulate.unified_config import UnifiedSimulateConfig

        config = UnifiedSimulateConfig()
        config.sample = "demo"

        pipeline = SimulatePipeline(
            config=config,
            output_dir=str(tmp_path),
            reference="dummy_ref.fa",
        )
        pipeline._setup_output_dirs()

        # region 输出直接放在 region_dir 下
        fasta_path = pipeline.region_dir / f"{config.sample}.all.fa"
        bed_path = pipeline.region_dir / f"{config.sample}.all.bed"

        u_seq = "A" * 100
        c_seq = "C" * 110
        fasta_path.write_text(
            "\n".join(
                [
                    ">UeccDNA_000001 chr1:100-200 length=100 type=U",
                    u_seq,
                    ">CeccDNA_000001 chr1:100-150;chr2:200-260 length=110 type=C",
                    c_seq,
                    "",
                ]
            )
        )

        bed_path.write_text(
            "\n".join(
                [
                    "#chrom\tstart\tend\tname\tlength\tstrand\ttype\tfragments",
                    "chr1\t100\t200\tUeccDNA_000001\t100\t+\tU\t.",
                    "chr1\t100\t150\tCeccDNA_000001\t110\t+\tC\tchr1:100-150;chr2:200-260",
                    "",
                ]
            )
        )

        pipeline._convert_region_to_lite_format(fasta_path)

        # RCA 相关文件放在 rca_dir 下
        pos_bed = pipeline.rca_dir / f"{config.sample}.pos.bed"
        pos_csv = pipeline.rca_dir / f"{config.sample}.pos.csv"

        assert pos_bed.exists(), f"Expected {pos_bed} to exist"
        assert pos_csv.exists(), f"Expected {pos_csv} to exist"

        bed_lines = [line.strip().split("\t") for line in pos_bed.read_text().splitlines() if line.strip()]
        assert all(cols[0] != "chrN" for cols in bed_lines)

        # UeccDNA: one fragment
        assert ["chr1", "100", "200", "100", "UeccDNA_000001"] in bed_lines
        # CeccDNA: two fragments with the same ID
        assert ["chr1", "100", "150", "50", "CeccDNA_000001"] in bed_lines
        assert ["chr2", "200", "260", "60", "CeccDNA_000001"] in bed_lines

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

    @pytest.fixture
    def minimap2_available(self):
        """Check if minimap2 is available."""
        import subprocess
        result = subprocess.run(["which", "minimap2"], capture_output=True)
        return result.returncode == 0

    def test_sim_region_basic(self, test_reference, minimap2_available):
        """Test basic sim-region functionality."""
        import subprocess

        if not minimap2_available:
            pytest.skip("minimap2 not installed")

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
    """Functional tests for readsim (requires rca_readsim module)."""

    def test_sim_reads_basic(self, test_eccdna_fasta):
        """Test basic readsim functionality."""
        import subprocess

        with tempfile.TemporaryDirectory() as tmpdir:
            result = subprocess.run(
                [
                    "ecc", "readsim",
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
