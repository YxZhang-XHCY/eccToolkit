"""Tests for background DNA functionality."""

import tempfile
from pathlib import Path

import pytest


class TestBackgroundImports:
    """Test that background module can be imported."""

    def test_import_background(self):
        """Test background module import."""
        from ecctoolkit.simulate.rca_readsim.background import (
            BackgroundDNAGenerator,
            calculate_background_count,
        )
        assert callable(calculate_background_count)

    def test_import_config_with_background(self):
        """Test config with background params."""
        from ecctoolkit.simulate.rca_readsim.config import SimConfig, BackgroundDNAParams
        config = SimConfig()
        assert hasattr(config, 'background')
        assert isinstance(config.background, BackgroundDNAParams)
        assert config.background.enabled == False
        assert config.background.ratio == 0.1


class TestBackgroundDNACalculation:
    """Test background DNA count calculation."""

    def test_calculate_count(self):
        """Test background count calculation."""
        from ecctoolkit.simulate.rca_readsim.background import calculate_background_count

        # 10% of 1000 = 100
        assert calculate_background_count(1000, 0.1) == 100

        # 30% of 10000 = 3000
        assert calculate_background_count(10000, 0.3) == 3000

        # 0% = 0
        assert calculate_background_count(1000, 0.0) == 0


class TestBackgroundDNAGenerator:
    """Test background DNA generator."""

    @pytest.fixture
    def test_reference(self):
        """Create a temporary test reference genome."""
        import random
        random.seed(42)

        tmpdir = tempfile.mkdtemp()
        ref_path = Path(tmpdir) / "test_ref.fa"

        # Create small reference
        with open(ref_path, 'w') as f:
            for chrom, length in [('chr1', 10000), ('chr2', 5000)]:
                f.write(f'>{chrom}\n')
                seq = ''.join(random.choices('ACGT', k=length))
                for i in range(0, length, 80):
                    f.write(seq[i:i+80] + '\n')

        yield ref_path

        # Cleanup
        import shutil
        shutil.rmtree(tmpdir)

    def test_generator_init(self, test_reference):
        """Test generator initialization."""
        from ecctoolkit.simulate.rca_readsim.background import BackgroundDNAGenerator
        from ecctoolkit.simulate.rca_readsim.config import BackgroundDNAParams

        params = BackgroundDNAParams(
            enabled=True,
            fasta_path=str(test_reference),
            ratio=0.1,
            min_length=100,
            max_length=1000
        )

        generator = BackgroundDNAGenerator(
            fasta_path=str(test_reference),
            params=params,
            seed=42
        )

        assert len(generator.chromosomes) == 2
        assert generator.total_genome_length == 15000

    def test_generate_molecules(self, test_reference):
        """Test molecule generation."""
        from ecctoolkit.simulate.rca_readsim.background import BackgroundDNAGenerator
        from ecctoolkit.simulate.rca_readsim.config import BackgroundDNAParams

        params = BackgroundDNAParams(
            enabled=True,
            fasta_path=str(test_reference),
            ratio=0.1,
            min_length=100,
            max_length=500,
            length_distribution="uniform"
        )

        generator = BackgroundDNAGenerator(
            fasta_path=str(test_reference),
            params=params,
            seed=42
        )

        ecc_db = {}
        molecules, updated_ecc_db = generator.generate_molecules(10, ecc_db)

        assert len(molecules) == 10
        assert len(updated_ecc_db) == 10

        # Check molecule properties
        for mol in molecules:
            assert mol.is_background == True
            assert mol.background_chrom in ['chr1', 'chr2']
            assert mol.background_start >= 0
            assert mol.background_end > mol.background_start
            assert 100 <= (mol.background_end - mol.background_start) <= 500


class TestModelsWithBackground:
    """Test models with background fields."""

    def test_linear_molecule_background_fields(self):
        """Test LinearMolecule has background fields."""
        from ecctoolkit.simulate.rca_readsim.models import (
            LinearMolecule,
            Segment,
            SegmentType,
            Strand,
        )

        # Create background molecule
        mol = LinearMolecule(
            molecule_id="bg_1",
            segments=[Segment(
                ecc_id="BG_chr1_100_500",
                ecc_offset=0,
                length=400,
                strand=Strand.FORWARD,
                segment_type=SegmentType.BACKGROUND
            )],
            source_graph_id="bg_1",
            is_background=True,
            background_chrom="chr1",
            background_start=100,
            background_end=500
        )

        assert mol.is_background == True
        assert mol.background_chrom == "chr1"
        assert mol.background_start == 100
        assert mol.background_end == 500

    def test_sequenced_read_background_fields(self):
        """Test SequencedRead has background fields."""
        from ecctoolkit.simulate.rca_readsim.models import SequencedRead

        read = SequencedRead(
            read_id="read_1",
            platform="NGS",
            sequence="ACGT" * 50,
            is_background=True,
            background_chrom="chr1",
            background_start=100,
            background_end=500
        )

        assert read.is_background == True

        truth = read.to_truth_dict()
        assert truth["source_type"] == "background"
        assert truth["background_chrom"] == "chr1"
        assert truth["background_start"] == 100
        assert truth["background_end"] == 500

    def test_sequenced_read_eccdna_fields(self):
        """Test SequencedRead eccDNA type."""
        from ecctoolkit.simulate.rca_readsim.models import SequencedRead

        read = SequencedRead(
            read_id="read_1",
            platform="NGS",
            sequence="ACGT" * 50,
            is_background=False,
            source_ecc_ids=["eccDNA_1"]
        )

        assert read.is_background == False

        truth = read.to_truth_dict()
        assert truth["source_type"] == "eccDNA"


class TestCLIBackgroundOptions:
    """Test CLI with background options."""

    def test_sim_reads_help_shows_background_options(self):
        """Test sim-reads --help shows background options."""
        import subprocess
        result = subprocess.run(
            ["ecc", "sim-reads", "--help"],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert "--ref" in result.stdout or "--reference" in result.stdout
        assert "--linear-ratio" in result.stdout


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
