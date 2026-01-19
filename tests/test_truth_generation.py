"""Tests for truth file generation and read ID parsing."""

import pytest
import tempfile
import os
from pathlib import Path

from ecctoolkit.simulate.lite import (
    parse_art_read_id,
    parse_pbsim_read_id,
    parse_hifi_simple_read_id,
    generate_truth_from_fastq,
)


class TestReadIdParsing:
    """Test read ID parsing functions."""

    def test_parse_art_read_id_basic(self):
        """Test ART read ID parsing with basic format."""
        result = parse_art_read_id("mol_123-45/1")
        assert result['mol_id'] == 'mol_123'
        assert result['read_index'] == 45
        assert result['read_number'] == 1

    def test_parse_art_read_id_r2(self):
        """Test ART read ID parsing for R2."""
        result = parse_art_read_id("ecc_1-100/2")
        assert result['mol_id'] == 'ecc_1'
        assert result['read_index'] == 100
        assert result['read_number'] == 2

    def test_parse_art_read_id_no_slash(self):
        """Test ART read ID parsing without /1 /2 suffix."""
        result = parse_art_read_id("mol_123-45")
        assert result['mol_id'] == 'mol_123'
        assert result['read_index'] == 45
        assert result['read_number'] == 1  # default

    def test_parse_art_read_id_with_annotation(self):
        """Test ART read ID parsing with annotation."""
        result = parse_art_read_id("mol_123-45/1 some annotation")
        assert result['mol_id'] == 'mol_123'
        assert result['read_index'] == 45
        assert result['read_number'] == 1

    def test_parse_pbsim_read_id_basic(self):
        """Test PBSIM2 read ID parsing."""
        result = parse_pbsim_read_id("mol_123_0000001")
        assert result['mol_id'] == 'mol_123'
        assert result['read_index'] == 1

    def test_parse_pbsim_read_id_complex_prefix(self):
        """Test PBSIM2 read ID parsing with complex prefix."""
        result = parse_pbsim_read_id("ecc_10_branch_5_0000025")
        assert result['mol_id'] == 'ecc_10_branch_5'
        assert result['read_index'] == 25

    def test_parse_hifi_simple_read_id_basic(self):
        """Test HiFi simple mode read ID parsing."""
        result = parse_hifi_simple_read_id("ecc_1_hifi_0")
        assert result['mol_id'] == 'ecc_1'
        assert result['read_index'] == 0

    def test_parse_hifi_simple_read_id_complex(self):
        """Test HiFi simple mode read ID parsing with complex ID."""
        result = parse_hifi_simple_read_id("mol_123_branch_hifi_15")
        assert result['mol_id'] == 'mol_123_branch'
        assert result['read_index'] == 15


class TestTruthGeneration:
    """Test truth file generation."""

    @pytest.fixture
    def temp_files(self):
        """Create temporary test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a mock pool CSV
            pool_csv = Path(tmpdir) / "pool.csv"
            pool_csv.write_text(
                "id\tseq\tlength\tecc_ids\trepeat_count\thas_chimera\tis_background\n"
                "mol_1\tACGT\t4\tecc_1\t3\tFalse\tFalse\n"
                "mol_2\tTGCA\t4\tecc_2\t5\tTrue\tFalse\n"
            )

            # Create a mock FASTQ (HiFi simple format)
            fastq = Path(tmpdir) / "reads.fastq"
            fastq.write_text(
                "@mol_1_hifi_0\n"
                "ACGT\n"
                "+\n"
                "IIII\n"
                "@mol_2_hifi_0\n"
                "TGCA\n"
                "+\n"
                "IIII\n"
            )

            yield {
                "tmpdir": tmpdir,
                "pool_csv": str(pool_csv),
                "fastq": str(fastq),
            }

    def test_generate_truth_from_fastq(self, temp_files):
        """Test truth file generation from FASTQ."""
        truth_path = os.path.join(temp_files["tmpdir"], "truth.tsv")

        count = generate_truth_from_fastq(
            fastq_path=temp_files["fastq"],
            pool_csv_path=temp_files["pool_csv"],
            output_truth_path=truth_path,
            platform="HiFi",
            read_id_parser="hifi_simple",
        )

        assert count == 2
        assert os.path.exists(truth_path)

        # Check content
        with open(truth_path) as f:
            lines = f.readlines()
            assert len(lines) == 3  # header + 2 records
            assert "read_id\tplatform" in lines[0]
            assert "mol_1_hifi_0" in lines[1]
            assert "mol_2_hifi_0" in lines[2]

    def test_generate_truth_nonexistent_fastq(self, temp_files):
        """Test truth generation with non-existent FASTQ."""
        truth_path = os.path.join(temp_files["tmpdir"], "truth.tsv")

        count = generate_truth_from_fastq(
            fastq_path="/nonexistent/path.fastq",
            pool_csv_path=temp_files["pool_csv"],
            output_truth_path=truth_path,
            platform="HiFi",
        )

        assert count == 0

    def test_generate_truth_nonexistent_pool(self, temp_files):
        """Test truth generation with non-existent pool CSV."""
        truth_path = os.path.join(temp_files["tmpdir"], "truth.tsv")

        count = generate_truth_from_fastq(
            fastq_path=temp_files["fastq"],
            pool_csv_path="/nonexistent/pool.csv",
            output_truth_path=truth_path,
            platform="HiFi",
        )

        assert count == 0


class TestConfigWithoutMode:
    """Test that mode field has been removed from configs."""

    def test_full_output_scale_config_no_mode(self):
        """Test FullOutputScaleConfig has no mode field."""
        from ecctoolkit.simulate.unified_config import FullOutputScaleConfig
        import dataclasses

        field_names = [f.name for f in dataclasses.fields(FullOutputScaleConfig)]
        assert 'mode' not in field_names
        assert 'cov_ngs' in field_names
        assert 'cov_hifi' in field_names
        assert 'cov_ont' in field_names

    def test_output_scale_params_no_mode(self):
        """Test OutputScaleParams has no mode field."""
        from ecctoolkit.simulate._wip_full.rca_readsim.config import OutputScaleParams
        import dataclasses

        field_names = [f.name for f in dataclasses.fields(OutputScaleParams)]
        assert 'mode' not in field_names
        assert 'Cov_NGS' in field_names
        assert 'Cov_HiFi' in field_names
        assert 'Cov_ONT' in field_names
