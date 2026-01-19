"""Tests for the process module."""

import os
import tempfile
from pathlib import Path

import pandas as pd
import pytest


class TestProcessImports:
    """Test that process module can be imported."""

    def test_import_merge(self):
        """Test merge module import."""
        from ecctoolkit.process.merge import merge_files
        assert callable(merge_files)

    def test_import_parse(self):
        """Test parse module import."""
        from ecctoolkit.process.parse import parse_eccdna
        assert callable(parse_eccdna)

    def test_import_filter(self):
        """Test filter module import."""
        from ecctoolkit.process.filter import filter_by_annotation
        assert callable(filter_by_annotation)

    def test_import_report(self):
        """Test report module import."""
        from ecctoolkit.process.report import generate_reports
        assert callable(generate_reports)


class TestCLICommands:
    """Test CLI command availability."""

    def test_merge_help(self):
        """Test merge --help."""
        import subprocess
        result = subprocess.run(
            ["ecc", "merge", "--help"],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert "Merge multiple CSV files" in result.stdout

    def test_filter_help(self):
        """Test filter --help."""
        import subprocess
        result = subprocess.run(
            ["ecc", "filter", "--help"],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert "Filter eccDNA" in result.stdout

    def test_parse_help(self):
        """Test parse --help."""
        import subprocess
        result = subprocess.run(
            ["ecc", "parse", "--help"],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert "Parse eccDNA CSV" in result.stdout


class TestMerge:
    """Tests for merge functionality."""

    def test_merge_basic(self):
        """Test basic merge functionality."""
        from ecctoolkit.process.merge import merge_files

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create test CSV files
            df1 = pd.DataFrame({
                'eChr': ['chr1', 'chr2'],
                'eStart': [100, 200],
                'eEnd': [500, 600],
            })
            df2 = pd.DataFrame({
                'eChr': ['chr3'],
                'eStart': [300],
                'eEnd': [700],
            })

            df1.to_csv(tmpdir / "sample1.csv", index=False)
            df2.to_csv(tmpdir / "sample2.csv", index=False)

            output_file = tmpdir / "merged.csv"
            merge_files(str(tmpdir), str(output_file))

            # Check result
            assert output_file.exists()
            result = pd.read_csv(output_file)
            assert len(result) == 3
            assert "sample" in result.columns
            assert set(result["sample"]) == {"sample1", "sample2"}


class TestParse:
    """Tests for parse functionality."""

    def test_parse_seqname(self):
        """Test parsing seqname to coordinates."""
        from ecctoolkit.process.parse import parse_eccdna

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create test CSV with seqname
            df = pd.DataFrame({
                'seqname': ['chr1:100-500', 'chr2:200-600'],
                'other_col': ['a', 'b'],
            })
            input_file = tmpdir / "input.csv"
            output_file = tmpdir / "output.csv"
            df.to_csv(input_file, index=False)

            parse_eccdna(str(input_file), str(output_file))

            # Check result
            result = pd.read_csv(output_file)
            assert "eChr" in result.columns
            assert "eStart" in result.columns
            assert "eEnd" in result.columns
            assert result["eChr"].tolist() == ["chr1", "chr2"]
            assert result["eStart"].tolist() == [100, 200]
            assert result["eEnd"].tolist() == [500, 600]


class TestFilter:
    """Tests for filter functionality."""

    def test_filter_by_percent(self):
        """Test filtering by annotation percentage."""
        from ecctoolkit.process.filter import filter_by_annotation

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create test CSV
            df = pd.DataFrame({
                'eChr': ['chr1', 'chr2', 'chr3', 'chr4'],
                'anno_Percent': [90, 50, 85, 70],
            })
            input_file = tmpdir / "input.csv"
            output_file = tmpdir / "output.csv"
            df.to_csv(input_file, index=False)

            filter_by_annotation(str(input_file), str(output_file), min_percent=80)

            # Check result
            result = pd.read_csv(output_file)
            assert len(result) == 2
            assert all(result["anno_Percent"] >= 80)


class TestReport:
    """Tests for report functionality."""

    def test_report_basic(self):
        """Test basic report generation."""
        from ecctoolkit.process.report import generate_reports

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create test CSV
            df = pd.DataFrame({
                'sample': ['s1', 's1', 's2'],
                'eLength': [500, 1000, 800],
            })
            input_file = tmpdir / "input.csv"
            df.to_csv(input_file, index=False)

            output_prefix = tmpdir / "report"
            generate_reports(str(input_file), str(output_prefix))

            # Check result
            assert (tmpdir / "report_sample_counts.csv").exists()
            assert (tmpdir / "report_length_stats.csv").exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
