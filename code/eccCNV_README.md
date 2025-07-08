# eccCNV.py

## Description

eccCNV is a comprehensive tool for analyzing the enrichment of eccDNA (extrachromosomal circular DNA) in Copy Number Variation (CNV) regions. The tool uses permutation testing to determine whether eccDNA elements are significantly enriched or depleted in CNV gain, loss, or neutral regions compared to random genomic distributions.

## Features

- **CNV Region Analysis**: Analyzes eccDNA enrichment in gain, loss, and neutral CNV regions
- **Permutation Testing**: Uses bedtools shuffle for robust statistical testing
- **Parallel Processing**: Multi-core support for efficient computation
- **Comprehensive Statistics**: Mann-Whitney U tests, Kruskal-Wallis tests, and correlation analysis
- **Quality Control**: Built-in QC plots and data validation
- **Rich Visualization**: Multiple plot types including violin plots, scatter plots, and heatmaps
- **Detailed Reporting**: Generates comprehensive markdown reports
- **Export Support**: Exports data for external visualization tools

## Installation Requirements

### System Dependencies
```bash
# Install bedtools (required)
conda install -c bioconda bedtools
# or
sudo apt-get install bedtools  # Ubuntu/Debian
```

### Python Dependencies
```bash
pip install pandas numpy matplotlib seaborn scipy statsmodels
```

### Required Software
- **bedtools** (version 2.30.0 or higher)
- **Python** 3.7 or higher

### Required Python Packages
- `pandas` - Data manipulation and analysis
- `numpy` - Numerical computing  
- `matplotlib` - Plotting and visualization
- `seaborn` - Statistical data visualization
- `scipy` - Statistical functions and tests
- `statsmodels` - Statistical modeling and multiple testing correction
- `pathlib` - Path handling (standard library)
- `multiprocessing` - Parallel processing (standard library)
- `subprocess` - System command execution (standard library)

## Usage

### Basic Usage
```bash
python eccCNV.py -i eccDNA_file.csv -c CNV_regions.bed
```

### Advanced Usage
```bash
python eccCNV.py \
    -i sample1.csv sample2.csv sample3.csv \
    -c CNV_regions.bed \
    -w ./analysis_output \
    -n 10000 \
    --cores 16 \
    --prefix HeLa_experiment \
    --keep-sex-chr
```

### Multiple File Analysis
```bash
python eccCNV.py \
    -i rep1_eccDNA.csv rep2_eccDNA.csv rep3_eccDNA.csv \
    -c CNV_data.bed \
    -w results/ \
    --prefix timecourse_24h
```

## Parameters

### Required Arguments
- `-i`, `--input`: Input eccDNA CSV file(s). Multiple files will be merged.
- `-c`, `--cnv`: CNV BED file with specific format (see Input Formats section)

### Optional Arguments
- `-g`, `--genome`: Genome file for bedtools shuffle (default: built-in hg38)
- `-w`, `--workdir`: Working directory for output (default: current directory)
- `-n`, `--nshuffle`: Number of permutations for enrichment testing (default: 1000)
- `--cores`: Number of CPU cores to use (default: auto-detect, max available - 2)
- `--prefix`: Prefix for all output files and directories
- `--keep-sex-chr`: Keep sex chromosomes (chrX, chrY) in analysis (default: filter out)
- `--version`: Show version information

## Input File Formats

### eccDNA CSV Files
Required columns:
- `eChr` - Chromosome name (e.g., "chr1", "chr2")
- `eStart` - Start position (0-based or 1-based)
- `eEnd` - End position
- `eName` - eccDNA element name/identifier

Optional columns:
- `eLength` - Length of eccDNA element (calculated if missing)
- Additional metadata columns are preserved

**Example:**
```csv
eChr,eStart,eEnd,eName,eLength
chr1,1000000,1005000,ecc1,5000
chr1,2000000,2003000,ecc2,3000
chr2,500000,510000,ecc3,10000
```

### CNV BED File
Tab-separated file with specific column order (no header):

**Columns:**
1. `chr` - Chromosome name
2. `start` - Start position (0-based)
3. `end` - End position (1-based)
4. `segment_mean` - Log2 ratio of copy number
5. `probe_count` - Number of probes in segment
6. `cnv_status` - CNV status: '+' (gain), '-' (loss), '0' (neutral)
7. `cell_line` - Cell line identifier
8. `id` - Unique segment ID
9. `extra` - Additional information (optional)

**Example:**
```
chr1	1000000	2000000	0.5	150	+	HeLa	seg1	info
chr1	3000000	4000000	-0.3	200	-	HeLa	seg2	info
chr2	5000000	6000000	0.1	120	0	HeLa	seg3	info
```

### Genome File (Optional)
If providing a custom genome file, use the format:
```
chr1	248956422
chr2	242193529
chr3	198295559
...
```

## Output Structure

### Directory Organization
```
output_directory/
├── [prefix_]CNV_results/          # Analysis results
├── [prefix_]CNV_figures/          # Generated plots
├── [prefix_]CNV_plot_data/        # Data for external plotting
└── temp/                          # Temporary files
```

### Results Files
1. **[prefix_]CNV_enrichment_results.csv** - Main enrichment statistics
2. **[prefix_]CNV_detailed_stats.csv** - Per-CNV region eccDNA density
3. **[prefix_]CNV_statistical_tests.csv** - Pairwise statistical comparisons
4. **[prefix_]CNV_analysis_report.md** - Comprehensive analysis report

### Figure Files
1. **[prefix_]quality_control.png** - eccDNA distribution and size plots
2. **[prefix_]CNV_analysis.png** - CNV-eccDNA relationship plots
3. **[prefix_]CNV_enrichment_heatmap.png** - Summary heatmap with significance
4. **[prefix_]null_dist_CNV_*.png** - Permutation test distributions

### Plot Data Files (for external visualization)
1. **[prefix_]chromosome_distribution.csv** - eccDNA counts per chromosome
2. **[prefix_]size_distribution.csv** - eccDNA size distribution data
3. **[prefix_]CNV_eccdna_density_data.csv** - Raw CNV-eccDNA density data
4. **[prefix_]CNV_summary_statistics.csv** - CNV group summary statistics
5. **[prefix_]CNV_enrichment_summary_data.csv** - Enrichment values for plotting
6. **[prefix_]CNV_null_distributions.csv** - Complete null distribution data

## Statistical Methods

### Permutation Testing
- Uses bedtools shuffle with `-chrom` flag to maintain chromosome distribution
- Generates null distribution by randomly shuffling eccDNA positions
- Calculates empirical p-values based on null distribution
- Default: 1000 permutations (adjustable with `-n` parameter)

### Multiple Testing Correction
- Applies Benjamini-Hochberg FDR correction
- Calculates q-values for all enrichment tests
- Controls false discovery rate across multiple CNV categories

### Statistical Tests
- **Mann-Whitney U test**: Pairwise comparisons between CNV groups
- **Kruskal-Wallis test**: Overall comparison of all three groups
- **Spearman correlation**: Correlation between segment mean and eccDNA density

### Enrichment Metrics
- **Fold enrichment**: Observed / Expected ratio
- **Z-score**: Standard deviations from null mean
- **Percentile**: Percentage of null values below observed
- **P-value**: Empirical probability from permutation test
- **Q-value**: FDR-corrected p-value

## Example Workflows

### Single Sample Analysis
```bash
# Basic analysis of one eccDNA dataset
python eccCNV.py \
    -i HeLa_eccDNA.csv \
    -c HeLa_CNV.bed \
    -w HeLa_analysis/
```

### Multi-replicate Analysis
```bash
# Combine multiple replicates
python eccCNV.py \
    -i HeLa_rep1.csv HeLa_rep2.csv HeLa_rep3.csv \
    -c HeLa_CNV.bed \
    -w HeLa_combined/ \
    --prefix HeLa_3rep
```

### High-resolution Analysis
```bash
# Use more permutations for higher precision
python eccCNV.py \
    -i eccDNA_data.csv \
    -c CNV_regions.bed \
    -n 10000 \
    --cores 32 \
    --prefix high_res
```

### Sex Chromosome Analysis
```bash
# Include sex chromosomes in analysis
python eccCNV.py \
    -i eccDNA_data.csv \
    -c CNV_regions.bed \
    --keep-sex-chr \
    --prefix full_genome
```

## Performance Considerations

### Memory Usage
- Scales with number of eccDNA elements and CNV regions
- Parallel processing increases memory usage
- Typical usage: 1-4 GB RAM for standard datasets

### Processing Time
- Depends on number of permutations and CPU cores
- 1000 permutations: ~5-15 minutes on 8 cores
- 10000 permutations: ~30-60 minutes on 16 cores

### Optimization Tips
- Use appropriate number of cores for your system
- Increase permutations for higher statistical precision
- Filter out unnecessary chromosomes if not needed

## Quality Control

### Input Validation
- Checks for required columns in input files
- Validates coordinate ranges and data types
- Reports data filtering steps and statistics

### Chromosome Filtering
- Automatically filters mitochondrial DNA (chrM)
- Optionally filters sex chromosomes (chrX, chrY)
- Removes non-standard chromosome contigs

### Data Quality Metrics
- eccDNA size distribution analysis
- Chromosome distribution validation
- CNV region statistics and validation

## Interpretation Guidelines

### P-value Significance Levels
- `***`: q < 0.001 (highly significant)
- `**`: q < 0.01 (very significant)  
- `*`: q < 0.05 (significant)
- No asterisk: q ≥ 0.05 (not significant)

### Fold Enrichment Interpretation
- `> 2.0`: Strong enrichment
- `1.5-2.0`: Moderate enrichment
- `1.2-1.5`: Weak enrichment
- `0.8-1.2`: No enrichment/depletion
- `< 0.8`: Depletion

### CNV Status Definitions
- **Gain**: Copy number gain regions (segment mean > 0, status '+')
- **Loss**: Copy number loss regions (segment mean < 0, status '-')
- **Neutral**: Copy number neutral regions (segment mean ≈ 0, status '0')

## Integration with Other Tools

### Upstream Analysis
- Compatible with standard CNV calling tools (DNAcopy, CBS, etc.)
- Accepts eccDNA calls from Circle-Map, FLED, ECCFinder
- Works with any tool producing standard coordinate files

### Downstream Analysis
- Exports data compatible with R/Bioconductor
- Provides data for custom visualization in Python/R
- Integrates with pathway analysis tools

### Visualization Integration
```r
# Example R integration
library(readr)
cnv_data <- read_csv("CNV_eccdna_density_data.csv")
library(ggplot2)
ggplot(cnv_data, aes(x=cnv_status, y=eccDNA_per_kb)) + 
  geom_violin() + geom_boxplot(width=0.1)
```

## Troubleshooting

### Common Issues
1. **bedtools not found**: Install bedtools using conda or package manager
2. **Memory errors**: Reduce number of cores or filter input data
3. **Empty results**: Check input file formats and coordinate systems
4. **Permission errors**: Ensure write access to output directory

### Error Messages
- **"CNV file not found"**: Check CNV file path and permissions
- **"Missing required columns"**: Verify input file column names
- **"Bedtools error"**: Check bedtools installation and file formats

### Performance Issues
- Large datasets: Consider filtering by chromosome or region
- Slow permutations: Reduce number of permutations or increase cores
- Memory issues: Use fewer cores or split analysis by chromosome

## Advanced Features

### Custom Genome Files
- Support for different genome builds (hg19, hg38, mm10, etc.)
- Custom chromosome sets for specialized analyses
- Integration with UCSC genome browser formats

### Batch Processing
- Process multiple datasets with consistent parameters
- Automated report generation for multiple samples
- Comparative analysis across experimental conditions

### Statistical Customization
- Adjustable permutation parameters
- Custom significance thresholds
- Alternative statistical tests

## Limitations and Notes

- Requires bedtools installation
- Assumes standard chromosome naming (chr1, chr2, etc.)
- CNV file format is strictly defined
- Large datasets may require significant computational resources
- Coordinate system assumptions may need validation

## Citation and References

If you use eccCNV in your research, please cite:
- **Bedtools**: Quinlan AR and Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26(6):841-842.
- **eccCNV**: CNV Region Enrichment Analysis Pipeline v1.0

## Version History

- **v1.0**: Initial release with comprehensive CNV enrichment analysis
- Features include permutation testing, parallel processing, and rich visualization