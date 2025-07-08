# eccDNA-MultiScale.py

## Description

eccDNA-MultiScale is a comprehensive tool for multi-scale window analysis of eccDNA (extrachromosomal circular DNA) data. The tool performs correlation analysis across different genomic window sizes to understand the spatial distribution patterns of eccDNA elements across multiple sample groups and replicates.

## Features

- **Multi-scale window analysis**: Analyzes eccDNA distribution across various genomic window sizes
- **Replicate group comparison**: Supports analysis of three replicate groups
- **Correlation analysis**: Computes both Pearson and Spearman correlations between groups
- **Comprehensive visualization**: Generates heatmaps and trend plots
- **Data filtering**: Automatically filters out windows with no eccDNA across all groups
- **Batch processing**: Processes multiple samples and window sizes efficiently

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
pip install pandas numpy scipy matplotlib seaborn
```

### Required Software
- **bedtools** (for window generation and intersection)
- **Python** 3.7 or higher

### Required Python Packages
- `pandas` - Data manipulation and analysis
- `numpy` - Numerical computing
- `scipy` - Statistical functions and correlation analysis
- `matplotlib` - Plotting and visualization
- `seaborn` - Statistical data visualization
- `argparse` - Command-line argument parsing (standard library)
- `subprocess` - System command execution (standard library)
- `os` - Operating system interface (standard library)

## Usage

### Basic Usage
```bash
python eccDNA-MultiScale.py \
    -1 sample1-rep1.csv sample2-rep1.csv \
    -2 sample1-rep2.csv sample2-rep2.csv \
    -3 sample1-rep3.csv sample2-rep3.csv \
    -f genome.fai \
    -o analysis_output
```

### Advanced Usage
```bash
python eccDNA-MultiScale.py \
    -1 group1_sample1.csv group1_sample2.csv group1_sample3.csv \
    -2 group2_sample1.csv group2_sample2.csv group2_sample3.csv \
    -3 group3_sample1.csv group3_sample2.csv group3_sample3.csv \
    -f hg38.fai \
    -o multiscale_analysis \
    --min_window 5000 \
    --max_window 200000 \
    --window_step 5000
```

## Parameters

### Required Arguments
- `-1`, `--group1`: CSV files for replicate group 1 (sample names will end with -1)
- `-2`, `--group2`: CSV files for replicate group 2 (sample names will end with -2)  
- `-3`, `--group3`: CSV files for replicate group 3 (sample names will end with -3)
- `-f`, `--fai`: FAI file of the reference genome
- `-o`, `--output_prefix`: Output prefix for all generated files

### Optional Arguments
- `--min_window`: Minimum window size in bp (default: 10000)
- `--max_window`: Maximum window size in bp (default: 100000)
- `--window_step`: Window size step in bp (default: 10000)

## Input File Formats

### eccDNA CSV Files
Each CSV file should contain eccDNA coordinates:

**Required columns:**
- `eChr` - Chromosome name (e.g., "chr1", "chr2")
- `eStart` - Start position
- `eEnd` - End position

**Example:**
```csv
eChr,eStart,eEnd,eName
chr1,1000000,1005000,ecc1
chr1,2000000,2003000,ecc2
chr2,500000,510000,ecc3
```

### Genome FAI File
Standard FASTA index file format:

```
chr1	248956422	52	80	81
chr2	242193529	253404903	80	81
chr3	198295559	500657651	80	81
...
```

You can generate this file using:
```bash
samtools faidx genome.fasta
```

## Analysis Workflow

### Step 1: Window Generation
- Uses bedtools makewindows to generate genomic windows of specified sizes
- Windows are created genome-wide based on the FAI file
- Each window size is processed separately

### Step 2: Intersection Analysis
- Converts eccDNA CSV files to BED format
- Uses bedtools intersect to count eccDNA overlaps in each window
- Generates count matrices for each window size

### Step 3: Data Integration
- Merges all samples into unified matrices
- Groups samples by replicate number (based on filename suffix)
- Fills missing values with zeros

### Step 4: Filtering
- Removes windows where all three groups have zero eccDNA counts
- Ensures meaningful correlation analysis

### Step 5: Correlation Analysis
- Calculates group means for each replicate
- Computes Pearson and Spearman correlations between groups
- Performs pairwise comparisons: Group1 vs Group2, Group1 vs Group3, Group2 vs Group3

### Step 6: Visualization and Export
- Creates correlation heatmaps across window sizes
- Generates trend plots showing correlation changes
- Exports all data and results

## Output Files

### Matrix Files
- `{prefix}_window_{size}kb_matrix.csv` - Count matrices for each window size
- Contains genomic coordinates and sample counts

### Visualization Files
- `{prefix}_correlation_heatmaps.png/pdf` - Correlation heatmaps
- `{prefix}_correlation_trends.png/pdf` - Correlation trend plots

### Data Export Files
- `{prefix}_pearson_heatmap_data.csv` - Pearson correlation matrix data
- `{prefix}_spearman_heatmap_data.csv` - Spearman correlation matrix data
- `{prefix}_correlation_trends_data.csv` - Trend analysis data
- `{prefix}_correlation_summary.csv` - Complete correlation statistics

## Sample Organization

### Replicate Groups
The tool expects three replicate groups, identified by filename suffixes:
- **Group 1**: Files ending with `-1`
- **Group 2**: Files ending with `-2`
- **Group 3**: Files ending with `-3`

### Example File Organization
```
experiment/
├── treatment_A_rep1.csv    → becomes treatment_A_rep1-1
├── treatment_A_rep2.csv    → becomes treatment_A_rep2-2
├── treatment_A_rep3.csv    → becomes treatment_A_rep3-3
├── treatment_B_rep1.csv    → becomes treatment_B_rep1-1
├── treatment_B_rep2.csv    → becomes treatment_B_rep2-2
└── treatment_B_rep3.csv    → becomes treatment_B_rep3-3
```

## Statistical Methods

### Correlation Analysis
- **Pearson correlation**: Measures linear relationships
- **Spearman correlation**: Measures monotonic relationships
- Calculated between group means across genomic windows

### Data Processing
- **Zero filtering**: Removes uninformative windows
- **Group averaging**: Reduces noise by averaging replicates
- **Multi-scale analysis**: Reveals scale-dependent patterns

## Example Workflows

### Time Course Analysis
```bash
# Analyze three time points with replicates
python eccDNA-MultiScale.py \
    -1 0h_rep1.csv 0h_rep2.csv 0h_rep3.csv \
    -2 24h_rep1.csv 24h_rep2.csv 24h_rep3.csv \
    -3 48h_rep1.csv 48h_rep2.csv 48h_rep3.csv \
    -f hg38.fai \
    -o timecourse_analysis
```

### Treatment Comparison
```bash
# Compare control vs two treatments
python eccDNA-MultiScale.py \
    -1 control_rep1.csv control_rep2.csv \
    -2 drug_A_rep1.csv drug_A_rep2.csv \
    -3 drug_B_rep1.csv drug_B_rep2.csv \
    -f genome.fai \
    -o treatment_comparison
```

### High-resolution Analysis
```bash
# Fine-scale analysis with small windows
python eccDNA-MultiScale.py \
    -1 sample1.csv sample2.csv \
    -2 sample3.csv sample4.csv \
    -3 sample5.csv sample6.csv \
    -f genome.fai \
    -o high_resolution \
    --min_window 1000 \
    --max_window 50000 \
    --window_step 1000
```

## Interpretation Guidelines

### Correlation Values
- **r > 0.8**: Strong positive correlation
- **0.5 < r < 0.8**: Moderate positive correlation
- **0.2 < r < 0.5**: Weak positive correlation
- **-0.2 < r < 0.2**: No correlation
- **r < -0.2**: Negative correlation

### Window Size Effects
- **Small windows (1-10kb)**: Capture local eccDNA clustering
- **Medium windows (10-50kb)**: Regional distribution patterns
- **Large windows (50-200kb)**: Broad genomic organization

### Visualization Interpretation
- **Heatmaps**: Show correlation patterns across scales
- **Trend plots**: Reveal optimal window sizes for correlation
- **Matrix data**: Enable custom analysis and visualization

## Performance Considerations

### Computational Requirements
- Memory usage scales with number of windows and samples
- Processing time increases with window resolution
- Temporary files are automatically cleaned up

### Optimization Tips
- Use reasonable window size ranges for your analysis
- Consider chromosome-specific analysis for large genomes
- Monitor disk space for large datasets

## Use Cases

### Comparative Genomics
- Compare eccDNA distribution between cell lines
- Analyze treatment effects on eccDNA organization
- Study developmental or disease-related changes

### Scale-dependent Analysis
- Identify optimal resolution for eccDNA analysis
- Understand hierarchical organization of eccDNA
- Optimize experimental design for future studies

### Quality Control
- Assess replicate reproducibility
- Identify batch effects or technical artifacts
- Validate experimental consistency

## Integration with Other Tools

### Upstream Analysis
- Compatible with Circle-Map, FLED, ECCFinder outputs
- Works with any tool producing coordinate files
- Standard CSV format ensures broad compatibility

### Downstream Analysis
```r
# Example R integration
library(readr)
correlation_data <- read_csv("analysis_correlation_summary.csv")
library(ggplot2)
ggplot(correlation_data, aes(x=Window_Size_kb, y=Pearson_R, color=Comparison)) +
  geom_line() + geom_point()
```

## Troubleshooting

### Common Issues
1. **bedtools not found**: Install bedtools using conda or package manager
2. **FAI file errors**: Ensure FAI file matches expected format
3. **Empty results**: Check eccDNA CSV file format and coordinate validity
4. **Memory errors**: Reduce window resolution or analyze subsets

### File Format Issues
- Verify CSV files have required columns (eChr, eStart, eEnd)
- Check coordinate systems are consistent
- Ensure chromosome names match between files and FAI

### Performance Issues
- Reduce window size range for faster processing
- Use appropriate step sizes for your analysis goals
- Monitor temporary file space usage

## Limitations and Notes

- Requires exactly three replicate groups
- Assumes standard chromosome naming conventions
- Window generation depends on bedtools availability
- Large-scale analysis may require significant computational resources
- Chinese language comments in code (functionality is language-independent)

## Advanced Features

### Custom Window Strategies
- Flexible window size parameterization
- Adaptive resolution based on data density
- Integration with genomic features

### Statistical Extensions
- Support for additional correlation methods
- Multiple testing correction options
- Confidence interval calculations

## Future Enhancements

Potential improvements could include:
- Support for variable numbers of groups
- Integration with genomic annotation
- Advanced statistical testing
- Interactive visualization options
- Parallel processing optimization