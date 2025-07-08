# CircExpress-Integrator.py

## Description

CircExpress-Integrator is a comprehensive tool for analyzing associations between eccDNA (extrachromosomal circular DNA) enrichment/depletion patterns and RNA-seq differential expression genes (DEGs). The tool performs per-sample association analysis using fold-change (FC) stratified gradient curves and generates detailed statistical reports with visualizations.

## Features

- **Multi-sample Analysis**: Processes multiple samples simultaneously
- **FC-stratified Analysis**: Analyzes associations across different fold-change thresholds
- **Statistical Testing**: Performs Fisher's exact tests with FDR correction
- **Comprehensive Visualization**: Generates gradient curves, heatmaps, and summary plots
- **Gene List Export**: Exports detailed gene lists for each analysis category
- **Enhanced Logging**: Provides detailed logging with timestamps and progress tracking

## Installation Requirements

### Dependencies
```bash
pip install pandas numpy matplotlib seaborn scipy statsmodels
```

### Required Python Packages
- `pandas` - Data manipulation and analysis
- `numpy` - Numerical computing
- `matplotlib` - Plotting and visualization
- `seaborn` - Statistical data visualization
- `scipy` - Statistical functions (Fisher's exact test)
- `statsmodels` - Statistical modeling (FDR correction)

## Usage

### Basic Usage
```bash
python CircExpress-Integrator.py --eccdna INPUT_ECCDNA.csv --degs INPUT_DEGS.tsv --outdir OUTPUT_DIR
```

### Advanced Usage
```bash
python CircExpress-Integrator.py \
    --eccdna merged_eccdna_results.csv \
    --degs deg_table.tsv \
    --outdir eccdna_deg_analysis \
    --fc "1,2,4,6,8" \
    --verbose
```

## Parameters

### Required Parameters
- `--eccdna`: Path to merged eccDNA CSV file containing enrichment/depletion results
- `--degs`: Path to DEG table in TSV format (e.g., from GEPIA2)

### Optional Parameters
- `--outdir`: Output directory for results (default: "eccdna_deg_results")
- `--fc`: Comma-separated log2FC thresholds for analysis (default: "1,2,4,6")
- `--verbose`, `-v`: Enable verbose logging for detailed progress tracking

## Input File Formats

### eccDNA File (CSV)
Required columns:
- `Sample`: Sample identifier
- `gene`: Gene symbol/identifier
- `direction`: Either "enrichment" or "depletion"
- `significant`: Boolean indicating statistical significance

Example:
```csv
Sample,gene,direction,significant
Sample1,TP53,enrichment,True
Sample1,MYC,depletion,False
Sample2,BRCA1,enrichment,True
```

### DEG File (TSV)
Tab-separated file with columns (no header):
1. `gene`: Gene symbol/identifier
2. `ensembl_id`: Ensembl gene ID
3. `tumor_exp`: Tumor expression value
4. `normal_exp`: Normal expression value
5. `log2fc`: Log2 fold change
6. `qvalue`: Adjusted p-value

Example:
```tsv
TP53	ENSG00000141510	5.2	3.1	0.74	0.001
MYC	ENSG00000136997	8.9	6.2	0.52	0.005
```

## Output Files

### Directory Structure
```
OUTPUT_DIR/
├── analysis_YYYYMMDD_HHMMSS.log
├── all_samples_gradient_analysis.csv
├── summary_heatmaps.png
├── overlap_summary_fc2.png
├── sample_SAMPLE1/
│   ├── SAMPLE1_gradient_results.csv
│   ├── SAMPLE1_gradient_summary.csv
│   ├── SAMPLE1_gradient_analysis.png
│   └── gene_lists/
│       ├── SAMPLE1_fc2_gene_lists.json
│       ├── SAMPLE1_fc2_enriched_up.txt
│       └── SAMPLE1_fc2_depleted_down.txt
```

### Output Descriptions

#### Main Results
- `all_samples_gradient_analysis.csv`: Combined results for all samples and FC thresholds
- `analysis_*.log`: Detailed log file with timestamps and analysis progress

#### Per-Sample Results
- `*_gradient_results.csv`: Detailed results including gene lists for each FC threshold
- `*_gradient_summary.csv`: Summary statistics without gene lists
- `*_gradient_analysis.png`: Four-panel visualization showing:
  - Overlap counts across FC thresholds
  - Odds ratios (log scale)
  - -log10 p-values
  - DEG counts

#### Visualizations
- `summary_heatmaps.png`: Four heatmaps showing odds ratios and p-values across samples
- `overlap_summary_fc2.png`: Bar plots of gene overlaps at FC threshold = 2

#### Gene Lists
- `*_gene_lists.json`: Complete gene lists in JSON format
- `*_enriched_up.txt`: Genes that are both eccDNA-enriched and up-regulated
- `*_depleted_down.txt`: Genes that are both eccDNA-depleted and down-regulated

## Statistical Methods

### Fisher's Exact Test
- One-tailed test (alternative='greater') for testing enrichment
- Tests associations between eccDNA enrichment/depletion and DEG up/down-regulation
- Calculates odds ratios and p-values

### Multiple Testing Correction
- FDR correction using Benjamini-Hochberg method
- Applied per sample across FC thresholds

### Association Types Tested
1. **Enriched-Up**: eccDNA enriched genes ∩ Up-regulated DEGs
2. **Depleted-Down**: eccDNA depleted genes ∩ Down-regulated DEGs

## Example Workflow

1. **Prepare Input Data**:
   - eccDNA enrichment/depletion results from analysis tools
   - DEG table from differential expression analysis

2. **Run Analysis**:
   ```bash
   python CircExpress-Integrator.py --eccdna data.csv --degs degs.tsv --outdir results
   ```

3. **Review Results**:
   - Check log file for analysis summary
   - Examine gradient plots for individual samples
   - Review summary heatmaps for cross-sample patterns
   - Extract significant gene lists for further analysis

## Limitations and Notes

- Requires gene symbols to be consistent between eccDNA and DEG datasets
- Fisher's exact test assumes independence of gene expression measurements
- Multiple testing correction is applied per sample, not globally
- Visualization quality depends on the number of samples and data distribution
- Large datasets may require significant memory for processing

## Output Interpretation

### Significant Associations
- **p < 0.05**: Statistically significant association
- **Odds Ratio > 1**: Positive association (more overlap than expected)
- **Odds Ratio < 1**: Negative association (less overlap than expected)

### Gradient Analysis
- **Increasing FC thresholds**: Generally fewer DEGs but potentially stronger associations
- **Stable patterns**: Consistent associations across FC thresholds suggest robust relationships
- **Sample-specific patterns**: Variations between samples may indicate biological heterogeneity