# degs_enrichment_analysis.py

## Description

degs_enrichment_analysis.py is a script for analyzing the relationship between differentially expressed genes (DEGs) and genes identified in eccDNA enrichment/depletion analyses. The tool performs overlap analysis across different log2 fold-change thresholds to identify genes that are both differentially expressed and associated with eccDNA changes.

## Features

- **Multi-threshold analysis**: Analyzes overlap at different log2FC thresholds (1, 2, 4, 6)
- **Sample-specific analysis**: Processes multiple sample combinations
- **Comprehensive overlap detection**: Identifies four types of gene overlaps
- **Detailed gene lists**: Exports specific gene lists for each overlap category
- **Statistical summaries**: Provides count statistics and detailed gene lists

## Installation Requirements

### Dependencies
```bash
pip install pandas numpy
```

### Required Python Packages
- `pandas` - Data manipulation and analysis
- `numpy` - Numerical computing
- `collections` - Data structures (standard library)

## Usage

### Basic Usage
```bash
python degs_enrichment_analysis.py
```

Note: The script currently has hardcoded file paths. You may need to modify the file paths in the `main()` function to match your input files.

### Input File Modification
Edit the following lines in the script to specify your input files:
```python
degs_file = "your_DEGs_file.txt"
enrichment_file = "your_enrichment_results.csv"
```

## Input File Formats

### DEG File (TSV format)
The script expects a tab-separated file with no header containing differential expression results:

**Columns:**
1. `gene` - Gene symbol/identifier
2. `ensembl` - Ensembl gene ID
3. `tumor_expr` - Tumor expression value
4. `normal_expr` - Normal expression value
5. `log2FC` - Log2 fold change
6. `qvalue` - Adjusted p-value

**Example:**
```tsv
TP53	ENSG00000141510	5.2	3.1	2.5	0.001
MYC	ENSG00000136997	8.9	6.2	-1.8	0.005
BRCA1	ENSG00000012048	3.2	4.1	-0.9	0.02
```

### Enrichment File (CSV format)
CSV file containing eccDNA enrichment/depletion analysis results:

**Required columns:**
- `Sample` - Sample identifier
- `Sample2` - Secondary sample identifier (for sample combinations)
- `gene` - Gene symbol/identifier
- `direction` - Either "enrichment" or "depletion"
- `significant` - Boolean indicating statistical significance

**Example:**
```csv
Sample,Sample2,gene,direction,significant
HeLa_24h,Control,TP53,enrichment,True
HeLa_24h,Control,MYC,depletion,False
HeLa_48h,Control,BRCA1,enrichment,True
```

## Analysis Logic

### Log2FC Thresholds
The script analyzes four predefined thresholds:
- **log2FC ≥ 1**: Up-regulated genes (≥2-fold increase)
- **log2FC ≥ 2**: Up-regulated genes (≥4-fold increase)
- **log2FC ≥ 4**: Up-regulated genes (≥16-fold increase)
- **log2FC ≥ 6**: Up-regulated genes (≥64-fold increase)

Corresponding down-regulated thresholds use negative values.

### Overlap Categories
For each sample and threshold combination, the script identifies:

1. **Up-regulated & Enriched**: Genes that are both up-regulated DEGs and eccDNA-enriched
2. **Up-regulated & Depleted**: Genes that are up-regulated DEGs but eccDNA-depleted
3. **Down-regulated & Enriched**: Genes that are down-regulated DEGs but eccDNA-enriched
4. **Down-regulated & Depleted**: Genes that are both down-regulated DEGs and eccDNA-depleted

## Output

### Console Output
The script provides detailed console output including:

1. **Statistical Summary Table**: Counts for each overlap category across all samples and thresholds
2. **Detailed Gene Lists**: Specific gene names for each category (shown for log2FC ≥ 2 threshold)

### Output File
- **degs_enrichment_overlap_analysis.csv**: Complete results including gene lists for all thresholds

### Output Columns
- `Sample` - Combined sample identifier (Sample_Sample2)
- `log2FC_threshold` - Log2FC threshold used
- `n_up_DEGs` - Number of up-regulated DEGs at this threshold
- `n_down_DEGs` - Number of down-regulated DEGs at this threshold
- `n_enriched` - Number of eccDNA-enriched genes in this sample
- `n_depleted` - Number of eccDNA-depleted genes in this sample
- `n_up_enriched` - Number of up-regulated & enriched genes
- `n_up_depleted` - Number of up-regulated & depleted genes
- `n_down_enriched` - Number of down-regulated & enriched genes
- `n_down_depleted` - Number of down-regulated & depleted genes
- `up_enriched_genes` - List of up-regulated & enriched gene names
- `up_depleted_genes` - List of up-regulated & depleted gene names
- `down_enriched_genes` - List of down-regulated & enriched gene names
- `down_depleted_genes` - List of down-regulated & depleted gene names

## Example Workflow

### 1. Prepare Input Files
```bash
# Ensure you have DEG results from tools like GEPIA2, DESeq2, etc.
# Example DEG file: DEGs_GBM.GEPIA2.Log2FC_1.qValue_0.01.txt

# Ensure you have eccDNA enrichment results
# Example: merged_gene_eccdna_enrichment_results.csv
```

### 2. Modify Script (if needed)
```python
# Edit file paths in main() function
degs_file = "path/to/your/DEGs_file.txt"
enrichment_file = "path/to/your/enrichment_results.csv"
```

### 3. Run Analysis
```bash
python degs_enrichment_analysis.py
```

### 4. Review Results
```bash
# Check console output for summary
# Open degs_enrichment_overlap_analysis.csv for detailed results
```

## Example Output Interpretation

### Console Output Example
```
差异基因与富集/贫化基因关系分析结果
================================================================================

统计汇总:
Sample               log2FC_threshold  n_up_DEGs  n_down_DEGs  n_enriched  n_depleted  n_up_enriched  n_up_depleted  n_down_enriched  n_down_depleted
HeLa_24h_Control     1                 150        120          45          38          12             8              5                15
HeLa_24h_Control     2                 95         85           45          38          8              5              3                12
...

详细基因列表 (log2FC >= 2):
================================================================================

样本: HeLa_24h_Control
上调且富集的基因 (8个): TP53, MYC, EGFR, KRAS, PIK3CA, AKT1, MTOR, CCND1
上调且贫化的基因 (5个): BRCA1, ATM, CHEK2, TP73, MDM2
下调且富集的基因 (3个): CDKN1A, RB1, PTEN
下调且贫化的基因 (12个): APAF1, BAX, BCL2L11, CASP3, CASP7, CASP9, FAS, FASLG, TNF, TRADD, TRAF2, XIAP
```

## Use Cases

### Cancer Research
- Identify oncogenes that are both overexpressed and eccDNA-enriched
- Find tumor suppressors that are downregulated and eccDNA-depleted
- Study therapeutic targets with coordinated expression and eccDNA changes

### Comparative Analysis
- Compare different treatment conditions or time points
- Identify condition-specific gene regulation patterns
- Track dynamic changes in gene expression and eccDNA association

### Biomarker Discovery
- Find genes with strong correlation between expression and eccDNA status
- Identify potential diagnostic or prognostic markers
- Validate findings across multiple datasets

## Statistical Considerations

### Threshold Selection
- Higher thresholds (log2FC ≥ 4, 6) identify more extreme changes
- Lower thresholds (log2FC ≥ 1, 2) capture moderate but potentially significant changes
- Consider biological relevance when interpreting results

### Multiple Testing
- No multiple testing correction is applied in this script
- Consider applying FDR correction for large-scale analyses
- Interpret results in context of experimental design

### Sample Size
- Results depend on the number of significant genes in each category
- Small overlap numbers may not be statistically meaningful
- Consider power analysis for study design

## Limitations and Notes

- File paths are currently hardcoded in the script
- No statistical significance testing for overlaps
- Assumes gene identifiers are consistent between input files
- Chinese language comments and output text
- Limited error handling for file format issues

## Troubleshooting

### Common Issues
1. **File not found**: Check file paths in the `main()` function
2. **Column name errors**: Verify input file formats match expected structure
3. **Empty results**: Check that gene identifiers match between files
4. **Encoding issues**: Ensure files are in UTF-8 encoding

### Data Quality Checks
- Verify gene identifier consistency
- Check for duplicate entries
- Validate threshold ranges
- Confirm significance criteria

## Future Enhancements

Potential improvements could include:
- Command-line argument parsing
- Statistical significance testing for overlaps
- Visualization of results
- Support for different file formats
- Multiple testing correction options
- Batch processing capabilities