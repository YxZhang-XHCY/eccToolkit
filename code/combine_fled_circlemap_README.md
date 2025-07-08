# combine_fled_circlemap.py

## Description

combine_fled_circlemap.py is a comprehensive tool for processing and combining eccDNA (extrachromosomal circular DNA) detection results from two popular analysis tools: FLED and Circle-Map. The script standardizes, filters, and merges outputs from both tools to create unified datasets for comparative analysis.

## Features

- **Dual tool support**: Processes both FLED and Circle-Map outputs
- **Automatic filtering**: Applies quality filters specific to each tool
- **Batch processing**: Handles multiple samples simultaneously
- **Data standardization**: Standardizes column formats and naming conventions
- **Combined output**: Creates unified datasets for comparative analysis
- **Comprehensive logging**: Detailed progress reporting and error handling

## Installation Requirements

### Dependencies
```bash
pip install pandas
```

### Required Python Packages
- `pandas` - Data manipulation and analysis
- `argparse` - Command-line argument parsing (standard library)
- `logging` - Logging functionality (standard library)
- `os` - Operating system interface (standard library)
- `glob` - File pattern matching (standard library)
- `typing` - Type hints (standard library)

## Usage

### Basic Usage
```bash
python combine_fled_circlemap.py \
    --fled-dir /path/to/fled/results \
    --circlemap-dir /path/to/circlemap/results \
    --output-dir /path/to/output
```

### Advanced Usage
```bash
python combine_fled_circlemap.py \
    --fled-dir FLED_results/ \
    --circlemap-dir CircleMap_results/ \
    --fled-pattern "*.DiGraph.OnesegJunction.out" \
    --circlemap-pattern "*_CM.bed" \
    --output-dir combined_results/
```

## Parameters

### Required Parameters
- `--fled-dir`: Directory containing FLED output files
- `--circlemap-dir`: Directory containing Circle-Map output files
- `--output-dir`: Directory to save the final result files

### Optional Parameters
- `--fled-pattern`: Filename pattern for FLED files (default: "*.DiGraph.OnesegJunction.out")
- `--circlemap-pattern`: Filename pattern for Circle-Map files (default: "*_CM.bed")

## Input File Formats

### FLED Output Files
Expected format: `*.DiGraph.OnesegJunction.out`

**Columns (tab-separated, no header):**
1. `chrom` - Chromosome name
2. `start` - Start position
3. `end` - End position
4. `strand` - Strand information
5. `label` - Label/identifier
6. `Tag` - Tag information
7. `Nfullpass` - Number of full pass reads
8. `Nbreakpoint` - Number of breakpoints
9. `L_Pvalue` - Left p-value
10. `R_Pvalue` - Right p-value
11. `L_covRatio` - Left coverage ratio
12. `R_covRatio` - Right coverage ratio
13. `readID` - Read identifiers

### Circle-Map Output Files
Expected format: `*_CM.bed`

**Columns (tab-separated, no header):**
1. `Chromosome` - Chromosome name
2. `Start` - Start position
3. `End` - End position
4. `Discordants` - Number of discordant reads
5. `Split_reads` - Number of split reads
6. `Circle_score` - Circle score
7. `Mean_coverage` - Mean coverage
8. `Standard_deviation` - Coverage standard deviation
9. `Coverage_increase_start` - Coverage increase at start
10. `Coverage_increase_end` - Coverage increase at end
11. `Coverage_continuity` - Coverage continuity

## Quality Filters

### FLED Filtering Criteria
- `Nfullpass > 0` - Excludes records with no full pass reads
- Cleans and standardizes readID format
- Rounds floating-point values to 2 decimal places

### Circle-Map Filtering Criteria
- `Circle_score > 50` - Minimum circle score threshold
- `Split_reads > 2` - Minimum split reads
- `Discordants > 2` - Minimum discordant reads
- `Coverage_increase_start > 0.33` - Minimum coverage increase at start
- `Coverage_increase_end > 0.33` - Minimum coverage increase at end
- `Length < 1e7` - Maximum length filter (10 Mb)

## Output Files

### Individual Tool Results
1. **FLED_filtered_results.csv** - Filtered FLED results with additional columns:
   - `name` - Region identifier (chrom-start-end)
   - `sample` - Sample name extracted from filename
   - `Length` - Calculated region length
   - All original FLED columns with cleaned data

2. **CircleMap_filtered_results.csv** - Filtered Circle-Map results with additional columns:
   - `name` - Region identifier (Chromosome-Start-End)
   - `sample` - Sample name extracted from filename
   - `Length` - Calculated region length
   - All original Circle-Map columns with rounded values

### Combined Results
3. **FLED_and_CircleMap_combined.csv** - Unified dataset containing:
   - `name` - Region identifier
   - `sample` - Sample name
   - `Length` - Region length
   - `Chromosome` - Chromosome name
   - `Start` - Start position
   - `End` - End position
   - `Source` - Tool source (FLED or Circle-Map)

## Sample Name Extraction

### FLED Files
- Removes suffix: `.DiGraph.OnesegJunction.out`
- Example: `Sample1.DiGraph.OnesegJunction.out` → `Sample1`

### Circle-Map Files
- Removes suffix: `_CM.bed`
- Example: `Sample1_CM.bed` → `Sample1`

## Data Processing Pipeline

### Step 1: File Discovery
- Scans directories for files matching specified patterns
- Reports found files and any missing data

### Step 2: Individual Processing
- Processes each file according to tool-specific logic
- Applies quality filters
- Standardizes data formats
- Extracts sample names

### Step 3: Data Combination
- Merges results from both tools
- Standardizes column names
- Sorts by chromosome and position
- Creates unified dataset

### Step 4: Output Generation
- Saves individual tool results
- Creates combined dataset
- Generates processing statistics

## Example Workflow

### Directory Structure
```
project/
├── FLED_results/
│   ├── Sample1.DiGraph.OnesegJunction.out
│   ├── Sample2.DiGraph.OnesegJunction.out
│   └── Sample3.DiGraph.OnesegJunction.out
├── CircleMap_results/
│   ├── Sample1_CM.bed
│   ├── Sample2_CM.bed
│   └── Sample3_CM.bed
└── combined_results/
    ├── FLED_filtered_results.csv
    ├── CircleMap_filtered_results.csv
    └── FLED_and_CircleMap_combined.csv
```

### Command Execution
```bash
python combine_fled_circlemap.py \
    --fled-dir FLED_results/ \
    --circlemap-dir CircleMap_results/ \
    --output-dir combined_results/
```

### Expected Output
```
[INFO] 2024-01-01 10:00:00 - Processing FLED file: FLED_results/Sample1.DiGraph.OnesegJunction.out
[INFO] 2024-01-01 10:00:01 -   Found 150 valid records in sample Sample1
[INFO] 2024-01-01 10:00:02 - Processing Circle-Map file: CircleMap_results/Sample1_CM.bed
[INFO] 2024-01-01 10:00:03 -   Found 75 valid records in sample Sample1
[INFO] 2024-01-01 10:00:04 - FLED results saved. Total records: 450
[INFO] 2024-01-01 10:00:05 - Circle-Map results saved. Total records: 225
[INFO] 2024-01-01 10:00:06 - Combined results saved. Total records: 675
[INFO] 2024-01-01 10:00:07 - ✔ Processing complete! All output files are saved in: combined_results/
```

## Performance Considerations

### Memory Usage
- Loads all files into memory simultaneously
- Memory usage scales with total dataset size
- Consider processing in batches for very large datasets

### Processing Speed
- Vectorized operations for improved performance
- Parallel processing could be added for very large datasets
- I/O operations are typically the bottleneck

## Error Handling

### File Processing
- Graceful handling of missing or malformed files
- Continues processing even if some files fail
- Detailed error logging for debugging

### Data Validation
- Checks for empty datasets
- Validates required columns
- Handles edge cases in data formatting

## Use Cases

### Comparative Analysis
- Compare eccDNA detection between FLED and Circle-Map
- Identify consensus regions detected by both tools
- Assess tool-specific sensitivity and specificity

### Quality Assessment
- Evaluate filtering effectiveness
- Compare detection statistics across samples
- Identify high-confidence eccDNA regions

### Downstream Analysis
- Prepare standardized datasets for visualization
- Create input for statistical analysis
- Generate data for publication figures

## Integration with Other Tools

### Visualization
- Compatible with genomic visualization tools
- Can be imported into R/Bioconductor for plotting
- Suitable for Circos plot generation

### Statistical Analysis
- Ready for statistical testing
- Compatible with machine learning workflows
- Suitable for comparative genomics studies

## Limitations and Notes

- Assumes specific file naming conventions
- Quality filters are hardcoded (tool-specific)
- No coordinate validation or genome build checking
- Limited to standard tool outputs
- Chinese comments in code (functionality is language-independent)

## Troubleshooting

### Common Issues
1. **No files found**: Check directory paths and file patterns
2. **Empty output**: Verify input file formats and quality filters
3. **Memory errors**: Process smaller batches or increase available memory
4. **Sample name extraction**: Ensure consistent file naming conventions

### Validation Steps
1. Check input file formats match expected specifications
2. Verify quality filter criteria are appropriate
3. Examine log output for processing statistics
4. Validate output file completeness and format