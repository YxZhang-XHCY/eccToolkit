# GenoLink_v2.py

## Description

GenoLink_v2 is a tool for generating genomic link relationship data suitable for visualization tools like Circos. It processes eccDNA (extrachromosomal circular DNA) data to create pairwise connections between genomic regions within the same element groups (eName). The tool supports chromosome-specific filtering and provides customizable color schemes for visualization.

## Features

- **Chromosome-specific filtering**: Filter for groups where all elements are on the same chromosome
- **Color schemes**: Two visualization modes - highlight largest groups or random colors
- **Parallel processing**: Multi-core processing for improved performance
- **Data validation**: Comprehensive input data validation
- **Statistical reporting**: Detailed statistics about groups and output
- **Flexible output**: TSV format compatible with Circos and other visualization tools

## Installation Requirements

### Dependencies
```bash
pip install pandas numpy tqdm
```

### Required Python Packages
- `pandas` - Data manipulation and analysis
- `numpy` - Numerical computing
- `tqdm` - Progress bar display
- `multiprocessing` - Built-in parallel processing (standard library)
- `pathlib` - Path handling (standard library)
- `argparse` - Command-line argument parsing (standard library)
- `logging` - Logging functionality (standard library)

## Usage

### Basic Usage
```bash
python GenoLink_v2.py input_file.csv
```

### Advanced Usage
```bash
python GenoLink_v2.py input_file.csv \
    --output output_links.tsv \
    --target-chr chr1 \
    --color-scheme highlight_largest \
    --processes 8 \
    --seed 42
```

### No Chromosome Filtering
```bash
python GenoLink_v2.py input_file.csv \
    --target-chr "" \
    --color-scheme random
```

## Parameters

### Required Parameters
- `input`: Path to input CSV file containing genomic region data

### Optional Parameters
- `--output`, `-o`: Output TSV file path (auto-generated if not specified)
- `--processes`, `-p`: Number of processes to use (default: CPU core count)
- `--target-chr`: Target chromosome for filtering (default: "chr1", empty string for no filtering)
- `--color-scheme`: Color scheme for visualization
  - `highlight_largest`: Highlight the largest group with bright colors, others in gray
  - `random`: Random colors for all groups
- `--seed`: Random seed for reproducible color generation

## Input File Format

### Required CSV Columns
- `eName`: Element name/identifier (groups with same eName will be linked)
- `eChr`: Chromosome name (e.g., "chr1", "chr2")
- `eStart`: Start position (integer)
- `eEnd`: End position (integer)

### Example Input
```csv
eName,eChr,eStart,eEnd
element1,chr1,1000,2000
element1,chr1,5000,6000
element1,chr1,10000,11000
element2,chr1,15000,16000
element2,chr1,20000,21000
element3,chr2,3000,4000
element3,chr2,8000,9000
```

## Output File Format

### TSV Format
The output file contains pairwise connections between regions within the same eName group:

```tsv
Chr1    Start1    End1    Chr2    Start2    End2    RGB
chr1    1000      2000    chr1    5000      6000    255,0,0
chr1    1000      2000    chr1    10000     11000   255,0,0
chr1    5000      6000    chr1    10000     11000   255,0,0
chr1    15000     16000   chr1    20000     21000   128,128,128
```

### Column Descriptions
- `Chr1`, `Start1`, `End1`: First genomic region
- `Chr2`, `Start2`, `End2`: Second genomic region
- `RGB`: Color specification in R,G,B format

## Processing Logic

### Chromosome Filtering
When `--target-chr` is specified (default: "chr1"):
1. Groups elements by `eName`
2. Only retains groups where ALL elements are on the target chromosome
3. Filters out mixed-chromosome groups

### Link Generation
1. Within each eName group, generates all pairwise combinations
2. Skips groups with only one element (no links possible)
3. Assigns colors based on the selected color scheme

### Color Schemes

#### highlight_largest
- Identifies the largest group (most elements)
- Assigns bright colors to the largest group
- Uses gray (128,128,128) for all other groups
- Helps highlight the most significant genomic element

#### random
- Assigns random RGB colors to each group
- Provides visual distinction between all groups
- Useful for exploring multiple groups simultaneously

## Example Workflows

### 1. Chromosome-specific Analysis
```bash
# Focus on chr1 elements only
python GenoLink_v2.py eccdna_data.csv \
    --target-chr chr1 \
    --color-scheme highlight_largest \
    --output chr1_links.tsv
```

### 2. Genome-wide Analysis
```bash
# Process all chromosomes
python GenoLink_v2.py eccdna_data.csv \
    --target-chr "" \
    --color-scheme random \
    --output genome_wide_links.tsv
```

### 3. High-performance Processing
```bash
# Use maximum CPU cores
python GenoLink_v2.py large_dataset.csv \
    --processes 16 \
    --seed 123
```

## Output Statistics

The tool provides comprehensive statistics:

### Input Statistics
- Total number of rows
- Number of unique eName groups
- Chromosome distribution
- Group size distribution

### Processing Statistics
- Number of groups after filtering
- Number of valid groups (>1 element)
- Top 5 largest groups with sizes
- Identification of the largest group

### Output Statistics
- Total number of links generated
- Number of unique colors used
- Color usage distribution
- File size and location

## Performance Considerations

### Memory Usage
- Processes data in memory using pandas
- Memory usage scales with input file size
- Large datasets may require sufficient RAM

### Processing Speed
- Parallel processing significantly improves performance
- Processing time depends on:
  - Number of groups
  - Group sizes (combinatorial explosion for large groups)
  - Number of CPU cores available

### Scalability
- Efficient for datasets with moderate group sizes
- Very large groups (>1000 elements) may require significant processing time
- Consider splitting large datasets by chromosome

## Data Validation

### Input Validation
- Checks for required columns
- Validates numeric data types for coordinates
- Warns about invalid coordinates (eStart >= eEnd)
- Reports missing or malformed data

### Processing Validation
- Handles empty results gracefully
- Reports groups with insufficient elements
- Provides warnings for edge cases

## Integration with Visualization Tools

### Circos Integration
The output TSV format is compatible with Circos link plots:
1. Use Chr1, Start1, End1 and Chr2, Start2, End2 as link coordinates
2. Use RGB column for color specification
3. Can be directly loaded into Circos configuration

### Other Visualization Tools
The TSV format can be adapted for:
- IGV (Integrative Genomics Viewer)
- Custom plotting with matplotlib/plotly
- Genomic visualization web tools

## Limitations and Notes

- Only processes elements within the same eName group
- Chromosome filtering is strict (excludes mixed-chromosome groups)
- Large groups generate quadratic numbers of links
- Color assignment is group-based, not link-based
- Requires sufficient memory for large datasets
- Chinese language comments in the code (functionality is language-independent)

## Error Handling

- Validates input file existence and format
- Handles empty or invalid data gracefully
- Provides informative error messages
- Logs processing progress and statistics
- Continues processing even with warnings

## Output File Naming

If no output file is specified, the tool automatically generates names:
- With chromosome filtering: `input.chr1.highlight_largest.tsv`
- Without filtering: `input.random.tsv`
- Reflects the processing parameters used