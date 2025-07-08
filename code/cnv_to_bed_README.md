# cnv_to_bed.py

## Description

cnv_to_bed.py is a utility script for converting CNV (Copy Number Variation) data from a single TSV file into separate BED files for each cell line. The script handles coordinate system conversion from 1-based to 0-based coordinates as required by the BED format specification.

## Features

- **Coordinate conversion**: Converts 1-based coordinates to 0-based BED format
- **Multi-cell line support**: Splits data by cell line into separate files
- **BED format compliance**: Generates proper BED format with standard coordinate system
- **Batch processing**: Processes all cell lines in a single run
- **Automatic output directory creation**: Creates output directory if it doesn't exist

## Installation Requirements

### Dependencies
```bash
pip install pandas
```

### Required Python Packages
- `pandas` - Data manipulation and CSV/TSV file handling
- `argparse` - Command-line argument parsing (standard library)
- `os` - Operating system interface (standard library)

## Usage

### Basic Usage
```bash
python cnv_to_bed.py -i input_cnv.tsv -o output_directory
```

### Example Usage
```bash
python cnv_to_bed.py \
    --input CNV_data.tsv \
    --outdir cnv_bed_files/
```

## Parameters

### Required Parameters
- `-i`, `--input`: Input CNV TSV file path
- `-o`, `--outdir`: Output directory for BED files

## Input File Format

### CNV TSV Format
The input file should be a tab-separated values file with the following columns:

- `Cell_Lines`: Cell line identifier (used for grouping and output file naming)
- `Chromosome`: Chromosome name (e.g., "chr1", "chr2")
- `Start`: Start position (1-based coordinate)
- `End`: End position (1-based coordinate)
- `SegmentMean`: Segment mean value
- `NumProbes`: Number of probes in the segment
- `Status`: CNV status (e.g., "gain", "loss", "normal")
- `ModelID`: Model identifier
- `ProfileID`: Profile identifier

### Example Input
```tsv
Cell_Lines	Chromosome	Start	End	SegmentMean	NumProbes	Status	ModelID	ProfileID
HeLa	chr1	1000000	2000000	0.5	150	loss	model1	profile1
HeLa	chr1	3000000	4000000	-0.3	200	loss	model1	profile1
A549	chr1	1500000	2500000	0.8	180	gain	model2	profile2
A549	chr2	5000000	6000000	0.2	120	normal	model2	profile2
```

## Output File Format

### BED Format
The script generates separate BED files for each cell line with the following format:

```
chromosome	start	end	SegmentMean	NumProbes	Status	ModelID	ProfileID
```

### Output Files
- Files are named using the pattern: `{Cell_Lines}.cnv.bed`
- Each file contains only the data for that specific cell line
- Coordinates are converted from 1-based to 0-based (BED standard)

### Example Output
For the input example above, two files would be created:

**HeLa.cnv.bed**:
```
chr1	999999	2000000	0.5	150	loss	model1	profile1
chr1	2999999	4000000	-0.3	200	loss	model1	profile1
```

**A549.cnv.bed**:
```
chr1	1499999	2500000	0.8	180	gain	model2	profile2
chr2	4999999	6000000	0.2	120	normal	model2	profile2
```

## Coordinate System Conversion

### 1-based to 0-based Conversion
- **Input (1-based)**: Start and end coordinates as provided in the CNV file
- **Output (0-based)**: Start coordinate is decremented by 1 (Start0 = Start - 1)
- **End coordinate**: Remains unchanged (BED format is 0-based start, 1-based end)

### BED Format Specification
- Column 1: Chromosome name
- Column 2: Start position (0-based, inclusive)
- Column 3: End position (1-based, exclusive)
- Columns 4-8: Additional CNV information

## Processing Logic

1. **Load Data**: Read the input TSV file using pandas
2. **Coordinate Conversion**: Create Start0 column with 0-based coordinates
3. **Group by Cell Line**: Split data by Cell_Lines column
4. **Generate BED Files**: Write each group to a separate BED file
5. **Output**: Save files with appropriate naming convention

## Use Cases

### CNV Analysis Workflows
- Preparing CNV data for visualization in genome browsers
- Converting data for input to other bioinformatics tools
- Splitting multi-cell line datasets for parallel processing

### Integration with Other Tools
- Input for BEDTools operations
- Visualization in IGV (Integrative Genomics Viewer)
- Input for custom genomic analysis scripts

## Example Workflow

1. **Obtain CNV Data**: Get CNV results from analysis tools like DNAcopy, CBS, or commercial platforms
2. **Format Data**: Ensure data includes required columns
3. **Convert to BED**: Use cnv_to_bed.py to generate per-cell line BED files
4. **Downstream Analysis**: Use BED files for visualization or further analysis

```bash
# Step 1: Convert CNV data
python cnv_to_bed.py -i CNV_results.tsv -o cnv_bed_output/

# Step 2: List generated files
ls cnv_bed_output/
# Output: HeLa.cnv.bed  A549.cnv.bed  U87.cnv.bed  ...

# Step 3: Use with other tools
bedtools intersect -a cnv_bed_output/HeLa.cnv.bed -b genes.bed > HeLa_cnv_genes.bed
```

## Performance Considerations

### Memory Usage
- Loads entire TSV file into memory
- Memory usage scales with input file size
- Suitable for typical CNV datasets (millions of segments)

### Processing Speed
- Fast processing using pandas vectorized operations
- Processing time depends on file size and number of cell lines
- I/O operations are the main bottleneck

## Limitations and Notes

- Requires all necessary columns to be present in the input file
- Assumes standard CNV file format with specific column names
- Does not validate coordinate ranges or chromosome names
- Output files overwrite existing files with the same names
- No error handling for malformed input data

## Error Handling

The script has minimal error handling:
- Pandas will raise errors for missing columns
- File I/O errors will be raised if permissions are insufficient
- Directory creation is handled with `exist_ok=True`

## Output Validation

To validate the output:
1. Check that all expected cell line files are created
2. Verify coordinate conversion (start coordinates should be 1 less than input)
3. Confirm that the number of lines matches the input for each cell line
4. Validate BED format compliance using external tools

## Integration Examples

### With BEDTools
```bash
# Find CNV overlaps with genes
bedtools intersect -a HeLa.cnv.bed -b genes.bed -wa -wb > HeLa_cnv_gene_overlaps.bed

# Merge adjacent CNV segments
bedtools merge -i HeLa.cnv.bed -d 1000 > HeLa_cnv_merged.bed
```

### With R/Bioconductor
```r
# Read BED file in R
cnv_data <- read.table("HeLa.cnv.bed", sep="\t", 
                       col.names=c("chr", "start", "end", "segmean", "nprobes", "status", "modelid", "profileid"))

# Convert to GRanges object
library(GenomicRanges)
cnv_gr <- makeGRangesFromDataFrame(cnv_data, keep.extra.columns=TRUE)
```