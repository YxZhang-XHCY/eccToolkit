# README for eccToolkit

Welcome to the `eccToolkit` repository. This repository contains the analytical code employed in the research paper titled "Dynamics of Extrachromosomal Circular DNA in Rice."

## Overview

eccToolkit is a specialized software suite designed for the analysis and classification of eccDNA (extrachromosomal circular DNA) elements. It comprises two main scripts: `eccDNA_Repeat_Classifier.py` and `eccDNA_Detected_eccDNA.py`, which together facilitate the detection and categorization of eccDNA sequences from genomic data.

## Features

- **eccDNA Detection**: Identifies eccDNA elements within genomic data.
- **Repeat Classification**: Classifies detected eccDNA elements based on their repeat characteristics.
- **Comprehensive Analysis**: Integrates various data processing steps to provide detailed insights into eccDNA profiles.

## Requirements

- Python environment with necessary dependencies installed (e.g., Pandas, NumPy).
- Genomic data files in appropriate formats (e.g., BED, FASTA).

## Installation

### Using Mamba (Recommended)

1. Install Mamba via Miniconda or Anaconda, if not already installed.

2. Create and activate the eccToolkit environment:

   ```shell
   mamba env create -f eccToolkit.yml
   conda activate eccToolkit
   ```

### Using Conda

1. Alternatively, use Conda for installation:

   ```shell
   conda env create -f eccToolkit.yml
   conda activate eccToolkit
   ```

## Script Usage and Parameters

### eccDNA_Detected_eccDNA.py

This script is responsible for the initial detection of eccDNA elements from genomic data.

#### Parameters:

- `-i`, `--input` (Required): Path to the input genomic data file.
- `-o`, `--output` (Required): Path for the output file where detected eccDNA elements will be listed.
- Additional parameters specific to the detection algorithm (e.g., sensitivity settings).

#### Example Usage:

```shell
python eccDNA_Detected_eccDNA.py -i input_genomic_data.bed -o detected_eccDNA.csv
```

#### Input/Output:

- **Input**: BED format genomic data.
- **Output**: CSV file listing detected eccDNA elements.

### eccDNA_Repeat_Classifier.py

This script classifies the detected eccDNA elements based on their repeat characteristics.

#### Parameters:

- `-i`, `--input` (Required): Path to the input file containing detected eccDNA elements.
- `-o`, `--output` (Required): Path for the output file where classified eccDNA elements will be listed.
- Additional parameters for classification criteria (e.g., repeat length, type).

#### Example Usage:

```shell
python eccDNA_Repeat_Classifier.py -i detected_eccDNA.csv -o classified_eccDNA.csv
```

#### Input/Output:

- **Input**: CSV file from `eccDNA_Detected_eccDNA.py`.
- **Output**: CSV file with eccDNA elements classified by repeat type.

## Comprehensive Workflow

To fully utilize eccToolkit, follow these steps:

1. Run `eccDNA_Detected_eccDNA.py` to detect eccDNA elements.
2. Use the output of the first script as the input for `eccDNA_Repeat_Classifier.py` to classify the detected elements.

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE). Kindly refer to the `LICENSE` file for detailed terms and conditions.
