# eccToolkit Cleanup and Renaming Plan

## Overview
This document outlines the recommended cleanup and renaming strategy for the eccToolkit codebase to reduce redundancy and improve maintainability.

## Conda Environment Setup

### Installation Instructions
```bash
# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate ecctoolkit

# Manual installations (if needed)
# Circle-Map: https://github.com/iprada/Circle-Map
# CircleSeeker: https://github.com/TreesLab/CircleSeeker
# FLED: https://github.com/FangLabTools/FLED
```

## Script Cleanup Strategy

### 1. Files to Remove (8 files - redundant/superseded)

#### Remove these files as they are superseded by newer versions:
```bash
# Version 1 files superseded by v2
rm code/te_motif_composition_analysis.py        # superseded by v2
rm code/GenoLink.py                             # superseded by v2

# Older eccDNA analysis versions
rm code/eccdna_analysis.py                      # superseded by eccdna_analyzer.py
rm code/Com_eccdna_analysis.py                  # duplicate of 222222.py

# Basic hotspot detection versions
rm code/eccdna_hotspot_detector.py              # superseded by parallel version
rm code/eccdna_hotspot_permtest.py              # superseded by parallel version

# TAD analysis duplicate
rm code/eccdna_tad_analysis.py                  # superseded by eccTAD.py

# Merge utility duplicate
rm code/merge_csvs.py                           # superseded by merge_csv_with_sample.py
```

### 2. Files to Rename (4 files - test/development files)

#### Rename poorly named files to descriptive names:
```bash
# Test files → Descriptive names
mv code/test2.py code/eccdna_deg_correlation_analysis_comprehensive.py
mv code/test3.py code/eccdna_deg_correlation_analysis_simple.py

# Numbered files → Descriptive names  
mv code/22222.py code/te_composition_ratio_analysis.py
mv code/222222.py code/eccdna_composite_te_analysis.py
```

### 3. Canonical Scripts to Keep (Primary Analysis Tools)

#### Core eccDNA Analysis
- `eccdna_analyzer.py` - **Primary eccDNA analyzer (v3.0)**
- `eccDNA-SHARP.py` - Hotspot analysis with refined precision
- `eccDNA-MultiScale.py` - Multi-scale bedtools-based analysis

#### Hotspot Detection
- `eccdna_hotspot_permtest_parallel.py` - **Primary hotspot detector**
- `eccdna_hotspot_refined.py` - Specialized boundary detection

#### TAD Analysis
- `eccTAD.py` - **Primary TAD boundary enrichment analysis**

#### Transposable Element Analysis
- `te_motif_composition_analysis_v2.py` - **Primary TE motif analysis**
- `te_percentage_distribution_analysis.py`
- `te_anno_percent_by_class_analysis.py`

#### Genomic Linkage
- `GenoLink_v2.py` - **Primary genomic linkage visualization**

#### Data Processing
- `merge_csv_with_sample.py` - **Primary CSV merge utility**

### 4. Execution Plan

#### Phase 1: Backup and Prepare
```bash
# Create backup of current code directory
cp -r code/ code_backup_$(date +%Y%m%d)/

# Create deprecated directory for removed files
mkdir code/deprecated/
```

#### Phase 2: Move Deprecated Files
```bash
# Move superseded files to deprecated directory
mv code/te_motif_composition_analysis.py code/deprecated/
mv code/GenoLink.py code/deprecated/
mv code/eccdna_analysis.py code/deprecated/
mv code/Com_eccdna_analysis.py code/deprecated/
mv code/eccdna_hotspot_detector.py code/deprecated/
mv code/eccdna_hotspot_permtest.py code/deprecated/
mv code/eccdna_tad_analysis.py code/deprecated/
mv code/merge_csvs.py code/deprecated/
```

#### Phase 3: Rename Files
```bash
# Rename test and numbered files
mv code/test2.py code/eccdna_deg_correlation_analysis_comprehensive.py
mv code/test3.py code/eccdna_deg_correlation_analysis_simple.py
mv code/22222.py code/te_composition_ratio_analysis.py
mv code/222222.py code/eccdna_composite_te_analysis.py
```

#### Phase 4: Update Documentation
- Update `eccToolkit.md` with new script names
- Create usage examples for canonical scripts
- Update any references in other documentation

### 5. Post-Cleanup Directory Structure

After cleanup, the main analysis scripts will be:

```
code/
├── Core Analysis
│   ├── eccdna_analyzer.py                           # Primary eccDNA analyzer
│   ├── eccDNA-SHARP.py                              # Hotspot analysis
│   ├── eccDNA-MultiScale.py                         # Multi-scale analysis
│   └── eccCNV.py                                    # CNV analysis
├── Hotspot Detection
│   ├── eccdna_hotspot_permtest_parallel.py         # Primary hotspot detector
│   └── eccdna_hotspot_refined.py                   # Specialized boundary detection
├── TAD Analysis
│   └── eccTAD.py                                    # Primary TAD analysis
├── Transposable Elements
│   ├── te_motif_composition_analysis_v2.py         # Primary TE motif analysis
│   ├── te_percentage_distribution_analysis.py
│   ├── te_anno_percent_by_class_analysis.py
│   └── te_processor.py
├── Expression Analysis
│   ├── eccdna_deg_correlation_analysis_comprehensive.py  # Renamed from test2.py
│   ├── eccdna_deg_correlation_analysis_simple.py        # Renamed from test3.py
│   ├── te_composition_ratio_analysis.py                 # Renamed from 22222.py
│   └── eccdna_composite_te_analysis.py                  # Renamed from 222222.py
├── Genomic Linkage
│   └── GenoLink_v2.py                               # Primary linkage visualization
├── Data Processing
│   ├── merge_csv_with_sample.py                     # Primary CSV merge
│   └── [other processing scripts]
└── deprecated/
    └── [8 deprecated files]
```

### 6. Benefits of This Cleanup

1. **Reduced Redundancy**: Eliminates 8 duplicate/superseded files
2. **Improved Naming**: Clear, descriptive names for all analysis scripts
3. **Better Maintainability**: Single canonical version for each analysis type
4. **Enhanced Documentation**: Clear purpose and usage for each script
5. **Streamlined Workflows**: Easier to identify the correct script for each task

### 7. Risk Mitigation

1. **Backup Strategy**: All original files preserved in `code_backup_YYYYMMDD/`
2. **Deprecation Directory**: Superseded files moved to `code/deprecated/` rather than deleted
3. **Incremental Changes**: Changes can be applied incrementally and tested
4. **Documentation Updates**: All changes documented with clear rationale

This cleanup plan will result in a more maintainable, clearly organized codebase with 37 scripts instead of the current 45, eliminating redundancy while preserving all functionality.