# eccToolkit Code Style Guide

## Overview
This document establishes the coding standards and style guidelines for the eccToolkit project. All Python scripts should follow these conventions to ensure consistency, maintainability, and professional quality.

## Table of Contents
1. [General Principles](#general-principles)
2. [Naming Conventions](#naming-conventions)
3. [Code Formatting](#code-formatting)
4. [Documentation Standards](#documentation-standards)
5. [Import Organization](#import-organization)
6. [Error Handling](#error-handling)
7. [Language Usage](#language-usage)
8. [Type Hints](#type-hints)
9. [File Structure](#file-structure)
10. [Tools and Automation](#tools-and-automation)

---

## General Principles

### 1. **Follow PEP 8**
- All Python code should adhere to [PEP 8](https://peps.python.org/pep-0008/) standards
- Use automated tools to enforce consistency
- Prioritize readability over cleverness

### 2. **Be Consistent**
- Use the same patterns throughout the codebase
- Follow established conventions within each file
- Maintain consistency in naming, formatting, and structure

### 3. **Write Self-Documenting Code**
- Use descriptive variable and function names
- Write clear, concise comments
- Structure code logically with proper organization

---

## Naming Conventions

### Variables and Functions
```python
# ✅ Good - snake_case
sample_count = 100
eccdna_file_path = "/path/to/file.csv"
deg_threshold = 2.0

def analyze_eccdna_overlap(sample_data, reference_data):
    """Function names should be descriptive and use snake_case."""
    pass

# ❌ Bad - mixed styles
sampleCount = 100  # camelCase
eccDNA_file_path = "/path/to/file.csv"  # mixed case
degThreshold = 2.0  # camelCase
```

### Classes
```python
# ✅ Good - PascalCase
class EccDNAAnalyzer:
    """Class names should use PascalCase."""
    pass

class TransposonProcessor:
    """Clear, descriptive class names."""
    pass

# ❌ Bad
class eccDNA_analyzer:  # mixed case
    pass
```

### Constants
```python
# ✅ Good - SCREAMING_SNAKE_CASE
DEFAULT_THRESHOLD = 1.0
MAX_ITERATIONS = 1000
BRIGHT_COLORS = ['#FF0000', '#00FF00', '#0000FF']

# ❌ Bad
default_threshold = 1.0  # lowercase
maxIterations = 1000  # camelCase
```

### File Names
```python
# ✅ Good - descriptive, snake_case
eccdna_hotspot_detector.py
te_composition_analysis.py
deg_correlation_analyzer.py

# ❌ Bad
test2.py  # non-descriptive
22222.py  # meaningless numbers
eccDNA-SHARP.py  # mixed case with hyphens
```

---

## Code Formatting

### Line Length
- **Maximum line length: 88 characters** (Black's default)
- Break long lines at logical points
- Use parentheses for multi-line expressions

```python
# ✅ Good
result = analyze_eccdna_overlap(
    sample_data=sample_df,
    reference_data=reference_df,
    threshold=2.0,
    output_dir="/path/to/output"
)

# ❌ Bad
result = analyze_eccdna_overlap(sample_data=sample_df, reference_data=reference_df, threshold=2.0, output_dir="/path/to/output")
```

### Indentation
- **4 spaces per indentation level**
- No tabs
- Consistent continuation line indentation

```python
# ✅ Good
if (condition_one and 
    condition_two and 
    condition_three):
    do_something()

# ❌ Bad
if (condition_one and 
condition_two and 
condition_three):
    do_something()
```

### Spacing
```python
# ✅ Good - proper spacing
result = function_call(a, b, c)
value = x + y * z
data[key] = value

# ❌ Bad - inconsistent spacing
result=function_call(a,b,c)
value = x+y*z
data[ key ]=value
```

### Blank Lines
```python
# ✅ Good - proper blank line usage
import os
import sys

import pandas as pd
import numpy as np


class EccDNAAnalyzer:
    """Class definition with proper spacing."""
    
    def __init__(self):
        self.data = None
    
    def load_data(self, file_path):
        """Method with proper spacing."""
        pass


def standalone_function():
    """Standalone function with proper spacing."""
    pass
```

---

## Documentation Standards

### Docstrings
Use **Google-style docstrings** for all functions, classes, and modules:

```python
def analyze_eccdna_overlap(sample_data, reference_data, threshold=1.0):
    """Analyze overlap between eccDNA and reference datasets.
    
    This function calculates the overlap between eccDNA regions and
    reference genomic features using the specified threshold.
    
    Args:
        sample_data (pd.DataFrame): DataFrame containing eccDNA data with
            columns ['chr', 'start', 'end', 'sample'].
        reference_data (pd.DataFrame): Reference genomic features DataFrame.
        threshold (float, optional): Overlap threshold ratio. Defaults to 1.0.
    
    Returns:
        dict: Dictionary containing overlap statistics including:
            - 'overlap_count': Number of overlapping regions
            - 'overlap_ratio': Ratio of overlapping regions
            - 'significant_overlaps': List of significant overlaps
    
    Raises:
        ValueError: If input DataFrames are empty or missing required columns.
        FileNotFoundError: If reference files are not found.
    
    Example:
        >>> sample_df = pd.read_csv('sample_data.csv')
        >>> reference_df = pd.read_csv('reference_data.csv')
        >>> results = analyze_eccdna_overlap(sample_df, reference_df, threshold=2.0)
        >>> print(f"Overlap count: {results['overlap_count']}")
    """
    pass
```

### Comments
```python
# ✅ Good - clear, descriptive comments
# Calculate the enrichment ratio using hypergeometric test
enrichment_ratio = calculate_enrichment(observed, expected)

# Process each sample separately to avoid memory issues
for sample_name in sample_list:
    process_sample(sample_name)

# ❌ Bad - obvious or unclear comments
# Add 1 to counter
counter += 1

# Do stuff
process_data()
```

### Module-Level Documentation
```python
#!/usr/bin/env python3
"""
eccDNA Hotspot Detection Module

This module provides functions for detecting eccDNA hotspots using
statistical analysis and permutation testing.

Main functions:
    - detect_hotspots: Main hotspot detection function
    - calculate_significance: Statistical significance testing
    - plot_hotspot_distribution: Visualization functions

Example:
    python eccdna_hotspot_detector.py -i input.csv -o output.csv -t 8

Author: eccToolkit Team
License: MIT
"""
```

---

## Import Organization

### Import Order (using isort)
```python
# 1. Standard library imports
import os
import sys
from pathlib import Path

# 2. Third-party imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# 3. Local imports
from .utils import helper_functions
from .analyzers import EccDNAAnalyzer
```

### Import Style
```python
# ✅ Good - specific imports
from pandas import DataFrame, read_csv
from scipy.stats import hypergeom, fisher_exact

# ✅ Good - aliasing for readability
import matplotlib.pyplot as plt
import seaborn as sns

# ❌ Bad - star imports
from pandas import *
from scipy.stats import *
```

---

## Error Handling

### Exception Handling Patterns
```python
# ✅ Good - specific exception handling
def load_data(file_path):
    """Load data with proper error handling."""
    try:
        data = pd.read_csv(file_path)
        return data
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        raise
    except pd.errors.EmptyDataError:
        logger.error(f"Empty data file: {file_path}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error loading {file_path}: {e}")
        raise

# ❌ Bad - bare except
def load_data(file_path):
    try:
        data = pd.read_csv(file_path)
        return data
    except:
        pass
```

### Input Validation
```python
def analyze_data(data, threshold=1.0):
    """Analyze data with proper validation."""
    # Input validation
    if data is None or data.empty:
        raise ValueError("Input data cannot be None or empty")
    
    if not isinstance(threshold, (int, float)):
        raise TypeError("Threshold must be a number")
    
    if threshold <= 0:
        raise ValueError("Threshold must be positive")
    
    # Required columns validation
    required_columns = ['chr', 'start', 'end']
    missing_columns = [col for col in required_columns if col not in data.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")
    
    # Process data
    return process_data(data, threshold)
```

---

## Language Usage

### English Only
```python
# ✅ Good - English only
def calculate_enrichment_ratio(observed_count, expected_count):
    """Calculate enrichment ratio for statistical analysis."""
    if expected_count == 0:
        return float('inf')
    return observed_count / expected_count

# Log messages in English
logger.info("Starting eccDNA analysis")
logger.info(f"Processing {len(samples)} samples")

# ❌ Bad - mixed languages
def 计算富集比例(observed_count, expected_count):
    """计算统计分析的富集比例"""  # Chinese docstring
    if expected_count == 0:
        return float('inf')
    return observed_count / expected_count

# Mixed language log messages
logger.info("Starting eccDNA analysis")
logger.info(f"正在处理 {len(samples)} 个样本")  # Chinese
```

### Variable Names
```python
# ✅ Good - descriptive English names
sample_count = len(samples)
enrichment_threshold = 2.0
significant_genes = []

# ❌ Bad - mixed languages
sample_count = len(samples)
富集阈值 = 2.0  # Chinese variable name
significant_genes = []
```

---

## Type Hints

### Function Type Hints
```python
from typing import List, Dict, Optional, Union, Tuple
import pandas as pd

def analyze_eccdna_overlap(
    sample_data: pd.DataFrame,
    reference_data: pd.DataFrame,
    threshold: float = 1.0,
    output_dir: Optional[str] = None
) -> Dict[str, Union[int, float, List[str]]]:
    """Function with comprehensive type hints."""
    pass

def process_samples(
    sample_list: List[str],
    config: Dict[str, Union[str, int, float]]
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Process multiple samples with type hints."""
    pass
```

### Class Type Hints
```python
class EccDNAAnalyzer:
    """EccDNA analysis class with type hints."""
    
    def __init__(self, config: Dict[str, Union[str, int, float]]) -> None:
        self.config = config
        self.data: Optional[pd.DataFrame] = None
        self.results: Dict[str, Union[int, float]] = {}
    
    def load_data(self, file_path: str) -> pd.DataFrame:
        """Load data with type hints."""
        pass
```

---

## File Structure

### Standard File Template
```python
#!/usr/bin/env python3
"""
Module Name: Brief description of the module

Detailed description of what this module does, its main functions,
and how it fits into the larger eccToolkit project.

Example:
    python module_name.py -i input.csv -o output.csv

Author: eccToolkit Team
License: MIT
"""

# Standard library imports
import os
import sys
from pathlib import Path

# Third-party imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Local imports
from .utils import helper_functions

# Module-level constants
DEFAULT_THRESHOLD = 1.0
MAX_ITERATIONS = 1000

# Module-level variables
logger = logging.getLogger(__name__)


class MainClass:
    """Main class for this module."""
    
    def __init__(self):
        pass


def main_function():
    """Main processing function."""
    pass


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Brief description of what this script does",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Usage examples:
    python module_name.py -i input.csv -o output.csv
    python module_name.py -i input.csv -o output.csv --threshold 2.0
        '''
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help='Input file path')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file path')
    
    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Input validation
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    # Main processing
    try:
        result = main_function()
        logger.info("Processing completed successfully")
    except Exception as e:
        logger.error(f"Processing failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
```

---

## Tools and Automation

### Recommended Tools
1. **Black** - Code formatting
2. **isort** - Import sorting
3. **flake8** - Style checking
4. **mypy** - Type checking
5. **pre-commit** - Git hooks

### Configuration Files

#### pyproject.toml
```toml
[tool.black]
line-length = 88
target-version = ['py38']
include = '\.pyi?$'
extend-exclude = '''
/(
  # directories
  \.eggs
  | \.git
  | \.hg
  | \.mypy_cache
  | \.tox
  | \.venv
  | build
  | dist
)/
'''

[tool.isort]
profile = "black"
multi_line_output = 3
line_length = 88
known_first_party = ["eccToolkit"]

[tool.mypy]
python_version = "3.8"
warn_return_any = true
warn_unused_configs = true
disallow_untyped_defs = true
```

#### .flake8
```ini
[flake8]
max-line-length = 88
extend-ignore = E203, W503
exclude = 
    .git,
    __pycache__,
    build,
    dist,
    *.egg-info,
    deprecated/
```

---

## Implementation Strategy

### Phase 1: Critical Issues (Week 1)
1. **Language standardization** - Convert all Chinese comments/variables to English
2. **Basic PEP 8 compliance** - Fix indentation, spacing, naming
3. **Import organization** - Standardize import order

### Phase 2: Documentation (Week 2)
1. **Docstring standardization** - Convert all to Google-style
2. **Type hints** - Add to all public functions
3. **Error handling** - Implement consistent patterns

### Phase 3: Advanced Features (Week 3)
1. **Pre-commit hooks** - Automate style checking
2. **Continuous integration** - Add style checks to CI/CD
3. **Code review guidelines** - Establish review process

### Phase 4: Long-term Maintenance
1. **Regular style audits** - Monthly code style reviews
2. **Tool updates** - Keep formatting tools updated
3. **Documentation updates** - Maintain style guide

---

## Enforcement

### Pre-commit Hooks
```yaml
# .pre-commit-config.yaml
repos:
  - repo: https://github.com/psf/black
    rev: 23.3.0
    hooks:
      - id: black
        language_version: python3.8
  
  - repo: https://github.com/pycqa/isort
    rev: 5.12.0
    hooks:
      - id: isort
        args: ["--profile", "black"]
  
  - repo: https://github.com/pycqa/flake8
    rev: 6.0.0
    hooks:
      - id: flake8
  
  - repo: https://github.com/pre-commit/mirrors-mypy
    rev: v1.3.0
    hooks:
      - id: mypy
```

### CI/CD Integration
```yaml
# .github/workflows/style-check.yml
name: Code Style Check
on: [push, pull_request]

jobs:
  style-check:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.8'
      - name: Install dependencies
        run: |
          pip install black isort flake8 mypy
      - name: Run style checks
        run: |
          black --check .
          isort --check-only .
          flake8 .
          mypy .
```

---

## Benefits of Standardization

### 1. **Improved Maintainability**
- Consistent code is easier to understand and modify
- Reduced cognitive load when switching between files
- Easier onboarding for new contributors

### 2. **Better Collaboration**
- Consistent style reduces merge conflicts
- Easier code reviews with standardized patterns
- Improved team productivity

### 3. **Professional Quality**
- Consistent codebase appears more professional
- Easier to package and distribute
- Better documentation and examples

### 4. **Automated Quality Control**
- Reduced manual code review overhead
- Consistent enforcement through tooling
- Prevention of style regressions

---

## Conclusion

This style guide provides a comprehensive foundation for maintaining consistent, high-quality code in the eccToolkit project. By following these guidelines and using the recommended tools, we can ensure that the codebase remains maintainable, readable, and professional as it continues to grow and evolve.

Regular adherence to these standards will significantly improve the project's overall quality and make it easier for both current and future contributors to work with the code effectively.