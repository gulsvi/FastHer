# FastHer: Fast Local Heritability Estimation from GWAS Summary Statistics

## Overview

FastHer is an R package for fast and accurate local heritability estimation using GWAS summary statistics and LD matrices. It implements the REML function optimized by eigen-decomposition of the LD matrix.

## Features

- Fast REML-based heritability estimation
- Works with GWAS summary data (Z-statistics)
- Uses pre-computed LD matrices
- Python integration for efficient eigen-decomposition
- Comprehensive output with standard errors

## Installation

Install from GitHub:
devtools::install_github("gulsvi/FastHer")

Load the package:
library(FastHer)

## Quick Start

### Basic Usage

Run example analysis with included data:
result <- FastHer(
    path.in = system.file("data", package = "FastHer"),
    file.Z = "ukb_200k_chr22_Z_BPIFC",
    file.POS = "ukb_200k_chr22_SNPorder_BPIFC", 
    file.LD = "ukb_200k_chr22_LD_BPIFC",
    path.out = "Results/",
    results_file = "results.txt"
)

print(result)

### Try Demos

List available demos:
demo(package = "FastHer")

Run BPIFC example:
demo("FastHer_BPIFC_example")

Run HIRA example:
demo("FastHer_HIRA_example")

## Documentation

Get help for main function:
?FastHer

Get help for eigen decomposition:
?eigen.Python

## Main Functions

### FastHer()

Main function for local heritability estimation.

Parameters:
- path.in: Path to input data files
- file.Z: Z-statistics file name (without extension)
- file.POS: SNPs info file name (without extension)
- file.LD: LD matrix file name (without extension)
- path.out: Output directory (default: "Results/")
- results_file: Output file name (default: "FastHer_results.txt")
- control_list: Optimization control parameters
- check_dependencies: Check for required packages (default: FALSE)

Returns: Data frame with heritability estimates and statistics.

### eigen.Python()

Eigen-decomposition using Python's Scipy.

Parameters:
- U: Symmetric correlation SNP-by-SNP matrix (LD matrix)
- only.values: Whether to return only eigenvalues (default: FALSE)

Returns: List with eigenvalues and eigenvectors.

## Data Format

### Input Files Required:

1. Z-statistics file (.txt): Contains Z-scores and sample sizes
   Format:
   Predictor    Z       n
   rs10000      0.123   200000
   rs10001      -0.456  200000

2. POS file (.bim): SNP information and order (PLINK BIM format)

3. LD matrix (.npz): Linkage disequilibrium matrix in numpy compressed format

## Example

### Custom Optimization Parameters

Example with HIRA data:
result <- FastHer(
    path.in = "data/",
    file.Z = "ukb_200k_chr22_Z_HIRA",
    file.POS = "ukb_200k_chr22_SNPorder_HIRA",
    file.LD = "ukb_200k_chr22_LD_HIRA",
    control_list = list(factr = 1e7, maxit = 5000)
)

## Output

The function returns a data frame with:
- initial_h2: Initial heritability estimate
- initial_se: Standard error of initial estimate
- m: Number of SNPs in the region
- rankU: Rank after eigenvalue pruning
- FH_h2: Final heritability estimate
- FH_se: Standard error of final estimate
- FH_d: Additional parameter estimate
- FH_sed: Standard error of additional parameter
- convergence: Optimization convergence code
- iterations: Number of iterations
- hard_settings: Whether hard optimization settings were used
- time_sec: Computation time in seconds

## Method

FastHer implements a REML-based approach for local heritability estimation that:
1. Performs eigen-decomposition of the LD matrix
2. Projects Z-statistics onto the eigenbasis
3. Uses optimization to estimate heritability parameters
4. Provides standard errors via Hessian matrix

## Dependencies

R Packages:
- reticulate - Python integration
- data.table - Efficient data handling
- Matrix - Matrix operations
- compiler - Byte compilation

Python Packages (automatically checked):
- scipy - Scientific computing
- numpy - Numerical operations

## Requirements
- R version >= 4.0.0
- Python with scipy and numpy

## Troubleshooting

Common Issues:
1. Python not found: Ensure Python is installed and accessible
2. Missing data files: Check file paths and extensions
3. Memory issues: Large LD matrices may require significant RAM

Check Dependencies:
result <- FastHer(..., check_dependencies = TRUE)

## Citation

If you use FastHer in your research, please cite:
Svishcheva, G. R. (2025). FastHer: Fast Local Heritability Estimation from GWAS Summary Statistics (Version 0.1.0) [Computer software]. Retrieved from https://github.com/gulsvi/FastHer

## Contact

- Author: Gulnara R. Svishcheva
- Email: gulsvi@mail.ru
- GitHub: gulsvi

## License

This package is licensed under the MIT License.