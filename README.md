
FastHer/
+-- README.md
+-- R/
¦   L-- fasther.R                      # R program
+-- Results/
¦   +-- README.md                      # readme
¦   +-- FH_h2.txt                      # h2 (in current analysis)
¦   +-- FH_se.txt                      # SE(h2) (in current analysis)
¦   +-- results_HIRA.txt               # results (HIRA gene)
¦   L-- results_BPIFC.txt              # results (BPIFC gene)
+-- data/
¦   +-- README.md                          # readme
¦   +-- ukb_200k_chr22_LD_HIRA.npz         # LD matrix
¦   +-- ukb_200k_chr22_Z_HIRA.txt          # Z-statistics
¦   +-- ukb_200k_chr22_SNPorder_HIRA.bim   # SNP info
¦   +-- ukb_200k_chr22_LD_BPIFC.npz        # LD matrix
¦   +-- ukb_200k_chr22_Z_BPIFC.txt         # Z-statistics
¦   L-- ukb_200k_chr22_SNPorder_BPIFC.bim  # SNP info
L-- demo/
    +-- Example_200k_HIRA_run.R         # demo with HIRA gene
    +-- Example_200k_BPIFC_run.R        # demo with BPIFC gene
    +-- Example_200k_HIRA_run.log       # log file
    +-- Example_200k_BPIFC_run.log      # log file
    +-- Example_200k_HIRA_run.sh        # script for run
    L-- Example_200k_BPIFC_run.sh       # script for run



# FastHer

Fast and Accurate Local Heritability Estimation 

## Description

FastHer is an R package for fast and accurate estimation of local heritability using GWAS summary statistics and LD matrices. 

## Features

- **Fast computation** using optimized likelihood function based on eigen-decomposition
- **Memory efficient** handling of large LD matrices
- **Robust optimization** with proper convergence checking
- **Comprehensive output** including standard errors and diagnostics

## Installation

```r
# Install required packages
install.packages(c("reticulate", "data.table", "compiler"))

# Install Python dependencies
reticulate::py_install(c("scipy", "numpy"))
```

## Usage

```r
# Load the function
source("FastHer.R")

# Run FastHer
result <- FastHer(
    path = "/path/to/your/data/",
    file.Z = "z_scores",
    file.MAF = "snp_info", 
    file.LD = "ld_matrix"
)

# View results
print(result)
```

## Input Files

1. **Z-score file** (`file.Z.txt`): 
   - Columns: `Predictor`, `Z`, `n`
   - Tab-separated, with header

2. **POS file** (`file.POS.bim`):
   - Standard PLINK .bim format
   - SNP information

3. **LD matrix** (`file.LD.npz`):
   - NumPy compressed format
   - LD matrix (SNP-by-SNP correlation matrix)

## Output

The function returns a data frame with:
- Heritability estimate (hÂ²) and standard error
- Additional parameter estimates
- Convergence information
- Computation time

## Citation

If you use FastHer in your research, please cite:
[Publication in preparation]

## License

MIT License