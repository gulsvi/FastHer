# FastHer: Fast Local Heritability Estimation from GWAS Summary Statistics
 
## Overview
 
FastHer estimates local SNP heritability (h2) from GWAS summary statistics (z-statistics) and a
linkage disequilibrium (LD) matrix by exact restricted maximum likelihood (REML). The LD matrix is
eigen-decomposed once; every REML iteration then costs O(r), where r is the rank of the matrix.
 
## Features
 
- Exact REML for h2 and for the model-data concordance parameter d, with analytical standard errors
- Two REML implementations: exact profile REML (default) and two-parameter REML, selected automatically
- Wald test of d = 1 (`wald_d`, `p_d`)
- SNP weights computed from allele frequencies (LDAK-type, `power = -0.25` by default)
- Regularization of reference LD matrices (adaptive eigenvalue truncation guided by the Wald test for d)
- In-memory function `fasther_fit()` and a file-based interface `FastHer()`
 
## Installation
 
~~~r
devtools::install_github("gulsvi/FastHer")
library(FastHer)
~~~
 
Eigen-decomposition of large matrices (more than 100 SNPs) uses Python (scipy) through `reticulate`.
If scipy is not available, R is used instead (slower, with a warning). To install the Python modules:
 
~~~r
reticulate::py_install(c("numpy", "scipy"))
~~~
 
## Quick start
 
### From data in memory
 
`z` (z-statistics), `n` (sample size, one number or one value per SNP) and `ld` (SNP-by-SNP LD matrix)
must be in the same SNP order.
 
~~~r
fit <- fasther_fit(z, n, ld)
fit
as.data.frame(fit)
~~~
 
With SNP weights (`maf`: allele frequencies strictly between 0 and 1):
 
~~~r
fit <- fasther_fit(z, n, ld, weights = fasther_weights(maf, power = -0.25))
~~~
 
With regularization of a reference LD matrix (e.g. from 1000 Genomes):
 
~~~r
fit <- fasther_fit(z, n, ref_ld, regularize = TRUE, reg_tau = 0.97)
~~~
 
### From files
 
~~~r
result <- FastHer(
    path.in = system.file("extdata", package = "FastHer"),
    file.Z = "ukb_200k_chr22_Z_BPIFC",
    file.POS = "ukb_200k_chr22_SNPorder_BPIFC",
    file.LD = "ukb_200k_chr22_LD_BPIFC",
    path.out = tempdir()
)
~~~
 
Add `use_weights = TRUE` for SNP weights (allele frequencies are read from a column MAF, AF, AF1, EAF
or FREQ of the z-statistics file, or from the column given by `maf_col`) and `regularize = TRUE` for a
reference LD matrix. Note that `path.out` must end with "/" when it is not the default.
 
### Demos
 
~~~r
demo(package = "FastHer")
demo("FastHer_BPIFC_example")
demo("FastHer_weights_example")
~~~
 
## Main functions
 
- `fasther_fit(z, n, ld, weights, reml, tau, control_list, eigen_method, regularize, reg_tau, reg_control)`:
  local heritability from data in memory. `reml` is "auto" (default), "profile" or "two-parameter";
  `tau` is the fraction of LD variance kept in the eigenvalue truncation (default 0.99999).
- `fasther_weights(maf, power = -0.25)`: SNP weights from allele frequencies.
- `FastHer(path.in, file.Z, file.POS, file.LD, ...)`: file-based interface to `fasther_fit()`; writes
  `FH_h2.txt`, `FH_se.txt` and a results table.
- `eigen.Python(U, only.values = FALSE)`: eigen-decomposition with scipy.
 
## Input files
 
1. Z-statistics file (.txt, tab-separated): columns `Predictor`, `Z`, `n`; optionally an allele-frequency
   column (`MAF`, `AF`, `AF1`, ...) for SNP weights.
2. SNP order file (.bim): SNP names in column 2, in the order of the LD matrix.
3. LD matrix (.npz): numpy file with the key `LD`.
 
## Output
 
`FastHer()` returns a one-row data frame (`as.data.frame()` of a `fasther` object gives the same):
 
- initial_h2, initial_se: one-parameter estimate used as starting value
- m, rankU: number of SNPs and rank after eigenvalue truncation
- weighted, regularized, reg_fraction, reg_converged: weighting and regularization information
- FH_h2, FH_se: heritability and its standard error
- FH_d, FH_sed: concordance parameter d and its standard error
- wald_d, p_d: Wald statistic and p-value for d = 1
- method, convergence, iterations, hard_settings, time_sec
 
## Method
 
1. The LD matrix (after optional SNP weighting) is eigen-decomposed once.
2. The z-statistics are projected on the eigenvectors.
3. h2 and d are estimated by REML using only the eigenvalues and the projections; standard errors come
   from the expected information matrix.
4. For a reference LD matrix, components with the smallest eigenvalues are removed one by one, starting
   from 97% of the LD variance, until the Wald test does not reject d = 1.
 
## Example data
 
`inst/extdata` contains a UK Biobank-based locus (BPIFC) and a simulated locus with allele frequencies
(`sim5k_chr22_1KG_*`); see `inst/extdata/README.md`.
 
## License
 
MIT. Author: Gulnara R. Svishcheva.
