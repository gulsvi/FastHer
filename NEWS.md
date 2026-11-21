# FastHer 1.1.0
 
## New features
 
* `fasther_fit()`: new main function working on data in memory (z-statistics, sample size and LD matrix).
* Two REML implementations, selected with `reml`: exact profile REML (default, via `"auto"`) and
  two-parameter REML. In `"auto"` mode the profile REML is used and the two-parameter REML takes over
  when the estimate is not finite or lies at the boundary.
* Wald test of the concordance parameter d = 1 (`wald_d`, `p_d`).
* SNP weights: `fasther_weights()` computes LDAK-type weights from allele frequencies; they are passed
  with `weights` (`fasther_fit()`) or `use_weights`, `power` and `maf_col` (`FastHer()`).
* Regularization of reference LD matrices (`regularize`, `reg_tau`, `reg_control`): components with the
  smallest eigenvalues are removed one by one, starting from 97% of the LD variance, until the Wald test
  does not reject d = 1.
* New example with allele frequencies (`sim5k_chr22_1KG_*`) and a demo for SNP weights.
* Tests (`testthat`), including regression tests against FastHer 1.0.0.
 
## Changes
 
* `FastHer()` is now a file-based wrapper of `fasther_fit()` with the same arguments and output files;
  the result table has additional columns (`weighted`, `regularized`, `reg_fraction`, `reg_converged`,
  `wald_d`, `p_d`, `method`).
* The two-parameter REML reproduces FastHer 1.0.0; the default profile REML agrees with it to about 1e-8
  in h2 on the example loci.
* Example data moved to `inst/extdata`; the large HIRA example was removed from the package.
* `FastHer()` no longer changes the number of threads or attaches packages; the unused dependencies
  `data.table`, `compiler` and `Matrix` were removed.
* The eigen-decomposition uses Python (scipy) for more than 100 SNPs and falls back to R, with a warning,
  when scipy is not available.
 
# FastHer 1.0.0
 
* Initial release.
