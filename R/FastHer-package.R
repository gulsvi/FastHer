#' FastHer: Fast and Accurate Local Heritability REML Estimator
#'
#' Estimates local SNP heritability (h2) and the model-data concordance
#' parameter d from GWAS summary statistics (z-statistics) and an LD matrix by
#' exact restricted maximum likelihood (REML). The LD matrix is
#' eigen-decomposed once; each REML iteration then costs O(r), where r is the
#' rank of the matrix.
#'
#' @author Gulnara R. Svishcheva
#' @keywords internal
"_PACKAGE"
