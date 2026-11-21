# Helpers for the tests: read an example locus shipped in inst/extdata
 
skip_if_no_python <- function() {
    testthat::skip_if_not_installed("reticulate")
    testthat::skip_if_not(reticulate::py_module_available("numpy"), "numpy is not available")
}
 
example_dir <- function() system.file("extdata", package = "FastHer")
 
read_example <- function(stem_z, stem_pos, stem_ld) {
    d <- example_dir()
    SNP <- utils::read.table(file.path(d, paste0(stem_z, ".txt")), header = TRUE, stringsAsFactors = FALSE)
    POS <- utils::read.table(file.path(d, paste0(stem_pos, ".bim")), header = FALSE, stringsAsFactors = FALSE)
    ld <- as.matrix(reticulate::import("numpy")$load(file.path(d, paste0(stem_ld, ".npz")))$f[["LD"]])
    z <- SNP$Z
    n <- SNP$n
    names(z) <- names(n) <- SNP$Predictor
    ord <- POS[, 2]
    maf <- if ("AF1" %in% names(SNP)) stats::setNames(SNP$AF1, SNP$Predictor)[ord] else NULL
    list(z = z[ord], n = n[ord], ld = ld, maf = maf)
}
 
read_bpifc <- function() {
    read_example("ukb_200k_chr22_Z_BPIFC", "ukb_200k_chr22_SNPorder_BPIFC", "ukb_200k_chr22_LD_BPIFC")
}
 
read_sim5k <- function() {
    read_example("sim5k_chr22_1KG_Z", "sim5k_chr22_1KG_SNPorder", "sim5k_chr22_1KG_LD")
}
 
# Small synthetic locus: AR(1) LD, z simulated from the FastHer model
make_locus <- function(m = 60, rho = 0.6, n = 1e5, h2 = 0.02, seed = 11) {
    set.seed(seed)
    U <- rho^abs(outer(seq_len(m), seq_len(m), "-"))
    ev <- eigen(U, symmetric = TRUE)
    lam <- pmax(ev$values, 0)
    z <- c(ev$vectors %*% (sqrt(n * h2 / m * lam^2 + (1 - h2) * lam) * rnorm(m)))
    list(z = z, n = n, ld = U)
}
