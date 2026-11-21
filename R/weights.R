###################################################################
###                        SNP weights                          ###
###################################################################
 
#' SNP weights from minor allele frequencies
#'
#' Computes weights for the weighted FastHer model. With the LDAK power
#' \code{power}, the expected per-SNP heritability is proportional to
#' \code{[MAF (1 - MAF)]^(1 + power)}, so the weight of a SNP is
#' \code{[MAF (1 - MAF)]^((1 + power) / 2)}. The weights are scaled so that
#' their mean square equals 1. With \code{power = -1} all weights are equal to 1
#' (no weighting).
#'
#' @param maf Numeric vector of minor (or alternative) allele frequencies,
#'   strictly between 0 and 1, in the same order as the z-statistics
#' @param power LDAK power (default -0.25)
#' @return Numeric vector of weights, to be passed to \code{\link{fasther_fit}}
#' @export
fasther_weights <- function(maf, power = -0.25) {
    if (!is.numeric(maf) || length(maf) == 0L || !all(is.finite(maf)) || any(maf <= 0 | maf >= 1)) {
        stop("'maf' must be a numeric vector with values strictly between 0 and 1.", call. = FALSE)
    }
    if (!is.numeric(power) || length(power) != 1L || !is.finite(power)) {
        stop("'power' must be a single number.", call. = FALSE)
    }
    w <- (maf * (1 - maf))^((1 + power) / 2)
    .check_weights(w, length(w))
}
 
# Checks weights and scales them to mean square 1 (needed so that the
# weighted LD matrix keeps trace m)
.check_weights <- function(w, m) {
    if (!is.numeric(w) || length(w) != m || !all(is.finite(w)) || any(w <= 0)) {
        stop("'weights' must be a vector of ", m, " finite positive numbers (one per SNP).",
             call. = FALSE)
    }
    w <- as.vector(w)
    w / sqrt(mean(w^2))
}
