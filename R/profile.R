###################################################################
###          Exact profile REML (one-parameter search)          ###
###################################################################
 
#' Expected Hessian of the REML log-likelihood with respect to (h2, d)
#'
#' @param dd_1 Vector (n/m) * Uval - 1
#' @param n Sample size
#' @param n_r n - r
#' @param h2 Heritability
#' @param d Concordance parameter
#' @return 2 x 2 matrix with dimnames h2, d
#' @keywords internal
.expected_hessian <- function(dd_1, n, n_r, h2, d) {
    Wr <- dd_1 * h2 + 1
    K <- dd_1 / Wr
    SD.h2 <- -0.5 * (sum(K^2) + n_r / (1 - h2)^2)
    SD.d <- -0.5 * n / d^2
    SD.h2.d <- (-0.5 / d) * sum(K) + 0.5 * n_r / (d * (1 - h2))
    matrix(c(SD.h2, SD.h2.d, SD.h2.d, SD.d), nrow = 2, ncol = 2,
           dimnames = list(c("h2", "d"), c("h2", "d")))
}
 
#' Profile REML estimator of h2 (and d)
#'
#' For fixed h2 the REML objective has a closed-form minimum over d,
#' d(h2) = (S(h2) + (n - zUz) / (1 - h2)) / n, with
#' S(h2) = sum(Vz2 / Uval / W), W = (n/m) * Uval * h2 + (1 - h2).
#' Substituting it back gives the one-parameter profile objective
#' Phi(h2) = sum(log W) + (n - r) * log(1 - h2) + n * log(S(h2) + (n - zUz) / (1 - h2)),
#' which is minimized over [0, 1) by a bounded scalar search (optimize()).
#' Standard errors are obtained from the expected Hessian in (h2, d)
#' evaluated at the estimate.
#'
#' @param Vz2 Squared projections of z-statistics on the eigenvectors
#' @param Uval Positive eigenvalues of the LD matrix
#' @param n Sample size
#' @param m Number of SNPs
#' @param r Rank (number of eigenvalues used)
#' @param eps Upper bound is 1 - eps
#' @param tol Tolerance of the scalar search
#' @return List with estimate, se, hessian, objective, code, message, iterations, type
#' @keywords internal
.reml_profile <- function(Vz2, Uval, n, m, r, eps = 1e-9, tol = 1e-12) {
    Vz2dUval <- Vz2 / Uval
    zUz <- sum(Vz2dUval)
    dd_1 <- Uval * (n / m) - 1
    n_r <- n - r
    n_zUz <- n - zUz
 
    Phi <- function(h2) {
        Wr <- dd_1 * h2 + 1
        a <- sum(Vz2dUval / Wr) + n_zUz / (1 - h2)
        if (!is.finite(a) || a <= 0) return(1e300)
        sum(log(Wr)) + n_r * log(1 - h2) + n * log(a)
    }
 
    opt <- stats::optimize(Phi, interval = c(0, 1 - eps), tol = tol)
    h2_hat <- opt$minimum
    obj <- opt$objective
    # the minimum may lie on the boundary h2 = 0
    if (Phi(0) <= obj) {
        h2_hat <- 0
        obj <- Phi(0)
    }
 
    Wr <- dd_1 * h2_hat + 1
    d_hat <- (sum(Vz2dUval / Wr) + n_zUz / (1 - h2_hat)) / n
 
    hessian <- .expected_hessian(dd_1, n, n_r, h2_hat, d_hat)
    se <- sqrt(diag(solve(-hessian, tol = 1e-20)))
 
    list(estimate = c(h2 = h2_hat, d = d_hat),
         se = se,
         hessian = hessian,
         objective = obj,
         code = 0L,
         message = "profile REML",
         iterations = NA_integer_,
         type = "profile")
}
