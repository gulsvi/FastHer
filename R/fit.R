###################################################################
###        Main function: local heritability from memory        ###
###################################################################
 
#' Local heritability by exact REML from summary statistics
#'
#' Estimates local heritability h2 and the concordance parameter d from
#' z-statistics and an LD matrix of a genomic region. The LD matrix is
#' eigen-decomposed once; the REML fit then uses only the eigenvalues and the
#' projections of z on the eigenvectors.
#'
#' @param z Numeric vector of z-statistics (one per SNP, in the order of the LD matrix)
#' @param n Sample size: one number or a vector (one value per SNP)
#' @param ld SNP-by-SNP LD (correlation) matrix of the region; rows and columns must be in
#'   the same SNP order as z (if both have names, the names are compared)
#' @param weights Optional SNP weights (one per SNP, e.g. from
#'   \code{\link{fasther_weights}}); they are scaled to mean square 1 and applied
#'   as z -> w z and LD -> W LD W before the eigen-decomposition
#' @param reml REML algorithm: "auto" (default; profile REML, switching to the
#'   two-parameter REML if the estimate is not finite or lies near the boundary
#'   of [0, 1]), "profile" or "two-parameter"
#' @param tau Fraction of the total LD variance retained in the eigenvalue
#'   truncation (default 0.99999, as in FastHer 0.1.0)
#' @param control_list Optimization control for the two-parameter REML (optional)
#' @param eigen_method "auto" (default: Python/scipy for more than 100 SNPs,
#'   otherwise R), "python" or "R"
#' @param regularize If TRUE, the (reference) LD matrix is regularized by
#'   truncating its spectrum: starting from the fraction \code{reg_tau} of the LD
#'   variance, components are dropped from the tail one at a time while the Wald
#'   statistic for d = 1 exceeds 1.96 (default FALSE). In this case \code{tau} is
#'   not used
#' @param reg_tau Starting fraction of the LD variance retained when
#'   \code{regularize = TRUE} (default 0.97)
#' @param reg_control Optional list with elements \code{min_fraction} (stop when the
#'   retained fraction is not larger than this value, default 0.95),
#'   \code{min_rank} (stop when the rank would not exceed this value, default 5)
#'   and \code{wald} (Wald cut-off, default 1.96)
#' @return An object of class "fasther": a list with h2, se_h2, d, se_d,
#'   wald_d, p_d, m, rank, n_eff, weighted, method, fallback, convergence, iterations,
#'   hard_settings, initial_h2, initial_se, time_sec. With regularization it
#'   also has regularized, reg_fraction (retained fraction of the LD variance),
#'   reg_converged (TRUE if the Wald statistic for d reached the cut-off) and
#'   reg_path (a data frame with the dropping path).
#' @export
fasther_fit <- function(z, n, ld, weights = NULL, reml = c("auto", "profile", "two-parameter"),
                        tau = 0.99999, control_list = NULL,
                        eigen_method = c("auto", "python", "R"),
                        regularize = FALSE, reg_tau = 0.97, reg_control = NULL) {
    start_time <- Sys.time()
    reml <- match.arg(reml)
    eigen_method <- match.arg(eigen_method)
 
    ## input checks
    if (!is.numeric(z) || !all(is.finite(z))) {
        stop("'z' must be a numeric vector without missing or infinite values.", call. = FALSE)
    }
    m <- length(z)
    ld <- as.matrix(ld)
    if (nrow(ld) != m || ncol(ld) != m) {
        stop("'ld' must be a square matrix with as many rows as 'z' has elements (", m, ").",
             call. = FALSE)
    }
    if (!is.numeric(n) || !all(is.finite(n)) || any(n <= 0) || !(length(n) %in% c(1L, m))) {
        stop("'n' must be a positive number or a vector with one value per SNP.", call. = FALSE)
    }
    if (length(n) == 1L) n <- rep(n, m)
    if (!is.numeric(tau) || length(tau) != 1L || tau <= 0 || tau >= 1) {
        stop("'tau' must be a number in (0, 1).", call. = FALSE)
    }
 
    .check_snp_names(z, ld, weights)

    ## harmonize sample sizes
    n_eff <- 1 / mean(1 / n)
    z <- z * sqrt(n_eff / n)

    ## SNP weights (optional)
    weighted <- !is.null(weights)
    if (weighted) {
        w <- .check_weights(weights, m)
        z <- z * w
        ld <- ld * tcrossprod(w)
    }
 
    ## eigen-decomposition (once) and truncation
    eigU <- .eigen_decomp(ld, method = eigen_method)
    Uval <- eigU$values
    Uvec <- eigU$vectors
    Uval[Uval < 0] <- 0
    cum <- cumsum(Uval) / sum(Uval)
    if (regularize) {
        ctl <- .reg_control(reg_control)
        if (!is.numeric(reg_tau) || length(reg_tau) != 1L || reg_tau <= 0 || reg_tau >= 1) {
            stop("'reg_tau' must be a number in (0, 1).", call. = FALSE)
        }
        keep <- seq_len(min(length(Uval), max(sum(cum < reg_tau), ctl$min_rank + 1L)))
    } else {
        threshold <- Uval[min(which(cum > tau))]
        keep <- Uval >= threshold
    }
    Uvec <- Uvec[, keep, drop = FALSE]
    Uval <- Uval[keep]
    rankU <- length(Uval)
 
    ## projections of z on the eigenvectors
    Vz2 <- c(crossprod(z, Uvec))^2
 
    ## REML
    reg <- NULL
    if (regularize) {
        reg <- .reml_regularized(Vz2, Uval, cum[keep], n_eff, m, reml, control_list, reg_tau, ctl)
        fit <- reg$fit
        rankU <- reg$rank
    } else {
        fit <- .reml_fit(Vz2, Uval, n_eff, m, rankU, reml = reml, control_list = control_list)
    }
    if (fit$code != 0) {
        warning("Optimization did not converge properly. Code: ", fit$code, call. = FALSE)
    }
 
    d <- fit$estimate[["d"]]
    se_d <- fit$se[["d"]]
    wald_d <- (d - 1) / se_d
    out <- list(h2 = fit$estimate[["h2"]], se_h2 = fit$se[["h2"]],
                d = d, se_d = se_d,
                wald_d = wald_d, p_d = 2 * stats::pnorm(-abs(wald_d)),
                m = m, rank = rankU, n_eff = n_eff, weighted = weighted,
                regularized = regularize,
                reg_fraction = if (regularize) reg$fraction else NA_real_,
                reg_converged = if (regularize) reg$converged else NA,
                reg_path = reg$path,
                method = fit$method, fallback = fit$fallback,
                convergence = fit$code, iterations = as.numeric(fit$iterations),
                hard_settings = fit$hard_settings,
                initial_h2 = fit$initial_h2, initial_se = fit$initial_se,
                time_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    class(out) <- "fasther"
    out
}
 
#' Print a FastHer result
#' @param x Object of class "fasther"
#' @param digits Number of significant digits
#' @param ... Not used
#' @export
print.fasther <- function(x, digits = 4, ...) {
    cat("FastHer: local heritability (", x$method, " REML)\n", sep = "")
    cat("  h2 = ", format(x$h2, digits = digits), " (SE ", format(x$se_h2, digits = digits), ")\n", sep = "")
    cat("  d  = ", format(x$d, digits = digits), " (SE ", format(x$se_d, digits = digits),
        "); Wald(d = 1) = ", format(x$wald_d, digits = 3), ", p = ", format(x$p_d, digits = 3), "\n", sep = "")
    cat("  m = ", x$m, " SNPs, rank = ", x$rank, ", effective n = ", format(x$n_eff, digits = 8), "\n", sep = "")
    if (isTRUE(x$regularized)) {
        cat("  regularized LD: retained fraction = ", format(x$reg_fraction, digits = 4),
            if (isTRUE(x$reg_converged)) "" else " (Wald cut-off not reached)", "\n", sep = "")
    }
    invisible(x)
}
 
#' Convert a FastHer result to a one-row data frame
#' @param x Object of class "fasther"
#' @param row.names,optional,... Not used
#' @export
as.data.frame.fasther <- function(x, row.names = NULL, optional = FALSE, ...) {
    data.frame(initial_h2 = x$initial_h2, initial_se = x$initial_se,
               m = x$m, rankU = x$rank, weighted = x$weighted,
               regularized = x$regularized, reg_fraction = x$reg_fraction,
               reg_converged = x$reg_converged,
               FH_h2 = x$h2, FH_se = x$se_h2, FH_d = x$d, FH_sed = x$se_d,
               wald_d = x$wald_d, p_d = x$p_d,
               method = x$method, convergence = x$convergence,
               iterations = x$iterations, hard_settings = x$hard_settings,
               time_sec = x$time_sec, stringsAsFactors = FALSE)
}
