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
#' @param ld SNP-by-SNP LD (correlation) matrix of the region
#' @param reml REML algorithm: "auto" (default; profile REML, switching to the
#'   two-parameter REML if the estimate is not finite or lies near the boundary
#'   of [0, 1]), "profile" or "two-parameter"
#' @param tau Fraction of the total LD variance retained in the eigenvalue
#'   truncation (default 0.99999, as in FastHer 0.1.0)
#' @param control_list Optimization control for the two-parameter REML (optional)
#' @param eigen_method "auto" (default: Python/scipy for more than 100 SNPs,
#'   otherwise R), "python" or "R"
#' @return An object of class "fasther": a list with h2, se_h2, d, se_d,
#'   wald_d, p_d, m, rank, n_eff, method, fallback, convergence, iterations,
#'   hard_settings, initial_h2, initial_se, time_sec.
#' @export
fasther_fit <- function(z, n, ld, reml = c("auto", "profile", "two-parameter"),
                        tau = 0.99999, control_list = NULL,
                        eigen_method = c("auto", "python", "R")) {
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
 
    ## harmonize sample sizes
    n_eff <- 1 / mean(1 / n)
    z <- z * sqrt(n_eff / n)
 
    ## eigen-decomposition (once) and truncation
    eigU <- .eigen_decomp(ld, method = eigen_method)
    Uval <- eigU$values
    Uvec <- eigU$vectors
    Uval[Uval < 0] <- 0
    cum <- cumsum(Uval) / sum(Uval)
    threshold <- Uval[min(which(cum > tau))]
    keep <- Uval >= threshold
    Uvec <- Uvec[, keep, drop = FALSE]
    Uval <- Uval[keep]
    rankU <- length(Uval)
 
    ## projections of z on the eigenvectors
    Vz2 <- c(crossprod(z, Uvec))^2
 
    ## REML
    fit <- .reml_fit(Vz2, Uval, n_eff, m, rankU, reml = reml, control_list = control_list)
    if (fit$code != 0) {
        warning("Optimization did not converge properly. Code: ", fit$code, call. = FALSE)
    }
 
    d <- fit$estimate[["d"]]
    se_d <- fit$se[["d"]]
    wald_d <- (d - 1) / se_d
    out <- list(h2 = fit$estimate[["h2"]], se_h2 = fit$se[["h2"]],
                d = d, se_d = se_d,
                wald_d = wald_d, p_d = 2 * stats::pnorm(-abs(wald_d)),
                m = m, rank = rankU, n_eff = n_eff,
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
    invisible(x)
}
 
#' Convert a FastHer result to a one-row data frame
#' @param x Object of class "fasther"
#' @param row.names,optional,... Not used
#' @export
as.data.frame.fasther <- function(x, row.names = NULL, optional = FALSE, ...) {
    data.frame(initial_h2 = x$initial_h2, initial_se = x$initial_se,
               m = x$m, rankU = x$rank,
               FH_h2 = x$h2, FH_se = x$se_h2, FH_d = x$d, FH_sed = x$se_d,
               wald_d = x$wald_d, p_d = x$p_d,
               method = x$method, convergence = x$convergence,
               iterations = x$iterations, hard_settings = x$hard_settings,
               time_sec = x$time_sec, stringsAsFactors = FALSE)
}
