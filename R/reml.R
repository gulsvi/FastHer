###################################################################
###        REML dispatcher: profile / two-parameter / auto      ###
###################################################################
 
#' Fit REML on an already decomposed locus
#'
#' The LD matrix is decomposed once (outside this function); here only the
#' squared projections Vz2 and the eigenvalues Uval are used.
#'
#' reml = "profile": one-parameter profile REML (see .reml_profile).
#' reml = "two-parameter": joint REML in (h2, d), L-BFGS-B started from the
#'   one-parameter estimate at d = 1 (Brent's method); stricter settings are
#'   used if |d - 1| >= 0.3 (as in FastHer 1.0.0).
#' reml = "auto" (default): profile REML; if the estimate is not finite or lies
#'   near the boundary of [0, 1] (closer than `boundary`), the two-parameter
#'   REML is used instead.
#'
#' @param Vz2 Squared projections of z-statistics on the eigenvectors
#' @param Uval Positive eigenvalues of the LD matrix
#' @param n Sample size
#' @param m Number of SNPs
#' @param r Rank (number of eigenvalues used)
#' @param reml "auto", "profile" or "two-parameter"
#' @param control_list Optimization control for the two-parameter REML
#' @param boundary Distance from 0 and 1 below which "auto" switches to two-parameter REML
#' @return List with estimate, se, hessian, code, iterations, method, fallback,
#'   hard_settings, initial_h2, initial_se
#' @keywords internal
.reml_fit <- function(Vz2, Uval, n, m, r, reml = c("auto", "profile", "two-parameter"),
                      control_list = NULL, boundary = 1e-3) {
    reml <- match.arg(reml)
    fallback <- FALSE
 
    if (reml != "two-parameter") {
        fit <- tryCatch(.reml_profile(Vz2, Uval, n, m, r), error = function(e) NULL)
        ok <- !is.null(fit) && all(is.finite(c(fit$estimate, fit$se)))
        near_boundary <- ok && (fit$estimate[["h2"]] < boundary ||
                                fit$estimate[["h2"]] > 1 - boundary)
        if (reml == "profile" || (ok && !near_boundary)) {
            if (!ok) {
                stop("Profile REML did not return finite estimates; ",
                     "try reml = 'two-parameter'.", call. = FALSE)
            }
            fit$method <- "profile"
            fit$fallback <- FALSE
            fit$hard_settings <- FALSE
            fit$initial_h2 <- NA_real_
            fit$initial_se <- NA_real_
            return(fit)
        }
        fallback <- TRUE
    }
 
    ini <- .LocalHeritability.ini(Vz2, Uval, n, m, r)
    fit <- .LocalHeritability(Vz2, Uval, n, m, r, ini$h2, control_list)
    hard <- FALSE
    if (abs(fit$estimate[2] - 1) >= 0.3) {
        hard_control <- list(factr = 1, pgtol = 1e-10, maxit = 10000, fnscale = -1, lmm = 20)
        fit <- .LocalHeritability(Vz2, Uval, n, m, r, ini$h2, hard_control)
        hard <- TRUE
    }
    fit$method <- "two-parameter"
    fit$fallback <- fallback
    fit$hard_settings <- hard
    fit$initial_h2 <- ini$h2
    fit$initial_se <- ini$se
    fit
}
