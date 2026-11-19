###################################################################
###                 Heritability Estimation Core                ###
###################################################################

#' Core Heritability Estimation Function (2 parameters)
#' 
#' @param Vz2 Squared projections of Z-statistics
#' @param Uval Eigenvalues of LD matrix
#' @param n Sample size
#' @param m Number of SNPs
#' @param r Rank after eigenvalue pruning
#' @param start Initial heritability value
#' @param control_list Optimization control parameters (optional)
#' @return List with estimation results
#' @keywords internal

.LocalHeritability <- function(Vz2, Uval, n, m, r, start, control_list = NULL) { 
    eps <- 1e-9
    Vz2dUval <- Vz2 / Uval
    zUz <- sum(Vz2dUval)
    dd <- Uval * (n / m)
    dd_1 <- dd - 1
    n_r <- n - r
    n_zUz <- n - zUz
    Dob <- n * log(2 * pi)
    
    # Gradient function
    Grad <- function(parm) {
        h2 <- parm[1]
        d <- parm[2]
        Wr <- dd_1 * h2 + 1
        K <- dd_1 / Wr
        common_expr <- n_r - n_zUz / (d * (1 - h2))
        Vz2dUvaldWr <- Vz2dUval / Wr
        part1_h2 <- -0.5 * (sum(K) - sum(K * Vz2dUvaldWr) / d)
        part2_h2 <- (0.5 / (1 - h2)) * common_expr
        part1_d <- (-0.5 / d) * (r - sum(Vz2dUvaldWr) / d)
        part2_d <- (-0.5 / d) * common_expr
        grad <- c(h2 = part1_h2 + part2_h2, d = part1_d + part2_d)
        return(grad)
    }
    
    # Hessian function (expected)
    Hess <- function(parm) {
        h2 <- parm[1]
        d <- parm[2]
        Wr <- dd_1 * h2 + 1
        K <- dd_1 / Wr
        SD.h2 <- (-0.5) * (sum(K^2) + n_r / ((1 - h2)^2))
        SD.d <- (-0.5 * n) / (d^2)
        SD.h2.d <- (-0.5 / d) * sum(K) + 0.5 * n_r / (d * (1 - h2))
        hess <- matrix(c(SD.h2, SD.h2.d, SD.h2.d, SD.d), nrow = 2, ncol = 2, dimnames = list(c("h2", "d"), c("h2", "d")))
        return(hess)
    }
    
    # Log-likelihood function
    Lh <- function(parm) {
        h2 <- parm[1]
        d <- parm[2]
        sigE <- (1 - h2) * d
        Wr <- dd * h2 + (1 - h2)
        delta <- n_r * log(1 - h2) + n * log(d) + (n_zUz) / sigE
        logLH <- sum(log(Wr)) + sum(Vz2dUval / Wr) / d + delta + Dob
        -0.5 * logLH  # Convert to maximization problem
    }
    
    # Validate control_list parameters
    control_list <- .validate_control_list(control_list)
    
    fit0 <- stats::optim(par = c(h2 = start, d = 1), fn = Lh, gr = Grad, method = "L-BFGS-B", 
                  lower = c(eps, eps), upper = c(1 - eps, 100), hessian = TRUE, control = control_list)
    
    # Calculate standard errors
    hessian <- Hess(fit0$par)
    se <- sqrt(diag(solve(-hessian)))
    
    fit <- list(
        maximum = fit0$value,
        estimate = fit0$par,
        se = se,
        gradient = Grad(fit0$par),
        hessian = hessian,
        code = fit0$convergence,
        message = fit0$message,
        type = "maximization",
        iterations = fit0$counts[1],
        control_used = control_list
    )
    return(fit)
}


###################################################################
###    Validate and Complete Control List for Optimization      ###
###################################################################

#' Validate and Complete Control List for Optimization
#' 
#' @param control_list User-provided control parameters
#' @return Validated control list with default values for missing parameters
#' @keywords internal

.validate_control_list <- function(control_list) {
    default_control <- list(factr = 1e7, pgtol = 1e-7, maxit = 1000,  fnscale = -1, lmm = 10)
    
    if (is.null(control_list)) {
        return(default_control)
    }
    
    # Replace only those parameters that are provided by the user
    for (param in names(control_list)) {
        if (param %in% names(default_control)) {
            default_control[[param]] <- control_list[[param]]
        }
    }
    
    return(default_control)
}

###################################################################
###                 Initial Heritability Estimation             ###
###################################################################

#' Initial Heritability Estimation (1 parameter)
#' 
#' @param Vz2 Squared projections of Z-statistics
#' @param Uval Eigenvalues of LD matrix
#' @param n Sample size
#' @param m Number of SNPs in genomic region
#' @param r Rank after eigenvalue pruning
#' @return List with initial h2 estimate and standard error
#' @keywords internal

.LocalHeritability.ini <- function(Vz2, Uval, n, m, r) {
    Vz2dUval <- Vz2 / Uval
    zUz <- sum(Vz2dUval)
    dd <- Uval * (n / m)
    dd_1 <- dd - 1
    n_r <- n - r
    n_zUz <- n - zUz
    Dob <- n * log(2 * pi)
    
    Lh <- function(h2) {
        sigE <- (1 - h2)
        Wval <- dd_1 * h2 + 1
        delta <- (n_r) * log(sigE) + (n_zUz) / sigE
        logLH <- sum(log(Wval)) + sum(Vz2dUval / Wval) + delta + Dob
        -0.5 * logLH
    }
    
    Hessi.exp <- function(h2) {
        Wval <- dd_1 * h2 + 1
        Part1 <- -0.5 * sum((dd_1 / Wval)^2)
        Part2 <- -0.5 * n_r / ((1 - h2)^2)
        sqrt(-1 / (Part1 + Part2))
    }
    
    B <- stats::optimize(f = Lh, interval = c(0, 1), maximum = TRUE)
    h2 <- B$maximum
    
    return(list(h2 = h2, se = Hessi.exp(h2)))
}

