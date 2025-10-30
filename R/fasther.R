#' FastHer: Fast and Accurate Local Heritability REML-Estimator using
#' GWAS summary statistics and LD matrix
#' 
#' @author Gulnara R. Svishcheva
#' @description Implements the FastHer method for local heritability estimation
#' using REML function optimized by eigen-decomposition of LD matrix

###################################################################
###                     Eigenvalue Decomposition                ###
###################################################################

#' Eigen-decomposition Using Python's Scipy
#' 
#' @param U A symmetric correlation SNP-by-SNP matrix (LD matrix)
#' @param only.values Logical, whether to return only eigenvalues
#' @return List with eigenvalues and eigenvectors
#' @importFrom reticulate import
#' @export

eigen.Python <- function(U, only.values = FALSE) {
    U <- as.matrix(U)
    if (!only.values) {
        scipy <- reticulate::import("scipy.linalg")
        eigen_result <- scipy$eigh(U)
        Uval <- eigen_result[[1]]
        Uvec <- eigen_result[[2]]
        Uval <- rev(Uval)
        Uvec <- Uvec[, length(Uval):1]
        return(list(values = Uval, vectors = Uvec))
    } else {
        scipy <- reticulate::import("scipy.linalg")
        eigen_result <- scipy$eigh(U, eigvals_only = TRUE)
        Uval <- rev(as.numeric(eigen_result))
        return(list(values = Uval, vectors = NA))
    }
}


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

LocalHeritability <- function(Vz2, Uval, n, m, r, start, control_list = NULL) { 
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
    control_list <- validate_control_list(control_list)
    
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

validate_control_list <- function(control_list) {
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

LocalHeritability.ini <- function(Vz2, Uval, n, m, r) {
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

###################################################################
###          Check and Install FastHer Dependencies            ###
###################################################################

#' Check and Install FastHer Dependencies
#' @keywords internal

check_fasther_dependencies <- function() {
    required_r_packages <- c("reticulate", "data.table", "compiler", "Matrix")
    missing_r <- required_r_packages[!required_r_packages %in% installed.packages()[,"Package"]]
    
    if (length(missing_r) > 0) {
        stop("Missing R packages: ", paste(missing_r, collapse = ", "),
             "\nPlease run: install.packages(c('", paste(missing_r, collapse = "', '"), "'))")
    }
    
    library(reticulate)
    use_virtualenv("r-reticulate", required = TRUE)
    cat("Python is available. Checking packages...\n")
    required_py_packages <- c("scipy", "numpy")
    for (pkg in required_py_packages) {
        cat("Checking", pkg, "...")
        if (!py_module_available(pkg)) {
            cat(" MISSING\n")
            stop("Python package '", pkg, "' is not available.",
                 "\nPlease run: reticulate::py_install('", pkg, "')")
        } else {
            cat(" OK\n")
        }
    }
    
    cat("All dependencies are available!\n")
    return(TRUE)
}

check_fasther_dependencies()

###################################################################
###                     Main FastHer Function                   ###
###################################################################

#' FastHer: Fast and Accurancy Local Heritability REML Estimation
#' 
#' @param path Path to data files
#' @param file.Z Z-statistics file name (without extension)
#' @param file.POS POS file name (without extension) 
#' @param file.LD LD matrix file name (without extension)
#' @param results_file Output file name (optional)
#' @param control_list Optimization control parameters (optional)
#' @param check_dependencies Logical, whether to check for required packages (default: TRUE)
#' @return Data frame with estimation results
#' @export
#' @examples
#' \dontrun{
#' # Example with default parameters:
#' result <- FastHer(path.in = "data/", 
#'                  file.Z = "z_statistics",
#'                  file.POS = "snp_info", 
#'                  file.LD = "ld_matrix")
#'
#' # Example with custom optimization parameters:
#' result <- FastHer(path.in = "data/",
#'                  file.Z = "z_statistics", 
#'                  file.POS = "snp_info",
#'                  file.LD = "ld_matrix",
#'                  control_list = list(factr = 1e7, maxit = 5000))
#' }

FastHer <- function(path.in, file.Z, file.POS, file.LD, path.out = NULL, results_file = NULL, 
                   control_list = NULL, check_dependencies = FALSE) {
    start_time <- Sys.time()
    
    # Input validation
    if (missing(path.in)) stop("Parameter 'path.in' must be specified")
	
    if (is.null(path.out)) path.out <- "Results/"
    if (!dir.exists(path.out)) {
        dir.create(path.out)
    }
    if (is.null(results_file)) results_file <- "FastHer_results.txt"
    
    # Check dependencies if requested
    if (check_dependencies) {
        check_fasther_dependencies()
    }
    
    library(reticulate)
    library(data.table)
    
    # Set up parallel processing
    max_cores <- parallel::detectCores()
    optimal_threads <- max(1, max_cores - 2)
    data.table::setDTthreads(optimal_threads)
    Sys.setenv(OMP_NUM_THREADS = optimal_threads)
    
    # Load data
    SNP <- read.table(paste0(path.in, file.Z, ".txt"), header = TRUE, stringsAsFactors = FALSE)
    POS <- read.table(paste0(path.in, file.POS, ".bim"), header = FALSE, stringsAsFactors = FALSE)
    np <- reticulate::import("numpy")
    data <- np$load(paste0(path.in, file.LD, ".npz"))
    
    # Data preparation
    z <- SNP$Z
    n_vec <- SNP$n
    names(z) <- names(n_vec) <- SNP$Predictor
    
    LD_array <- data$f[["LD"]]
    U <- as.matrix(LD_array)
    m <- ncol(U)
    
    order_SNP <- POS[, 2]
    n_vec <- n_vec[order_SNP]
    z <- z[order_SNP]
    
    colnames(U) <- rownames(U) <- order_SNP
    
    # Harmonize sample sizes
    n <- 1 / mean(1 / n_vec)
    z <- z * sqrt(n / n_vec)
    
    # Eigen-decomposition and pruning
    if (m > 100) {
        eigU <- eigen.Python(U)
    } else {
        eigU <- eigen(U, symmetric = TRUE)
    }
    
    Uval <- eigU$values
    Uvec <- eigU$vectors
    Uval[Uval < 0] <- 0
    
    # Keep components explaining 99.999% of variance
    cum <- cumsum(Uval) / sum(Uval)
    threshold <- Uval[min(which(cum > 0.99999))]
    Uvec <- Uvec[, Uval >= threshold]
    Uval <- Uval[Uval >= threshold]
    rankU <- length(Uval)
    
    # Project Z-statistics
    Vz <- c(t(crossprod(z, Uvec)))
    Vz2 <- Vz^2
    
    # Initial estimation
    res <- LocalHeritability.ini(Vz2, Uval, n, m, rankU)
    initial_h2 <- res$h2
    initial_se <- res$se
    
    # Main estimation
    myres <- LocalHeritability(Vz2, Uval, n, m, rankU, initial_h2, control_list)
    hard_settings <- FALSE

    if(abs(myres$estimate[2] - 1) >= 0.3) {
        hard.control_list <- list(factr = 1, pgtol = 1e-10, maxit = 10000, fnscale = -1, lmm = 20)
        myres <- LocalHeritability(Vz2, Uval, n, m, rankU, initial_h2, hard.control_list)
        hard_settings <- TRUE
     }
    
    # Extract results
    FH_h2 <- myres$estimate[1]
    FH_se <- myres$se[1]
    FH_d <- myres$estimate[2]
    FH_sed <- myres$se[2]
    convergence <- myres$code
    iterations <- myres$iterations
    
    # Print convergence information
    if (myres$code == 0) {
        cat("Optimization converged successfully in", iterations, "iterations.\n")
    } else {
        warning("Optimization did not converge properly. Code: ", myres$code)
    }
    
    # Prepare results
    FH_Tab <- data.frame(
        initial_h2 = initial_h2,
        initial_se = initial_se,
        m = round(m, 1),
        rankU = round(rankU, 1),
        FH_h2 = FH_h2,
        FH_se = FH_se,
        FH_d = FH_d,
        FH_sed = FH_sed,
        convergence = myres$code,
        iterations = myres$iterations,
        hard_settings = hard_settings,
        time_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    )
    
    write.table(FH_h2, paste0(path.out, "FH_h2.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(FH_se, paste0(path.out, "FH_se.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(FH_Tab, paste0(path.out, results_file), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    return(FH_Tab)
}