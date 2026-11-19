###################################################################
###                     Main FastHer Function                   ###
###################################################################
#' FastHer: Fast and Accurate Local Heritability REML Estimation
#' 
#' @param path.in Path to input data files
#' @param file.Z Z-statistics file name (without extension)
#' @param file.POS POS file name (without extension) 
#' @param file.LD LD matrix file name (without extension)
#' @param path.out Path to output directory (optional, default: "Results/")
#' @param results_file Output file name (optional, default: "FastHer_results.txt")
#' @param control_list Optimization control parameters (optional)
#' @param check_dependencies Logical, whether to check for required packages (default: FALSE)
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
        .check_fasther_dependencies()
    }
    
    library(reticulate)
    library(data.table)
    
    # Set up parallel processing
    max_cores <- parallel::detectCores()
    optimal_threads <- max(1, max_cores - 2)
    data.table::setDTthreads(optimal_threads)
    Sys.setenv(OMP_NUM_THREADS = optimal_threads)
    
    # Load data
    SNP <- read.table(paste0(path.in,"/", file.Z, ".txt"), header = TRUE, stringsAsFactors = FALSE)
    POS <- read.table(paste0(path.in,"/", file.POS, ".bim"), header = FALSE, stringsAsFactors = FALSE)
    np <- reticulate::import("numpy")
    data <- np$load(paste0(path.in,"/", file.LD, ".npz"))
    
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
    eigU <- .eigen_decomp(U)
    
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
    res <- .LocalHeritability.ini(Vz2, Uval, n, m, rankU)
    initial_h2 <- res$h2
    initial_se <- res$se
    
    # Main estimation
    myres <- .LocalHeritability(Vz2, Uval, n, m, rankU, initial_h2, control_list)
    hard_settings <- FALSE

    if(abs(myres$estimate[2] - 1) >= 0.3) {
        hard.control_list <- list(factr = 1, pgtol = 1e-10, maxit = 10000, fnscale = -1, lmm = 20)
        myres <- .LocalHeritability(Vz2, Uval, n, m, rankU, initial_h2, hard.control_list)
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