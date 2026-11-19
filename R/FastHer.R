###################################################################
###        File-based interface (backward compatible)           ###
###################################################################
 
# Reads allele frequencies for SNP weighting, in the order of the .bim file.
# Looks for a column named MAF, AF, AF1, EAF or FREQ (case-insensitive) in the
# z-statistics file, or for the column 'maf_col' given by the user. If there is
# none, the 3rd column of the .bim file is used when it holds frequencies in (0, 1).
.read_maf <- function(SNP, POS, order_SNP, maf_col = NULL) {
    nm <- names(SNP)
    col <- maf_col
    if (is.null(col)) {
        hit <- match(c("maf", "af", "af1", "eaf", "freq"), tolower(nm), nomatch = 0L)
        hit <- hit[hit > 0L]
        if (length(hit) > 0L) col <- nm[hit[1]]
    }
    if (!is.null(col)) {
        if (!col %in% nm) {
            stop("Column '", col, "' not found in the z-statistics file.", call. = FALSE)
        }
        maf <- SNP[[col]]
        names(maf) <- SNP$Predictor
        return(maf[order_SNP])
    }
    maf <- suppressWarnings(as.numeric(POS[, 3]))
    if (length(maf) == nrow(POS) && all(is.finite(maf)) && all(maf > 0 & maf < 1)) {
        return(maf)
    }
    stop("SNP weights need allele frequencies: add a column MAF (or AF) to the z-statistics ",
         "file, name the column with 'maf_col', or put the frequencies in the 3rd column ",
         "of the .bim file.", call. = FALSE)
}
 
#' FastHer: Fast and Accurate Local Heritability REML Estimation
#'
#' File-based interface of \code{\link{fasther_fit}}: reads the z-statistics,
#' the SNP order file and the LD matrix of a locus, estimates local
#' heritability h2 and the concordance parameter d, and writes the results.
#'
#' @param path.in Path to input data files
#' @param file.Z Z-statistics file name (without extension; columns Predictor, Z, n)
#' @param file.POS SNP order file name (without extension; .bim format, SNP names in column 2)
#' @param file.LD LD matrix file name (without extension; .npz with key "LD")
#' @param path.out Path to output directory (optional, default: "Results/")
#' @param results_file Output file name (optional, default: "FastHer_results.txt")
#' @param control_list Optimization control parameters for the two-parameter REML (optional)
#' @param check_dependencies Logical, whether to check for required packages (default: FALSE)
#' @param reml REML algorithm: "auto" (default), "profile" or "two-parameter"
#'   (see \code{\link{fasther_fit}})
#' @param tau Fraction of the LD variance retained in the eigenvalue truncation
#'   (default 0.99999)
#' @param use_weights Logical, whether to use SNP weights computed from allele
#'   frequencies (default FALSE); see \code{\link{fasther_weights}}
#' @param power LDAK power for the SNP weights (default -0.25)
#' @param maf_col Name of the allele-frequency column in the z-statistics file
#'   (optional; by default a column MAF, AF, AF1, EAF or FREQ, or the 3rd
#'   column of the .bim file, is used)
#' @param regularize Logical, whether to regularize a reference LD matrix
#'   (default FALSE); see \code{\link{fasther_fit}}
#' @param reg_tau Starting fraction of the LD variance for the regularization
#'   (default 0.97)
#' @param reg_control Optional list of regularization settings (see
#'   \code{\link{fasther_fit}})
#' @return Data frame with estimation results (one row)
#' @export
#' @examples
#' \dontrun{
#' # Example with default parameters:
#' result <- FastHer(path.in = "data/",
#'                  file.Z = "z_statistics",
#'                  file.POS = "snp_info",
#'                  file.LD = "ld_matrix")
#'
#' # With SNP weights (allele frequencies in a column MAF of the z-statistics file)
#' # and regularization of a reference LD matrix:
#' result <- FastHer(path.in = "data/",
#'                  file.Z = "z_statistics",
#'                  file.POS = "snp_info",
#'                  file.LD = "ld_matrix",
#'                  use_weights = TRUE, regularize = TRUE)
#' }
FastHer <- function(path.in, file.Z, file.POS, file.LD, path.out = NULL, results_file = NULL,
                    control_list = NULL, check_dependencies = FALSE,
                    reml = c("auto", "profile", "two-parameter"), tau = 0.99999,
                    use_weights = FALSE, power = -0.25, maf_col = NULL,
                    regularize = FALSE, reg_tau = 0.97, reg_control = NULL) {
    start_time <- Sys.time()
    reml <- match.arg(reml)
 
    if (missing(path.in)) stop("Parameter 'path.in' must be specified")
    if (!is.logical(use_weights) || length(use_weights) != 1L || is.na(use_weights)) {
        stop("'use_weights' must be TRUE or FALSE.", call. = FALSE)
    }
    if (is.null(path.out)) path.out <- "Results/"
    if (!dir.exists(path.out)) dir.create(path.out)
    if (is.null(results_file)) results_file <- "FastHer_results.txt"
 
    if (check_dependencies) .check_fasther_dependencies()
 
    ## load data
    SNP <- utils::read.table(paste0(path.in, "/", file.Z, ".txt"), header = TRUE, stringsAsFactors = FALSE)
    POS <- utils::read.table(paste0(path.in, "/", file.POS, ".bim"), header = FALSE, stringsAsFactors = FALSE)
    np <- reticulate::import("numpy")
    data <- np$load(paste0(path.in, "/", file.LD, ".npz"))
 
    z <- SNP$Z
    n_vec <- SNP$n
    names(z) <- names(n_vec) <- SNP$Predictor
    U <- as.matrix(data$f[["LD"]])
 
    ## consistency of the SNPs in the three files
    if (anyDuplicated(POS[, 2]) > 0L) stop("The .bim file has duplicated SNP names.", call. = FALSE)
    if (anyDuplicated(SNP$Predictor) > 0L) {
        stop("The z-statistics file has duplicated SNP names (column Predictor).", call. = FALSE)
    }
    miss <- setdiff(POS[, 2], SNP$Predictor)
    if (length(miss) > 0L) {
        stop(length(miss), " SNPs of the .bim file are not in the z-statistics file (e.g. ",
             paste(utils::head(miss, 3), collapse = ", "), ").", call. = FALSE)
    }
    if (nrow(U) != nrow(POS) || ncol(U) != nrow(POS)) {
        stop("The LD matrix is ", nrow(U), " x ", ncol(U), " but the .bim file has ", nrow(POS),
             " SNPs.", call. = FALSE)
    }
 
    order_SNP <- POS[, 2]
    n_vec <- n_vec[order_SNP]
    z <- z[order_SNP]
 
    weights <- NULL
    if (use_weights) {
        weights <- fasther_weights(.read_maf(SNP, POS, order_SNP, maf_col), power = power)
    }
 
    res <- fasther_fit(z, n_vec, U, weights = weights, reml = reml, tau = tau,
                       control_list = control_list, regularize = regularize,
                       reg_tau = reg_tau, reg_control = reg_control)
    message("FastHer (", res$method, " REML): h2 = ", format(res$h2, digits = 4),
            ", SE = ", format(res$se_h2, digits = 4))
 
    FH_Tab <- as.data.frame(res)
    FH_Tab$time_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
 
    utils::write.table(FH_Tab$FH_h2, paste0(path.out, "FH_h2.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    utils::write.table(FH_Tab$FH_se, paste0(path.out, "FH_se.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    utils::write.table(FH_Tab, paste0(path.out, results_file), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
 
    return(FH_Tab)
}
