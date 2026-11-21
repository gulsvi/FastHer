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
###        Eigen-decomposition: Python for m > 100, else R      ###
###################################################################
 
#' Check that Python module scipy is available
#' @keywords internal
.python_scipy_available <- function() {
    requireNamespace("reticulate", quietly = TRUE) &&
        isTRUE(tryCatch(reticulate::py_module_available("scipy"),
                        error = function(e) FALSE))
}
 
#' Eigen-decomposition of an LD matrix
#'
#' Default rule (as in FastHer 1.0.0): Python/scipy for regions with more than
#' 100 SNPs, base R otherwise. If scipy is not available, base R is used with a
#' warning (the result is the same, only the speed may differ).
#'
#' @param U Symmetric SNP-by-SNP correlation (LD) matrix
#' @param only.values Logical, return only eigenvalues
#' @param method "auto" (default), "python" or "R"
#' @return List with `values` (decreasing) and `vectors`
#' @keywords internal
.eigen_decomp <- function(U, only.values = FALSE, method = c("auto", "python", "R")) {
    method <- match.arg(method)
    U <- as.matrix(U)
    use_python <- switch(method, auto = ncol(U) > 100, python = TRUE, R = FALSE)
    if (use_python && !.python_scipy_available()) {
        warning("Python module 'scipy' was not found: the eigen-decomposition is done in R ",
                "and may be slower for large regions. To use Python: reticulate::py_install('scipy').",
                call. = FALSE)
        use_python <- FALSE
    }
    if (use_python) {
        res <- eigen.Python(U, only.values = only.values)
        return(list(values = as.numeric(res$values),
                    vectors = if (only.values) NA else as.matrix(res$vectors)))
    }
    e <- eigen(U, symmetric = TRUE, only.values = only.values)
    list(values = e$values, vectors = if (only.values) NA else e$vectors)
}
