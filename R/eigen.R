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


