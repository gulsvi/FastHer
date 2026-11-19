###################################################################
###          Check and Install FastHer Dependencies            ###
###################################################################

#' Check and Install FastHer Dependencies
#' @keywords internal

.check_fasther_dependencies <- function() {
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

