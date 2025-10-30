#!/usr/bin/env Rscript
# Example script for FastHer
# Load FastHer

setwd("..")
source("R/fasther.R")

# Example usage
example_analysis <- function() {
    cat("Running FastHer example...\n")
    
    # Replace with your actual file paths
    result <- FastHer(
        path.in = "data/",
        file.Z = "ukb_200k_chr22_Z_BPIFC",
        file.POS = "ukb_200k_chr22_SNPorder_BPIFC", 
        file.LD = "ukb_200k_chr22_LD_BPIFC",
        path.out = "Results/",
        results_file = "results_BPIFC.txt"
    )
    
    cat("\n=== FastHer Results ===\n")
#   print(result)
    
    return(result)
}

# Run example if called directly
if (sys.nframe() == 0) {
    example_analysis()
}