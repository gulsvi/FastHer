#' FastHer Example Analysis
#' 
#' This demo shows how to use FastHer for local heritability estimation
#' using GWAS summary statistics and LD matrices.

library(FastHer)

# Example usage of FastHer
example_analysis <- function() {
    cat("Running FastHer example analysis...\n")
    
    # Example file paths (replace with your actual paths)
    result <- FastHer(
        path.in = "data/",
        file.Z = "ukb_200k_chr22_Z_BPIFC",
        file.POS = "ukb_200k_chr22_SNPorder_BPIFC", 
        file.LD = "ukb_200k_chr22_LD_BPIFC",
        path.out = "Results/",
        results_file = "results_BPIFC.txt"
    )
    
    cat("\n=== FastHer Results ===\n")
    print(result)
    
    return(result)
}

# Run the example
if (sys.nframe() == 0) {
    example_analysis()
}