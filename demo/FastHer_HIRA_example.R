library(FastHer)

cat("=== FastHer Package Demo ===\n")
cat("Running real analysis with package data...\n\n")

# Get the path to package data
pkg_path <- system.file(package = "FastHer")
data_path <- file.path(pkg_path, "data")

cat("Package path:", pkg_path, "\n")
cat("Data path:", data_path, "\n")

# Check if data directory exists
if (!dir.exists(data_path)) {
    cat("Data directory not found at expected location.\n")
    cat("Trying alternative locations...\n")
    # Try current directory for development
    if (dir.exists("data")) {
        data_path <- "data"
        cat("Using local data directory\n")
    } else {
        cat("No data directory found. Demo cannot run.\n")
        return()
    }
}

# List available data files
data_files <- list.files(data_path)
cat("Available data files in", data_path, ":\n")
cat(paste("-", data_files, collapse = "\n"), "\n\n")

# Check if required files exist
required_files <- c(
    "ukb_200k_chr22_Z_HIRA.txt",
    "ukb_200k_chr22_SNPorder_HIRA.bim", 
    "ukb_200k_chr22_LD_HIRA.npz"
)

missing_files <- required_files[!file.exists(file.path(data_path, required_files))]
if (length(missing_files) > 0) {
    cat("Missing required files:", paste(missing_files, collapse = ", "), "\n")
    cat("Available files:", paste(list.files(data_path), collapse = ", "), "\n")
} else {
    cat("All required data files found! Running analysis...\n\n")
    
    # Create results directory
    if (!dir.exists("Results")) {
        dir.create("Results")
    }
    
    # Run analysis
    result <- FastHer(
        path.in = data_path,
        file.Z = "ukb_200k_chr22_Z_HIRA",
        file.POS = "ukb_200k_chr22_SNPorder_HIRA", 
        file.LD = "ukb_200k_chr22_LD_HIRA",
        path.out = "Results/",
        results_file = "demo_results_HIRA.txt",
        check_dependencies = TRUE
    )
    
    cat("\n=== FastHer Results ===\n")
    print(result)
    cat("\n Analysis completed successfully!\n")
    cat(" Results saved to: Results/demo_results_HIRA.txt\n")
}

cat("\nDemo completed!\n")
