library(FastHer)
 
cat("=== FastHer demo: SNP weights ===\n")
cat("Simulated locus (1000G chr22, N = 5000, LDAK model, power -0.25, true h2 = 0.1)\n\n")
 
data_path <- system.file("extdata", package = "FastHer")
if (data_path == "" && dir.exists("inst/extdata")) data_path <- "inst/extdata"
files <- paste0("sim5k_chr22_1KG_", c("Z", "SNPorder", "LD"))
out <- paste0(file.path(tempdir(), "FastHer_results"), "/")
dir.create(out, showWarnings = FALSE)
 
# Without weights
res0 <- FastHer(data_path, files[1], files[2], files[3], path.out = out,
                results_file = "unweighted.txt")
# With SNP weights computed from the allele frequencies (column AF1 of the z file)
res1 <- FastHer(data_path, files[1], files[2], files[3], path.out = out,
                results_file = "weighted.txt", use_weights = TRUE)
 
cols <- c("m", "rankU", "weighted", "FH_h2", "FH_se", "FH_d", "FH_sed", "wald_d", "p_d")
print(rbind(unweighted = res0, weighted = res1)[, cols])
cat("\nDemo completed. Results were written to:", out, "\n")
