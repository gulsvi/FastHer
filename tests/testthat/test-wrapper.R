test_that("FastHer() writes the result files and returns the table", {
    skip_if_no_python()
    out <- paste0(tempfile("fh_"), "/")
    dir.create(out)
    res <- suppressMessages(FastHer(example_dir(), "ukb_200k_chr22_Z_BPIFC",
                                    "ukb_200k_chr22_SNPorder_BPIFC", "ukb_200k_chr22_LD_BPIFC",
                                    path.out = out))
    expect_true(all(c("FastHer_results.txt", "FH_h2.txt", "FH_se.txt") %in% list.files(out)))
    expect_equal(res$FH_h2, 0.02355058, tolerance = 1e-4)
    expect_false(res$weighted)
})
 
test_that("FastHer() takes allele frequencies from the z-statistics file", {
    skip_if_no_python()
    out <- paste0(tempfile("fh_"), "/")
    dir.create(out)
    res <- suppressMessages(FastHer(example_dir(), "sim5k_chr22_1KG_Z", "sim5k_chr22_1KG_SNPorder",
                                    "sim5k_chr22_1KG_LD", path.out = out, use_weights = TRUE))
    expect_true(res$weighted)
    expect_equal(res$FH_h2, 0.09179231, tolerance = 1e-4)
})
 
test_that("FastHer() explains missing allele frequencies", {
    skip_if_no_python()
    out <- paste0(tempfile("fh_"), "/")
    dir.create(out)
    expect_error(suppressMessages(FastHer(example_dir(), "ukb_200k_chr22_Z_BPIFC",
                                          "ukb_200k_chr22_SNPorder_BPIFC", "ukb_200k_chr22_LD_BPIFC",
                                          path.out = out, use_weights = TRUE)),
                 "allele frequencies")
})
