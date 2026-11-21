test_that("SNP names of z, ld and weights are compared when they are available", {
    a <- make_locus()
    nm <- paste0("rs", 1:60)
    z <- stats::setNames(a$z, nm)
    ld <- a$ld
    dimnames(ld) <- list(nm, nm)
    f <- fasther_fit(z, a$n, ld, eigen_method = "R")
    expect_s3_class(f, "fasther")
    expect_error(fasther_fit(z, a$n, ld[rev(nm), rev(nm)], eigen_method = "R"), "do not match")
    ld2 <- ld
    colnames(ld2) <- rev(nm)
    expect_error(fasther_fit(z, a$n, ld2, eigen_method = "R"), "Row names and column names")
    expect_error(fasther_fit(stats::setNames(a$z, rep("rs1", 60)), a$n, ld, eigen_method = "R"), "duplicated")
    w <- stats::setNames(rep(1, 60), rev(nm))
    expect_error(fasther_fit(z, a$n, ld, weights = w, eigen_method = "R"), "weights")
})
 
test_that("FastHer() checks that the files describe the same SNPs", {
    skip_if_no_python()
    stems <- c("ukb_200k_chr22_Z_BPIFC", "ukb_200k_chr22_SNPorder_BPIFC", "ukb_200k_chr22_LD_BPIFC")
    d <- example_dir()
    out <- paste0(tempfile("fh_"), "/")
    dir.create(out)
    tmp <- tempfile("fh_in")
    dir.create(tmp)
    file.copy(file.path(d, paste0(stems[3], ".npz")), tmp)
    SNP <- utils::read.table(file.path(d, paste0(stems[1], ".txt")), header = TRUE, stringsAsFactors = FALSE)
    POS <- utils::read.table(file.path(d, paste0(stems[2], ".bim")), header = FALSE, stringsAsFactors = FALSE)
    ## a SNP of the .bim file is missing in the z-statistics file
    utils::write.table(SNP[-1, ], file.path(tmp, paste0(stems[1], ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(POS, file.path(tmp, paste0(stems[2], ".bim")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    expect_error(suppressMessages(FastHer(tmp, stems[1], stems[2], stems[3], path.out = out)),
                 "not in the z-statistics file")
    ## the .bim file has fewer SNPs than the LD matrix
    utils::write.table(SNP, file.path(tmp, paste0(stems[1], ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(POS[-1, ], file.path(tmp, paste0(stems[2], ".bim")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    expect_error(suppressMessages(FastHer(tmp, stems[1], stems[2], stems[3], path.out = out)),
                 "but the .bim file has")
})
