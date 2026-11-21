test_that("profile and two-parameter REML agree", {
    a <- make_locus()
    p <- fasther_fit(a$z, a$n, a$ld, reml = "profile", eigen_method = "R")
    t2 <- fasther_fit(a$z, a$n, a$ld, reml = "two-parameter", eigen_method = "R")
    expect_s3_class(p, "fasther")
    expect_equal(p$method, "profile")
    expect_equal(t2$method, "two-parameter")
    expect_lt(abs(p$h2 - t2$h2), 1e-3)
    expect_lt(abs(p$d - t2$d), 1e-2)
    expect_equal(p$wald_d, (p$d - 1) / p$se_d)
    expect_true(p$p_d >= 0 && p$p_d <= 1)
})
 
test_that("auto mode returns a valid method", {
    a <- make_locus()
    f <- fasther_fit(a$z, a$n, a$ld, eigen_method = "R")
    expect_true(f$method %in% c("profile", "two-parameter"))
})
 
test_that("sample sizes are harmonized", {
    a <- make_locus()
    n <- c(rep(1e5, 30), rep(2e5, 30))
    f <- fasther_fit(a$z, n, a$ld, eigen_method = "R")
    expect_equal(f$n_eff, 1 / mean(1 / n))
})
 
test_that("unit weights give the unweighted result", {
    a <- make_locus()
    f0 <- fasther_fit(a$z, a$n, a$ld, eigen_method = "R")
    f1 <- fasther_fit(a$z, a$n, a$ld, weights = rep(1, 60), eigen_method = "R")
    expect_false(f0$weighted)
    expect_true(f1$weighted)
    expect_equal(f1$h2, f0$h2, tolerance = 1e-8)
})
 
test_that("input errors are reported clearly", {
    a <- make_locus()
    expect_error(fasther_fit(a$z[-1], a$n, a$ld), "square matrix")
    expect_error(fasther_fit(replace(a$z, 1, NA), a$n, a$ld), "without missing")
    expect_error(fasther_fit(a$z, -1, a$ld), "positive")
    expect_error(fasther_fit(a$z, a$n, a$ld, tau = 1.5), "tau")
    expect_error(fasther_fit(a$z, a$n, a$ld, weights = rep(1, 10), eigen_method = "R"), "finite positive")
    expect_error(fasther_fit(a$z, a$n, a$ld, regularize = TRUE, reg_tau = 1.2), "reg_tau")
    expect_error(fasther_fit(a$z, a$n, a$ld, regularize = TRUE, reg_control = list(foo = 1)), "reg_control")
})
 
test_that("regularization works on a small locus", {
    a <- make_locus()
    f <- fasther_fit(a$z, a$n, a$ld, regularize = TRUE, eigen_method = "R")
    expect_true(f$regularized)
    expect_gte(f$rank, 6)
    expect_true(is.data.frame(f$reg_path))
})
 
test_that("print and as.data.frame work", {
    a <- make_locus()
    f <- fasther_fit(a$z, a$n, a$ld, eigen_method = "R")
    expect_output(print(f), "FastHer: local heritability")
    df <- as.data.frame(f)
    expect_equal(nrow(df), 1)
    expect_true(all(c("FH_h2", "FH_se", "FH_d", "FH_sed", "wald_d", "p_d", "method") %in% names(df)))
})
