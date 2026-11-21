test_that("Python and R eigen-decompositions agree", {
    skip_if_not_installed("reticulate")
    skip_if_not(reticulate::py_module_available("scipy"), "scipy is not available")
    U <- 0.5^abs(outer(1:30, 1:30, "-"))
    e1 <- .eigen_decomp(U, method = "python")
    e2 <- .eigen_decomp(U, method = "R")
    expect_equal(e1$values, e2$values, tolerance = 1e-8)
    expect_equal(e1$vectors %*% (e1$values * t(e1$vectors)), U, tolerance = 1e-8)
})
