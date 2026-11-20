test_that("BPIFC locus reproduces FastHer 1.0.0", {
    skip_if_no_python()
    a <- read_bpifc()
    f <- fasther_fit(a$z, a$n, a$ld, reml = "two-parameter")
    expect_equal(f$m, 1032)
    expect_equal(f$rank, 951)
    expect_equal(f$h2, 0.02355058, tolerance = 1e-4)
    expect_equal(f$se_h2, 0.0014569, tolerance = 1e-4)
    expect_equal(f$d, 1.001212, tolerance = 1e-4)
    expect_equal(f$se_d, 0.003420417, tolerance = 1e-4)
    p <- fasther_fit(a$z, a$n, a$ld)
    expect_equal(p$method, "profile")
    expect_equal(p$h2, f$h2, tolerance = 1e-6)
})
 
test_that("simulated locus with allele frequencies: weighted and unweighted REML", {
    skip_if_no_python()
    a <- read_sim5k()
    f0 <- fasther_fit(a$z, a$n, a$ld)
    f1 <- fasther_fit(a$z, a$n, a$ld, weights = fasther_weights(a$maf))
    expect_equal(f0$m, 617)
    expect_equal(f0$rank, 499)
    expect_equal(f0$h2, 0.08183056, tolerance = 1e-4)
    expect_equal(f0$se_h2, 0.01319207, tolerance = 1e-4)
    expect_equal(f1$h2, 0.09179231, tolerance = 1e-4)
    expect_equal(f1$se_h2, 0.01533379, tolerance = 1e-4)
    expect_equal(f1$d, 0.9819051, tolerance = 1e-4)
    expect_gt(f1$p_d, 0.05)
})
 
test_that("regularization drops components until the Wald test passes", {
    skip_if_no_python()
    a <- read_bpifc()
    m <- length(a$z)
    n <- 2e5
    h <- 0.025
    set.seed(2619)
    ev <- eigen(a$ld, symmetric = TRUE)
    lam <- pmax(ev$values, 0)
    z <- c(ev$vectors %*% (sqrt(n * h / m * lam^2 + (1 - h) * lam) * rnorm(m)))
    E <- matrix(rnorm(m * m, sd = 0.03), m)
    E <- (E + t(E)) / 2
    diag(E) <- 0
    ref <- a$ld + E
    f <- fasther_fit(z, n, ref, regularize = TRUE, reg_tau = 0.999, reml = "two-parameter")
    expect_true(f$regularized)
    expect_equal(nrow(f$reg_path), 10)
    expect_equal(f$rank, 747)
    expect_equal(f$h2, 0.029663, tolerance = 1e-4)
    expect_true(f$reg_converged)
    expect_lte(abs(f$wald_d), 1.96)
    expect_gt(abs(f$reg_path$wald[1]), 1.96)
    expect_true(all(diff(f$reg_path$rank) == -1))
})
