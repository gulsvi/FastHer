test_that("fasther_weights scales to mean square 1 and follows the LDAK exponent", {
    maf <- c(0.05, 0.1, 0.2, 0.3, 0.5)
    w <- fasther_weights(maf)
    expect_equal(mean(w^2), 1)
    raw <- (maf * (1 - maf))^0.375
    expect_equal(w, raw / sqrt(mean(raw^2)))
})
 
test_that("power = -1 gives equal weights", {
    expect_equal(fasther_weights(c(0.1, 0.3, 0.45), power = -1), rep(1, 3))
})
 
test_that("fasther_weights rejects invalid input", {
    expect_error(fasther_weights(c(0.2, 0, 0.3)), "strictly between")
    expect_error(fasther_weights(c(0.2, NA)), "strictly between")
    expect_error(fasther_weights(0.2, power = c(1, 2)), "single number")
})
