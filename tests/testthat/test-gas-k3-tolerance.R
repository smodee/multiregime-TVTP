test_that("GAS model K=3 filtering does not crash with tolerance error", {
  # Reproduces the bug from GitHub issue #9: K=3 GAS filtering crashed with
  # "X_pred must be valid probabilities that sum to 1" due to floating-point
  # drift exceeding the strict 1e-10 tolerance.

  set.seed(42)

  par <- c(-2, 0, 2, 1, 1, 1, 0.9, 0.85, 0.88,
           rep(0.05, 3), rep(0.9, 3))
  par <- set_parameter_attributes(par, K = 3, model_type = "gas",
                                   diag_probs = TRUE, equal_variances = FALSE)

  y <- c(rnorm(200, -2), rnorm(200, 0), rnorm(200, 2))

  expect_no_error(
    Rfiltering_GAS(par, y, B_burnin = 5, C = 5, use_fallback = FALSE)
  )
})

test_that("K=3 filtered probabilities are valid after re-normalization", {
  set.seed(42)

  par <- c(-2, 0, 2, 1, 1, 1, 0.9, 0.85, 0.88,
           rep(0.05, 3), rep(0.9, 3))
  par <- set_parameter_attributes(par, K = 3, model_type = "gas",
                                   diag_probs = TRUE, equal_variances = FALSE)

  y <- c(rnorm(200, -2), rnorm(200, 0), rnorm(200, 2))

  result <- Rfiltering_GAS(par, y, B_burnin = 5, C = 5,
                            use_fallback = FALSE, diagnostics = TRUE)
  X_t <- attr(result, "X.t")

  # All filtered probabilities should be positive

  expect_true(all(X_t > 0))

  # Each column should sum to 1 (within reasonable tolerance)
  col_sums <- colSums(X_t)
  expect_true(all(abs(col_sums - 1) < 1e-10))
})

test_that("GAS model K=3 simulation does not crash", {
  set.seed(42)

  par <- c(-2, 0, 2, 1, 1, 1, 0.9, 0.85, 0.88,
           rep(0.05, 3), rep(0.9, 3))
  par <- set_parameter_attributes(par, K = 3, model_type = "gas",
                                   diag_probs = TRUE, equal_variances = FALSE)

  expect_no_error(
    dataGASCD(M = 2, N = 500, par = par)
  )
})
