# Tests for transform_parameters / untransform_parameters

test_that("diagonal transform round-trip preserves parameters", {
  par <- c(-1, 1, 0.5, 0.6, 0.8, 0.9)
  par <- set_parameter_attributes(par, K = 2, model_type = "constant",
                                  diag_probs = TRUE, equal_variances = FALSE)

  par_t <- transform_parameters(par)
  par_back <- untransform_parameters(par_t)

  expect_equal(as.numeric(par_back), as.numeric(par), tolerance = 1e-10)
})

test_that("off-diagonal transform is identity for transition params", {
  par <- c(-1, 1, 0.5, 0.6, -0.5, 0.3)
  par <- set_parameter_attributes(par, K = 2, model_type = "constant",
                                  diag_probs = FALSE, equal_variances = FALSE)

  par_t <- transform_parameters(par)

  # mu unchanged
  expect_equal(par_t[1:2], par[1:2])
  # sigma2 log-transformed
  expect_equal(par_t[3:4], log(par[3:4]), tolerance = 1e-12)
  # transition params unchanged (identity)
  expect_equal(par_t[5:6], par[5:6])
})

test_that("off-diagonal transform round-trip preserves parameters", {
  par <- c(-1, 1, 0.5, 0.6, -0.5, 0.3)
  par <- set_parameter_attributes(par, K = 2, model_type = "constant",
                                  diag_probs = FALSE, equal_variances = FALSE)

  par_t <- transform_parameters(par)
  par_back <- untransform_parameters(par_t)

  expect_equal(as.numeric(par_back), as.numeric(par), tolerance = 1e-10)
})

test_that("off-diagonal transform round-trip for TVP model", {
  # mu1, mu2, sigma2_1, sigma2_2, x12, x21, A1, A2
  par <- c(-1, 1, 0.5, 0.8, -0.3, 0.2, 0.1, -0.1)
  par <- set_parameter_attributes(par, K = 2, model_type = "tvp",
                                  diag_probs = FALSE, equal_variances = FALSE)

  par_t <- transform_parameters(par)
  par_back <- untransform_parameters(par_t)

  expect_equal(as.numeric(par_back), as.numeric(par), tolerance = 1e-10)
})

test_that("off-diagonal transform round-trip for GAS model", {
  # mu1, mu2, sigma2_1, sigma2_2, x12, x21, A1, A2, B1, B2
  par <- c(-1, 1, 0.5, 0.8, -0.3, 0.2, 0.1, -0.1, 0.5, 0.6)
  par <- set_parameter_attributes(par, K = 2, model_type = "gas",
                                  diag_probs = FALSE, equal_variances = FALSE)

  par_t <- transform_parameters(par)
  par_back <- untransform_parameters(par_t)

  expect_equal(as.numeric(par_back), as.numeric(par), tolerance = 1e-10)
})

test_that("diagonal transform rejects out-of-bound probabilities", {
  par <- c(0, 1, 0.5, 0.6, 0.0, 0.9)  # p11 = 0 is invalid

  par <- set_parameter_attributes(par, K = 2, model_type = "constant",
                                  diag_probs = TRUE, equal_variances = FALSE)
  expect_error(transform_parameters(par), "strictly between 0 and 1")
})

test_that("off-diagonal transform accepts any real values for transition params", {
  par <- c(-1, 1, 0.5, 0.6, -100, 100)
  par <- set_parameter_attributes(par, K = 2, model_type = "constant",
                                  diag_probs = FALSE, equal_variances = FALSE)
  # Should not error - softmax params are unconstrained
  par_t <- transform_parameters(par)
  expect_true(all(is.finite(par_t)))
})
