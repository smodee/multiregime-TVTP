# End-to-end model estimation tests
skip_on_cran()

# --- Constant model ---

test_that("constant model estimation works with off-diagonal parameterization", {
  set.seed(100)
  y <- c(rnorm(200, -1, 0.5), rnorm(200, 1, 0.7))

  result <- estimate_constant_model(y, K = 2, diag_probs = FALSE,
                                    n_starts = 5, verbose = 0)
  expect_equal(result$diagnostics$convergence_code, 0)
  expect_true(is.finite(result$diagnostics$neg_log_likelihood))
  # Regime means should be roughly -1 and 1 (order may differ)
  sorted_mu <- sort(result$mu_est)
  expect_true(abs(sorted_mu[1] - (-1)) < 0.3)
  expect_true(abs(sorted_mu[2] - 1) < 0.3)
})

test_that("constant model estimation still works with diagonal parameterization", {
  set.seed(100)
  y <- c(rnorm(200, -1, 0.5), rnorm(200, 1, 0.7))

  result <- estimate_constant_model(y, K = 2, diag_probs = TRUE,
                                    n_starts = 5, verbose = 0)
  expect_equal(result$diagnostics$convergence_code, 0)
  expect_true(is.finite(result$diagnostics$neg_log_likelihood))
})

# --- Data simulation ---

test_that("dataTVPCD generates valid data with off-diagonal params", {
  set.seed(42)
  # mu1, mu2, sigma2, x12, x21, A1, A2
  par <- c(-1, 1, 0.5, -0.5, 0.3, 0.05, -0.05)
  par <- set_parameter_attributes(par, K = 2, model_type = "tvp",
                                  diag_probs = FALSE, equal_variances = TRUE)

  data <- dataTVPCD(5, 500, par)
  expect_true(is.matrix(data))
  expect_equal(nrow(data), 5)
  expect_equal(ncol(data), 500)
  expect_true(all(is.finite(data)))
})

test_that("dataTVPCD generates valid data with diagonal params", {
  set.seed(42)
  # mu1, mu2, sigma2, p11, p22, A1, A2
  par <- c(-1, 1, 0.5, 0.8, 0.9, 0.05, -0.05)
  par <- set_parameter_attributes(par, K = 2, model_type = "tvp",
                                  diag_probs = TRUE, equal_variances = TRUE)

  data <- dataTVPCD(5, 500, par)
  expect_true(is.matrix(data))
  expect_equal(nrow(data), 5)
  expect_equal(ncol(data), 500)
  expect_true(all(is.finite(data)))
})

test_that("dataTVPXExoCD generates valid data with off-diagonal params", {
  set.seed(42)
  # mu1, mu2, sigma2, x12, x21, A1, A2
  par <- c(-1, 1, 0.5, -0.3, 0.2, 0.1, -0.1)
  par <- set_parameter_attributes(par, K = 2, model_type = "exogenous",
                                  diag_probs = FALSE, equal_variances = TRUE)
  X_Exo <- rnorm(1200)

  data <- dataTVPXExoCD(5, 500, par, X_Exo)
  expect_true(is.matrix(data))
  expect_equal(nrow(data), 5)
  expect_equal(ncol(data), 500)
  expect_true(all(is.finite(data)))
})

test_that("dataGASCD generates valid data with off-diagonal params", {
  set.seed(42)
  # mu1, mu2, sigma2_1, sigma2_2, x12, x21, A1, A2, B1, B2
  par <- c(-1, 1, 0.5, 0.8, -0.3, 0.2, 0.05, -0.05, 0.5, 0.6)
  par <- set_parameter_attributes(par, K = 2, model_type = "gas",
                                  diag_probs = FALSE, equal_variances = FALSE)

  data <- dataGASCD(5, 500, par)
  expect_true(is.matrix(data))
  expect_equal(nrow(data), 5)
  expect_equal(ncol(data), 500)
  expect_true(all(is.finite(data)))
})

# --- Filtering ---

test_that("Rfiltering_TVP returns finite likelihood with off-diagonal params", {
  set.seed(42)
  y <- c(rnorm(200, -1, sqrt(0.5)), rnorm(200, 1, sqrt(0.5)))

  par <- c(-1, 1, 0.5, -0.3, 0.2, 0.05, -0.05)
  par <- set_parameter_attributes(par, K = 2, model_type = "tvp",
                                  diag_probs = FALSE, equal_variances = TRUE)

  nll <- Rfiltering_TVP(par, y, 100, 50)
  expect_true(is.finite(nll))
  expect_true(nll > 0)
})

test_that("Rfiltering_TVPXExo returns finite likelihood with off-diagonal params", {
  set.seed(42)
  y <- c(rnorm(200, -1, sqrt(0.5)), rnorm(200, 1, sqrt(0.5)))
  X_Exo <- rnorm(400)

  par <- c(-1, 1, 0.5, -0.3, 0.2, 0.1, -0.1)
  par <- set_parameter_attributes(par, K = 2, model_type = "exogenous",
                                  diag_probs = FALSE, equal_variances = TRUE)

  nll <- Rfiltering_TVPXExo(par, X_Exo, y, 100, 50)
  expect_true(is.finite(nll))
  expect_true(nll > 0)
})

test_that("Rfiltering_GAS returns finite likelihood with off-diagonal params", {
  set.seed(42)
  y <- c(rnorm(200, -1, sqrt(0.5)), rnorm(200, 1, sqrt(0.8)))

  # Need separate variances for GAS to set up quadrature properly
  par <- c(-1, 1, 0.5, 0.8, -0.3, 0.2, 0.05, -0.05, 0.5, 0.6)
  par <- set_parameter_attributes(par, K = 2, model_type = "gas",
                                  diag_probs = FALSE, equal_variances = FALSE)

  nll <- Rfiltering_GAS(par, y, 100, 50)
  expect_true(is.finite(nll))
  expect_true(nll > 0)
})

test_that("Rfiltering_Const returns finite likelihood with off-diagonal params", {
  set.seed(42)
  y <- c(rnorm(200, -1, sqrt(0.5)), rnorm(200, 1, sqrt(0.5)))

  par <- c(-1, 1, 0.5, -0.3, 0.2)
  par <- set_parameter_attributes(par, K = 2, model_type = "constant",
                                  diag_probs = FALSE, equal_variances = TRUE)

  nll <- Rfiltering_Const(par, y, 100, 50)
  expect_true(is.finite(nll))
  expect_true(nll > 0)
})

# --- K=3 tests ---

test_that("constant model estimation works with K=3 off-diagonal", {
  set.seed(200)
  y <- c(rnorm(150, -2, 0.5), rnorm(150, 0, 0.5), rnorm(150, 2, 0.5))

  result <- estimate_constant_model(y, K = 3, diag_probs = FALSE,
                                    n_starts = 5, verbose = 0)
  expect_equal(result$diagnostics$convergence_code, 0)
  expect_true(is.finite(result$diagnostics$neg_log_likelihood))
})

test_that("Rfiltering_Const works with K=3 off-diagonal", {
  set.seed(42)
  y <- c(rnorm(100, -2, 0.5), rnorm(100, 0, 0.5), rnorm(100, 2, 0.5))

  # K=3: mu1, mu2, mu3, sigma2, x12, x13, x21, x23, x31, x32
  par <- c(-2, 0, 2, 0.5, -1, -1, -1, -1, -1, -1)
  par <- set_parameter_attributes(par, K = 3, model_type = "constant",
                                  diag_probs = FALSE, equal_variances = TRUE)

  nll <- Rfiltering_Const(par, y, 100, 50)
  expect_true(is.finite(nll))
  expect_true(nll > 0)
})
