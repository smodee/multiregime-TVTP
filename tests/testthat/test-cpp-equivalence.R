# Test numerical equivalence between R and C filtering implementations
# Covers: {constant, tvp, exo, gas} x {K=2, K=3} x {diag} x {equal_var, separate_var}

test_that("C constant model matches R for K=2 diag equal_var", {
  skip_if_not(cpp_available(), "C backend not available")
  set.seed(42)
  par <- set_parameter_attributes(c(-1, 1, 1, 0.95, 0.90), K = 2, model_type = "constant",
                                   diag_probs = TRUE, equal_variances = TRUE)
  y <- dataConstCD(200, 1, par)[, 1]
  nll_r <- Rfiltering_Const(par, y, 5, 0, use_cpp = FALSE)
  nll_c <- Rfiltering_Const(par, y, 5, 0, use_cpp = TRUE)
  expect_equal(nll_c, nll_r, tolerance = 1e-10)
})

test_that("C constant model matches R for K=2 diag separate_var", {
  skip_if_not(cpp_available(), "C backend not available")
  set.seed(42)
  par <- set_parameter_attributes(c(-1, 1, 0.8, 1.2, 0.95, 0.90), K = 2, model_type = "constant",
                                   diag_probs = TRUE, equal_variances = FALSE)
  y <- dataConstCD(200, 1, par)[, 1]
  nll_r <- Rfiltering_Const(par, y, 5, 0, use_cpp = FALSE)
  nll_c <- Rfiltering_Const(par, y, 5, 0, use_cpp = TRUE)
  expect_equal(nll_c, nll_r, tolerance = 1e-10)
})

test_that("C constant model matches R for K=3 diag", {
  skip_if_not(cpp_available(), "C backend not available")
  set.seed(42)
  par <- set_parameter_attributes(c(-2, 0, 2, 1, 1, 1, 0.90, 0.85, 0.80), K = 3,
                                   model_type = "constant", diag_probs = TRUE, equal_variances = FALSE)
  y <- dataConstCD(200, 1, par)[, 1]
  nll_r <- Rfiltering_Const(par, y, 5, 0, use_cpp = FALSE)
  nll_c <- Rfiltering_Const(par, y, 5, 0, use_cpp = TRUE)
  expect_equal(nll_c, nll_r, tolerance = 1e-10)
})

test_that("C TVP model matches R for K=2 diag", {
  skip_if_not(cpp_available(), "C backend not available")
  set.seed(42)
  par <- set_parameter_attributes(c(-1, 1, 1, 0.95, 0.90, 0.05, 0.05), K = 2,
                                   model_type = "tvp", diag_probs = TRUE, equal_variances = TRUE)
  y <- dataTVPCD(200, 1, par)[, 1]
  nll_r <- Rfiltering_TVP(par, y, 5, 0, use_cpp = FALSE)
  nll_c <- Rfiltering_TVP(par, y, 5, 0, use_cpp = TRUE)
  expect_equal(nll_c, nll_r, tolerance = 1e-10)
})

test_that("C exogenous model matches R for K=2 diag", {
  skip_if_not(cpp_available(), "C backend not available")
  set.seed(42)
  par <- set_parameter_attributes(c(-1, 1, 1, 0.95, 0.90, 0.05, 0.05), K = 2,
                                   model_type = "exogenous", diag_probs = TRUE, equal_variances = TRUE)
  X_Exo <- rnorm(200)
  y <- dataTVPXExoCD(200, 1, par, X_Exo = X_Exo)[, 1]
  nll_r <- Rfiltering_TVPXExo(par, X_Exo, y, 5, 0, use_cpp = FALSE)
  nll_c <- Rfiltering_TVPXExo(par, X_Exo, y, 5, 0, use_cpp = TRUE)
  expect_equal(nll_c, nll_r, tolerance = 1e-10)
})

test_that("C GAS model matches R for K=2 diag (same GH seed)", {
  skip_if_not(cpp_available(), "C backend not available")
  set.seed(42)
  par <- set_parameter_attributes(c(-1, 1, 1, 0.95, 0.90, 0.3, 0.3, 0.5, 0.5), K = 2,
                                   model_type = "gas", diag_probs = TRUE, equal_variances = TRUE)
  y <- dataGASCD(200, 1, par)[, 1]
  set.seed(123)
  nll_r <- Rfiltering_GAS(par, y, 5, 0, use_cpp = FALSE)
  set.seed(123)
  nll_c <- Rfiltering_GAS(par, y, 5, 0, use_cpp = TRUE)
  expect_equal(nll_c, nll_r, tolerance = 1e-8)
})

test_that("C GAS model fallback to constant works", {
  skip_if_not(cpp_available(), "C backend not available")
  set.seed(42)
  # A ~ 0 triggers fallback to constant
  par <- set_parameter_attributes(c(-1, 1, 1, 0.95, 0.90, 1e-6, 1e-6, 0.5, 0.5), K = 2,
                                   model_type = "gas", diag_probs = TRUE, equal_variances = TRUE)
  y <- rnorm(100)
  set.seed(123)
  nll_r <- Rfiltering_GAS(par, y, 5, 0, use_cpp = FALSE)
  set.seed(123)
  nll_c <- Rfiltering_GAS(par, y, 5, 0, use_cpp = TRUE)
  expect_equal(nll_c, nll_r, tolerance = 1e-10)
})

test_that("C constant model works with short series (T=10)", {
  skip_if_not(cpp_available(), "C backend not available")
  set.seed(42)
  par <- set_parameter_attributes(c(-1, 1, 1, 0.95, 0.90), K = 2, model_type = "constant",
                                   diag_probs = TRUE, equal_variances = TRUE)
  y <- rnorm(10)
  nll_r <- Rfiltering_Const(par, y, 0, 0, use_cpp = FALSE)
  nll_c <- Rfiltering_Const(par, y, 0, 0, use_cpp = TRUE)
  expect_equal(nll_c, nll_r, tolerance = 1e-10)
})
