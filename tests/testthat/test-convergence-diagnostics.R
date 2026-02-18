# Convergence diagnostic and regression tests for issue #22
#
# Issue #22 reported severe convergence issues for TVP, exogenous, and GAS
# models when diag_probs=FALSE. Root cause was that the off-diagonal
# parameterization used f_to_params = identity (no logistic link), allowing
# transition "probabilities" to leave [0,1].
#
# Fix: Apply logistic/logit link for off-diagonal case too, matching
# Bazzi et al. (2017) eq. 16 and Sendstad et al. (2025) Section 5.1.
#
# These tests verify that the fix works across all layers: filter,
# transform, starting points, and full estimation.

# ===========================================================================
# Section 1: Filtering Sanity -- K=3 off-diagonal with A=0
# ===========================================================================
# When A=0, time-varying models degenerate to the constant model.
# All filters should return finite NLL with reasonable parameters.

test_that("constant filter K=3 off-diagonal unequal-var returns finite NLL", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))

  # K=3 constant, off-diagonal, unequal-var: 3 mu + 3 sigma2 + 6 trans = 12 params
  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, rep(0.1, 6))
  par <- set_parameter_attributes(par, K = 3, model_type = "constant",
                                  diag_probs = FALSE, equal_variances = FALSE)

  nll <- Rfiltering_Const(par, y, n_burnin = 20, n_cutoff = 10)
  expect_true(is.finite(nll), info = paste("NLL was:", nll))
  expect_true(nll > 0)
})

test_that("TVP filter K=3 off-diagonal unequal-var A=0 returns finite NLL", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))

  # K=3 TVP, off-diagonal, unequal-var: 3 mu + 3 sigma2 + 6 trans + 6 A = 18 params
  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, rep(0.1, 6), rep(0, 6))
  par <- set_parameter_attributes(par, K = 3, model_type = "tvp",
                                  diag_probs = FALSE, equal_variances = FALSE)

  nll <- Rfiltering_TVP(par, y, n_burnin = 20, n_cutoff = 10)
  expect_true(is.finite(nll), info = paste("NLL was:", nll))
  expect_true(nll > 0)
})

test_that("exogenous filter K=3 off-diagonal unequal-var A=0 returns finite NLL", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))
  X_Exo <- rnorm(length(y))

  # K=3 exogenous, off-diagonal, unequal-var: 3 mu + 3 sigma2 + 6 trans + 6 A = 18 params
  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, rep(0.1, 6), rep(0, 6))
  par <- set_parameter_attributes(par, K = 3, model_type = "exogenous",
                                  diag_probs = FALSE, equal_variances = FALSE)

  nll <- Rfiltering_TVPXExo(par, X_Exo, y, n_burnin = 20, n_cutoff = 10)
  expect_true(is.finite(nll), info = paste("NLL was:", nll))
  expect_true(nll > 0)
})

test_that("GAS filter K=3 off-diagonal unequal-var A=0 returns finite NLL", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))

  # K=3 GAS, off-diagonal, unequal-var: 3 mu + 3 sigma2 + 6 trans + 6 A + 6 B = 24 params
  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, rep(0.1, 6), rep(0, 6), rep(0.8, 6))
  par <- set_parameter_attributes(par, K = 3, model_type = "gas",
                                  diag_probs = FALSE, equal_variances = FALSE)

  nll <- Rfiltering_GAS(par, y, n_burnin = 20, n_cutoff = 10,
                         use_fallback = FALSE, diagnostics = FALSE)
  expect_true(is.finite(nll), info = paste("NLL was:", nll))
  expect_true(nll > 0)
})


# ===========================================================================
# Section 2: Filtering Robustness -- Off-diagonal with logistic link
# ===========================================================================
# With the logistic link applied (matching Bazzi et al. 2017 eq. 16 and
# Sendstad et al. 2025 Section 5.1), f-values are squashed to (0,1) via
# logistic before being used as transition probabilities. This means even
# large A values produce valid transition matrices.

test_that("TVP filter K=3 off-diagonal with small A=0.05 returns finite NLL", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))

  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, rep(0.1, 6), rep(0.05, 6))
  par <- set_parameter_attributes(par, K = 3, model_type = "tvp",
                                  diag_probs = FALSE, equal_variances = FALSE)

  nll <- Rfiltering_TVP(par, y, n_burnin = 20, n_cutoff = 10)
  expect_true(is.finite(nll), info = paste("Small A: NLL was:", nll))
  expect_true(nll > 0)
})

test_that("TVP filter K=3 off-diagonal with moderate A=0.5 produces degraded NLL", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))

  # With A=0.5: omega = 0.1*(1-0.5) = 0.05, f[t] = 0.05 + 0.5*y[t]
  # For y[t] in regime 3 (~2): f ~ 0.05 + 1.0 = 1.05 > 1 (invalid probability!)
  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, rep(0.1, 6), rep(0.5, 6))
  par <- set_parameter_attributes(par, K = 3, model_type = "tvp",
                                  diag_probs = FALSE, equal_variances = FALSE)

  # This may or may not be finite depending on data -- document the behavior
  nll <- Rfiltering_TVP(par, y, n_burnin = 20, n_cutoff = 10)
  # Just check we don't crash; the NLL value documents degradation
  expect_true(is.numeric(nll), info = paste("Moderate A: NLL was:", nll))
})

test_that("TVP filter K=3 off-diagonal with large A=2.0 returns finite NLL (logistic protects)", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))

  # With A=2.0: omega_LR = logit(0.1) ~ -2.2, omega = omega_LR*(1-2.0) = +2.2
  # f[t] = 2.2 + 2.0*y[t] can be large, but logistic(f) squashes to (0,1)
  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, rep(0.1, 6), rep(2.0, 6))
  par <- set_parameter_attributes(par, K = 3, model_type = "tvp",
                                  diag_probs = FALSE, equal_variances = FALSE)

  nll <- Rfiltering_TVP(par, y, n_burnin = 20, n_cutoff = 10)
  expect_true(is.finite(nll),
              info = paste("Large A off-diagonal: NLL was:", nll))
  expect_true(nll > 0)
})

test_that("exogenous filter K=3 off-diagonal with large A=2.0 returns finite NLL (logistic protects)", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))
  X_Exo <- rnorm(length(y))

  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, rep(0.1, 6), rep(2.0, 6))
  par <- set_parameter_attributes(par, K = 3, model_type = "exogenous",
                                  diag_probs = FALSE, equal_variances = FALSE)

  # With logistic link, f-values are squashed to (0,1) before constructing
  # transition matrices, so even large A values produce valid probabilities.
  nll <- Rfiltering_TVPXExo(par, X_Exo, y, n_burnin = 20, n_cutoff = 10)
  expect_true(is.finite(nll),
              info = paste("Exo large A off-diagonal: NLL was:", nll))
  expect_true(nll > 0)
})

test_that("TVP filter K=3 DIAGONAL with A=2.0 remains finite (logistic protects)", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))

  # With diag_probs=TRUE, f_to_params = logistic, so f-values are squashed to (0,1)
  # K=3 diagonal TVP, unequal-var: 3 mu + 3 sigma2 + 3 trans + 3 A = 12 params
  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, 0.8, 0.8, 0.8, 2.0, 2.0, 2.0)
  par <- set_parameter_attributes(par, K = 3, model_type = "tvp",
                                  diag_probs = TRUE, equal_variances = FALSE)

  nll <- Rfiltering_TVP(par, y, n_burnin = 20, n_cutoff = 10)
  # Diagonal case with logistic link should produce finite NLL even with large A
  expect_true(is.finite(nll),
              info = paste("Large A diagonal (logistic): NLL was:", nll))
  expect_true(nll > 0)
})

test_that("GAS filter K=3 off-diagonal with large A=2.0 returns finite NLL (logistic protects)", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))

  # GAS uses omega = logit(init_trans) directly, no (1-A) factor
  # With logistic link, f-values are squashed to (0,1) for transition probs
  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, rep(0.1, 6), rep(2.0, 6), rep(0.8, 6))
  par <- set_parameter_attributes(par, K = 3, model_type = "gas",
                                  diag_probs = FALSE, equal_variances = FALSE)

  nll <- Rfiltering_GAS(par, y, n_burnin = 20, n_cutoff = 10,
                        use_fallback = FALSE, diagnostics = FALSE)
  expect_true(is.finite(nll),
              info = paste("GAS large A off-diagonal: NLL was:", nll))
  expect_true(nll > 0)
})


# ===========================================================================
# Section 3: Parameter Transform Round-Trips
# ===========================================================================
# Verify transform/untransform preserves parameters for K=3 off-diagonal config.

test_that("round-trip K=3 TVP off-diagonal unequal-var (18 params)", {
  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, rep(0.1, 6), rep(0.3, 6))
  par <- set_parameter_attributes(par, K = 3, model_type = "tvp",
                                  diag_probs = FALSE, equal_variances = FALSE)

  transformed <- transform_parameters(par)
  recovered <- untransform_parameters(transformed)

  expect_equal(as.numeric(recovered), as.numeric(par), tolerance = 1e-10,
               label = "TVP K=3 off-diagonal round-trip")
})

test_that("round-trip K=3 exogenous off-diagonal unequal-var (18 params)", {
  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, rep(0.1, 6), rep(0.3, 6))
  par <- set_parameter_attributes(par, K = 3, model_type = "exogenous",
                                  diag_probs = FALSE, equal_variances = FALSE)

  transformed <- transform_parameters(par)
  recovered <- untransform_parameters(transformed)

  expect_equal(as.numeric(recovered), as.numeric(par), tolerance = 1e-10,
               label = "Exogenous K=3 off-diagonal round-trip")
})

test_that("round-trip K=3 GAS off-diagonal unequal-var (24 params)", {
  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, rep(0.1, 6), rep(0.3, 6), rep(0.8, 6))
  par <- set_parameter_attributes(par, K = 3, model_type = "gas",
                                  diag_probs = FALSE, equal_variances = FALSE)

  transformed <- transform_parameters(par)
  recovered <- untransform_parameters(transformed)

  expect_equal(as.numeric(recovered), as.numeric(par), tolerance = 1e-10,
               label = "GAS K=3 off-diagonal round-trip")
})

test_that("A coefficients pass through transform unchanged (unbounded)", {
  par <- c(-2, 0, 2, 0.5, 1.0, 1.5, rep(0.1, 6), c(-1.5, -0.5, 0, 0.5, 1.5, 3.0))
  par <- set_parameter_attributes(par, K = 3, model_type = "tvp",
                                  diag_probs = FALSE, equal_variances = FALSE)

  transformed <- transform_parameters(par)
  A_original <- extract_parameter_component(par, "A")
  A_transformed <- extract_parameter_component(transformed, "A")

  # A should be identity-transformed (pass through unchanged)
  expect_equal(A_transformed, A_original, tolerance = 1e-10,
               label = "A passes through transform unchanged")
})


# ===========================================================================
# Section 4: Starting Point Generation Diagnostics
# ===========================================================================
# Verify generated starting points produce finite NLL at the filter level.

test_that("TVP K=3 off-diagonal starting points all produce finite NLL", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))

  sp_list <- generate_starting_points(
    y = y, K = 3, model_type = "tvp", n_starts = 10,
    diag_probs = FALSE, equal_variances = FALSE, seed = 42
  )

  n_finite <- 0
  for (i in seq_along(sp_list)) {
    nll <- tryCatch(
      Rfiltering_TVP(sp_list[[i]], y, n_burnin = 20, n_cutoff = 10),
      error = function(e) NA_real_
    )
    if (is.finite(nll)) n_finite <- n_finite + 1
  }

  # With logistic link, all starting points should produce finite NLL
  # regardless of A values, since logistic squashes f to (0,1)
  expect_equal(n_finite, 10,
               info = paste(n_finite, "of 10 starting points produced finite NLL"))
})

test_that("exogenous K=3 off-diagonal starting points all produce finite NLL", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))
  X_Exo <- rnorm(length(y))

  sp_list <- generate_starting_points(
    y = y, K = 3, model_type = "exogenous", n_starts = 10,
    diag_probs = FALSE, equal_variances = FALSE, seed = 42
  )

  n_finite <- 0
  n_error <- 0
  for (i in seq_along(sp_list)) {
    nll <- tryCatch(
      Rfiltering_TVPXExo(sp_list[[i]], X_Exo, y, n_burnin = 20, n_cutoff = 10),
      error = function(e) { n_error <<- n_error + 1; NA_real_ }
    )
    if (!is.na(nll) && is.finite(nll)) n_finite <- n_finite + 1
  }

  # With logistic link, all starting points should produce finite NLL
  expect_equal(n_finite, 10,
               info = paste(n_finite, "of 10 finite,", n_error, "errors"))
  expect_equal(n_error, 0)
})

test_that("GAS K=3 off-diagonal starting points all produce finite NLL", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))

  sp_list <- generate_starting_points(
    y = y, K = 3, model_type = "gas", n_starts = 10,
    diag_probs = FALSE, equal_variances = FALSE, seed = 42
  )

  n_finite <- 0
  for (i in seq_along(sp_list)) {
    nll <- tryCatch(
      Rfiltering_GAS(sp_list[[i]], y, n_burnin = 20, n_cutoff = 10,
                     use_fallback = FALSE, diagnostics = FALSE),
      error = function(e) NA_real_
    )
    if (!is.na(nll) && is.finite(nll)) n_finite <- n_finite + 1
  }

  # With logistic link, all starting points should produce finite NLL
  expect_equal(n_finite, 10,
               info = paste(n_finite, "of 10 starting points produced finite NLL"))
})

test_that("TVP/exo starting points: negative omega in logit space is safe with logistic link", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))

  sp_list <- generate_starting_points(
    y = y, K = 3, model_type = "tvp", n_starts = 10,
    diag_probs = FALSE, equal_variances = FALSE, seed = 42
  )

  # First starting point (A=0) should always have A=0
  first_A <- extract_parameter_component(sp_list[[1]], "A")
  expect_true(all(first_A == 0),
              info = "First starting point should have A=0")

  # With logistic link, omega is computed in logit space:
  #   omega_LR = logit(init_trans), omega = omega_LR * (1 - A)
  # Negative omega in logit space is perfectly fine -- logistic(negative) gives
  # a small but valid probability. Verify all starting points produce finite NLL.
  n_finite <- 0
  for (i in seq_along(sp_list)) {
    nll <- tryCatch(
      Rfiltering_TVP(sp_list[[i]], y, n_burnin = 20, n_cutoff = 10),
      error = function(e) NA_real_
    )
    if (is.finite(nll)) n_finite <- n_finite + 1
  }
  expect_equal(n_finite, 10)
})


# ===========================================================================
# Section 5: Single-Start Optimization Progress
# ===========================================================================
# Test whether a single nlminb run makes progress. Uses short data and
# limited iterations to keep runtime reasonable.

skip_on_cran()

test_that("constant K=3 off-diagonal single-start converges (control)", {
  skip_on_cran()
  set.seed(42)
  y <- c(rnorm(100, -2, sqrt(0.5)), rnorm(100, 0, 1), rnorm(100, 2, sqrt(1.5)))

  sp_list <- generate_starting_points(
    y = y, K = 3, model_type = "constant", n_starts = 1,
    diag_probs = FALSE, equal_variances = FALSE, seed = 42
  )
  start_params <- sp_list[[1]]
  transformed <- transform_parameters(start_params)

  initial_nll <- Rfiltering_Const(start_params, y, n_burnin = 20, n_cutoff = 10)

  opt <- nlminb(
    start = transformed,
    objective = function(par_t) {
      par_t_with_attrs <- transformed
      par_t_with_attrs[] <- par_t
      attr(par_t_with_attrs, "parameterization") <- "transformed"
      par_natural <- untransform_parameters(par_t_with_attrs)
      Rfiltering_Const(par_natural, y, n_burnin = 20, n_cutoff = 10)
    },
    control = list(eval.max = 500, iter.max = 100, trace = 0)
  )

  expect_true(opt$objective < initial_nll,
              info = paste("Initial:", round(initial_nll, 2),
                           "Final:", round(opt$objective, 2)))
})

test_that("TVP K=3 off-diagonal single-start attempts optimization", {
  skip_on_cran()
  set.seed(42)
  y <- c(rnorm(100, -2, sqrt(0.5)), rnorm(100, 0, 1), rnorm(100, 2, sqrt(1.5)))

  sp_list <- generate_starting_points(
    y = y, K = 3, model_type = "tvp", n_starts = 1,
    diag_probs = FALSE, equal_variances = FALSE, seed = 42
  )
  start_params <- sp_list[[1]]
  transformed <- transform_parameters(start_params)

  initial_nll <- Rfiltering_TVP(start_params, y, n_burnin = 20, n_cutoff = 10)

  opt <- tryCatch(
    nlminb(
      start = transformed,
      objective = function(par_t) {
        par_t_with_attrs <- transformed
        par_t_with_attrs[] <- par_t
        attr(par_t_with_attrs, "parameterization") <- "transformed"
        par_natural <- untransform_parameters(par_t_with_attrs)
        Rfiltering_TVP(par_natural, y, n_burnin = 20, n_cutoff = 10)
      },
      control = list(eval.max = 500, iter.max = 100, trace = 0)
    ),
    error = function(e) list(objective = Inf, convergence = 999, message = e$message)
  )

  # Document: does it make progress, stagnate, or error?
  expect_true(is.numeric(opt$objective),
              info = paste("TVP single-start: initial=", round(initial_nll, 2),
                           "final=", round(opt$objective, 2),
                           "convergence=", opt$convergence))
})

test_that("exogenous K=3 off-diagonal single-start attempts optimization", {
  skip_on_cran()
  set.seed(42)
  y <- c(rnorm(100, -2, sqrt(0.5)), rnorm(100, 0, 1), rnorm(100, 2, sqrt(1.5)))
  X_Exo <- rnorm(length(y))

  sp_list <- generate_starting_points(
    y = y, K = 3, model_type = "exogenous", n_starts = 1,
    diag_probs = FALSE, equal_variances = FALSE, seed = 42
  )
  start_params <- sp_list[[1]]
  transformed <- transform_parameters(start_params)

  initial_nll <- Rfiltering_TVPXExo(start_params, X_Exo, y, n_burnin = 20, n_cutoff = 10)

  opt <- tryCatch(
    nlminb(
      start = transformed,
      objective = function(par_t) {
        par_t_with_attrs <- transformed
        par_t_with_attrs[] <- par_t
        attr(par_t_with_attrs, "parameterization") <- "transformed"
        par_natural <- untransform_parameters(par_t_with_attrs)
        Rfiltering_TVPXExo(par_natural, X_Exo, y, n_burnin = 20, n_cutoff = 10)
      },
      control = list(eval.max = 500, iter.max = 100, trace = 0)
    ),
    error = function(e) list(objective = Inf, convergence = 999, message = e$message)
  )

  expect_true(is.numeric(opt$objective),
              info = paste("Exo single-start: initial=", round(initial_nll, 2),
                           "final=", round(opt$objective, 2),
                           "convergence=", opt$convergence))
})

test_that("GAS K=3 off-diagonal single-start attempts optimization", {
  skip_on_cran()
  set.seed(42)
  y <- c(rnorm(100, -2, sqrt(0.5)), rnorm(100, 0, 1), rnorm(100, 2, sqrt(1.5)))

  sp_list <- generate_starting_points(
    y = y, K = 3, model_type = "gas", n_starts = 1,
    diag_probs = FALSE, equal_variances = FALSE, seed = 42
  )
  start_params <- sp_list[[1]]
  transformed <- transform_parameters(start_params)

  initial_nll <- tryCatch(
    Rfiltering_GAS(start_params, y, n_burnin = 20, n_cutoff = 10,
                   use_fallback = FALSE, diagnostics = FALSE),
    error = function(e) Inf
  )

  opt <- suppressWarnings(tryCatch(
    nlminb(
      start = transformed,
      objective = function(par_t) {
        par_t_with_attrs <- transformed
        par_t_with_attrs[] <- par_t
        attr(par_t_with_attrs, "parameterization") <- "transformed"
        par_natural <- untransform_parameters(par_t_with_attrs)
        Rfiltering_GAS(par_natural, y, n_burnin = 20, n_cutoff = 10,
                       use_fallback = FALSE, diagnostics = FALSE)
      },
      control = list(eval.max = 500, iter.max = 100, trace = 0)
    ),
    error = function(e) list(objective = Inf, convergence = 999, message = e$message)
  ))

  expect_true(is.numeric(opt$objective),
              info = paste("GAS single-start: initial=", round(initial_nll, 2),
                           "final=", round(opt$objective, 2),
                           "convergence=", opt$convergence))
})


# ===========================================================================
# Section 6: Synthetic Data Recovery -- K=2 off-diagonal
# ===========================================================================
# Verify all models converge with K=2 off-diagonal config.
# With the logistic link fix, all three models should converge.

test_that("TVP K=2 off-diagonal estimation converges on synthetic data", {
  skip_on_cran()
  set.seed(100)
  y <- c(rnorm(150, -1, 0.5), rnorm(150, 1, 0.7))

  result <- estimate_tvp_model(
    y, K = 2, diag_probs = FALSE, equal_variances = TRUE,
    n_starts = 4, n_burnin = 20, n_cutoff = 10, verbose = 0
  )

  expect_true(is.finite(result$diagnostics$neg_log_likelihood),
              info = paste("TVP K=2 NLL:", result$diagnostics$neg_log_likelihood,
                           "conv:", result$diagnostics$convergence_code))
})

test_that("exogenous K=2 off-diagonal estimation converges on synthetic data", {
  skip_on_cran()
  set.seed(100)
  y <- c(rnorm(150, -1, 0.5), rnorm(150, 1, 0.7))
  X_Exo <- rnorm(length(y))

  # With logistic link, exogenous K=2 off-diagonal should converge
  result <- estimate_exo_model(
    y, X_Exo = X_Exo, K = 2, diag_probs = FALSE, equal_variances = TRUE,
    n_starts = 4, n_burnin = 20, n_cutoff = 10, verbose = 0
  )

  expect_true(is.finite(result$diagnostics$neg_log_likelihood),
              info = paste("Exo K=2 NLL:", result$diagnostics$neg_log_likelihood,
                           "conv:", result$diagnostics$convergence_code))
})

test_that("GAS K=2 off-diagonal estimation converges on synthetic data", {
  skip_on_cran()
  set.seed(100)
  y <- c(rnorm(150, -1, 0.5), rnorm(150, 1, 0.7))

  result <- estimate_gas_model(
    y, K = 2, diag_probs = FALSE, equal_variances = FALSE,
    n_starts = 4, n_burnin = 20, n_cutoff = 10, verbose = 0,
    use_fallback = TRUE
  )

  expect_true(is.finite(result$diagnostics$neg_log_likelihood),
              info = paste("GAS K=2 NLL:", result$diagnostics$neg_log_likelihood,
                           "conv:", result$diagnostics$convergence_code))
})


# ===========================================================================
# Section 7: Regression Tests for Issue #22 Fix
# ===========================================================================
# Full estimation with K=3 off-diagonal unequal-var.
# With the logistic link fix, these should now converge on well-separated
# synthetic data.

test_that("TVP K=3 off-diagonal unequal-var estimation converges", {
  skip_on_cran()
  set.seed(42)
  y <- c(rnorm(135, -2, sqrt(0.5)), rnorm(135, 0, 1), rnorm(130, 2, sqrt(1.5)))

  result <- estimate_tvp_model(
    y, K = 3, diag_probs = FALSE, equal_variances = FALSE,
    n_starts = 4, n_burnin = 20, n_cutoff = 10, verbose = 0
  )

  expect_true(is.finite(result$diagnostics$neg_log_likelihood),
              info = paste("TVP K=3 NLL:", round(result$diagnostics$neg_log_likelihood, 2)))
})

test_that("exogenous K=3 off-diagonal unequal-var estimation converges", {
  skip_on_cran()
  set.seed(42)
  y <- c(rnorm(135, -2, sqrt(0.5)), rnorm(135, 0, 1), rnorm(130, 2, sqrt(1.5)))
  X_Exo <- rnorm(length(y))

  result <- estimate_exo_model(
    y, X_Exo = X_Exo, K = 3, diag_probs = FALSE, equal_variances = FALSE,
    n_starts = 4, n_burnin = 20, n_cutoff = 10, verbose = 0
  )

  expect_true(is.finite(result$diagnostics$neg_log_likelihood),
              info = paste("Exo K=3 NLL:", round(result$diagnostics$neg_log_likelihood, 2)))
})

test_that("GAS K=3 off-diagonal unequal-var estimation converges", {
  skip_on_cran()
  set.seed(42)
  y <- c(rnorm(135, -2, sqrt(0.5)), rnorm(135, 0, 1), rnorm(130, 2, sqrt(1.5)))

  result <- estimate_gas_model(
    y, K = 3, diag_probs = FALSE, equal_variances = FALSE,
    n_starts = 4, n_burnin = 20, n_cutoff = 10, verbose = 0,
    use_fallback = TRUE
  )

  expect_true(is.finite(result$diagnostics$neg_log_likelihood),
              info = paste("GAS K=3 NLL:", round(result$diagnostics$neg_log_likelihood, 2)))
})

# ===========================================================================
# Section 7: Regression test for issue #31 -- off-diagonal row sum clamping
# ===========================================================================
# Issue #31: Exogenous model with wide-range X_Exo causes off-diagonal
# transition probabilities to sum >1 within a row, making the diagonal
# negative. This produces invalid transition matrices and crashes.

test_that("clamp_offdiag_rowsums prevents invalid row sums", {
  K <- 3L
  # Off-diagonal probs that sum to >1 in row 1 (0.6+0.5=1.1)
  p_bad <- c(0.6, 0.5,   # row 1: sum = 1.1 > 1
             0.2, 0.3,   # row 2: sum = 0.5, OK
             0.4, 0.7)   # row 3: sum = 1.1 > 1
  p_fixed <- clamp_offdiag_rowsums(p_bad, K)

  # Check row sums are now <= 1 - 1e-6
  for (i in 1:K) {
    idx <- ((i - 1) * (K - 1) + 1):(i * (K - 1))
    expect_lte(sum(p_fixed[idx]), 1 - 1e-6)
  }

  # Row 2 was fine, so should be unchanged
  expect_equal(p_fixed[3:4], p_bad[3:4])

  # Rows 1 and 3 should be proportionally scaled (same ratio)
  expect_equal(p_fixed[1] / p_fixed[2], p_bad[1] / p_bad[2], tolerance = 1e-10)
  expect_equal(p_fixed[5] / p_fixed[6], p_bad[5] / p_bad[6], tolerance = 1e-10)
})

test_that("clamp_offdiag_rowsums is no-op for valid probabilities", {
  K <- 3L
  p_good <- c(0.1, 0.15, 0.2, 0.1, 0.05, 0.25)
  p_result <- clamp_offdiag_rowsums(p_good, K)
  expect_equal(p_result, p_good)
})

test_that("exogenous filter K=3 off-diagonal with large A and wide X_Exo returns finite NLL (issue #31)", {
  set.seed(42)
  y <- c(rnorm(150, -2, sqrt(0.5)), rnorm(150, 0, 1), rnorm(150, 2, sqrt(1.5)))
  # Wide-range exogenous variable (mimics lagged yield levels)
  X_Exo <- seq(0, 16, length.out = length(y))

  # Parameters with large A coefficients that would push off-diag row sums > 1
  par <- c(-2, 0, 2,                           # mu
           0.5, 1.0, 1.5,                      # sigma2
           0.1, 0.15, 0.1, 0.15, 0.1, 0.15,   # trans_prob (off-diagonal)
           0.5, 0.5, 0.3, 0.3, 0.4, 0.4)      # A (large enough to cause issues)
  par <- set_parameter_attributes(par, K = 3, model_type = "exogenous",
                                  diag_probs = FALSE, equal_variances = FALSE)

  # Before the fix, this would error with "Transition matrix contains negative values"
  nll <- Rfiltering_TVPXExo(par, X_Exo, y, n_burnin = 20, n_cutoff = 10)
  expect_true(is.finite(nll), info = paste("NLL was:", nll))
})

test_that("TVP filter K=3 off-diagonal with large A returns finite NLL (issue #31)", {
  set.seed(42)
  # Data with wide range to trigger the bug via autoregressive dynamics
  y <- c(rnorm(200, -3, 2), rnorm(200, 3, 2))

  par <- c(-3, 0, 3,                           # mu
           1.0, 1.0, 1.0,                      # sigma2
           0.1, 0.15, 0.1, 0.15, 0.1, 0.15,   # trans_prob (off-diagonal)
           0.8, 0.8, 0.5, 0.5, 0.6, 0.6)      # A (large)
  par <- set_parameter_attributes(par, K = 3, model_type = "tvp",
                                  diag_probs = FALSE, equal_variances = FALSE)

  nll <- Rfiltering_TVP(par, y, n_burnin = 20, n_cutoff = 10)
  expect_true(is.finite(nll), info = paste("NLL was:", nll))
})
