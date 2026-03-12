# tests/testthat/test-nll-improvement.R
#
# Issue #28: GAS model correctness test coverage gaps
#
# Gap 2: Verify that nlminb improves the NLL from starting parameters.
# Gap 1: Verify that filtered probabilities identify the correct regime
#         when using true (known) parameters.

skip_on_cran()

# =============================================================================
# SHARED HELPERS
# =============================================================================

n_burnin <- 20L
n_cutoff <- 10L

# Helper: extract starting params for the best optimization run
best_start_params <- function(result) {
  objectives <- sapply(result$all_results, function(x) x$objective)
  result$all_results[[which.min(objectives)]]$start_params
}

# =============================================================================
# GAP 2: NLL IMPROVEMENT FROM STARTING POINT
# =============================================================================

test_that("optimizer improves NLL from starting parameters (constant model)", {
  set.seed(100)
  y <- c(rnorm(150, -1, 0.5), rnorm(150, 1, 0.7))

  result <- estimate_constant_model(
    y, K = 2, diag_probs = TRUE, equal_variances = FALSE,
    n_starts = 4, n_burnin = n_burnin, n_cutoff = n_cutoff, verbose = 0
  )

  start_params <- best_start_params(result)
  initial_nll <- Rfiltering_Const(start_params, y, n_burnin, n_cutoff)
  final_nll <- result$diagnostics$neg_log_likelihood

  expect_true(is.finite(initial_nll))
  expect_true(is.finite(final_nll))
  expect_true(initial_nll > 0)
  expect_true(final_nll > 0)
  expect_lt(final_nll, initial_nll)
})

test_that("optimizer improves NLL from starting parameters (TVP model)", {
  set.seed(100)
  y <- c(rnorm(150, -1, 0.5), rnorm(150, 1, 0.7))

  result <- estimate_tvp_model(
    y, K = 2, diag_probs = FALSE, equal_variances = TRUE,
    n_starts = 4, n_burnin = n_burnin, n_cutoff = n_cutoff, verbose = 0
  )

  start_params <- best_start_params(result)
  initial_nll <- Rfiltering_TVP(start_params, y, n_burnin, n_cutoff)
  final_nll <- result$diagnostics$neg_log_likelihood

  expect_true(is.finite(initial_nll))
  expect_true(is.finite(final_nll))
  expect_true(initial_nll > 0)
  expect_true(final_nll > 0)
  expect_lt(final_nll, initial_nll)
})

test_that("optimizer improves NLL from starting parameters (exogenous model)", {
  set.seed(100)
  y <- c(rnorm(150, -1, 0.5), rnorm(150, 1, 0.7))
  X_Exo <- rnorm(length(y))

  result <- estimate_exo_model(
    y, X_Exo, K = 2, diag_probs = FALSE, equal_variances = TRUE,
    n_starts = 4, n_burnin = n_burnin, n_cutoff = n_cutoff, verbose = 0
  )

  start_params <- best_start_params(result)
  initial_nll <- Rfiltering_TVPXExo(start_params, X_Exo, y, n_burnin, n_cutoff)
  final_nll <- result$diagnostics$neg_log_likelihood

  expect_true(is.finite(initial_nll))
  expect_true(is.finite(final_nll))
  expect_true(initial_nll > 0)
  expect_true(final_nll > 0)
  expect_lt(final_nll, initial_nll)
})

test_that("optimizer improves NLL from starting parameters (GAS model)", {
  set.seed(100)
  y <- c(rnorm(150, -1, 0.5), rnorm(150, 1, 0.7))

  result <- estimate_gas_model(
    y, K = 2, diag_probs = TRUE, equal_variances = FALSE,
    n_starts = 4, n_burnin = n_burnin, n_cutoff = n_cutoff, verbose = 0
  )

  start_params <- best_start_params(result)
  initial_nll <- Rfiltering_GAS(start_params, y, n_burnin, n_cutoff)
  final_nll <- result$diagnostics$neg_log_likelihood

  expect_true(is.finite(initial_nll))
  expect_true(is.finite(final_nll))
  expect_true(initial_nll > 0)
  expect_true(final_nll > 0)
  expect_lt(final_nll, initial_nll)
})

# =============================================================================
# GAP 1: FILTERED PROBABILITY REGIME VERIFICATION (TRUE PARAMS)
# =============================================================================

test_that("GAS filter identifies correct regime when given true parameters", {
  # Well-separated regimes: means at -2 and 2, small variances
  # K=2 diagonal, unequal variances, with moderate GAS dynamics
  par <- c(-2, 2, 0.3, 0.3, 0.9, 0.9, 0.05, -0.05, 0.5, 0.5)
  par <- set_parameter_attributes(par, K = 2, model_type = "gas",
                                  diag_probs = TRUE, equal_variances = FALSE)

  set.seed(123)
  data <- dataGASCD(1, 500, par)
  y <- as.numeric(data[1, ])
  true_states <- attr(data, "simulation_info")$true_states[1, ]

  # Run filter with the true parameters
  filter_burnin <- 50L
  filter_cutoff <- 0L
  nll <- Rfiltering_GAS(par, y, filter_burnin, filter_cutoff, diagnostics = TRUE)

  X_t <- attr(nll, "X.t")
  valid_idx <- attr(nll, "valid_indices")

  # Determine assigned regime at each valid time point
  assigned <- apply(X_t[, valid_idx, drop = FALSE], 2, which.max)
  agreement <- mean(assigned == true_states[valid_idx])

  expect_gt(agreement, 0.80,
            label = sprintf("Regime agreement rate %.1f%% should exceed 80%%",
                            agreement * 100))
})
