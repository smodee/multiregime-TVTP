# Tests for early stopping mechanism in multi-start optimization

# --- Unit tests for create_early_stop_objective ---

test_that("early stopping triggers on stagnation", {
  # Objective that never improves after the first call

  fake_obj <- function(par) 100.0
  config <- create_early_stop_config(enabled = TRUE, patience = 10, max_evals = 10000)
  es <- create_early_stop_objective(fake_obj, config)

  # First call sets best_value, then 10 more with no improvement
  vals <- numeric(12)
  for (i in 1:12) {
    vals[i] <- es$objective(c(0, 0))
  }

  expect_true(es$tracker$early_stopped)
  expect_equal(es$tracker$stop_reason, "stagnation")
  # The call that exceeds patience should return the penalty
  expect_equal(vals[12], 1e20)
  # Earlier calls should return the real value
  expect_equal(vals[1], 100.0)
})

test_that("early stopping triggers on objective explosion", {
  fake_obj <- function(par) 1e11  # Above default max_objective
  config <- create_early_stop_config(enabled = TRUE, max_objective = 1e10)
  es <- create_early_stop_objective(fake_obj, config)

  val <- es$objective(c(0, 0))
  expect_true(es$tracker$early_stopped)
  expect_equal(es$tracker$stop_reason, "objective_explosion")
  expect_equal(val, 1e20)
})

test_that("early stopping triggers on max evaluations exceeded", {
  call_count <- 0
  fake_obj <- function(par) {
    call_count <<- call_count + 1
    100 - call_count * 0.001  # Slowly improving, never stagnates
  }
  config <- create_early_stop_config(enabled = TRUE, patience = 100000,
                                     max_evals = 50)
  es <- create_early_stop_objective(fake_obj, config)

  vals <- numeric(55)
  for (i in 1:55) {
    vals[i] <- es$objective(c(0, 0))
  }
  expect_true(es$tracker$early_stopped)
  expect_equal(es$tracker$stop_reason, "max_evals_exceeded")
  # The 51st call (eval_count > 50) should trigger
  expect_equal(vals[55], 1e20)
})

test_that("early stopping does NOT trigger on healthy convergence", {
  i <- 0
  fake_obj <- function(par) {
    i <<- i + 1
    100 * exp(-i / 10)  # Exponentially decreasing
  }
  config <- create_early_stop_config(enabled = TRUE, patience = 500,
                                     max_evals = 50000)
  es <- create_early_stop_objective(fake_obj, config)

  for (j in 1:200) {
    es$objective(c(0, 0))
  }
  expect_false(es$tracker$early_stopped)
  expect_equal(es$tracker$stop_reason, "")
})

test_that("subsequent calls after early stop return penalty immediately", {
  fake_obj <- function(par) 1e11
  config <- create_early_stop_config(enabled = TRUE, max_objective = 1e10)
  es <- create_early_stop_objective(fake_obj, config)

  es$objective(c(0, 0))  # Triggers early stop
  expect_equal(es$tracker$eval_count, 1L)

  # Subsequent calls should NOT increment eval_count
  val2 <- es$objective(c(0, 0))
  expect_equal(es$tracker$eval_count, 1L)
  expect_equal(val2, 1e20)
})

test_that("non-finite objective values return penalty without triggering early stop", {
  call_count <- 0
  fake_obj <- function(par) {
    call_count <<- call_count + 1
    if (call_count == 3) return(NaN)
    50.0
  }
  config <- create_early_stop_config(enabled = TRUE, patience = 100,
                                     max_evals = 50000)
  es <- create_early_stop_objective(fake_obj, config)

  val1 <- es$objective(c(0, 0))
  val2 <- es$objective(c(0, 0))
  val3 <- es$objective(c(0, 0))  # NaN call

  expect_equal(val1, 50.0)
  expect_equal(val2, 50.0)
  expect_equal(val3, 1e20)  # NaN returns penalty
  expect_false(es$tracker$early_stopped)  # But not early-stopped
  expect_equal(es$tracker$eval_count, 3L)
})

test_that("create_early_stop_config returns correct defaults", {
  config <- create_early_stop_config()
  expect_false(config$enabled)
  expect_equal(config$patience, 3000L)
  expect_equal(config$rel_tol, 1e-8)
  expect_equal(config$max_objective, 1e10)
  expect_equal(config$max_evals, 50000L)
})

test_that("create_early_stop_config respects custom values", {
  config <- create_early_stop_config(enabled = TRUE, patience = 200,
                                     rel_tol = 1e-6, max_objective = 1e8,
                                     max_evals = 10000)
  expect_true(config$enabled)
  expect_equal(config$patience, 200L)
  expect_equal(config$rel_tol, 1e-6)
  expect_equal(config$max_objective, 1e8)
  expect_equal(config$max_evals, 10000L)
})

test_that("stagnation respects rel_tol for tiny improvements", {
  # Improvements just barely above rel_tol should reset the counter
  call_count <- 0
  fake_obj <- function(par) {
    call_count <<- call_count + 1
    # Decrease by a tiny fraction each call (larger than rel_tol=1e-6)
    100 * (1 - call_count * 1e-5)
  }
  config <- create_early_stop_config(enabled = TRUE, patience = 20,
                                     rel_tol = 1e-6, max_evals = 100000)
  es <- create_early_stop_objective(fake_obj, config)

  for (j in 1:50) {
    es$objective(c(0, 0))
  }
  # Should NOT stagnate because improvements exceed rel_tol
  expect_false(es$tracker$early_stopped)
})

# --- Integration tests with model estimation ---

test_that("constant model estimation works with early stopping enabled", {
  skip_on_cran()
  set.seed(100)
  y <- c(rnorm(200, -1, 0.5), rnorm(200, 1, 0.7))

  # Use generous thresholds to avoid prematurely stopping all starts
  result <- estimate_constant_model(y, K = 2, diag_probs = FALSE,
                                    n_starts = 3, n_burnin = 20, n_cutoff = 10, verbose = 0,
                                    parallel = FALSE,
                                    early_stopping = TRUE,
                                    early_stop_patience = 5000,
                                    early_stop_max_evals = 500000)

  # Should produce valid results
  expect_true(is.finite(result$diagnostics$neg_log_likelihood))
  expect_true(result$diagnostics$early_stopping_enabled)
  expect_true("n_early_stopped" %in% names(result$diagnostics))
  expect_true(is.numeric(result$diagnostics$n_early_stopped))
})

test_that("early_stopping=FALSE preserves backward compatibility", {
  skip_on_cran()
  set.seed(100)
  y <- c(rnorm(200, -1, 0.5), rnorm(200, 1, 0.7))

  result <- estimate_constant_model(y, K = 2, diag_probs = FALSE,
                                    n_starts = 2, n_burnin = 20, n_cutoff = 10, verbose = 0,
                                    parallel = FALSE,
                                    early_stopping = FALSE)

  # New fields should exist but with zero/FALSE values
  expect_equal(result$diagnostics$n_early_stopped, 0)
  expect_false(result$diagnostics$early_stopping_enabled)

  # Each individual result should have the new fields
  for (r in result$all_results) {
    expect_false(isTRUE(r$early_stopped))
  }
})

test_that("early-stopped runs never win best-result selection", {
  skip_on_cran()
  set.seed(100)
  y <- c(rnorm(200, -1, 0.5), rnorm(200, 1, 0.7))

  # Use enough starts and generous enough thresholds that at least one converges
  result <- estimate_constant_model(y, K = 2, diag_probs = FALSE,
                                    n_starts = 5, n_burnin = 20, n_cutoff = 10, verbose = 0,
                                    parallel = FALSE,
                                    early_stopping = TRUE,
                                    early_stop_patience = 5000,
                                    early_stop_max_evals = 500000)

  # The best result should have a finite objective (not Inf from early stopping)
  expect_true(is.finite(result$diagnostics$neg_log_likelihood))

  # If any starts were early-stopped, verify the selected best is not one of them
  if (result$diagnostics$n_early_stopped > 0) {
    # Find the selected best start
    best_objectives <- sapply(result$all_results, function(x) x$objective)
    best_idx <- which.min(best_objectives)
    expect_false(isTRUE(result$all_results[[best_idx]]$early_stopped))
  }
})

test_that("TVP model estimation works with early stopping", {
  skip_on_cran()
  set.seed(42)
  y <- c(rnorm(200, -1, 0.5), rnorm(200, 1, 0.7))

  result <- estimate_tvp_model(y, K = 2, diag_probs = FALSE,
                               n_starts = 2, n_burnin = 20, n_cutoff = 10, verbose = 0,
                               parallel = FALSE,
                               early_stopping = TRUE,
                               early_stop_patience = 5000,
                               early_stop_max_evals = 500000)

  expect_true(is.finite(result$diagnostics$neg_log_likelihood))
  expect_true(result$diagnostics$early_stopping_enabled)
  expect_true("n_early_stopped" %in% names(result$diagnostics))
})
