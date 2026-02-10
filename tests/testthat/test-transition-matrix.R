# Tests for transition matrix construction and related functions

# --- transition_matrix_offdiagonal ---

test_that("off-diagonal K=2 uniform matrix with zero params", {
  P <- transition_matrix(c(0, 0), diag_probs = FALSE)
  expect_equal(P, matrix(0.5, 2, 2), tolerance = 1e-12)
})

test_that("off-diagonal K=3 uniform matrix with zero params", {
  P <- transition_matrix(rep(0, 6), diag_probs = FALSE)
  expect_equal(P, matrix(1/3, 3, 3), tolerance = 1e-12)
})

test_that("off-diagonal transition matrix has valid row sums", {
  set.seed(123)
  for (K in 2:4) {
    n_params <- K * (K - 1)
    for (i in 1:10) {
      params <- rnorm(n_params, sd = 2)
      P <- transition_matrix(params, diag_probs = FALSE)
      expect_equal(rowSums(P), rep(1, K), tolerance = 1e-10,
                   label = sprintf("K=%d, trial %d row sums", K, i))
    }
  }
})

test_that("off-diagonal transition matrix entries are all positive", {
  set.seed(456)
  for (K in 2:4) {
    n_params <- K * (K - 1)
    for (i in 1:10) {
      params <- rnorm(n_params, sd = 3)
      P <- transition_matrix(params, diag_probs = FALSE)
      expect_true(all(P > 0),
                  label = sprintf("K=%d, trial %d all positive", K, i))
    }
  }
})

test_that("off-diagonal negative params produce diagonal-dominant matrix", {
  P <- transition_matrix(rep(-5, 6), diag_probs = FALSE)
  for (i in 1:3) {
    expect_true(P[i, i] > 0.9,
                label = sprintf("row %d diagonal dominant", i))
  }
})

test_that("off-diagonal transition matrix rejects invalid param length", {
  expect_error(transition_matrix(c(1, 2, 3), diag_probs = FALSE),
               "Invalid length")
})

# --- transition_matrix_diagonal (unchanged) ---

test_that("diagonal parameterization still works correctly for K=2", {
  P <- transition_matrix(c(0.8, 0.9), diag_probs = TRUE)
  expect_equal(P[1, 1], 0.8, tolerance = 1e-12)
  expect_equal(P[2, 2], 0.9, tolerance = 1e-12)
  expect_equal(P[1, 2], 0.2, tolerance = 1e-12)
  expect_equal(P[2, 1], 0.1, tolerance = 1e-12)
})

test_that("diagonal parameterization still works for K=3", {
  P <- transition_matrix(c(0.8, 0.7, 0.9), diag_probs = TRUE)
  expect_equal(diag(P), c(0.8, 0.7, 0.9), tolerance = 1e-12)
  expect_equal(rowSums(P), rep(1, 3), tolerance = 1e-12)
})

# --- validate_transition_matrix ---

test_that("validation passes for valid off-diagonal matrices", {
  P <- transition_matrix(rnorm(6), diag_probs = FALSE)
  expect_true(validate_transition_matrix(P))
})

test_that("validation passes for valid diagonal matrices", {
  P <- transition_matrix(c(0.8, 0.9), diag_probs = TRUE)
  expect_true(validate_transition_matrix(P))
})
