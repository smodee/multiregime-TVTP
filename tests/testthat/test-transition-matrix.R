# Tests for transition matrix construction and related functions

# --- transition_matrix_offdiagonal ---

test_that("off-diagonal K=2 matrix with known probabilities", {
  # p12=0.3, p21=0.2 -> P = [[0.7, 0.3], [0.2, 0.8]]
  P <- transition_matrix(c(0.3, 0.2), diag_probs = FALSE)
  expect_equal(P[1, 1], 0.7, tolerance = 1e-12)
  expect_equal(P[1, 2], 0.3, tolerance = 1e-12)
  expect_equal(P[2, 1], 0.2, tolerance = 1e-12)
  expect_equal(P[2, 2], 0.8, tolerance = 1e-12)
})

test_that("off-diagonal K=3 matrix with equal off-diagonal probs", {
  # Each row has off-diagonal probs = 0.1, so diagonal = 0.8
  P <- transition_matrix(rep(0.1, 6), diag_probs = FALSE)
  for (i in 1:3) {
    expect_equal(P[i, i], 0.8, tolerance = 1e-12)
    for (j in 1:3) {
      if (i != j) expect_equal(P[i, j], 0.1, tolerance = 1e-12)
    }
  }
})

test_that("off-diagonal transition matrix has valid row sums", {
  set.seed(123)
  for (K in 2:4) {
    n_params <- K * (K - 1)
    for (i in 1:10) {
      # Generate valid off-diagonal probabilities that sum to < 1 per row
      max_per_entry <- 0.9 / (K - 1)
      params <- runif(n_params, 0.01, max_per_entry)
      P <- transition_matrix(params, diag_probs = FALSE)
      expect_equal(rowSums(P), rep(1, K), tolerance = 1e-10,
                   label = sprintf("K=%d, trial %d row sums", K, i))
    }
  }
})

test_that("off-diagonal transition matrix entries are all non-negative", {
  set.seed(456)
  for (K in 2:4) {
    n_params <- K * (K - 1)
    for (i in 1:10) {
      max_per_entry <- 0.8 / (K - 1)
      params <- runif(n_params, 0.01, max_per_entry)
      P <- transition_matrix(params, diag_probs = FALSE)
      expect_true(all(P >= 0),
                  label = sprintf("K=%d, trial %d all non-negative", K, i))
    }
  }
})

test_that("off-diagonal small probs produce diagonal-dominant matrix", {
  P <- transition_matrix(rep(0.02, 6), diag_probs = FALSE)
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
  P <- transition_matrix(c(0.2, 0.1), diag_probs = FALSE)
  expect_true(validate_transition_matrix(P))
})

test_that("validation passes for valid diagonal matrices", {
  P <- transition_matrix(c(0.8, 0.9), diag_probs = TRUE)
  expect_true(validate_transition_matrix(P))
})
