# Tests for convert_parameterization() and inverse softmax

test_that("convert diag to off-diag round-trip for K=2", {
  diag_params <- c(0.8, 0.9)
  off_diag <- multiregimeTVTP:::convert_parameterization(diag_params, from_diag = TRUE)

  # Reconstruct matrix from softmax params
  P <- transition_matrix(off_diag, diag_probs = FALSE)
  expect_equal(P[1, 1], 0.8, tolerance = 1e-10)
  expect_equal(P[2, 2], 0.9, tolerance = 1e-10)
  expect_equal(P[1, 2], 0.2, tolerance = 1e-10)
  expect_equal(P[2, 1], 0.1, tolerance = 1e-10)
})

test_that("convert diag to off-diag round-trip for K=3", {
  diag_params <- c(0.7, 0.8, 0.9)
  off_diag <- multiregimeTVTP:::convert_parameterization(diag_params, from_diag = TRUE)
  P <- transition_matrix(off_diag, diag_probs = FALSE)

  # Diagonal entries match
  expect_equal(diag(P), diag_params, tolerance = 1e-10)
  # Row sums = 1
  expect_equal(rowSums(P), rep(1, 3), tolerance = 1e-10)
})

test_that("convert off-diag to diag extracts correct diagonal", {
  params <- c(-1, -2)  # K=2: x_12, x_21
  P <- transition_matrix(params, diag_probs = FALSE)
  diag_params <- multiregimeTVTP:::convert_parameterization(params, from_diag = FALSE)

  expect_equal(diag_params, diag(P), tolerance = 1e-10)
})

test_that("inverse softmax is correct: x_ij = log(p_ij / p_ii)", {
  # Start from known probabilities
  diag_params <- c(0.6, 0.85)
  off_diag <- multiregimeTVTP:::convert_parameterization(diag_params, from_diag = TRUE)

  # For K=2: off_diag[1] = log(p12/p11), off_diag[2] = log(p21/p22)
  expect_equal(off_diag[1], log(0.4 / 0.6), tolerance = 1e-10)
  expect_equal(off_diag[2], log(0.15 / 0.85), tolerance = 1e-10)
})

test_that("full round-trip: diag -> off-diag -> matrix -> diag", {
  original <- c(0.75, 0.85, 0.90)
  off_diag <- multiregimeTVTP:::convert_parameterization(original, from_diag = TRUE)
  recovered <- multiregimeTVTP:::convert_parameterization(off_diag, from_diag = FALSE)
  expect_equal(recovered, original, tolerance = 1e-10)
})

test_that("convert_parameterization rejects invalid input", {
  expect_error(multiregimeTVTP:::convert_parameterization(numeric(0)),
               "non-empty numeric")
  expect_error(multiregimeTVTP:::convert_parameterization(0.5, from_diag = TRUE),
               "at least 2")
})
