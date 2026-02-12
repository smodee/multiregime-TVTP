# Tests for GAS score functions (logit-based derivatives)

test_that("raw score vector has correct length", {
  K <- 3
  eta <- c(0.1, 0.8, 0.1)
  X_t_prev <- c(0.3, 0.4, 0.3)
  p_trans <- rep(0.1, 6)  # Valid off-diagonal probabilities
  tot_lik <- sum(eta * X_t_prev)

  score <- multiregimeTVTP:::calculate_raw_score_vector(eta, tot_lik, X_t_prev, p_trans, K)
  expect_length(score, K * (K - 1))
})

test_that("raw score is zero when all regimes have equal likelihoods", {
  K <- 3
  eta <- c(1, 1, 1)  # Equal likelihoods
  X_t_prev <- c(1/3, 1/3, 1/3)
  p_trans <- rep(0.1, 6)  # Valid probabilities
  tot_lik <- sum(eta * X_t_prev)

  score <- multiregimeTVTP:::calculate_raw_score_vector(eta, tot_lik, X_t_prev, p_trans, K)

  # When all eta are equal, (eta[i] - eta[j]) = 0
  expect_equal(as.numeric(score), rep(0, 6), tolerance = 1e-12)
})

test_that("raw score has correct sign structure", {
  K <- 2
  # Regime 1 has much higher likelihood
  eta <- c(0.99, 0.01)
  X_t_prev <- c(0.5, 0.5)
  p_trans <- c(0.2, 0.2)  # Valid probabilities
  tot_lik <- sum(eta * X_t_prev)

  score <- multiregimeTVTP:::calculate_raw_score_vector(eta, tot_lik, X_t_prev, p_trans, K)
  # p_trans[1] is p_{1,2}: transition from regime 1 to regime 2
  # (eta[1] - eta[2]) > 0, so score[1] should be positive
  expect_true(score[1] > 0, label = "Score for p_12 positive when regime 1 is favored")
})

test_that("raw score rejects invalid inputs", {
  K <- 2
  eta <- c(0.5, 0.5)
  X_t_prev <- c(0.5, 0.5)
  p_trans <- c(0.2, 0.2)
  tot_lik <- 0.5

  expect_error(
    multiregimeTVTP:::calculate_raw_score_vector(eta, tot_lik, c(0.3, 0.3, 0.4), p_trans, K),
    "X_t_prev must be a numeric vector of length K"
  )
  expect_error(
    multiregimeTVTP:::calculate_raw_score_vector(eta, -1, X_t_prev, p_trans, K),
    "tot_lik must be a positive scalar"
  )
})

test_that("raw score rejects out-of-range probabilities", {
  K <- 2
  eta <- c(0.3, 0.7)
  X_t_prev <- c(0.5, 0.5)
  tot_lik <- 0.5

  # Negative probability
  expect_error(
    multiregimeTVTP:::calculate_raw_score_vector(eta, tot_lik, X_t_prev, c(-0.1, 0.2), K),
    "p_trans must contain valid probabilities"
  )
  # Probability > 1
  expect_error(
    multiregimeTVTP:::calculate_raw_score_vector(eta, tot_lik, X_t_prev, c(1.1, 0.2), K),
    "p_trans must contain valid probabilities"
  )
})

# --- Fisher Information ---

test_that("Fisher Information is positive and finite", {
  K <- 2
  mu <- c(-1, 1)
  sigma2 <- c(0.5, 1)
  X_t_prev <- c(0.5, 0.5)
  p_trans <- c(0.2, 0.1)  # Valid off-diagonal probabilities

  y <- rnorm(500)
  gh_setup <- setup_gauss_hermite_quadrature(y)

  fi <- multiregimeTVTP:::calculate_fisher_information(mu, sigma2, X_t_prev, p_trans, gh_setup, K)
  expect_true(is.finite(fi))
  expect_true(fi > 0)
})

test_that("Fisher Information depends on off-diagonal transition probs", {
  # With non-uniform X_t_prev, the transition matrix affects X_pred,
  # which in turn affects Fisher info. Before the fix, convert_to_valid_probs
  # misinterpreted off-diagonal probs as diagonal probs (wrong matrix).
  K <- 2
  mu <- c(-2, 2)
  sigma2 <- c(0.5, 0.5)
  X_t_prev <- c(0.8, 0.2)  # Non-uniform so P actually affects X_pred

  y <- rnorm(500)
  gh_setup <- setup_gauss_hermite_quadrature(y)

  fi_small <- multiregimeTVTP:::calculate_fisher_information(
    mu, sigma2, X_t_prev, c(0.05, 0.05), gh_setup, K
  )
  fi_large <- multiregimeTVTP:::calculate_fisher_information(
    mu, sigma2, X_t_prev, c(0.45, 0.45), gh_setup, K
  )

  # The key assertion: different transition probs produce different Fisher info
  # (with uniform X_t_prev they'd be equal since P drops out)
  expect_false(isTRUE(all.equal(fi_small, fi_large)),
               label = "Fisher Info differs for different transition probs")
})

test_that("Fisher Information works for K=3 off-diagonal", {
  K <- 3
  mu <- c(-2, 0, 2)
  sigma2 <- c(0.5, 0.5, 0.5)
  X_t_prev <- c(0.5, 0.3, 0.2)
  p_trans <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)  # 6 off-diagonal probs

  y <- rnorm(500)
  gh_setup <- setup_gauss_hermite_quadrature(y)

  fi <- multiregimeTVTP:::calculate_fisher_information(
    mu, sigma2, X_t_prev, p_trans, gh_setup, K
  )
  expect_true(is.finite(fi))
  expect_true(fi > 0)
})

test_that("Fisher Information is larger when regimes are well-separated", {
  K <- 2
  sigma2 <- c(0.5, 0.5)
  X_t_prev <- c(0.5, 0.5)
  p_trans <- c(0.2, 0.1)

  y <- rnorm(500)
  gh_setup <- setup_gauss_hermite_quadrature(y)

  fi_close <- multiregimeTVTP:::calculate_fisher_information(
    c(-0.1, 0.1), sigma2, X_t_prev, p_trans, gh_setup, K
  )
  fi_far <- multiregimeTVTP:::calculate_fisher_information(
    c(-3, 3), sigma2, X_t_prev, p_trans, gh_setup, K
  )

  expect_true(fi_far > fi_close,
              label = "Fisher Info larger for well-separated regimes")
})

# --- calculate_gas_score (integration) ---

test_that("calculate_gas_score returns complete result for off-diagonal", {
  K <- 2
  mu <- c(-1, 1)
  sigma2 <- c(0.5, 1)
  X_pred <- c(0.6, 0.4)
  trans_prob <- c(0.2, 0.1)  # Valid off-diagonal probabilities

  y <- rnorm(500)
  gh_setup <- setup_gauss_hermite_quadrature(y)

  result <- calculate_gas_score(
    y_obs = 0.5, mu = mu, sigma2 = sigma2,
    trans_prob = trans_prob, diag_probs = FALSE,
    X_pred = X_pred, gh_setup = gh_setup
  )

  expect_true(is.list(result))
  expect_true("scaled_score" %in% names(result))
  expect_true("fisher_info" %in% names(result))
  expect_true("raw_score" %in% names(result))
  expect_length(result$scaled_score, K * (K - 1))
  expect_length(result$raw_score, K * (K - 1))
  expect_true(all(is.finite(result$scaled_score)))
})

test_that("calculate_gas_score still works for diagonal parameterization", {
  K <- 2
  mu <- c(-1, 1)
  sigma2 <- c(0.5, 1)
  X_pred <- c(0.6, 0.4)
  trans_prob <- c(0.8, 0.9)

  y <- rnorm(500)
  gh_setup <- setup_gauss_hermite_quadrature(y)

  result <- calculate_gas_score(
    y_obs = 0.5, mu = mu, sigma2 = sigma2,
    trans_prob = trans_prob, diag_probs = TRUE,
    X_pred = X_pred, gh_setup = gh_setup
  )

  expect_true(is.list(result))
  expect_length(result$scaled_score, K)
  expect_true(all(is.finite(result$scaled_score)))
})
