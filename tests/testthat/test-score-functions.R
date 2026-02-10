# Tests for GAS score functions (softmax Jacobian)

test_that("raw score vector has correct length", {
  K <- 3
  eta <- c(0.1, 0.8, 0.1)
  X_t_prev <- c(0.3, 0.4, 0.3)
  p_trans <- rep(0, 6)
  tot_lik <- sum(eta * X_t_prev)

  score <- multiregimeTVTP:::calculate_raw_score_vector(eta, tot_lik, X_t_prev, p_trans, K)
  expect_length(score, K * (K - 1))
})

test_that("raw score is zero when all regimes have equal likelihoods", {
  K <- 3
  eta <- c(1, 1, 1)  # Equal likelihoods
  X_t_prev <- c(1/3, 1/3, 1/3)
  p_trans <- rep(0, 6)  # Uniform transition matrix
  tot_lik <- sum(eta * X_t_prev)

  score <- multiregimeTVTP:::calculate_raw_score_vector(eta, tot_lik, X_t_prev, p_trans, K)

  # When all eta are equal, eta[j] - row_lik_i = eta[j] - eta[j] = 0 for uniform P
  expect_equal(as.numeric(score), rep(0, 6), tolerance = 1e-12)
})

test_that("raw score has correct sign structure", {
  K <- 2
  # Regime 2 has much higher likelihood
  eta <- c(0.01, 0.99)
  X_t_prev <- c(0.5, 0.5)
  p_trans <- c(0, 0)  # Equal transition probs
  tot_lik <- sum(eta * X_t_prev)

  score <- multiregimeTVTP:::calculate_raw_score_vector(eta, tot_lik, X_t_prev, p_trans, K)
  # p_trans[1] is x_{1,2}: transition from regime 1 to regime 2
  # eta[2] - row_lik_1 should be positive (regime 2 has higher likelihood)
  # So score[1] > 0 (pushing towards more transition to regime 2)
  expect_true(score[1] > 0, label = "Score for x_12 positive when regime 2 is favored")
})

test_that("raw score rejects invalid inputs", {
  K <- 2
  eta <- c(0.5, 0.5)
  X_t_prev <- c(0.5, 0.5)
  p_trans <- c(0, 0)
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

test_that("raw score accepts unconstrained softmax params (negative, zero, positive)", {
  K <- 2
  eta <- c(0.3, 0.7)
  X_t_prev <- c(0.5, 0.5)
  tot_lik <- 0.5

  # Large negative params (high persistence)
  score1 <- multiregimeTVTP:::calculate_raw_score_vector(eta, tot_lik, X_t_prev, c(-10, -10), K)
  expect_true(all(is.finite(score1)))

  # Large positive params (high switching)
  score2 <- multiregimeTVTP:::calculate_raw_score_vector(eta, tot_lik, X_t_prev, c(10, 10), K)
  expect_true(all(is.finite(score2)))
})

# --- Fisher Information ---

test_that("Fisher Information is positive and finite", {
  K <- 2
  mu <- c(-1, 1)
  sigma2 <- c(0.5, 1)
  X_t_prev <- c(0.5, 0.5)
  p_trans <- c(0, 0)

  y <- rnorm(500)
  gh_setup <- setup_gauss_hermite_quadrature(y)

  fi <- multiregimeTVTP:::calculate_fisher_information(mu, sigma2, X_t_prev, p_trans, gh_setup, K)
  expect_true(is.finite(fi))
  expect_true(fi > 0)
})

test_that("Fisher Information is larger when regimes are well-separated", {
  K <- 2
  sigma2 <- c(0.5, 0.5)
  X_t_prev <- c(0.5, 0.5)
  p_trans <- c(0, 0)

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
  trans_prob <- c(-0.5, 0.3)

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
