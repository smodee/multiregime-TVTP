# Tests for generate_starting_points with softmax parameters

test_that("starting points have correct count for off-diagonal K=2", {
  set.seed(42)
  sp_list <- generate_starting_points(y = rnorm(100), K = 2, model_type = "constant",
                                      n_starts = 3,
                                      diag_probs = FALSE, equal_variances = TRUE)
  expect_true(is.list(sp_list))
  expect_length(sp_list, 3)
  # Each starting point should be a numeric vector
  expect_true(is.numeric(sp_list[[1]]))
})

test_that("off-diagonal starting points produce valid transition matrices", {
  set.seed(42)
  sp_list <- generate_starting_points(y = rnorm(200), K = 3, model_type = "constant",
                                      n_starts = 10,
                                      diag_probs = FALSE, equal_variances = TRUE)
  for (i in seq_along(sp_list)) {
    sp <- sp_list[[i]]
    # For K=3 constant, equal_var: 3 mu + 1 sigma2 + 6 trans = 10
    # Set attributes to extract trans_prob
    sp_attr <- set_parameter_attributes(sp, K = 3, model_type = "constant",
                                        diag_probs = FALSE, equal_variances = TRUE)
    trans_prob <- extract_parameter_component(sp_attr, "trans_prob")
    P <- transition_matrix(trans_prob, diag_probs = FALSE)

    expect_equal(rowSums(P), rep(1, 3), tolerance = 1e-10,
                 label = sprintf("trial %d row sums", i))
    expect_true(all(P > 0), label = sprintf("trial %d positive", i))
  }
})

test_that("diagonal starting points still work correctly", {
  set.seed(42)
  sp_list <- generate_starting_points(y = rnorm(100), K = 2, model_type = "constant",
                                      n_starts = 3,
                                      diag_probs = TRUE, equal_variances = TRUE)
  expect_true(is.list(sp_list))
  # Each should be a numeric vector
  for (sp in sp_list) {
    sp_attr <- set_parameter_attributes(sp, K = 2, model_type = "constant",
                                        diag_probs = TRUE, equal_variances = TRUE)
    trans_prob <- extract_parameter_component(sp_attr, "trans_prob")
    expect_true(all(trans_prob > 0 & trans_prob < 1))
  }
})

test_that("GAS model starting points have correct parameter count", {
  set.seed(42)
  sp_list <- generate_starting_points(y = rnorm(100), K = 2, model_type = "gas",
                                      n_starts = 3,
                                      diag_probs = FALSE, equal_variances = TRUE)
  # K=2 GAS off-diag equal_var: 2 mu + 1 sigma2 + 2 trans + 2 A + 2 B = 9
  expect_length(sp_list[[1]], 9)
})
