# Tests for softmax_row() utility function

test_that("softmax_row returns equal probabilities for zero parameters", {
  result <- multiregimeTVTP:::softmax_row(c(0, 0))
  expect_equal(result$offdiag, c(1/3, 1/3), tolerance = 1e-12)
  expect_equal(result$diag, 1/3, tolerance = 1e-12)
})

test_that("softmax_row probabilities sum to 1", {
  for (k in 2:5) {
    x <- rnorm(k - 1, sd = 2)
    result <- multiregimeTVTP:::softmax_row(x)
    total <- sum(result$offdiag) + result$diag
    expect_equal(total, 1, tolerance = 1e-12,
                 label = sprintf("K=%d, sum of probs", k))
  }
})

test_that("softmax_row produces all positive probabilities", {
  set.seed(42)
  for (i in 1:50) {
    x <- rnorm(sample(1:5, 1), sd = 3)
    result <- multiregimeTVTP:::softmax_row(x)
    expect_true(all(result$offdiag > 0), label = "offdiag > 0")
    expect_true(result$diag > 0, label = "diag > 0")
  }
})

test_that("softmax_row handles single parameter (K=2)", {
  result <- multiregimeTVTP:::softmax_row(0)
  expect_equal(result$offdiag, 0.5, tolerance = 1e-12)
  expect_equal(result$diag, 0.5, tolerance = 1e-12)
})

test_that("softmax_row gives high diagonal for large negative params", {
  result <- multiregimeTVTP:::softmax_row(c(-10, -10))
  expect_true(result$diag > 0.99)
  expect_true(all(result$offdiag < 0.005))
})

test_that("softmax_row gives low diagonal for large positive params", {
  result <- multiregimeTVTP:::softmax_row(c(10, 10))
  expect_true(result$diag < 0.001)
  expect_true(all(result$offdiag > 0.4))
})

test_that("softmax_row is numerically stable for extreme values", {
  # Very large positive
  result_large <- multiregimeTVTP:::softmax_row(c(500, 500))
  expect_true(all(is.finite(c(result_large$offdiag, result_large$diag))))
  expect_equal(sum(result_large$offdiag) + result_large$diag, 1, tolerance = 1e-10)

  # Very large negative
  result_small <- multiregimeTVTP:::softmax_row(c(-500, -500))
  expect_true(all(is.finite(c(result_small$offdiag, result_small$diag))))
  expect_equal(sum(result_small$offdiag) + result_small$diag, 1, tolerance = 1e-10)

  # Mixed extreme
  result_mixed <- multiregimeTVTP:::softmax_row(c(500, -500))
  expect_true(all(is.finite(c(result_mixed$offdiag, result_mixed$diag))))
  expect_equal(sum(result_mixed$offdiag) + result_mixed$diag, 1, tolerance = 1e-10)
})
