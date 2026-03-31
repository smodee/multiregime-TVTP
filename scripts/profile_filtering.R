#!/usr/bin/env Rscript
# Profile filtering functions to establish baseline performance
# Usage: Rscript scripts/profile_filtering.R
#
# This script measures single NLL call times for each model type
# before C++ implementation, providing baselines for speedup comparison.

library(multiregimeTVTP)

cat("=== Filtering Function Profiling ===\n")
cat("Date:", format(Sys.time()), "\n")
cat("R version:", R.version.string, "\n\n")

# Check for microbenchmark
if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  stop("Install microbenchmark: install.packages('microbenchmark')")
}
library(microbenchmark)

set.seed(42)

# --- Generate test data for each model type ---

# K=2, T=500 baseline
K <- 2
TT <- 500

cat("Generating K=2, T=500 test data...\n")

# Constant model data
const_data <- dataConstCD(T = TT, K = K, mu = c(-1, 1), sigma2 = c(1, 1),
                          P_diag = c(0.95, 0.90))
const_par <- set_parameter_attributes(
  c(-1, 1, 1, 1, 0.95, 0.90), K = K, model_type = "constant",
  diag_probs = TRUE, equal_variances = TRUE
)

# TVP model data
tvp_data <- dataTVPCD(T = TT, K = K, mu = c(-1, 1), sigma2 = c(1, 1),
                      P_diag = c(0.95, 0.90), A = c(0.05, 0.05))
tvp_par <- set_parameter_attributes(
  c(-1, 1, 1, 1, 0.95, 0.90, 0.05, 0.05), K = K, model_type = "tvp",
  diag_probs = TRUE, equal_variances = TRUE
)

# Exogenous model data
X_Exo <- rnorm(TT)
exo_data <- dataTVPXExoCD(T = TT, K = K, mu = c(-1, 1), sigma2 = c(1, 1),
                           P_diag = c(0.95, 0.90), A = c(0.05, 0.05),
                           X_Exo = X_Exo)
exo_par <- set_parameter_attributes(
  c(-1, 1, 1, 1, 0.95, 0.90, 0.05, 0.05), K = K, model_type = "exogenous",
  diag_probs = TRUE, equal_variances = TRUE
)

# GAS model data
gas_data <- dataGASCD(T = TT, K = K, mu = c(-1, 1), sigma2 = c(1, 1),
                      P_diag = c(0.95, 0.90), A = c(0.3, 0.3), B = c(0.5, 0.5))
gas_par <- set_parameter_attributes(
  c(-1, 1, 1, 1, 0.95, 0.90, 0.3, 0.3, 0.5, 0.5), K = K, model_type = "gas",
  diag_probs = TRUE, equal_variances = TRUE
)

# --- Benchmark single NLL calls ---

n_burnin <- 5
n_cutoff <- 0
n_eval <- 200  # Number of benchmark iterations

cat("\nBenchmarking single NLL calls (", n_eval, " iterations each)...\n\n")

results <- list()

cat("1. Constant model...\n")
results$constant <- microbenchmark(
  Rfiltering_Const(const_par, const_data$y, n_burnin, n_cutoff),
  times = n_eval
)
cat("   Median:", round(median(results$constant$time) / 1e6, 3), "ms\n")

cat("2. TVP model...\n")
results$tvp <- microbenchmark(
  Rfiltering_TVP(tvp_par, tvp_data$y, n_burnin, n_cutoff),
  times = n_eval
)
cat("   Median:", round(median(results$tvp$time) / 1e6, 3), "ms\n")

cat("3. Exogenous model...\n")
results$exo <- microbenchmark(
  Rfiltering_TVPXExo(exo_par, exo_data$y, X_Exo, n_burnin, n_cutoff),
  times = n_eval
)
cat("   Median:", round(median(results$exo$time) / 1e6, 3), "ms\n")

cat("4. GAS model...\n")
results$gas <- microbenchmark(
  Rfiltering_GAS(gas_par, gas_data$y, n_burnin, n_cutoff),
  times = n_eval
)
cat("   Median:", round(median(results$gas$time) / 1e6, 3), "ms\n")

# --- Summary table ---

cat("\n=== Baseline Summary (K=2, T=500) ===\n")
cat(sprintf("%-12s %10s %10s %10s\n", "Model", "Median(ms)", "Mean(ms)", "SD(ms)"))
cat(paste(rep("-", 44), collapse = ""), "\n")
for (name in names(results)) {
  times_ms <- results[[name]]$time / 1e6
  cat(sprintf("%-12s %10.3f %10.3f %10.3f\n",
              name, median(times_ms), mean(times_ms), sd(times_ms)))
}

cat("\nTypical NLL evaluations per estimation: 50,000 - 600,000\n")
cat("Estimated full estimation time (at 100k evals):\n")
for (name in names(results)) {
  med_ms <- median(results[[name]]$time) / 1e6
  est_time <- med_ms * 100000 / 1000  # seconds
  cat(sprintf("  %-12s %8.1f seconds (%4.1f minutes)\n", name, est_time, est_time / 60))
}

cat("\nDone.\n")
