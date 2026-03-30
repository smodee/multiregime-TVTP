#!/usr/bin/env Rscript
# Benchmark R vs C filtering implementations
# Usage: Rscript scripts/benchmark_cpp.R

library(multiregimeTVTP)

if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  stop("Install microbenchmark: install.packages('microbenchmark')")
}
library(microbenchmark)

cat("=== R vs C Filtering Benchmark ===\n")
cat("Date:", format(Sys.time()), "\n")
cat("cpp_available:", cpp_available(), "\n\n")

if (!cpp_available()) stop("C backend not available")

set.seed(42)
n_eval <- 200

results <- list()

for (TT in c(200, 500, 1000)) {
  cat(sprintf("--- T = %d ---\n", TT))

  # Constant K=2
  par_c <- set_parameter_attributes(c(-1, 1, 1, 0.95, 0.90), K = 2,
                                     model_type = "constant", diag_probs = TRUE, equal_variances = TRUE)
  y <- dataConstCD(TT, 1, par_c)[, 1]

  mb <- microbenchmark(
    R = Rfiltering_Const(par_c, y, 5, 0, use_cpp = FALSE),
    C = Rfiltering_Const(par_c, y, 5, 0, use_cpp = TRUE),
    times = n_eval
  )
  r_ms <- median(mb$time[mb$expr == "R"]) / 1e6
  c_ms <- median(mb$time[mb$expr == "C"]) / 1e6
  cat(sprintf("  Const K=2:  R=%7.3f ms  C=%7.3f ms  speedup=%.1fx\n", r_ms, c_ms, r_ms / c_ms))

  # TVP K=2
  par_t <- set_parameter_attributes(c(-1, 1, 1, 0.95, 0.90, 0.05, 0.05), K = 2,
                                     model_type = "tvp", diag_probs = TRUE, equal_variances = TRUE)
  y_t <- dataTVPCD(TT, 1, par_t)[, 1]

  mb_t <- microbenchmark(
    R = Rfiltering_TVP(par_t, y_t, 5, 0, use_cpp = FALSE),
    C = Rfiltering_TVP(par_t, y_t, 5, 0, use_cpp = TRUE),
    times = n_eval
  )
  r_ms <- median(mb_t$time[mb_t$expr == "R"]) / 1e6
  c_ms <- median(mb_t$time[mb_t$expr == "C"]) / 1e6
  cat(sprintf("  TVP K=2:    R=%7.3f ms  C=%7.3f ms  speedup=%.1fx\n", r_ms, c_ms, r_ms / c_ms))

  # GAS K=2
  par_g <- set_parameter_attributes(c(-1, 1, 1, 0.95, 0.90, 0.3, 0.3, 0.5, 0.5), K = 2,
                                     model_type = "gas", diag_probs = TRUE, equal_variances = TRUE)
  y_g <- dataGASCD(TT, 1, par_g)[, 1]

  mb_g <- microbenchmark(
    R = { set.seed(1); Rfiltering_GAS(par_g, y_g, 5, 0, use_cpp = FALSE) },
    C = { set.seed(1); Rfiltering_GAS(par_g, y_g, 5, 0, use_cpp = TRUE) },
    times = min(n_eval, 50)
  )
  r_ms <- median(mb_g$time[mb_g$expr == "R"]) / 1e6
  c_ms <- median(mb_g$time[mb_g$expr == "C"]) / 1e6
  cat(sprintf("  GAS K=2:    R=%7.3f ms  C=%7.3f ms  speedup=%.1fx\n", r_ms, c_ms, r_ms / c_ms))
}

cat("\nDone.\n")
