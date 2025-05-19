#' Test Script for GAS (Score-Driven) Regime-Switching Model
#'
#' This script tests the implementation of the GAS regime-switching model with
#' time-varying transition probabilities driven by the score of the predictive likelihood.
#'
#' Reference: Sendstad, Chronopoulos, & Li (2025) - The Value of Turning-Point 
#' Detection for Optimal Investment
#' Based on: Bazzi et al. (2017) - Time-Varying Transition Probabilities for Markov Regime
#' Switching Models

# Clear workspace
rm(list = ls())

# Set working directory to the project root if needed
# setwd("path/to/your/project")

# Load required source files
source("models/model_GAS.R")
source("simulation/sim_GAS.R")

# Set random seed for reproducibility
set.seed(42)

#-------------------------------------------------------------------------------
# Test 1: Simple 2-Regime Model Parameter Recovery
#-------------------------------------------------------------------------------
cat("\n========== Test 1: Simple 2-Regime Model Parameter Recovery ==========\n")

# Define true parameters for a simple 2-regime model
mu_true <- c(-1, 1)
sigma2_true <- c(0.5, 1.0)
init_trans_true <- c(0.2, 0.2)  # Off-diagonal elements for a 2x2 matrix
A_true <- c(0.1, 0.1)          # Score scaling parameters
B_true <- c(0.9, 0.9)          # Persistence parameters

# Generate data
N <- 1000
cat("Generating data with", N, "observations and 2 regimes...\n")
data_sim <- dataGASCD(1, N, mu_true, sigma2_true, init_trans_true, A_true, B_true)

# Estimate model
cat("Estimating model parameters...\n")
model_est <- estimate_gas_model(
  y = data_sim[1,],
  K = 2,
  B_burnin = 100,
  C = 50,
  verbose = TRUE
)

# Compare true vs. estimated parameters
cat("\nParameter Comparison:\n")
cat("Means:   True:", mu_true, "Estimated:", round(model_est$parameters$mu, 3), "\n")
cat("Var:     True:", sigma2_true, "Estimated:", round(model_est$parameters$sigma2, 3), "\n")
cat("Trans:   True:", init_trans_true, "Estimated:", round(model_est$parameters$init_trans, 3), "\n")
cat("A:       True:", A_true, "Estimated:", round(model_est$parameters$A, 3), "\n")
cat("B:       True:", B_true, "Estimated:", round(model_est$parameters$B, 3), "\n")

# Calculate RMSEs
mu_rmse <- sqrt(mean((model_est$parameters$mu - mu_true)^2))
sigma2_rmse <- sqrt(mean((model_est$parameters$sigma2 - sigma2_true)^2))
init_trans_rmse <- sqrt(mean((model_est$parameters$init_trans - init_trans_true)^2))
A_rmse <- sqrt(mean((model_est$parameters$A - A_true)^2))
B_rmse <- sqrt(mean((model_est$parameters$B - B_true)^2))

cat("\nRMSEs:\n")
cat("Mean RMSE:", mu_rmse, "\n")
cat("Variance RMSE:", sigma2_rmse, "\n")
cat("Trans Prob RMSE:", init_trans_rmse, "\n")
cat("A RMSE:", A_rmse, "\n")
cat("B RMSE:", B_rmse, "\n")

# Check if RMSEs are within acceptable bounds
tol <- 0.3  # Tolerance for parameter recovery
test1_pass <- all(c(mu_rmse, sigma2_rmse, init_trans_rmse, A_rmse, B_rmse) < tol)
cat("\nTest 1 Result:", ifelse(test1_pass, "PASS", "FAIL"), "\n")

# Plot filtered probabilities
if (test1_pass) {
  pdf("results/test_gas_model_2regime.pdf", width = 10, height = 6)
  par(mfrow = c(2, 1))
  
  # Plot the data
  plot(data_sim[1,], type = "l", main = "Simulated Data (2 Regimes)", 
       xlab = "Time", ylab = "Value")
  
  # Plot filtered probabilities
  matplot(model_est$filtered_probabilities, type = "l", main = "Filtered Regime Probabilities",
          xlab = "Time", ylab = "Probability", col = c("blue", "red"), lty = 1)
  legend("topright", legend = c("Regime 1", "Regime 2"), col = c("blue", "red"), lty = 1)
  
  dev.off()
  cat("Plots saved to results/test_gas_model_2regime.pdf\n")
}

#-------------------------------------------------------------------------------
# Test 2: 3-Regime Model Parameter Recovery
#-------------------------------------------------------------------------------
cat("\n========== Test 2: 3-Regime Model Parameter Recovery ==========\n")

# Define true parameters for a 3-regime model
mu_true <- c(-2, 0, 2)
sigma2_true <- c(0.2, 0.5, 1.0)
init_trans_true <- rep(0.2, 6)  # Off-diagonal elements for a 3x3 matrix
A_true <- rep(0.1, 6)          # Score scaling parameters
B_true <- rep(0.9, 6)          # Persistence parameters

# Generate data
N <- 2000  # Longer data for more complex model
cat("Generating data with", N, "observations and 3 regimes...\n")
data_sim <- dataGASCD(1, N, mu_true, sigma2_true, init_trans_true, A_true, B_true)

# Estimate model
cat("Estimating model parameters...\n")
model_est <- estimate_gas_model(
  y = data_sim[1,],
  K = 3,
  B_burnin = 100,
  C = 50,
  verbose = TRUE
)

# Compare true vs. estimated parameters
cat("\nParameter Comparison:\n")
cat("Means:   True:", mu_true, "Estimated:", round(model_est$parameters$mu, 3), "\n")
cat("Var:     True:", sigma2_true, "Estimated:", round(model_est$parameters$sigma2, 3), "\n")
cat("Trans:   True:", init_trans_true, "Estimated:", round(model_est$parameters$init_trans, 3), "\n")
cat("A (mean): True:", mean(A_true), "Estimated:", round(mean(model_est$parameters$A), 3), "\n")
cat("B (mean): True:", mean(B_true), "Estimated:", round(mean(model_est$parameters$B), 3), "\n")

# Calculate RMSEs
mu_rmse <- sqrt(mean((model_est$parameters$mu - mu_true)^2))
sigma2_rmse <- sqrt(mean((model_est$parameters$sigma2 - sigma2_true)^2))
init_trans_rmse <- sqrt(mean((model_est$parameters$init_trans - init_trans_true)^2))
A_rmse <- sqrt(mean((model_est$parameters$A - A_true)^2))
B_rmse <- sqrt(mean((model_est$parameters$B - B_true)^2))

cat("\nRMSEs:\n")
cat("Mean RMSE:", mu_rmse, "\n")
cat("Variance RMSE:", sigma2_rmse, "\n")
cat("Trans Prob RMSE:", init_trans_rmse, "\n")
cat("A RMSE:", A_rmse, "\n")
cat("B RMSE:", B_rmse, "\n")

# Check if RMSEs are within acceptable bounds
tol <- 0.4  # Slightly higher tolerance for more complex model
test2_pass <- all(c(mu_rmse, sigma2_rmse, init_trans_rmse, A_rmse, B_rmse) < tol)
cat("\nTest 2 Result:", ifelse(test2_pass, "PASS", "FAIL"), "\n")

# Plot filtered probabilities
if (test2_pass) {
  pdf("results/test_gas_model_3regime.pdf", width = 10, height = 8)
  par(mfrow = c(2, 1))
  
  # Plot the data
  plot(data_sim[1,], type = "l", main = "Simulated Data (3 Regimes)",
       xlab = "Time", ylab = "Value")
  
  # Plot filtered probabilities
  matplot(model_est$filtered_probabilities, type = "l", main = "Filtered Regime Probabilities",
          xlab = "Time", ylab = "Probability", col = c("blue", "green", "red"), lty = 1)
  legend("topright", legend = c("Regime 1", "Regime 2", "Regime 3"), 
         col = c("blue", "green", "red"), lty = 1)
  
  dev.off()
  cat("Plots saved to results/test_gas_model_3regime.pdf\n")
}

#-------------------------------------------------------------------------------
# Test 3: Varying A and B Impact
#-------------------------------------------------------------------------------
cat("\n========== Test 3: Impact of A and B Parameters ==========\n")

# Generate data with 2 regimes
N <- 1000
mu <- c(-1, 1)
sigma2 <- c(0.5, 1.0)
init_trans <- c(0.2, 0.2)

# Test different A and B combinations
cat("Testing impact of different A and B values...\n")

# Define grid of A and B values to test
A_values <- c(0, 0.1, 0.5)
B_values <- c(0.7, 0.9, 0.99)

# Create grid
param_grid <- expand.grid(A = A_values, B = B_values)
n_combinations <- nrow(param_grid)

# Initialize results
transition_rates <- numeric(n_combinations)
avg_durations <- numeric(n_combinations)

for (i in 1:n_combinations) {
  A_val <- param_grid$A[i]
  B_val <- param_grid$B[i]
  
  # Generate data with specific A and B values
  A_test <- rep(A_val, 2)
  B_test <- rep(B_val, 2)
  data_sim <- dataGASCD(1, N, mu, sigma2, init_trans, A_test, B_test)
  
  # Estimate model
  model_est <- estimate_gas_model(
    y = data_sim[1,],
    K = 2,
    B_burnin = 100,
    C = 50,
    verbose = FALSE
  )
  
  # Calculate persistence metrics
  persistence <- calculate_persistence(model_est$filtered_probabilities)
  
  # Store results
  transition_rates[i] <- persistence$transition_rate
  avg_durations[i] <- mean(unlist(persistence$avg_durations))
  
  cat(sprintf("A = %.2f, B = %.2f: Transition Rate = %.3f, Avg Duration = %.1f\n", 
              A_val, B_val, transition_rates[i], avg_durations[i]))
}

# Create result data frame
results <- data.frame(
  A = param_grid$A,
  B = param_grid$B,
  TransitionRate = transition_rates,
  AvgDuration = avg_durations
)

cat("\nResults summary:\n")
print(results)

# Plot results
pdf("results/test_gas_model_AB_impact.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

# Plot impact on transition rate
plot(results$A, results$TransitionRate, type = "n", 
     main = "Impact on Transition Rate",
     xlab = "A (Sensitivity)", ylab = "Transition Rate")

for (b in unique(results$B)) {
  subset_idx <- results$B == b
  points(results$A[subset_idx], results$TransitionRate[subset_idx], 
         col = which(unique(results$B) == b), pch = 16)
  lines(results$A[subset_idx], results$TransitionRate[subset_idx], 
        col = which(unique(results$B) == b))
}
legend("topleft", legend = paste("B =", unique(results$B)), 
       col = 1:length(unique(results$B)), lty = 1, pch = 16)

# Plot impact on average duration
plot(results$A, results$AvgDuration, type = "n", 
     main = "Impact on Average Duration",
     xlab = "A (Sensitivity)", ylab = "Average Duration")

for (b in unique(results$B)) {
  subset_idx <- results$B == b
  points(results$A[subset_idx], results$AvgDuration[subset_idx], 
         col = which(unique(results$B) == b), pch = 16)
  lines(results$A[subset_idx], results$AvgDuration[subset_idx], 
        col = which(unique(results$B) == b))
}
legend("topleft", legend = paste("B =", unique(results$B)), 
       col = 1:length(unique(results$B)), lty = 1, pch = 16)

dev.off()
cat("Plots saved to results/test_gas_model_AB_impact.pdf\n")

# Check if test passes (at least we should see A and B impact dynamics)
test3_pass <- length(unique(round(transition_rates, 2))) > 1 && 
              length(unique(round(avg_durations, 1))) > 1
cat("\nTest 3 Result:", ifelse(test3_pass, "PASS", "FAIL"), "\n")

#-------------------------------------------------------------------------------
# Test 4: Compare GAS with TVP and Constant
#-------------------------------------------------------------------------------
cat("\n========== Test 4: Model Comparison ==========\n")

# Load TVP and Constant models
source("models/model_TVP.R")
source("models/model_constant.R")

# Generate data with time-varying transition probabilities from GAS model
cat("Generating data with 3 regimes from GAS model...\n")
mu <- c(-2, 0, 2)
sigma2 <- c(0.2, 0.5, 1.0)
init_trans <- rep(0.2, 6)
A <- rep(0.2, 6)  # More sensitive to changes
B <- rep(0.9, 6)
N <- 1000

data_sim <- dataGASCD(1, N, mu, sigma2, init_trans, A, B)
y <- data_sim[1,]

# Estimate all three models
cat("Estimating all three models...\n")

gas_model <- estimate_gas_model(y, K = 3, verbose = FALSE)
tvp_model <- estimate_tvp_model(y, K = 3, verbose = FALSE)
const_model <- estimate_const_model(y, K = 3, verbose = FALSE)

# Compare model fit
gas_aic <- gas_model$diagnostics$aic
tvp_aic <- tvp_model$diagnostics$aic
const_aic <- const_model$diagnostics$aic

cat("\nModel Fit Comparison (AIC, lower is better):\n")
cat("GAS:     ", gas_aic, "\n")
cat("TVP:     ", tvp_aic, "\n")
cat("Constant:", const_aic, "\n")

# Check if GAS model performs better on GAS-generated data
test4_pass <- gas_aic < tvp_aic && gas_aic < const_aic
cat("\nTest 4 Result:", ifelse(test4_pass, "PASS", "FAIL"), "\n")

# Plot comparison
pdf("results/test_gas_model_comparison.pdf", width = 12, height = 8)
plot_model_comparison(
  empirical_data = y,
  gas_model = gas_model,
  tvp_model = tvp_model,
  constant_model = const_model
)
dev.off()
cat("Plots saved to results/test_gas_model_comparison.pdf\n")

#-------------------------------------------------------------------------------
# Test 5: Numeric Stability
#-------------------------------------------------------------------------------
cat("\n========== Test 5: Numeric Stability ==========\n")

# Test with different starting values
cat("Testing numeric stability with different starting values...\n")

# Generate data
mu_true <- c(-1, 1)
sigma2_true <- c(0.5, 1.0)
init_trans_true <- c(0.2, 0.2)
A_true <- c(0.1, 0.1)
B_true <- c(0.9, 0.9)
N <- 1000

data_sim <- dataGASCD(1, N, mu_true, sigma2_true, init_trans_true, A_true, B_true)
y <- data_sim[1,]

# Try different starting values
starting_points <- list(
  default = NULL,  # Default initialization
  close = c(mu_true * 0.9, sigma2_true * 1.1, init_trans_true * 0.9, A_true * 0.9, B_true * 0.9),
  far = c(c(0, 2), c(0.2, 1.5), c(0.1, 0.1), c(0.2, 0.2), c(0.8, 0.8))
)

results <- list()

for (name in names(starting_points)) {
  cat("Testing starting point:", name, "\n")
  
  # Try to estimate model with these starting values
  tryCatch({
    model <- estimate_gas_model(
      y = y,
      K = 2,
      B_burnin = 100,
      C = 50,
      initial_params = starting_points[[name]],
      verbose = FALSE
    )
    
    # Store results
    results[[name]] <- list(
      success = TRUE,
      loglik = model$diagnostics$loglik,
      aic = model$diagnostics$aic,
      mu = model$parameters$mu,
      sigma2 = model$parameters$sigma2,
      init_trans = model$parameters$init_trans,
      A = model$parameters$A,
      B = model$parameters$B
    )
    
    cat("  Successful. LogLik:", round(model$diagnostics$loglik, 2), "\n")
    
  }, error = function(e) {
    cat("  Failed:", e$message, "\n")
    results[[name]] <- list(
      success = FALSE,
      error = e$message
    )
  })
}

# Check if all succeeded
all_succeeded <- all(sapply(results, function(x) x$success))

# If successful, check if results are similar
if (all_succeeded) {
  # Compare log-likelihoods
  logliks <- sapply(results, function(x) x$loglik)
  loglik_range <- max(logliks) - min(logliks)
  
  # Compare parameter estimates (using mu as example)
  mu_rmse <- sapply(names(results), function(name) {
    if (name == "default") return(0)
    sqrt(mean((results[[name]]$mu - results$default$mu)^2))
  })
  
  cat("\nResults comparison:\n")
  cat("LogLik range:", loglik_range, "\n")
  cat("Mu RMSE vs default:", mu_rmse, "\n")
  
  # Pass if log-likelihoods are similar
  test5_pass <- loglik_range < 1.0  # Tolerance for loglik differences
} else {
  test5_pass <- FALSE
}

cat("\nTest 5 Result:", ifelse(test5_pass, "PASS", "FAIL"), "\n")

#-------------------------------------------------------------------------------
# Test 6: Performance on Long Time Series
#-------------------------------------------------------------------------------
cat("\n========== Test 6: Performance on Long Time Series ==========\n")

# Define true parameters
mu <- c(-1, 1)
sigma2 <- c(0.5, 1.0)
init_trans <- c(0.2, 0.2)
A <- c(0.1, 0.1)
B <- c(0.9, 0.9)

# Generate longer time series
N <- 5000
cat("Generating data with", N, "observations...\n")
data_sim <- dataGASCD(1, N, mu, sigma2, init_trans, A, B)
y <- data_sim[1,]

# Measure time to estimate model
cat("Estimating model on long time series...\n")
start_time <- Sys.time()
model_est <- estimate_gas_model(
  y = y,
  K = 2,
  B_burnin = 100,
  C = 50,
  verbose = FALSE
)
end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("Estimation completed in", round(elapsed, 2), "seconds\n")

# Check performance
test6_pass <- !is.na(model_est$diagnostics$loglik) && elapsed < 60  # Should finish in under 60 seconds
cat("\nTest 6 Result:", ifelse(test6_pass, "PASS", "FAIL"), "\n")

#-------------------------------------------------------------------------------
# Summarize all test results
#-------------------------------------------------------------------------------
cat("\n========== Summary of All Tests ==========\n")
all_tests <- c(test1_pass, test2_pass, test3_pass, test4_pass, test5_pass, test6_pass)
test_names <- c(
  "Simple 2-Regime Parameter Recovery",
  "3-Regime Parameter Recovery",
  "A and B Parameter Impact",
  "Model Comparison",
  "Numeric Stability",
  "Performance on Long Time Series"
)

results_df <- data.frame(
  Test = test_names,
  Result = ifelse(all_tests, "PASS", "FAIL")
)

print(results_df)

all_pass <- all(all_tests)
cat("\nOverall Result:", ifelse(all_pass, "ALL TESTS PASSED", "SOME TESTS FAILED"), "\n")

# Create a directory for test results if it doesn't exist
dir.create("results", showWarnings = FALSE)

# Save summary to file
write.csv(results_df, "results/gas_model_test_summary.csv", row.names = FALSE)
cat("Summary saved to results/gas_model_test_summary.csv\n")
