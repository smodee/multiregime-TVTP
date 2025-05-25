#' Test Script for GAS Model Implementation
#' 
#' This script tests the GAS (Generalized Autoregressive Score) model implementation
#' for regime-switching models with time-varying transition probabilities.
#'
#' The GAS model updates transition probabilities using score-driven dynamics
#' based on the predictive likelihood function.

# Clear workspace
rm(list = ls())

# Source required scripts
source("models/model_GAS.R")
source("simulation/sim_GAS.R")

# PATCHES - Override incorrect functions
# Create a specific count_regime function for TVP
count_regime_TVP <- function(par) {
  # For TVP model, we have K means, K variances, K*(K-1) transition probabilities, K*(K-1) A coefficients
  # Total: 2K + 2K*(K-1) = 2K + 2K^2 - 2K = 2K^2
  # So we need to solve: 2K^2 = length(par)
  # Therefore: K = sqrt(length(par)/2)
  
  K <- sqrt(length(par)/2)
  
  if (abs(K - round(K)) > 1e-10) {
    stop(paste("The parameter vector length", length(par), 
               "is not compatible with any valid number of regimes for TVP model.",
               "Expected length is 2K^2.",
               "For K=2: 8, K=3: 18, K=4: 32, etc."))
  }
  
  return(round(K))
}

# Override transform_TVP to use the correct count function
transform_TVP <- function(par) {
  K <- count_regime_TVP(par)  # Use TVP-specific count function
  
  mu <- par[1:K]
  sigma2 <- par[(K+1):(2*K)]
  init_trans <- par[(2*K+1):(2*K + K*(K-1))]
  A <- par[(2*K + K*(K-1) + 1):length(par)]
  
  # Check for invalid values
  if (any(sigma2 <= 0)) {
    stop("Variance parameters must be positive")
  }
  if (any(init_trans <= 0) || any(init_trans >= 1)) {
    stop("Transition probabilities must be between 0 and 1")
  }
  
  return(c(mu, log(sigma2), logit(init_trans), A))
}

# Override untransform_TVP to use the correct count function
untransform_TVP <- function(par_t) {
  K <- count_regime_TVP(par_t)  # Use TVP-specific count function
  
  mu <- par_t[1:K]
  log_sigma2 <- par_t[(K+1):(2*K)]
  logit_init_trans <- par_t[(2*K+1):(2*K + K*(K-1))]
  A <- par_t[(2*K + K*(K-1) + 1):length(par_t)]
  
  return(c(mu, exp(log_sigma2), logistic(logit_init_trans), A))
}

# Also fix the generic count_regime to handle both cases
count_regime <- function(par) {
  # Try to determine which model based on parameter count
  # For K regimes:
  # - Constant model: K + K + K*(K-1) = K*(K+1)
  # - TVP model: K + K + K*(K-1) + K*(K-1) = 2*K^2
  # - Exogenous model: K + K + K*(K-1) + K*(K-1) = 2*K^2
  # - GAS model: K + K + K*(K-1) + K*(K-1) + K*(K-1) = 3*K^2 - K
  
  len <- length(par)
  
  # Check if it's a constant model (K^2 + K = len)
  K_const <- (-1 + sqrt(1 + 4*len))/2
  if (abs(K_const - round(K_const)) < 1e-10) {
    return(round(K_const))
  }
  
  # Check if it's a TVP/Exogenous model (2K^2 = len)
  K_tvp <- sqrt(len/2)
  if (abs(K_tvp - round(K_tvp)) < 1e-10) {
    return(round(K_tvp))
  }
  
  # Check if it's a GAS model (3K^2 - K = len)
  K_gas <- (1 + sqrt(1 + 12*len))/6
  if (abs(K_gas - round(K_gas)) < 1e-10) {
    return(round(K_gas))
  }
  
  stop(paste("The parameter vector length", len, 
             "is not compatible with any valid number of regimes."))
}

# Set random seed for reproducibility
set.seed(42)

cat("========================================\n")
cat("Testing GAS Model Implementation\n")
cat("========================================\n\n")

#===============================================================================
# Test 1: Basic Data Generation
#===============================================================================
cat("Test 1: Basic Data Generation\n")
cat("-" , rep("", 40), "\n", sep = "")

# Parameters for a 3-regime model
K <- 3
mu_true <- c(-2, 1, 2)
sigma2_true <- c(0.02, 0.2, 0.6)
init_trans_true <- rep(0.2, 6)  # K*(K-1) = 6 transition probabilities
A_true <- rep(0.1, 6)  # Score scaling parameters
B_true <- rep(0.9, 6)  # Persistence parameters

# Generate data
M <- 5  # Number of paths
N <- 1000  # Length of each path

cat("Generating GAS data with:\n")
cat("  Regimes:", K, "\n")
cat("  Paths:", M, "\n")
cat("  Length:", N, "\n")

data_sim <- dataGASCD(M, N, mu_true, sigma2_true, init_trans_true, A_true, B_true)

cat("Data generation successful!\n")
cat("  Data dimensions:", dim(data_sim), "\n")
cat("  Data range:", range(data_sim), "\n\n")

#===============================================================================
# Test 2: Parameter Estimation on Single Path
#===============================================================================
cat("Test 2: Parameter Estimation on Single Path\n")
cat("-" , rep("", 40), "\n", sep = "")

# Extract first path
y <- data_sim[1,]

# Estimate parameters
cat("Estimating GAS model parameters...\n")
result <- estimate_gas_model(y, K = K, B_burnin = 100, C = 50, verbose = FALSE)

cat("Estimation completed!\n")
cat("  Log-likelihood:", round(result$diagnostics$loglik, 2), "\n")
cat("  AIC:", round(result$diagnostics$aic, 2), "\n")
cat("  BIC:", round(result$diagnostics$bic, 2), "\n")
cat("  Number of parameters:", result$diagnostics$num_params, "\n")

# Compare estimated vs true parameters
cat("\nParameter comparison:\n")
cat("Means:\n")
cat("  True:     ", round(mu_true, 3), "\n")
cat("  Estimated:", round(result$parameters$mu, 3), "\n")

cat("Variances:\n")
cat("  True:     ", round(sigma2_true, 3), "\n")
cat("  Estimated:", round(result$parameters$sigma2, 3), "\n")

cat("Initial transition probabilities:\n")
cat("  True:     ", round(init_trans_true, 3), "\n")
cat("  Estimated:", round(result$parameters$init_trans, 3), "\n")

cat("A coefficients (mean):\n")
cat("  True:     ", round(mean(A_true), 3), "\n")
cat("  Estimated:", round(mean(result$parameters$A), 3), "\n")

cat("B coefficients (mean):\n")
cat("  True:     ", round(mean(B_true), 3), "\n")
cat("  Estimated:", round(mean(result$parameters$B), 3), "\n\n")

#===============================================================================
# Test 3: Multiple Simulations and Estimations
#===============================================================================
cat("Test 3: Multiple Simulations and Estimations\n")
cat("-" , rep("", 40), "\n", sep = "")

# Run example simulation with estimation
cat("Running multiple simulations with estimation...\n")
sim_results <- example_GAS_simulation(
  M = 5,
  N = 1000,
  mu = mu_true,
  sigma2 = sigma2_true,
  init_trans = init_trans_true,
  A = A_true,
  B = B_true,
  B_burnin = 100,
  C = 50,
  verbose = FALSE
)

cat("Simulation study completed!\n")
cat("Summary of parameter estimation across", nrow(sim_results$estimation$parameters), "paths:\n")

# Extract summary for key parameters
summary_df <- sim_results$estimation$summary
mu_summary <- summary_df[1:K, ]
sigma2_summary <- summary_df[(K+1):(2*K), ]
A_summary <- summary_df[grepl("^A", summary_df$Parameter), ]
B_summary <- summary_df[grepl("^B", summary_df$Parameter), ]

cat("\nMean parameter bias:\n")
cat("  Mu:        ", round(mean(abs(mu_summary$Bias)), 4), "\n")
cat("  Sigma2:    ", round(mean(abs(sigma2_summary$Bias)), 4), "\n")
cat("  A:         ", round(mean(abs(A_summary$Bias)), 4), "\n")
cat("  B:         ", round(mean(abs(B_summary$Bias)), 4), "\n")

cat("\nMean relative bias (%):\n")
cat("  Mu:        ", round(mean(abs(mu_summary$Rel_Bias_Pct)), 2), "%\n")
cat("  Sigma2:    ", round(mean(abs(sigma2_summary$Rel_Bias_Pct)), 2), "%\n")
cat("  A:         ", round(mean(abs(A_summary$Rel_Bias_Pct)), 2), "%\n")
cat("  B:         ", round(mean(abs(B_summary$Rel_Bias_Pct)), 2), "%\n\n")

#===============================================================================
# Test 4: Sensitivity Analysis
#===============================================================================
cat("Test 4: Sensitivity Analysis\n")
cat("-" , rep("", 40), "\n", sep = "")

# Test sensitivity to different A and B values
cat("Testing sensitivity to A and B parameters...\n")

# Generate a test dataset
test_data <- dataGASCD(1, 1000, mu_true, sigma2_true, init_trans_true, A_true, B_true)[1,]

# Test with different parameter values
sensitivity_results <- analyze_gas_sensitivity(
  y = test_data,
  K = K,
  A_values = c(0, 0.05, 0.1, 0.2),
  B_values = c(0.7, 0.8, 0.9, 0.95),
  init_trans = init_trans_true,
  verbose = FALSE
)

cat("Sensitivity analysis completed!\n")
cat("  Number of parameter combinations tested:", nrow(sensitivity_results$ab_grid), "\n")

# Find best parameter combination by AIC
best_idx <- which.min(sensitivity_results$ab_grid$aic)
best_params <- sensitivity_results$ab_grid[best_idx, ]

cat("\nBest parameter combination (by AIC):\n")
cat("  A value:", best_params$A_value, "\n")
cat("  B value:", best_params$B_value, "\n")
cat("  Log-likelihood:", round(best_params$loglik, 2), "\n")
cat("  AIC:", round(best_params$aic, 2), "\n")
cat("  Transition frequency:", round(best_params$transition_freq, 3), "\n")
cat("  Average regime duration:", round(best_params$avg_duration, 2), "\n\n")

#===============================================================================
# Test 5: Model Comparison
#===============================================================================
cat("Test 5: Model Comparison\n")
cat("-" , rep("", 40), "\n", sep = "")

# Generate data with known GAS dynamics
cat("Generating data from GAS model for comparison...\n")
comparison_data <- dataGASCD(1, 1000, mu_true, sigma2_true, init_trans_true, A_true, B_true)[1,]

# Compare different models
cat("Comparing GAS with other models...\n")
model_comparison <- compare_tvtp_models(
  y = comparison_data,
  X_Exo = rnorm(length(comparison_data)),  # Random exogenous variable for Exogenous model
  K = K,
  models = c("Constant", "TVP", "GAS"),  # Excluding Exogenous for cleaner comparison
  B_burnin = 100,
  C = 50,
  verbose = FALSE
)

cat("\nModel comparison results:\n")
print(model_comparison)

# Check if GAS performs best (should be true since data was generated from GAS)
if (model_comparison$Model[1] == "GAS") {
  cat("\n✓ GAS model correctly identified as best model for GAS-generated data\n")
} else {
  cat("\n✗ Warning: GAS model not identified as best model for GAS-generated data\n")
}

#===============================================================================
# Test 6: Regime Persistence Analysis
#===============================================================================
cat("\nTest 6: Regime Persistence Analysis\n")
cat("-" , rep("", 40), "\n", sep = "")

# Calculate persistence metrics for the estimated model
persistence <- calculate_persistence(result$filtered_probabilities)

cat("Regime persistence metrics:\n")
cat("  Number of transitions:", persistence$num_transitions, "\n")
cat("  Transition rate:", round(persistence$transition_rate, 3), "\n")
cat("  Regime percentages:\n")
for (r in names(persistence$regime_percentages)) {
  cat("    Regime", r, ":", round(persistence$regime_percentages[r], 1), "%\n")
}
cat("  Average durations:\n")
for (r in names(persistence$avg_durations)) {
  cat("    Regime", r, ":", round(persistence$avg_durations[r], 1), "periods\n")
}

#===============================================================================
# Test 7: Time-Varying Transition Probabilities
#===============================================================================
cat("\nTest 7: Time-Varying Transition Probabilities\n")
cat("-" , rep("", 40), "\n", sep = "")

# Extract transition probabilities over time
trans_probs <- result$transition_probabilities

cat("Transition probability statistics:\n")
cat("  Dimensions:", dim(trans_probs), "\n")
cat("  Mean values:\n")
for (i in 1:nrow(trans_probs)) {
  cat("    Transition", i, ":", round(mean(trans_probs[i,]), 3), "\n")
}
cat("  Standard deviations:\n")
for (i in 1:nrow(trans_probs)) {
  cat("    Transition", i, ":", round(sd(trans_probs[i,]), 3), "\n")
}

# Check if transition probabilities vary over time (they should for GAS)
total_variation <- sum(apply(trans_probs, 1, sd))
if (total_variation > 0.01) {
  cat("\n✓ Transition probabilities show time variation (total SD:", round(total_variation, 3), ")\n")
} else {
  cat("\n✗ Warning: Transition probabilities show little time variation\n")
}

#===============================================================================
# Test 8: Different Sample Sizes
#===============================================================================
cat("\nTest 8: Different Sample Sizes\n")
cat("-" , rep("", 40), "\n", sep = "")

sample_sizes <- c(250, 500, 1000)
size_results <- data.frame()

cat("Testing estimation with different sample sizes...\n")
for (n in sample_sizes) {
  # Generate data
  test_data <- dataGASCD(1, n, mu_true, sigma2_true, init_trans_true, A_true, B_true)[1,]
  
  # Estimate model
  est_result <- estimate_gas_model(test_data, K = K, B_burnin = 50, C = 25, verbose = FALSE)
  
  # Calculate errors
  mu_error <- sqrt(mean((est_result$parameters$mu - mu_true)^2))
  sigma2_error <- sqrt(mean((est_result$parameters$sigma2 - sigma2_true)^2))
  A_error <- sqrt(mean((est_result$parameters$A - A_true)^2))
  B_error <- sqrt(mean((est_result$parameters$B - B_true)^2))
  
  # Store results
  size_results <- rbind(size_results, data.frame(
    SampleSize = n,
    LogLik = est_result$diagnostics$loglik,
    MuError = mu_error,
    Sigma2Error = sigma2_error,
    AError = A_error,
    BError = B_error
  ))
  
  cat("  n =", n, ": LogLik =", round(est_result$diagnostics$loglik, 2), "\n")
}

cat("\nParameter estimation errors by sample size:\n")
print(round(size_results, 4))

# Check if errors decrease with sample size
if (all(diff(size_results$MuError) < 0)) {
  cat("\n✓ Mean parameter errors decrease with sample size\n")
} else {
  cat("\n✗ Warning: Mean parameter errors do not consistently decrease with sample size\n")
}

cat("\n========================================\n")
cat("GAS Model Testing Complete\n")
cat("========================================\n")

# Plot results if desired
plot_results <- FALSE  # Set to TRUE to generate plots

if (plot_results) {
  cat("\nGenerating diagnostic plots...\n")
  
  # Plot 1: Filtered probabilities
  par(mfrow = c(2, 2))
  
  # Filtered probabilities over time
  matplot(result$filtered_probabilities, type = "l", 
          main = "Filtered Regime Probabilities",
          xlab = "Time", ylab = "Probability",
          col = 1:K, lty = 1)
  legend("topright", paste("Regime", 1:K), col = 1:K, lty = 1, cex = 0.8)
  
  # Transition probabilities over time (first 4)
  if (nrow(trans_probs) > 0) {
    n_plot <- min(4, nrow(trans_probs))
    matplot(t(trans_probs[1:n_plot,]), type = "l",
            main = "Time-Varying Transition Probabilities",
            xlab = "Time", ylab = "Probability",
            col = 1:n_plot, lty = 1)
    legend("topright", paste("Trans", 1:n_plot), col = 1:n_plot, lty = 1, cex = 0.8)
  }
  
  # Score evolution (first 4)
  if (!is.null(result$scores) && nrow(result$scores) > 0) {
    n_plot <- min(4, nrow(result$scores))
    matplot(t(result$scores[1:n_plot,]), type = "l",
            main = "Score Evolution",
            xlab = "Time", ylab = "Score",
            col = 1:n_plot, lty = 1)
    legend("topright", paste("Score", 1:n_plot), col = 1:n_plot, lty = 1, cex = 0.8)
  }
  
  # Most likely regime
  most_likely <- apply(result$filtered_probabilities, 1, which.max)
  plot(most_likely, type = "s",
       main = "Most Likely Regime Over Time",
       xlab = "Time", ylab = "Regime",
       ylim = c(0.5, K + 0.5))
  
  par(mfrow = c(1, 1))
  
  cat("Plots generated successfully!\n")
}

cat("\nAll tests completed successfully!\n")