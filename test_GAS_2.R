# test_gas_model.R
#
# This script tests the GAS (Generalized Autoregressive Score) regime-switching model
# with simulation data.
#
# Author: Samuel Modée
# Date: May 25, 2025

# Clear environment
rm(list = ls())

# Load required libraries and model implementation
source("helpers/utility_functions.R")
source("helpers/transition_helpers.R")
source("helpers/parameter_transforms.R")
source("models/model_GAS.R")

# Create a function to handle logging
log_message <- function(section = NULL, message, quiet = FALSE) {
  if (!quiet) {
    if (!is.null(section)) {
      cat("\n===", section, "===\n")
    }
    cat(message, "\n")
  }
}

# Set random seed for reproducibility
set.seed(42)

# Set up parameters
M <- 10  # Number of simulation paths
N <- 1000  # Length of each simulation path
B_burnin <- 200  # Burn-in period
C <- 50  # Cut-off period

# Set true parameter values
mu_true <- c(-2, -0.3, 2)  # Different means for each regime
sigma2_true <- c(0.23, 0.17, 0.6)  # Different variances for each regime
init_trans_true <- c(0.15, 0.19, 0.1, 0.2, 0.3, 0.05)  # Initial transition probabilities
A_true <- c(0.1, 0.15, 0.05, 0.08, 0.12, 0.06)  # Score scaling parameters (sensitivity)
B_true <- c(0.85, 0.9, 0.8, 0.88, 0.92, 0.87)  # Persistence parameters (memory)
K <- length(mu_true)  # Number of regimes

log_message("CONFIGURATION", sprintf(
  "K=%d regimes, M=%d paths, N=%d observations
   Means: %s
   Variances: %s
   Initial transition probabilities: %s
   A coefficients (sensitivity): %s
   B coefficients (persistence): %s",
  K, M, N,
  paste(round(mu_true, 2), collapse=", "),
  paste(round(sigma2_true, 2), collapse=", "),
  paste(round(init_trans_true, 2), collapse=", "),
  paste(round(A_true, 2), collapse=", "),
  paste(round(B_true, 2), collapse=", ")
))

# Generate simulation data
log_message("DATA GENERATION", "Generating simulation data...")
sim_start_time <- Sys.time()

# Generate data with burn-in and cut-off periods
sim_data <- dataGASCD(M, N+B_burnin+C, mu_true, sigma2_true, init_trans_true, A_true, B_true)

sim_end_time <- Sys.time()
log_message(NULL, sprintf("Data generation completed in %s", 
                          format(difftime(sim_end_time, sim_start_time), digits=4)))

# Storage for estimation results
n_transition <- K*(K-1)
param_estimates <- matrix(0, M, 2*K + 3*n_transition)
colnames(param_estimates) <- c(paste0("mu", 1:K), 
                               paste0("sigma2", 1:K), 
                               paste0("trans", 1:n_transition), 
                               paste0("A", 1:n_transition),
                               paste0("B", 1:n_transition))
diagnostics <- matrix(0, M, 3)
colnames(diagnostics) <- c("loglik", "aic", "bic")

# Estimate parameters for each simulation path
log_message("PARAMETER ESTIMATION", sprintf("Starting parameter estimation for %d simulation paths...", M))
total_start_time <- Sys.time()

# Progress bar function
progress_bar <- function(current, total, width=50) {
  percent <- round(current / total * 100)
  filled <- round(width * current / total)
  bar <- paste0("[", paste(rep("=", filled), collapse=""), 
                paste(rep(" ", width - filled), collapse=""), "]")
  return(sprintf("\r%s %3d%% (%d/%d)", bar, percent, current, total))
}

for (i in 1:M) {
  cat(progress_bar(i, M))
  
  # Extract the current simulation path
  current_data <- sim_data[i, ]
  
  # Initial parameter guesses
  mu_guess <- mu_true * 0.9 + rnorm(K, 0, 0.1*abs(mu_true))
  sigma2_guess <- sigma2_true * 0.9 + rnorm(K, 0, 0.1*sigma2_true)
  init_trans_guess <- pmax(pmin(init_trans_true + rnorm(n_transition, 0, 0.05), 0.95), 0.05)
  A_guess <- rep(0.1, n_transition)  # Start with moderate sensitivity
  B_guess <- rep(0.8, n_transition)  # Start with high persistence
  
  par_guess <- c(mu_guess, sigma2_guess, init_trans_guess, A_guess, B_guess)
  
  # Set parameter bounds
  lower_bounds <- c(rep(-Inf, K),               # No bounds on means
                    rep(-Inf, K),               # Variance >= 0 (log-transformed)
                    rep(-Inf, n_transition),    # Probabilities >= 0 (logit-transformed)
                    rep(-Inf, n_transition),    # A-coefficients >= 0 (logit-transformed)
                    rep(-Inf, n_transition))    # B-coefficients >= 0 (logit-transformed)
  
  upper_bounds <- c(rep(Inf, K),                # No bounds on means
                    rep(Inf, K),                # Variance unbounded
                    rep(Inf, n_transition),     # Probabilities <= 1 (logit-transformed)
                    rep(Inf, n_transition),     # A-coefficients <= 1 (logit-transformed)
                    rep(Inf, n_transition))     # B-coefficients <= 1 (logit-transformed)
  
  bounds <- list(lower = lower_bounds, upper = upper_bounds)
  
  # Estimate the model
  tryCatch({
    estimate <- estimate_gas_model(
      y = current_data,
      K = K,
      B_burnin = B_burnin,
      C = C,
      initial_params = par_guess,
      bounds = bounds,
      verbose = FALSE
    )
    
    # Extract estimated parameters
    param_estimates[i,] <- c(
      estimate$parameters$mu,
      estimate$parameters$sigma2,
      estimate$parameters$init_trans,
      estimate$parameters$A,
      estimate$parameters$B
    )
    
    # Store diagnostics
    diagnostics[i, 1] <- estimate$diagnostics$loglik
    diagnostics[i, 2] <- estimate$diagnostics$aic
    diagnostics[i, 3] <- estimate$diagnostics$bic
    
  }, error = function(e) {
    log_message(NULL, paste("ERROR in path", i, ":", e$message))
    param_estimates[i,] <- NA
    diagnostics[i,] <- NA
  })
}

# End progress bar with a newline
cat("\n")

total_end_time <- Sys.time()
total_time <- difftime(total_end_time, total_start_time, units = "mins")
log_message(NULL, sprintf("Parameter estimation completed in %s minutes", 
                          format(total_time, digits=4)))

# Calculate summary statistics for parameter estimates
estimates_summary <- data.frame(
  Parameter = colnames(param_estimates),
  True_Value = c(mu_true, sigma2_true, init_trans_true, A_true, B_true),
  Mean_Estimate = colMeans(param_estimates, na.rm = TRUE),
  SD_Estimate = apply(param_estimates, 2, sd, na.rm = TRUE),
  Bias = colMeans(param_estimates, na.rm = TRUE) - c(mu_true, sigma2_true, init_trans_true, A_true, B_true),
  Rel_Bias_Pct = 100 * (colMeans(param_estimates, na.rm = TRUE) - c(mu_true, sigma2_true, init_trans_true, A_true, B_true)) / 
    ifelse(abs(c(mu_true, sigma2_true, init_trans_true, A_true, B_true)) > 1e-10, 
           c(mu_true, sigma2_true, init_trans_true, A_true, B_true), 
           1e-10)
)

# Print summary statistics
log_message("RESULTS", "Parameter Estimation Summary:")
print(estimates_summary)

log_message(NULL, sprintf("Mean Log-Likelihood: %.2f", mean(diagnostics[,1], na.rm = TRUE)))
log_message(NULL, sprintf("Mean AIC: %.2f", mean(diagnostics[,2], na.rm = TRUE)))
log_message(NULL, sprintf("Mean BIC: %.2f", mean(diagnostics[,3], na.rm = TRUE)))

# Additional GAS-specific diagnostics
successful_runs <- !is.na(diagnostics[,1])
if (sum(successful_runs) > 0) {
  # Calculate regime persistence metrics for the first successful run
  first_success_idx <- which(successful_runs)[1]
  
  tryCatch({
    # Re-estimate the first successful model to get detailed output
    test_data <- sim_data[first_success_idx, ]
    test_estimate <- estimate_gas_model(
      y = test_data,
      K = K,
      B_burnin = B_burnin,
      C = C,
      verbose = FALSE
    )
    
    # Calculate persistence metrics
    persistence <- calculate_persistence(test_estimate$filtered_probabilities)
    
    log_message("GAS MODEL DYNAMICS", sprintf(
      "Regime dynamics for path %d:
       Number of transitions: %d
       Transition rate: %.3f
       Average regime durations: %s
       Regime occupancy percentages: %s",
      first_success_idx,
      persistence$num_transitions,
      persistence$transition_rate,
      paste(round(persistence$avg_durations, 2), collapse=", "),
      paste(round(persistence$regime_percentages, 1), collapse="%, ") 
    ))
    
  }, error = function(e) {
    log_message(NULL, paste("Could not calculate persistence metrics:", e$message))
  })
}

# Save results
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

# Create a timestamp ID in format YYYY-MM-DD-HHMMSS
timestamp_id <- format(Sys.time(), "%Y-%m-%d-%H%M%S")
results_filename <- paste0("gas_model_test_results_", timestamp_id, ".rds")

saveRDS(list(
  parameters = list(
    mu = mu_true,
    sigma2 = sigma2_true,
    init_trans = init_trans_true,
    A = A_true,
    B = B_true
  ),
  estimation = list(
    parameters = param_estimates,
    diagnostics = diagnostics,
    summary = estimates_summary
  ),
  settings = list(
    M = M,
    N = N,
    K = K,
    B_burnin = B_burnin,
    C = C
  )
), file = file.path(results_dir, results_filename))

# Create plots
plot_filename <- paste0("gas_model_test_plots_", timestamp_id, ".pdf")
pdf(file.path(results_dir, plot_filename), width = 12, height = 10)

# Plot 1: Sample path and regime classification
par(mfrow = c(2, 1))
plot(sim_data[1,], type = "l", main = "Sample Path from GAS Model",
     xlab = "Time", ylab = "Value", col = "darkblue")

# Add regime classification if available
if (exists("test_estimate") && !is.null(test_estimate$filtered_probabilities)) {
  # Get the most likely regime at each point
  regime_probs <- test_estimate$filtered_probabilities
  regime_classifications <- apply(regime_probs, 1, which.max)
  
  # Plot regime classification
  plot(regime_classifications, type = "l", main = "Estimated Regime Classification (GAS Model)",
       xlab = "Time", ylab = "Regime", ylim = c(1, K), yaxt = "n", col = "darkred", lwd = 2)
  axis(2, at = 1:K, labels = 1:K)
  
  # Add colored background for regimes
  for (r in 1:K) {
    regime_times <- which(regime_classifications == r)
    if (length(regime_times) > 0) {
      points(regime_times, rep(r, length(regime_times)), 
             col = rainbow(K)[r], pch = 16, cex = 0.3)
    }
  }
}

# Plot 2: Parameter recovery
par(mfrow = c(2, 3))

# Means
boxplot(param_estimates[,1:K], main = "Estimated Means",
        names = paste0("μ", 1:K), col = "lightblue",
        ylim = range(c(param_estimates[,1:K], mu_true), na.rm = TRUE))
points(1:K, mu_true, pch = 16, col = "red", cex = 1.5)
legend("topright", legend = "True value", pch = 16, col = "red")

# Variances
boxplot(param_estimates[,(K+1):(2*K)], main = "Estimated Variances",
        names = paste0("σ²", 1:K), col = "lightgreen",
        ylim = range(c(param_estimates[,(K+1):(2*K)], sigma2_true), na.rm = TRUE))
points(1:K, sigma2_true, pch = 16, col = "red", cex = 1.5)

# Transition probabilities
boxplot(param_estimates[,(2*K+1):(2*K+n_transition)], 
        main = "Estimated Transition Probabilities",
        names = paste0("p", 1:n_transition), col = "lightyellow",
        ylim = range(c(param_estimates[,(2*K+1):(2*K+n_transition)], init_trans_true), na.rm = TRUE),
        las = 2)
points(1:n_transition, init_trans_true, pch = 16, col = "red", cex = 1.5)

# A coefficients (sensitivity)
boxplot(param_estimates[,(2*K+n_transition+1):(2*K+2*n_transition)], 
        main = "Estimated A Coefficients (Sensitivity)",
        names = paste0("A", 1:n_transition), col = "lightpink",
        ylim = range(c(param_estimates[,(2*K+n_transition+1):(2*K+2*n_transition)], A_true), na.rm = TRUE),
        las = 2)
points(1:n_transition, A_true, pch = 16, col = "red", cex = 1.5)

# B coefficients (persistence)
boxplot(param_estimates[,(2*K+2*n_transition+1):(2*K+3*n_transition)], 
        main = "Estimated B Coefficients (Persistence)",
        names = paste0("B", 1:n_transition), col = "lightcyan",
        ylim = range(c(param_estimates[,(2*K+2*n_transition+1):(2*K+3*n_transition)], B_true), na.rm = TRUE),
        las = 2)
points(1:n_transition, B_true, pch = 16, col = "red", cex = 1.5)

# Model diagnostics
plot(diagnostics[,1], main = "Log-Likelihood Distribution",
     xlab = "Simulation Run", ylab = "Log-Likelihood", 
     col = "darkgreen", pch = 16)
abline(h = mean(diagnostics[,1], na.rm = TRUE), col = "red", lty = 2)

# Plot 3: GAS-specific analysis
if (exists("test_estimate")) {
  par(mfrow = c(2, 2))
  
  # Filtered probabilities over time
  if (!is.null(test_estimate$filtered_probabilities)) {
    matplot(test_estimate$filtered_probabilities, type = "l", 
            main = "Filtered Regime Probabilities Over Time",
            xlab = "Time", ylab = "Probability",
            col = rainbow(K), lty = 1, lwd = 2)
    legend("topright", legend = paste("Regime", 1:K), 
           col = rainbow(K), lty = 1, lwd = 2)
  }
  
  # Transition probabilities over time (show a subset)
  if (!is.null(test_estimate$transition_probabilities)) {
    # Show first few transition probabilities to avoid clutter
    n_show <- min(4, ncol(test_estimate$transition_probabilities))
    matplot(t(test_estimate$transition_probabilities[1:n_show, ]), type = "l",
            main = "Time-Varying Transition Probabilities",
            xlab = "Time", ylab = "Probability",
            col = 1:n_show, lty = 1, lwd = 1.5)
    legend("topright", legend = paste("Trans", 1:n_show), 
           col = 1:n_show, lty = 1, lwd = 1.5)
  }
  
  # Score evolution (if available)
  if (!is.null(test_estimate$scores)) {
    # Show scaled scores for first few transitions
    n_show <- min(4, nrow(test_estimate$scores))
    matplot(t(test_estimate$scaled_scores[1:n_show, ]), type = "l",
            main = "Scaled Score Evolution",
            xlab = "Time", ylab = "Scaled Score",
            col = 1:n_show, lty = 1, lwd = 1)
    legend("topright", legend = paste("Score", 1:n_show), 
           col = 1:n_show, lty = 1, lwd = 1)
    abline(h = 0, col = "gray", lty = 2)
  }
  
  # Parameter bias visualization
  bias_data <- estimates_summary$Bias[!is.na(estimates_summary$Bias)]
  bias_names <- estimates_summary$Parameter[!is.na(estimates_summary$Bias)]
  
  barplot(bias_data, names.arg = bias_names, 
          main = "Parameter Estimation Bias",
          ylab = "Bias (Estimate - True)", 
          col = c(rep("lightblue", K), rep("lightgreen", K), 
                  rep("lightyellow", n_transition), rep("lightpink", n_transition),
                  rep("lightcyan", n_transition)),
          las = 2)
  abline(h = 0, col = "red", lty = 2)
}

dev.off()

log_message(NULL, sprintf("Plots saved to %s", file.path(results_dir, plot_filename)))
log_message(NULL, sprintf("Results saved to %s", file.path(results_dir, results_filename)))

# Print final summary
log_message("FINAL SUMMARY", sprintf(
  "GAS Model Test Completed Successfully!
   
   Successful estimations: %d out of %d (%.1f%%)
   Average estimation bias:
   - Means: %.3f
   - Variances: %.3f  
   - Transition probabilities: %.3f
   - A coefficients: %.3f
   - B coefficients: %.3f
   
   Model fit:
   - Mean Log-Likelihood: %.2f
   - Mean AIC: %.2f
   - Mean BIC: %.2f",
  sum(successful_runs), M, 100 * sum(successful_runs) / M,
  mean(abs(estimates_summary$Bias[1:K]), na.rm = TRUE),
  mean(abs(estimates_summary$Bias[(K+1):(2*K)]), na.rm = TRUE),
  mean(abs(estimates_summary$Bias[(2*K+1):(2*K+n_transition)]), na.rm = TRUE),
  mean(abs(estimates_summary$Bias[(2*K+n_transition+1):(2*K+2*n_transition)]), na.rm = TRUE),
  mean(abs(estimates_summary$Bias[(2*K+2*n_transition+1):(2*K+3*n_transition)]), na.rm = TRUE),
  mean(diagnostics[,1], na.rm = TRUE),
  mean(diagnostics[,2], na.rm = TRUE),
  mean(diagnostics[,3], na.rm = TRUE)
))

log_message("COMPLETE", "Test completed successfully!")