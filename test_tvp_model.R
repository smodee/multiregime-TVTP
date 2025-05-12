# test_tvp_model.R
#
# This script tests the TVP (Time-Varying Probabilities) regime-switching model
# with simulation data.
#
# Author: Samuel Modée
# Date: May 12, 2025

# Clear environment
rm(list = ls())

# Load required libraries and model implementation
source("helpers/utility_functions.R")
source("helpers/transition_helpers.R")
source("helpers/parameter_transforms.R")
source("models/model_TVP.R")

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
set.seed(123)

# Set up parameters
K <- 3  # Number of regimes
M <- 20  # Number of simulation paths
N <- 1500  # Length of each simulation path
B <- 200  # Burn-in period
C <- 50  # Cut-off period

# Set true parameter values
mu_true <- c(-2, 1, 2)  # Different means for each regime
sigma2_true <- c(0.05, 0.2, 0.6)  # Different variances for each regime
#init_trans_true <- rep(0.2, K*(K-1))  
init_trans_true <- c(0.05,0.3,0.2,0.1,0.19,0.15) # Initial transition probabilities
A_true <- c(0.1, -0.1, 0.05, -0.05, 0.2, -0.2)  # Autoregressive factor weights

log_message("CONFIGURATION", sprintf(
  "K=%d regimes, M=%d paths, N=%d observations
   Means: %s
   Variances: %s
   Initial transition probabilities: %s
   A coefficients: %s",
  K, M, N,
  paste(round(mu_true, 2), collapse=", "),
  paste(round(sigma2_true, 2), collapse=", "),
  paste(round(init_trans_true, 2), collapse=", "),
  paste(round(A_true, 2), collapse=", ")
))

# Generate simulation data
log_message("DATA GENERATION", "Generating simulation data...")
sim_start_time <- Sys.time()

# Generate data with burn-in and cut-off periods
sim_data <- dataTVPCD(M, N+B+C, mu_true, sigma2_true, init_trans_true, A_true)

sim_end_time <- Sys.time()
log_message(NULL, sprintf("Data generation completed in %s", 
                          format(difftime(sim_end_time, sim_start_time), digits=4)))

# Storage for estimation results
param_estimates <- matrix(0, M, 2*K^2)
colnames(param_estimates) <- c(paste0("mu", 1:K), 
                               paste0("sigma2", 1:K), 
                               paste0("trans", 1:(K*(K-1))), 
                               paste0("A", 1:(K*(K-1))))
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
  init_trans_guess <- pmax(pmin(init_trans_true + rnorm(K*(K-1), 0, 0.05), 0.95), 0.05)
  A_guess <- rep(0, K*(K-1))  # Start with no effect
  
  par_guess <- c(mu_guess, sigma2_guess, init_trans_guess, A_guess)
  
  # Set parameter bounds
  lower_bounds <- c(rep(-Inf, K),         # No bounds on means
                    rep(-Inf, K),         # Variance >= 0 (log-transformed)
                    rep(-Inf, K*(K-1)),   # Probabilities >= 0 (logit-transformed)
                    rep(-1, K*(K-1)))     # A-coefficients bounded
  
  upper_bounds <- c(rep(Inf, K),          # No bounds on means
                    rep(Inf, K),          # Variance unbounded
                    rep(Inf, K*(K-1)),    # Probabilities <= 1 (logit-transformed)
                    rep(1, K*(K-1)))      # A-coefficients bounded
  
  bounds <- list(lower = lower_bounds, upper = upper_bounds)
  
  # Estimate the model
  tryCatch({
    estimate <- estimate_tvp_model(
      y = current_data,
      K = K,
      B = B,
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
      estimate$parameters$A
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
  True_Value = c(mu_true, sigma2_true, init_trans_true, A_true),
  Mean_Estimate = colMeans(param_estimates, na.rm = TRUE),
  SD_Estimate = apply(param_estimates, 2, sd, na.rm = TRUE),
  Bias = colMeans(param_estimates, na.rm = TRUE) - c(mu_true, sigma2_true, init_trans_true, A_true),
  Rel_Bias_Pct = 100 * (colMeans(param_estimates, na.rm = TRUE) - c(mu_true, sigma2_true, init_trans_true, A_true)) / 
    ifelse(abs(c(mu_true, sigma2_true, init_trans_true, A_true)) > 1e-10, 
           c(mu_true, sigma2_true, init_trans_true, A_true), 
           1e-10)
)

# Print summary statistics
log_message("RESULTS", "Parameter Estimation Summary:")
print(estimates_summary)

log_message(NULL, sprintf("Mean Log-Likelihood: %.2f", mean(diagnostics[,1], na.rm = TRUE)))
log_message(NULL, sprintf("Mean AIC: %.2f", mean(diagnostics[,2], na.rm = TRUE)))
log_message(NULL, sprintf("Mean BIC: %.2f", mean(diagnostics[,3], na.rm = TRUE)))

# Save results
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

# Create a timestamp ID in format YYYY-MM-DD-HHMMSS
timestamp_id <- format(Sys.time(), "%Y-%m-%d-%H%M%S")
results_filename <- paste0("tvp_model_test_results_", timestamp_id, ".rds")

saveRDS(list(
  parameters = list(
    mu = mu_true,
    sigma2 = sigma2_true,
    init_trans = init_trans_true,
    A = A_true
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
    B = B,
    C = C
  )
), file = file.path(results_dir, results_filename))

# Create plots
plot_filename <- paste0("tvp_model_test_plots_", timestamp_id, ".pdf")
pdf(file.path(results_dir, plot_filename), width = 10, height = 8)

# Plot 1: Sample path
par(mfrow = c(2, 1))
plot(sim_data[1,], type = "l", main = "Sample Path from TVP Model",
     xlab = "Time", ylab = "Value")

# Add a visualized regime classification if available
if (!is.null(estimate$filtered_probabilities)) {
  # Get the most likely regime at each point
  regime_probs <- estimate$filtered_probabilities
  regime_classifications <- apply(regime_probs, 1, which.max)
  
  # Plot regime classification
  plot(regime_classifications, type = "l", main = "Estimated Regime Classification",
       xlab = "Time", ylab = "Regime", ylim = c(1, K), yaxt = "n")
  axis(2, at = 1:K, labels = 1:K)
}

# Plot 2: Parameter recovery
par(mfrow = c(2, 2))

# Means
boxplot(param_estimates[,1:K], main = "Estimated Means",
        names = paste0("mu", 1:K), col = "lightblue",
        ylim = range(c(param_estimates[,1:K], mu_true), na.rm = TRUE))
points(1:K, mu_true, pch = 16, col = "red")
legend("topright", legend = "True value", pch = 16, col = "red")

# Variances
boxplot(param_estimates[,(K+1):(2*K)], main = "Estimated Variances",
        names = paste0("sigma2", 1:K), col = "lightgreen",
        ylim = range(c(param_estimates[,(K+1):(2*K)], sigma2_true), na.rm = TRUE))
points(1:K, sigma2_true, pch = 16, col = "red")

# Transition probabilities
n_trans <- K*(K-1)
boxplot(param_estimates[,(2*K+1):(2*K+n_trans)], 
        main = "Estimated Transition Probabilities",
        names = paste0("p", 1:n_trans), col = "lightyellow",
        ylim = range(c(param_estimates[,(2*K+1):(2*K+n_trans)], init_trans_true), na.rm = TRUE))
points(1:n_trans, init_trans_true, pch = 16, col = "red")

# A coefficients
boxplot(param_estimates[,(2*K+n_trans+1):(2*K+2*n_trans)], 
        main = "Estimated A Coefficients",
        names = paste0("A", 1:n_trans), col = "lightpink",
        ylim = range(c(param_estimates[,(2*K+n_trans+1):(2*K+2*n_trans)], A_true), na.rm = TRUE))
points(1:n_trans, A_true, pch = 16, col = "red")

dev.off()

log_message(NULL, sprintf("Plots saved to %s", file.path(results_dir, plot_filename)))
log_message(NULL, sprintf("Results saved to %s", file.path(results_dir, results_filename)))
log_message("COMPLETE", "Test completed successfully!")