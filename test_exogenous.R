# test_exogenous_model.R
#
# This script tests the exogenous regime-switching model with time-varying 
# transition probabilities using simulated data.
#
# Author: Samuel Modée
# Date: May 4, 2025

# Clear environment
rm(list = ls())

# Load required libraries and model implementation
source("helpers/utility_functions.R")
source("helpers/transition_helpers.R")
source("helpers/parameter_transforms.R")
source("models/model_exogenous.R")

# Set random seed for reproducibility
set.seed(42)

# Set up parameters
K <- 3  # Number of regimes
M <- 10  # Number of simulation paths
N <- 700  # Length of each simulation path
B <- 100  # Burn-in period
C <- 50  # Cut-off period

# Set true parameter values
mu_true <- c(-2, 0.2, 1.5)  # Different means for each regime
sigma2_true <- c(0.8, 0.2, 0.6)  # Different variances for each regime
#init_trans_true <- rep(1/K, K*(K-1))  # Initial transition probabilities
init_trans_true <- c(0.05,0.3,0.2,0.1,0.6,0.15)
A_true <- c(0.4, -0.1, 0.05, -0.75, 0.8, -0.2)  # Exogenous factor weights

cat("True parameter values:\n")
cat("Means (mu):", mu_true, "\n")
cat("Variances (sigma2):", sigma2_true, "\n")
cat("Initial transition probabilities:", init_trans_true, "\n")
cat("Exogenous factor weights (A):", A_true, "\n\n")

# Generate exogenous process
cat("Generating exogenous process...\n")
#X_Exo <- sin(2 * pi * (1:(N+B+C)) / 100)  # Sine wave with period of 100 time steps
X_Exo=rnorm(N+B+C) # Gaussian white noise

# Generate simulation data
cat("Generating simulation data...\n")
sim_start_time <- Sys.time()
sim_data <- dataexoCD(M, N+B+C, mu_true, sigma2_true, init_trans_true, A_true, X_Exo)
sim_end_time <- Sys.time()
cat("Data generation completed in", 
    format(difftime(sim_end_time, sim_start_time), digits = 4), "\n\n")

# Storage for estimation results
param_estimates <- matrix(0, M, 2*K^2)
colnames(param_estimates) <- c(paste0("mu", 1:K), 
                               paste0("sigma2", 1:K), 
                               paste0("trans", 1:(K*(K-1))), 
                               paste0("A", 1:(K*(K-1))))
diagnostics <- matrix(0, M, 3)
colnames(diagnostics) <- c("loglik", "aic", "bic")

# Estimate parameters for each simulation path
cat("Starting parameter estimation for", M, "simulation paths...\n")
total_start_time <- Sys.time()

for (i in 1:M) {
  path_start_time <- Sys.time()
  cat("Processing path", i, "of", M, "... ")
  
  # Extract the current simulation path
  current_data <- sim_data[i,]
  
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
    estimate <- estimate_exo_model(
      y = current_data,
      X_Exo = X_Exo,
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
    
    path_end_time <- Sys.time()
    path_time <- difftime(path_end_time, path_start_time, units = "secs")
    cat("completed in", round(path_time, 2), "seconds\n")
    
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    param_estimates[i,] <- NA
    diagnostics[i,] <- NA
  })
  
  # Progress update
  if (i > 1) {
    elapsed <- difftime(Sys.time(), total_start_time, units = "secs")
    time_per_path <- as.numeric(elapsed) / i
    remaining_paths <- M - i
    est_remaining <- time_per_path * remaining_paths
    
    cat("Progress:", i, "/", M, "paths -", 
        round(i/M*100), "% complete. Est. remaining time:", 
        round(est_remaining/60, 1), "minutes\n")
  }
}

total_end_time <- Sys.time()
total_time <- difftime(total_end_time, total_start_time, units = "mins")
cat("\nParameter estimation completed in", 
    format(total_time, digits = 4), "minutes\n")

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
cat("\nParameter Estimation Summary:\n")
print(estimates_summary)

cat("\nMean Log-Likelihood:", mean(diagnostics[,1], na.rm = TRUE), "\n")
cat("Mean AIC:", mean(diagnostics[,2], na.rm = TRUE), "\n")
cat("Mean BIC:", mean(diagnostics[,3], na.rm = TRUE), "\n")

# Save results
results_dir <- "results"

# Create a timestamp ID in format YYYYMMDDHHMMSS
timestamp_id <- format(Sys.time(), "%Y%m%d%H%M%S")
results_filename <- paste0("exo_model_test_results_", timestamp_id, ".rds")
dir.create(results_dir, showWarnings = FALSE)
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
plot_filename <- paste0("exo_model_test_plots_", timestamp_id, ".pdf")
pdf(file.path(results_dir, plot_filename), width = 10, height = 8)

# Plot 1: Sample path
par(mfrow = c(3, 1))
plot(sim_data[1,], type = "l", main = "Sample Path from Exogenous Model",
     xlab = "Time", ylab = "Value")

# Plot 2: Exogenous process
plot(X_Exo[1:N], type = "l", main = "Exogenous Process",
     xlab = "Time", ylab = "Value", col = "blue")

# Plot 3: Parameter recovery
par(mfrow = c(2, 2))

# Means
boxplot(param_estimates[,1:K], main = "Estimated Means",
        names = paste0("mu", 1:K), col = "lightblue")
points(1:K, mu_true, pch = 16, col = "red")
legend("topright", legend = "True value", pch = 16, col = "red")

# Variances
boxplot(param_estimates[,(K+1):(2*K)], main = "Estimated Variances",
        names = paste0("sigma2", 1:K), col = "lightgreen")
points(1:K, sigma2_true, pch = 16, col = "red")

# Transition probabilities
n_trans <- K*(K-1)
boxplot(param_estimates[,(2*K+1):(2*K+n_trans)], 
        main = "Estimated Transition Probabilities",
        names = paste0("p", 1:n_trans), col = "lightyellow")
points(1:n_trans, init_trans_true, pch = 16, col = "red")

# A coefficients
boxplot(param_estimates[,(2*K+n_trans+1):(2*K+2*n_trans)], 
        main = "Estimated A Coefficients",
        names = paste0("A", 1:n_trans), col = "lightpink")
points(1:n_trans, A_true, pch = 16, col = "red")

dev.off()

cat("\nPlots saved to", file.path(results_dir, plot_filename), "\n")
cat("Results saved to", file.path(results_dir, results_filename), "\n")
cat("Test completed successfully!\n")