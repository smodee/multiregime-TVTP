#' Main Simulation Script for Markov Regime-Switching Models
#' 
#' This script serves as the main entry point for running simulations and analyses
#' with Markov regime-switching models featuring time-varying transition probabilities.
#'
#' Author: Samuel Modée
#' Date: 2025-04-29

# Clear workspace
rm(list = ls())

# Set working directory to the project root if needed
# setwd("path/to/your/project")

# Source helper functions
source("helpers/utility_functions.R")
source("helpers/transition_helpers.R")
source("helpers/parameter_transforms.R")

# Source model implementations
source("models/model_constant.R")
source("models/model_TVP.R")
source("models/model_exogenous.R")

# Source simulation scripts
source("simulation/sim_TVP.R")
source("simulation/sim_exogenous.R")

# Set random seed for reproducibility
set.seed(42)

#-------------------------------------------------------------------------------
# Configuration Settings
#-------------------------------------------------------------------------------

# Define common parameters
K <- 3  # Number of regimes
N <- 1000  # Sample size
burn_in <- 100  # Burn-in period
cut_off <- 50  # Cut-off period

# Define model parameters
mu <- c(-2, 1, 2)
sigma2 <- c(0.05, 0.2, 0.6)
init_trans <- rep(0.2, K*(K-1))
A <- c(0.1, -0.1, 0.05, -0.05, 0.2, -0.2)

# Create output directory for results
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Configure which simulations to run
run_config <- list(
  "single_path_sim" = TRUE,
  "multiple_model_comparison" = TRUE,
  "different_A_values" = TRUE,
  "exo_with_different_processes" = TRUE,
  "performance_comparison" = FALSE  # Set to TRUE for more extensive comparisons
)

# Set output options
verbose <- TRUE  # Print detailed output
save_results <- TRUE  # Save results to files

#-------------------------------------------------------------------------------
# 1. Basic Model Simulations
#-------------------------------------------------------------------------------

if (run_config$single_path_sim) {
  cat("\n==== Running Single Path Simulations ====\n")
  
  # Generate data from TVP model
  cat("Generating TVP data...\n")
  tvp_data <- dataTVPCD(1, N + burn_in, mu, sigma2, init_trans, A, burn_in = burn_in)
  
  # Generate exogenous process
  X_Exo <- rnorm(N + burn_in)
  
  # Generate data from Exogenous model
  cat("Generating Exogenous model data...\n")
  exo_data <- dataexoCD(1, N + burn_in, mu, sigma2, init_trans, A, X_Exo, burn_in = burn_in)
  
  # Generate data from Constant model
  cat("Generating Constant model data...\n")
  const_data <- dataConstCD(1, N, mu, sigma2, init_trans)
  
  # Plot the generated data
  pdf(file.path(results_dir, "simulated_data.pdf"), width = 10, height = 7)
  par(mfrow = c(3, 1))
  plot(tvp_data[1,], type = "l", main = "TVP Model Data", 
       xlab = "Time", ylab = "Value")
  plot(exo_data[1,], type = "l", main = "Exogenous Model Data", 
       xlab = "Time", ylab = "Value")
  plot(const_data[1,], type = "l", main = "Constant Model Data", 
       xlab = "Time", ylab = "Value")
  dev.off()
  
  cat("Single path simulations completed. Plots saved to", 
      file.path(results_dir, "simulated_data.pdf"), "\n")
}

#-------------------------------------------------------------------------------
# 2. Model Comparison for a Single Dataset
#-------------------------------------------------------------------------------

if (run_config$multiple_model_comparison) {
  cat("\n==== Comparing Multiple Models on a Single Dataset ====\n")
  
  # Generate data with known properties for testing
  X_Exo <- sin(2 * pi * (1:N) / 100)  # Sine wave exogenous process
  true_data <- dataexoCD(1, N + burn_in, mu, sigma2, init_trans, A, 
                         X_Exo, burn_in = burn_in)
  
  # Extract the data series
  y <- true_data[1,]
  
  # Compare different models
  model_comparison <- compare_models(
    y = y, 
    K = K,
    models = c("TVP", "Constant", "AR", "GARCH"),
    verbose = verbose
  )
  
  # Print results
  cat("\nModel comparison results:\n")
  print(model_comparison)
  
  # Save results
  if (save_results) {
    write.csv(model_comparison, 
              file.path(results_dir, "model_comparison_results.csv"), 
              row.names = FALSE)
  }
}

#-------------------------------------------------------------------------------
# 3. TVP Model with Different A Values
#-------------------------------------------------------------------------------

if (run_config$different_A_values) {
  cat("\n==== Analyzing TVP Model with Different A Values ====\n")
  
  # Generate data 
  tvp_data <- dataTVPCD(1, N + burn_in, mu, sigma2, init_trans, A, burn_in = burn_in)
  y <- tvp_data[1,]
  
  # Test different A values
  A_values <- c(0, 0.1, 0.3, 0.5, 0.7)
  
  sensitivity_results <- analyze_tvp_sensitivity(
    y = y,
    K = K,
    A_values = A_values,
    init_trans = init_trans
  )
  
  # Print summary
  cat("\nTVP sensitivity results:\n")
  print(sensitivity_results$regime_dynamics)
  
  # Save results
  if (save_results) {
    saveRDS(sensitivity_results, 
            file.path(results_dir, "tvp_sensitivity_results.rds"))
    
    # Create a plot of key metrics
    pdf(file.path(results_dir, "tvp_sensitivity_plot.pdf"), width = 10, height = 8)
    
    # Set up plotting area with 2x2 panels
    par(mfrow = c(2, 2))
    
    # Plot 1: Log-likelihood
    plot(sensitivity_results$regime_dynamics$A_value, 
         sensitivity_results$regime_dynamics$loglik,
         type = "b", xlab = "A Value", ylab = "Log-Likelihood",
         main = "Model Fit vs. A Value")
    
    # Plot 2: Average regime duration
    plot(sensitivity_results$regime_dynamics$A_value, 
         sensitivity_results$regime_dynamics$avg_duration,
         type = "b", xlab = "A Value", ylab = "Average Duration",
         main = "Regime Persistence vs. A Value")
    
    # Plot 3: Transition frequency
    plot(sensitivity_results$regime_dynamics$A_value, 
         sensitivity_results$regime_dynamics$transition_freq,
         type = "b", xlab = "A Value", ylab = "Transition Frequency",
         main = "Regime Switching Rate vs. A Value")
    
    # Plot 4: Entropy
    plot(sensitivity_results$regime_dynamics$A_value, 
         sensitivity_results$regime_dynamics$avg_entropy,
         type = "b", xlab = "A Value", ylab = "Average Entropy",
         main = "Regime Uncertainty vs. A Value")
    
    dev.off()
    
    cat("TVP sensitivity plots saved to", 
        file.path(results_dir, "tvp_sensitivity_plot.pdf"), "\n")
  }
}

#-------------------------------------------------------------------------------
# 4. Exogenous Model with Different Processes
#-------------------------------------------------------------------------------

if (run_config$exo_with_different_processes) {
  cat("\n==== Testing Exogenous Model with Different Processes ====\n")
  
  # Generate test data
  test_data <- generate_test_datasets(N)
  y <- test_data$regime_switching  # Use regime switching data
  
  # Generate different exogenous processes
  exo_processes <- generate_test_processes(N)
  
  # Compare different processes with fixed data
  exo_comparison <- compare_exo_processes(
    y = y,
    exo_processes = exo_processes,
    K = K
  )
  
  # Print results
  cat("\nExogenous process comparison results:\n")
  print(exo_comparison)
  
  # Save results
  if (save_results) {
    write.csv(exo_comparison, 
              file.path(results_dir, "exo_comparison_results.csv"), 
              row.names = FALSE)
  }
  
  # Test transformations of the best process
  best_process_name <- exo_comparison$ExoProcess[1]
  best_process <- exo_processes[[best_process_name]]
  
  transformation_results <- test_exo_transformations(
    base_process = best_process,
    y = y,
    K = K
  )
  
  # Print results
  cat("\nTransformation comparison results for process", best_process_name, ":\n")
  print(transformation_results)
  
  # Save results
  if (save_results) {
    write.csv(transformation_results, 
              file.path(results_dir, "transformation_results.csv"), 
              row.names = FALSE)
  }
}

#-------------------------------------------------------------------------------
# 5. Performance Comparison (Optional - More Computationally Intensive)
#-------------------------------------------------------------------------------

if (run_config$performance_comparison) {
  cat("\n==== Running Performance Comparison (This may take a while) ====\n")
  
  # Run TVP model simulation study
  tvp_results <- run_tvp_simulation_study(
    num_repetitions = 5,
    sample_sizes = c(500, 1000),
    K = K,
    param_settings = list(
      mu = mu,
      sigma2 = sigma2,
      init_trans = init_trans,
      A = A
    ),
    seed = 123,
    output_dir = file.path(results_dir, "tvp_simulation")
  )
  
  # Create exogenous process generator
  sine_generator <- function(N) {
    return(sin(2 * pi * (1:N) / 100))
  }
  
  # Run exogenous model simulation study
  exo_results <- run_exo_simulation_study(
    num_repetitions = 5,
    sample_sizes = c(500, 1000),
    K = K,
    param_settings = list(
      mu = mu,
      sigma2 = sigma2,
      init_trans = init_trans,
      A = A
    ),
    exo_process = sine_generator,
    seed = 123,
    output_dir = file.path(results_dir, "exo_simulation")
  )
  
  # Save results
  if (save_results) {
    write.csv(tvp_results, 
              file.path(results_dir, "tvp_performance_results.csv"), 
              row.names = FALSE)
    write.csv(exo_results, 
              file.path(results_dir, "exo_performance_results.csv"), 
              row.names = FALSE)
  }
}

cat("\n==== All simulations completed successfully ====\n")
