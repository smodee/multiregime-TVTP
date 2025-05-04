#' Simulation Scripts for the Exogenous Model
#' 
#' This file provides simulation scripts and comparison tools for the exogenous 
#' time-varying transition probability model.

# Load required model implementations
source("models/model_exogenous.R")

#' Run a comprehensive simulation study for the exogenous model
#'
#' @param num_repetitions Number of simulation repetitions
#' @param sample_sizes Vector of sample sizes to test
#' @param K Number of regimes
#' @param param_settings List with parameter settings
#' @param exo_process Function to generate exogenous process
#' @param seed Random seed for reproducibility
#' @param output_dir Directory to save detailed results (NULL for no saving)
#' @return Data frame with simulation results
#' @details
#' Conducts a comprehensive simulation study for the exogenous model with different
#' sample sizes and measures estimation accuracy and computational performance.
#'
#' @examples
#' # Define parameter settings
#' param_settings <- list(
#'   mu = c(-2, 1, 2),
#'   sigma2 = c(0.02, 0.2, 0.6),
#'   init_trans = rep(0.2, 6),
#'   A = c(0.1, -0.1, 0.05, -0.05, 0.2, -0.2)
#' )
#' 
#' # Define exogenous process generator
#' exo_gen <- function(N) rnorm(N)
#' 
#' # Run simulation study
#' results <- run_exo_simulation_study(
#'   num_repetitions = 10,
#'   sample_sizes = c(500, 1000),
#'   K = 3,
#'   param_settings = param_settings,
#'   exo_process = exo_gen
#' )
run_exo_simulation_study <- function(num_repetitions = 100, 
                                     sample_sizes = c(250, 500, 1000), 
                                     K = 3,
                                     param_settings = NULL,
                                     exo_process = function(N) rnorm(N),
                                     seed = 123,
                                     output_dir = NULL) {
  # Set random seed for reproducibility
  set.seed(seed)
  
  if (is.null(param_settings)) {
    # Default parameter settings for K regimes
    param_settings <- list(
      mu = seq(-3, 3, length.out = K),
      sigma2 = seq(0.1, 1, length.out = K),
      init_trans = rep(0.2, K*(K-1)),
      A = rep(0.1, K*(K-1))
    )
  }
  
  # Ensure we have the right parameters for K regimes
  if (length(param_settings$mu) != K) {
    stop("Parameter 'mu' must have length K")
  }
  if (length(param_settings$sigma2) != K) {
    stop("Parameter 'sigma2' must have length K")
  }
  if (length(param_settings$init_trans) != K*(K-1)) {
    stop("Parameter 'init_trans' must have length K*(K-1)")
  }
  if (length(param_settings$A) != K*(K-1)) {
    stop("Parameter 'A' must have length K*(K-1)")
  }
  
  # Extract parameters
  mu <- param_settings$mu
  sigma2 <- param_settings$sigma2
  init_trans <- param_settings$init_trans
  A <- param_settings$A
  
  # Prepare results storage
  all_results <- data.frame()
  
  # Set burn-in and cut-off
  B <- 100
  C <- 50
  
  # Loop through sample sizes
  for (N in sample_sizes) {
    cat("\n----- Sample size:", N, "-----\n")
    
    # Loop through repetitions
    for (rep in 1:num_repetitions) {
      cat("Repetition", rep, "of", num_repetitions, "\n")
      
      # Generate exogenous process
      X_Exo <- exo_process(N + B + C)
      
      # Time the data generation
      start_time <- Sys.time()
      
      # Generate a single path
      data_sim <- dataexoCD(1, N + B + C, mu, sigma2, init_trans, A, X_Exo, burn_in = 0)
      
      # Extract the data
      y <- data_sim[1,]
      
      gen_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Time the estimation
      start_time <- Sys.time()
      
      # Estimate the model
      tryCatch({
        estimate <- estimate_exo_model(
          y = y,
          X_Exo = X_Exo,
          K = K,
          B = B,
          C = C,
          verbose = FALSE
        )
        
        est_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        
        # Extract estimated parameters
        mu_est <- estimate$parameters$mu
        sigma2_est <- estimate$parameters$sigma2
        init_trans_est <- estimate$parameters$init_trans
        A_est <- estimate$parameters$A
        
        # Calculate estimation errors
        mu_error <- sqrt(mean((mu_est - mu)^2))
        sigma2_error <- sqrt(mean((sigma2_est - sigma2)^2))
        init_trans_error <- sqrt(mean((init_trans_est - init_trans)^2))
        A_error <- sqrt(mean((A_est - A)^2))
        
        # Store results
        result <- data.frame(
          SampleSize = N,
          Repetition = rep,
          K = K,
          GenerationTime = gen_time,
          EstimationTime = est_time,
          LogLik = estimate$diagnostics$loglik,
          AIC = estimate$diagnostics$aic,
          BIC = estimate$diagnostics$bic,
          MuError = mu_error,
          Sigma2Error = sigma2_error,
          InitTransError = init_trans_error,
          AError = A_error
        )
        
        all_results <- rbind(all_results, result)
        
        cat("Estimation completed. LogLik:", round(estimate$diagnostics$loglik, 2), "\n")
        
        # Save detailed results if output directory is provided
        if (!is.null(output_dir)) {
          dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
          
          # Save the simulation parameters and estimate results
          sim_details <- list(
            parameters = param_settings,
            estimates = estimate,
            y = y,
            X_Exo = X_Exo,
            simulation_info = list(
              N = N,
              rep = rep,
              K = K,
              seed = seed,
              gen_time = gen_time,
              est_time = est_time
            )
          )
          
          save_file <- file.path(output_dir, 
                                 sprintf("exo_sim_N%d_rep%d_K%d.rds", N, rep, K))
          saveRDS(sim_details, save_file)
        }
      }, error = function(e) {
        cat("Error in estimation:", e$message, "\n")
        
        # Store error result
        result <- data.frame(
          SampleSize = N,
          Repetition = rep,
          K = K,
          GenerationTime = gen_time,
          EstimationTime = NA,
          LogLik = NA,
          AIC = NA,
          BIC = NA,
          MuError = NA,
          Sigma2Error = NA,
          InitTransError = NA,
          AError = NA,
          ErrorMessage = e$message
        )
        
        all_results <- rbind(all_results, result)
      })
    }
    
    # Summarize results for current sample size
    summary_stats <- all_results[all_results$SampleSize == N, ]
    
    cat("\nSummary for sample size", N, ":\n")
    cat("Mean estimation time:", mean(summary_stats$EstimationTime, na.rm = TRUE), "seconds\n")
    cat("Mean parameter errors:\n")
    cat("  Mu:", mean(summary_stats$MuError, na.rm = TRUE), "\n")
    cat("  Sigma2:", mean(summary_stats$Sigma2Error, na.rm = TRUE), "\n")
    cat("  InitTrans:", mean(summary_stats$InitTransError, na.rm = TRUE), "\n")
    cat("  A:", mean(summary_stats$AError, na.rm = TRUE), "\n")
  }
  
  return(all_results)
}

#' Test various exogenous process specifications
#'
#' @param y Observed data
#' @param exo_processes List of exogenous processes to test
#' @param K Number of regimes
#' @param B Burn-in period
#' @param C Cut-off period
#' @return Data frame comparing model fit for different exogenous processes
#' @details
#' Tests different exogenous processes as drivers for the transition probabilities
#' and compares their model fit using AIC, BIC and likelihood.
#'
#' @examples
#' # Generate some sample data
#' y <- rnorm(1000)
#' 
#' # Create different exogenous processes
#' exo_processes <- list(
#'   random = rnorm(1000),
#'   trend = 1:1000/1000,
#'   sine = sin(1:1000 * 2 * pi / 100)
#' )
#' 
#' # Compare the processes
#' comparison <- compare_exo_processes(y, exo_processes)
compare_exo_processes <- function(y, exo_processes, K = 3, B = 100, C = 50) {
  # Ensure y is a numeric vector
  if (!is.numeric(y)) {
    stop("Input 'y' must be a numeric vector")
  }
  
  # Check that exo_processes is a list
  if (!is.list(exo_processes)) {
    stop("Input 'exo_processes' must be a list of exogenous processes")
  }
  
  # Check that all exogenous processes have the correct length
  for (name in names(exo_processes)) {
    if (length(exo_processes[[name]]) != length(y)) {
      stop(paste("Exogenous process", name, "must have the same length as y"))
    }
  }
  
  # Prepare results storage
  results <- data.frame()
  
  # Loop through each exogenous process
  for (name in names(exo_processes)) {
    cat("Testing exogenous process:", name, "\n")
    
    X_Exo <- exo_processes[[name]]
    
    # Estimate the model
    tryCatch({
      estimate <- estimate_exo_model(
        y = y,
        X_Exo = X_Exo,
        K = K,
        B = B,
        C = C,
        verbose = FALSE
      )
      
      # Store results
      result <- data.frame(
        ExoProcess = name,
        LogLik = estimate$diagnostics$loglik,
        AIC = estimate$diagnostics$aic,
        BIC = estimate$diagnostics$bic,
        NumParams = estimate$diagnostics$num_params
      )
      
      results <- rbind(results, result)
      
      cat("  LogLik:", round(estimate$diagnostics$loglik, 2), 
          "AIC:", round(estimate$diagnostics$aic, 2), 
          "BIC:", round(estimate$diagnostics$bic, 2), "\n")
      
    }, error = function(e) {
      cat("  Error in estimation:", e$message, "\n")
      
      # Store error result
      result <- data.frame(
        ExoProcess = name,
        LogLik = NA,
        AIC = NA,
        BIC = NA,
        NumParams = NA,
        ErrorMessage = e$message
      )
      
      results <- rbind(results, result)
    })
  }
  
  # Sort results by AIC
  results <- results[order(results$AIC), ]
  
  return(results)
}

#' Create and test various transformations of an exogenous process
#'
#' @param base_process Base exogenous process
#' @param y Observed data
#' @param K Number of regimes
#' @param transformations List of transformation functions
#' @return Data frame comparing model fit for different transformations
#' @details
#' Applies different transformations to a base exogenous process and compares
#' their effectiveness as drivers for transition probabilities using model fit criteria.
#'
#' @examples
#' # Generate data and base process
#' y <- rnorm(1000)
#' base_process <- rnorm(1000)
#' 
#' # Define transformations
#' transformations <- list(
#'   identity = function(x) x,
#'   abs = function(x) abs(x),
#'   squared = function(x) x^2,
#'   lagged = function(x) c(NA, x[-length(x)])
#' )
#' 
#' # Compare transformations
#' result <- test_exo_transformations(base_process, y, transformations = transformations)
test_exo_transformations <- function(base_process, y, K = 3, 
                                     transformations = list(
                                       identity = function(x) x,
                                       abs = function(x) abs(x),
                                       squared = function(x) x^2,
                                       exp = function(x) exp(pmin(x, 5)),
                                       lagged = function(x) c(NA, x[-length(x)])
                                     )) {
  # Ensure base_process and y have the same length
  if (length(base_process) != length(y)) {
    stop("base_process and y must have the same length")
  }
  
  # Apply transformations to create different exogenous processes
  exo_processes <- list()
  
  for (name in names(transformations)) {
    exo_processes[[name]] <- transformations[[name]](base_process)
    
    # Replace NAs with 0 for any transformation that creates NAs
    exo_processes[[name]][is.na(exo_processes[[name]])] <- 0
  }
  
  # Compare the transformed processes
  results <- compare_exo_processes(y, exo_processes, K = K)
  
  return(results)
}

#' Generate a variety of exogenous processes for testing
#'
#' @param N Length of the processes
#' @param seed Random seed for reproducibility
#' @return List of different exogenous processes
#' @details
#' Generates a variety of exogenous processes with different characteristics
#' that can be used for model testing and comparison.
#'
#' @examples
#' # Generate processes of length 1000
#' processes <- generate_test_processes(1000)
#' 
#' # Plot the first few processes
#' par(mfrow=c(3,2))
#' for(i in 1:6) {
#'   name <- names(processes)[i]
#'   plot(processes[[name]], type="l", main=name)
#' }
generate_test_processes <- function(N, seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Initialize list to store processes
  processes <- list()
  
  # 1. Random walk
  processes$random_walk <- cumsum(rnorm(N))
  
  # 2. AR(1) process with high persistence
  ar1 <- numeric(N)
  ar1[1] <- rnorm(1)
  for (i in 2:N) {
    ar1[i] <- 0.9 * ar1[i-1] + rnorm(1, 0, 0.2)
  }
  processes$ar1 <- ar1
  
  # 3. Slow sine wave
  processes$slow_sine <- sin(1:N * 2 * pi / (N/5))
  
  # 4. Fast sine wave
  processes$fast_sine <- sin(1:N * 2 * pi / (N/20))
  
  # 5. Linear trend
  processes$trend <- seq(-1, 1, length.out = N)
  
  # 6. Step function (structural break)
  processes$step <- c(rep(-0.5, N/2), rep(0.5, N/2))
  
  # 7. Volatility clustering (GARCH-like)
  vol <- numeric(N)
  vol[1] <- abs(rnorm(1))
  for (i in 2:N) {
    if (i %% 100 < 50) {
      # High volatility periods
      vol[i] <- 0.9 * vol[i-1] + abs(rnorm(1, 0, 0.5))
    } else {
      # Low volatility periods
      vol[i] <- 0.9 * vol[i-1] + abs(rnorm(1, 0, 0.1))
    }
  }
  processes$volatility <- vol
  
  # 8. Noisy quadratic trend
  processes$quadratic <- (1:N - N/2)^2 / (N/2)^2 + rnorm(N, 0, 0.1)
  
  # 9. Exponential growth
  processes$exponential <- exp(seq(0, 1, length.out = N)) - 1
  
  # 10. Cyclical with increasing amplitude
  t <- 1:N
  processes$increasing_cycles <- (1 + 0.5 * t/N) * sin(2 * pi * t / (N/10))
  
  # 11. White noise
  processes$white_noise <- rnorm(N)
  
  # 12. Heavy-tailed noise
  processes$heavy_tailed <- rt(N, df = 3)
  
  return(processes)
}

#' Analyze how exogenous process affects transition dynamics
#'
#' @param data_gen_params Parameters for data generation
#' @param exo_process Exogenous process to analyze
#' @param K Number of regimes
#' @param N Sample size
#' @param A_values Different A values to test
#' @return List with analysis results
#' @details
#' Analyzes how different strengths of the exogenous process effect (via A parameter)
#' change the dynamics of the regime transitions.
#'
#' @examples
#' # Generate a sine wave exogenous process
#' exo <- sin(1:1000 * 2 * pi / 100)
#' 
#' # Define data generation parameters
#' params <- list(
#'   mu = c(-1, 1, 2),
#'   sigma2 = c(0.1, 0.2, 0.5),
#'   init_trans = rep(0.2, 6)
#' )
#' 
#' # Analyze with different A values
#' analysis <- analyze_exo_impact(params, exo, A_values = c(0, 0.2, 0.5))
analyze_exo_impact <- function(data_gen_params, exo_process, K = 3, N = 1000,
                               A_values = c(0, 0.1, 0.2, 0.5, 1.0)) {
  # Ensure exo_process has the right length
  if (length(exo_process) < N) {
    stop("Exogenous process must have length at least N")
  }
  
  # Extract parameters
  mu <- data_gen_params$mu
  sigma2 <- data_gen_params$sigma2
  init_trans <- data_gen_params$init_trans
  
  # Ensure we have the right parameters for K regimes
  if (length(mu) != K) {
    stop("Parameter 'mu' must have length K")
  }
  if (length(sigma2) != K) {
    stop("Parameter 'sigma2' must have length K")
  }
  if (length(init_trans) != K*(K-1)) {
    stop("Parameter 'init_trans' must have length K*(K-1)")
  }
  
  # Initialize result list
  results <- list(
    parameters = data_gen_params,
    exo_process = exo_process[1:N],
    A_values = A_values,
    N = N,
    K = K,
    simulations = list(),
    transition_dynamics = list(),
    regime_distributions = data.frame()
  )
  
  # Analyze for each A value
  for (a_idx in seq_along(A_values)) {
    a_val <- A_values[a_idx]
    cat("Analyzing with A =", a_val, "\n")
    
    # Set A coefficients to current value
    A <- rep(a_val, K*(K-1))
    
    # Generate data
    sim_data <- dataexoCD(1, N, mu, sigma2, init_trans, A, exo_process[1:N], burn_in = 0)
    y <- sim_data[1,]
    
    # Estimate model
    estimate <- estimate_exo_model(
      y = y,
      X_Exo = exo_process[1:N],
      K = K,
      B = 0,
      C = 0,
      verbose = FALSE
    )
    
    # Store simulation results
    results$simulations[[a_idx]] <- list(
      data = y,
      estimate = estimate,
      A_value = a_val
    )
    
    # Calculate regime probabilities over time
    regime_probs <- estimate$filtered_probabilities
    
    # Calculate transition matrices over time
    trans_matrices <- array(NA, dim = c(K, K, N))
    for (t in 1:N) {
      # Extract raw transition probabilities at time t
      p_trans_t <- estimate$transition_probabilities[,t]
      
      # Convert to transition matrix
      trans_matrices[,,t] <- transition_matrix(p_trans_t)
    }
    
    # Store transition dynamics
    results$transition_dynamics[[a_idx]] <- list(
      transition_matrices = trans_matrices,
      A_value = a_val
    )
    
    # Calculate regime distribution statistics
    avg_regime_probs <- colMeans(regime_probs)
    regime_switching_rate <- sum(diff(apply(regime_probs, 1, which.max)) != 0) / (N-1)
    
    # Add to regime distribution data frame
    results$regime_distributions <- rbind(
      results$regime_distributions,
      data.frame(
        A_value = a_val,
        regime_switching_rate = regime_switching_rate,
        t(avg_regime_probs)
      )
    )
  }
  
  # Set column names for regime distributions
  colnames(results$regime_distributions)[3:(K+2)] <- paste0("avg_regime", 0:(K-1))
  
  return(results)
}

#' Generate a function to create a custom exogenous process
#'
#' @param type Type of process ("ar1", "sine", "step", etc.)
#' @param params Parameters for the process
#' @return A function that generates the specified process
#' @details
#' Creates a function to generate custom exogenous processes with specific properties.
#' This is useful for testing how different types of exogenous processes affect
#' the transition dynamics.
#'
#' @examples
#' # Create an AR(1) process generator
#' ar1_gen <- create_exo_process_generator("ar1", list(phi = 0.9, sigma = 0.2))
#' 
#' # Generate a process of length 1000
#' exo_process <- ar1_gen(1000)
#' 
#' # Plot the generated process
#' plot(exo_process, type = "l")
create_exo_process_generator <- function(type, params = list()) {
  # Define the generator function based on type
  if (type == "ar1") {
    # Default AR(1) parameters
    if (is.null(params$phi)) params$phi <- 0.9
    if (is.null(params$sigma)) params$sigma <- 0.2
    if (is.null(params$mu)) params$mu <- 0
    
    # Create generator function
    generator <- function(N, seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
      
      x <- numeric(N)
      x[1] <- rnorm(1, params$mu, params$sigma)
      
      for (t in 2:N) {
        x[t] <- params$mu + params$phi * (x[t-1] - params$mu) + rnorm(1, 0, params$sigma)
      }
      
      return(x)
    }
  } else if (type == "sine") {
    # Default sine wave parameters
    if (is.null(params$period)) params$period <- 100
    if (is.null(params$amplitude)) params$amplitude <- 1
    if (is.null(params$phase)) params$phase <- 0
    
    # Create generator function
    generator <- function(N, seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
      
      t <- 1:N
      x <- params$amplitude * sin(2 * pi * t / params$period + params$phase)
      
      return(x)
    }
  } else if (type == "step") {
    # Default step function parameters
    if (is.null(params$step_points)) params$step_points <- c(0.5)
    if (is.null(params$levels)) params$levels <- c(-1, 1)
    
    # Create generator function
    generator <- function(N, seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
      
      # Initialize output
      x <- numeric(N)
      
      # Insert breaks
      break_points <- floor(params$step_points * N)
      break_points <- c(0, break_points, N)
      
      # Fill in levels
      for (i in 1:length(params$levels)) {
        start_idx <- break_points[i] + 1
        end_idx <- break_points[i+1]
        x[start_idx:end_idx] <- params$levels[i]
      }
      
      return(x)
    }
  } else if (type == "random_walk") {
    # Default random walk parameters
    if (is.null(params$sigma)) params$sigma <- 0.1
    
    # Create generator function
    generator <- function(N, seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
      
      x <- cumsum(rnorm(N, 0, params$sigma))
      
      return(x)
    }
  } else if (type == "garch") {
    # Default GARCH-like parameters
    if (is.null(params$omega)) params$omega <- 0.1
    if (is.null(params$alpha)) params$alpha <- 0.1
    if (is.null(params$beta)) params$beta <- 0.8
    
    # Create generator function
    generator <- function(N, seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
      
      # Initialize volatility and returns
      sigma2 <- numeric(N)
      x <- numeric(N)
      
      # Set initial volatility
      sigma2[1] <- params$omega / (1 - params$alpha - params$beta)
      x[1] <- rnorm(1, 0, sqrt(sigma2[1]))
      
      # Generate GARCH process
      for (t in 2:N) {
        sigma2[t] <- params$omega + params$alpha * x[t-1]^2 + params$beta * sigma2[t-1]
        x[t] <- rnorm(1, 0, sqrt(sigma2[t]))
      }
      
      return(x)
    }
  } else if (type == "composite") {
    # Composite process - combine multiple processes
    if (is.null(params$components)) {
      stop("For composite process, you must provide 'components' parameter")
    }
    
    # Create generator function
    generator <- function(N, seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
      
      # Initialize output
      x <- numeric(N)
      
      # Add each component
      for (component in params$components) {
        # Get the component generator
        comp_gen <- create_exo_process_generator(component$type, component$params)
        
        # Generate the component
        comp_x <- comp_gen(N)
        
        # Add to total with specified weight
        x <- x + component$weight * comp_x
      }
      
      return(x)
    }
  } else {
    stop("Unknown process type: ", type)
  }
  
  return(generator)
}
