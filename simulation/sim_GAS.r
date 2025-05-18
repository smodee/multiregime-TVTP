#' Simulation Scripts for the GAS Model
#' 
#' This file provides simulation scripts and comparison tools for the score-driven
#' time-varying transition probability model (GAS model).

# Load required model implementations
source("models/model_GAS.R")

#' Run a comprehensive simulation study for the GAS model
#'
#' @param num_repetitions Number of simulation repetitions
#' @param sample_sizes Vector of sample sizes to test
#' @param K Number of regimes
#' @param param_settings List with parameter settings
#' @param seed Random seed for reproducibility
#' @param output_dir Directory to save detailed results (NULL for no saving)
#' @return Data frame with simulation results
#' @details
#' Conducts a comprehensive simulation study for the GAS model with different
#' sample sizes and measures estimation accuracy and computational performance.
#'
#' @examples
#' # Define parameter settings
#' param_settings <- list(
#'   mu = c(-2, 1, 2),
#'   sigma2 = c(0.02, 0.2, 0.6),
#'   init_trans = rep(0.2, 6),
#'   A = rep(0.1, 6),
#'   B = rep(0.9, 6)
#' )
#' 
#' # Run simulation study
#' results <- run_gas_simulation_study(
#'   num_repetitions = 10,
#'   sample_sizes = c(500, 1000),
#'   K = 3,
#'   param_settings = param_settings
#' )
run_gas_simulation_study <- function(num_repetitions = 100, 
                                     sample_sizes = c(250, 500, 1000), 
                                     K = 3,
                                     param_settings = NULL,
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
      A = rep(0.1, K*(K-1)),
      B = rep(0.9, K*(K-1))
    )
  }
  
  # Ensure we have the right parameters for K regimes
  if (length(param_settings$mu) != K) {
    stop("Parameter 'mu' must have length K")
  }
  if (length(param_settings$sigma2) != K) {
    stop("Parameter 'sigma2' must have length K")
  }
  n_transition <- K*(K-1)
  if (length(param_settings$init_trans) != n_transition) {
    stop("Parameter 'init_trans' must have length K*(K-1)")
  }
  if (length(param_settings$A) != n_transition) {
    stop("Parameter 'A' must have length K*(K-1)")
  }
  if (length(param_settings$B) != n_transition) {
    stop("Parameter 'B' must have length K*(K-1)")
  }
  
  # Extract parameters
  mu <- param_settings$mu
  sigma2 <- param_settings$sigma2
  init_trans <- param_settings$init_trans
  A <- param_settings$A
  B <- param_settings$B
  
  # Prepare results storage
  all_results <- data.frame()
  
  # Set burn-in and cut-off
  B_burnin <- 100
  C <- 50
  
  # Loop through sample sizes
  for (N in sample_sizes) {
    cat("\n----- Sample size:", N, "-----\n")
    
    # Loop through repetitions
    for (rep in 1:num_repetitions) {
      cat("Repetition", rep, "of", num_repetitions, "\n")
      
      # Time the data generation
      start_time <- Sys.time()
      
      # Generate a single path
      data_sim <- dataGASCD(1, N + B_burnin + C, mu, sigma2, init_trans, A, B, burn_in = 0)
      
      # Extract the data
      y <- data_sim[1,]
      
      gen_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Time the estimation
      start_time <- Sys.time()
      
      # Estimate the model
      tryCatch({
        estimate <- estimate_gas_model(
          y = y,
          K = K,
          B_burnin = B_burnin,
          C = C,
          verbose = FALSE
        )
        
        est_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        
        # Extract estimated parameters
        mu_est <- estimate$parameters$mu
        sigma2_est <- estimate$parameters$sigma2
        init_trans_est <- estimate$parameters$init_trans
        A_est <- estimate$parameters$A
        B_est <- estimate$parameters$B
        
        # Calculate estimation errors
        mu_error <- sqrt(mean((mu_est - mu)^2))
        sigma2_error <- sqrt(mean((sigma2_est - sigma2)^2))
        init_trans_error <- sqrt(mean((init_trans_est - init_trans)^2))
        A_error <- sqrt(mean((A_est - A)^2))
        B_error <- sqrt(mean((B_est - B)^2))
        
        # Calculate regime persistence metrics
        persistence <- calculate_persistence(estimate$filtered_probabilities)
        
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
          AError = A_error,
          BError = B_error,
          TransitionRate = persistence$transition_rate,
          AvgDuration = mean(unlist(persistence$avg_durations))
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
                                 sprintf("gas_sim_N%d_rep%d_K%d.rds", N, rep, K))
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
          BError = NA,
          TransitionRate = NA,
          AvgDuration = NA,
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
    cat("  B:", mean(summary_stats$BError, na.rm = TRUE), "\n")
  }
  
  return(all_results)
}

#' Create plots to analyze GAS model simulation results
#' 
#' @param results Results data frame from run_gas_simulation_study
#' @param output_file Path to save the output PDF (NULL for no saving)
#' @param width Width of the output PDF
#' @param height Height of the output PDF
#' @return NULL (creates plots as a side effect)
#' @details
#' Creates various diagnostic plots to visualize the performance of the GAS model
#' across different sample sizes and parameter values.
plot_gas_simulation_results <- function(results, output_file = NULL, width = 10, height = 8) {
  # Check if results is empty
  if (nrow(results) == 0) {
    stop("Results data frame is empty.")
  }
  
  # Start plotting device if output_file is provided
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
    on.exit(dev.off())
  }
  
  # Set up plotting area with 2x2 panels
  par(mfrow = c(2, 2))
  
  # Get unique sample sizes for grouping
  sample_sizes <- sort(unique(results$SampleSize))
  colors <- rainbow(length(sample_sizes))
  
  # Plot 1: Estimation error by sample size
  boxplot(MuError ~ SampleSize, data = results, 
          main = "Mean Parameter Estimation Error",
          xlab = "Sample Size", ylab = "RMSE",
          col = "lightblue")
  
  # Plot 2: Estimation time by sample size
  boxplot(EstimationTime ~ SampleSize, data = results, 
          main = "Estimation Time",
          xlab = "Sample Size", ylab = "Time (seconds)",
          col = "lightgreen")
  
  # Plot 3: AIC by sample size
  boxplot(AIC ~ SampleSize, data = results, 
          main = "Model Fit (AIC)",
          xlab = "Sample Size", ylab = "AIC",
          col = "lightyellow")
  
  # Plot 4: Transition rate by sample size
  boxplot(TransitionRate ~ SampleSize, data = results, 
          main = "Detected Regime Switching Rate",
          xlab = "Sample Size", ylab = "Transition Rate",
          col = "lightpink")
  
  # Reset plotting area for the next set of plots
  par(mfrow = c(2, 2))
  
  # Plot 5: Error in A parameters
  boxplot(AError ~ SampleSize, data = results, 
          main = "A Parameter Estimation Error",
          xlab = "Sample Size", ylab = "RMSE",
          col = "lightblue")
  
  # Plot 6: Error in B parameters
  boxplot(BError ~ SampleSize, data = results, 
          main = "B Parameter Estimation Error",
          xlab = "Sample Size", ylab = "RMSE",
          col = "lightgreen")
  
  # Plot 7: Average regime duration
  boxplot(AvgDuration ~ SampleSize, data = results, 
          main = "Average Regime Duration",
          xlab = "Sample Size", ylab = "Duration",
          col = "lightyellow")
  
  # Plot 8: LogLik by sample size
  boxplot(LogLik ~ SampleSize, data = results, 
          main = "Log-Likelihood",
          xlab = "Sample Size", ylab = "LogLik",
          col = "lightpink")
  
  # If not saving to file, reset the plot layout
  if (is.null(output_file)) {
    par(mfrow = c(1, 1))
  }
}

#' Analyze sensitivity of GAS model to different A and B parameters
#'
#' @param y Time series to analyze
#' @param K Number of regimes
#' @param A_values Vector of A values to test
#' @param B_values Vector of B values to test
#' @param init_trans Initial transition probabilities
#' @param verbose Whether to print progress
#' @return List with analysis results
#' @details
#' Explores how different values of the sensitivity (A) and persistence (B)
#' parameters affect the model fit and regime dynamics.
#'
#' @examples
#' # Generate data and analyze with different parameter values
#' set.seed(123)
#' y <- rnorm(1000)
#' results <- analyze_gas_sensitivity(
#'   y, K = 2, 
#'   A_values = c(0, 0.05, 0.1, 0.2, 0.5),
#'   B_values = c(0.7, 0.8, 0.9, 0.95, 0.99)
#' )
analyze_gas_sensitivity <- function(y, K = 3, 
                                    A_values = seq(0, 0.5, by = 0.1),
                                    B_values = c(0.7, 0.8, 0.9, 0.95, 0.99),
                                    init_trans = NULL,
                                    verbose = TRUE) {
  # Initialize output
  n_transition <- K*(K-1)
  
  # Set default initial transition probabilities if not provided
  if (is.null(init_trans)) {
    init_trans <- rep(0.2, n_transition)
  }
  
  # Check that init_trans has the right length
  if (length(init_trans) != n_transition) {
    stop(sprintf("init_trans must have length K*(K-1) = %d", n_transition))
  }
  
  # Initialize results storage
  results <- list()
  
  # Grid of parameter combinations
  param_grid <- expand.grid(A_value = A_values, B_value = B_values)
  n_combinations <- nrow(param_grid)
  
  # Initialize result dataframes
  ab_results <- data.frame(
    A_value = numeric(n_combinations),
    B_value = numeric(n_combinations),
    loglik = numeric(n_combinations),
    aic = numeric(n_combinations),
    bic = numeric(n_combinations),
    transition_freq = numeric(n_combinations),
    avg_duration = numeric(n_combinations),
    avg_entropy = numeric(n_combinations)
  )
  
  # Store all model estimates for detailed analysis
  all_models <- list()
  
  # Create a progress bar
  if (verbose) {
    cat("Analyzing GAS model sensitivity with", n_combinations, "parameter combinations\n")
    pb <- txtProgressBar(min = 0, max = n_combinations, style = 3)
  }
  
  # Loop through parameter combinations
  for (i in 1:n_combinations) {
    A_val <- param_grid$A_value[i]
    B_val <- param_grid$B_value[i]
    
    # Use the same A and B values for all transitions
    A <- rep(A_val, n_transition)
    B <- rep(B_val, n_transition)
    
    # Create initial parameters with these A and B values
    initial_params <- c(
      # Mean and variance will be estimated from data
      rep(NA, 2*K),
      # Use provided initial transition probabilities
      init_trans,
      # Set A and B to specified values
      A, B
    )
    
    tryCatch({
      # Estimate the model with fixed A and B values
      model <- estimate_gas_model(
        y = y,
        K = K,
        B_burnin = 100,
        C = 50,
        initial_params = initial_params,
        verbose = FALSE
      )
      
      # Extract filtered probabilities
      probs <- model$filtered_probabilities
      
      # Calculate persistence metrics
      persistence <- calculate_persistence(probs)
      
      # Calculate entropy (measure of uncertainty about regime)
      entropy <- apply(probs, 1, function(p) {
        # Avoid log(0) issues
        p_clean <- pmax(p, 1e-10)
        -sum(p_clean * log(p_clean))
      })
      
      # Store results
      ab_results$A_value[i] <- A_val
      ab_results$B_value[i] <- B_val
      ab_results$loglik[i] <- model$diagnostics$loglik
      ab_results$aic[i] <- model$diagnostics$aic
      ab_results$bic[i] <- model$diagnostics$bic
      ab_results$transition_freq[i] <- persistence$transition_rate
      ab_results$avg_duration[i] <- mean(unlist(persistence$avg_durations))
      ab_results$avg_entropy[i] <- mean(entropy)
      
      # Store full model for further analysis
      all_models[[i]] <- model
      
    }, error = function(e) {
      # On error, set results to NA
      ab_results$A_value[i] <- A_val
      ab_results$B_value[i] <- B_val
      ab_results$loglik[i] <- NA
      ab_results$aic[i] <- NA
      ab_results$bic[i] <- NA
      ab_results$transition_freq[i] <- NA
      ab_results$avg_duration[i] <- NA
      ab_results$avg_entropy[i] <- NA
      
      if (verbose) {
        cat("\nError for A =", A_val, "B =", B_val, ":", e$message, "\n")
      }
    })
    
    # Update progress bar
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }
  
  if (verbose) {
    close(pb)
    cat("\nAnalysis completed.\n")
  }
  
  # Prepare return value
  results <- list(
    ab_grid = ab_results,
    models = all_models,
    settings = list(
      A_values = A_values,
      B_values = B_values,
      K = K
    )
  )
  
  return(results)
}

#' Plot sensitivity analysis results
#'
#' @param sensitivity_results Results from analyze_gas_sensitivity
#' @param output_file Path to save the output PDF (NULL for no saving)
#' @param width Width of the output PDF
#' @param height Height of the output PDF
#' @return NULL (creates plots as a side effect)
plot_gas_sensitivity <- function(sensitivity_results, output_file = NULL, width = 10, height = 10) {
  # Extract results
  results <- sensitivity_results$ab_grid
  A_values <- sensitivity_results$settings$A_values
  B_values <- sensitivity_results$settings$B_values
  
  # Check if results is empty
  if (nrow(results) == 0) {
    stop("Results data frame is empty.")
  }
  
  # Start plotting device if output_file is provided
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
    on.exit(dev.off())
  }
  
  # Create a 2x2 grid for plots
  par(mfrow = c(2, 2))
  
  # Create a matrix of unique A and B values for the heatmap
  A_unique <- sort(unique(results$A_value))
  B_unique <- sort(unique(results$B_value))
  n_A <- length(A_unique)
  n_B <- length(B_unique)
  
  # Function to create a heatmap matrix
  create_heatmap_matrix <- function(data_vec) {
    mat <- matrix(NA, nrow = n_B, ncol = n_A)
    for (i in 1:length(data_vec)) {
      a_idx <- match(results$A_value[i], A_unique)
      b_idx <- match(results$B_value[i], B_unique)
      mat[b_idx, a_idx] <- data_vec[i]
    }
    return(mat)
  }
  
  # Create heatmaps for different metrics
  loglik_mat <- create_heatmap_matrix(results$loglik)
  transition_mat <- create_heatmap_matrix(results$transition_freq)
  duration_mat <- create_heatmap_matrix(results$avg_duration)
  entropy_mat <- create_heatmap_matrix(results$avg_entropy)
  
  # Plot 1: Log-likelihood heatmap
  image(A_unique, B_unique, loglik_mat, 
        xlab = "A (Sensitivity)", ylab = "B (Persistence)",
        main = "Log-Likelihood",
        col = heat.colors(20, rev = TRUE))
  contour(A_unique, B_unique, loglik_mat, add = TRUE)
  
  # Plot 2: Transition frequency heatmap
  image(A_unique, B_unique, transition_mat, 
        xlab = "A (Sensitivity)", ylab = "B (Persistence)",
        main = "Regime Transition Frequency",
        col = terrain.colors(20))
  contour(A_unique, B_unique, transition_mat, add = TRUE)
  
  # Plot 3: Average duration heatmap
  image(A_unique, B_unique, duration_mat, 
        xlab = "A (Sensitivity)", ylab = "B (Persistence)",
        main = "Average Regime Duration",
        col = cm.colors(20))
  contour(A_unique, B_unique, duration_mat, add = TRUE)
  
  # Plot 4: Entropy heatmap
  image(A_unique, B_unique, entropy_mat, 
        xlab = "A (Sensitivity)", ylab = "B (Persistence)",
        main = "Average Regime Entropy",
        col = topo.colors(20))
  contour(A_unique, B_unique, entropy_mat, add = TRUE)
  
  # Reset plotting parameters if not saving to file
  if (is.null(output_file)) {
    par(mfrow = c(1, 1))
  }
}

#' Compare model performance for different types of time series
#'
#' @param sample_size Length of each generated time series
#' @param time_series_types List of time series generation functions
#' @param K Number of regimes
#' @param models Vector of model types to compare
#' @param seed Random seed for reproducibility
#' @return Data frame with model comparison results
#' @details
#' Generates different types of time series and compares the performance
#' of various models (GAS, TVP, Exogenous, Constant) on each type.
#'
#' @examples
#' # Define time series generators
#' generators <- list(
#'   stationary = function(n) rnorm(n),
#'   trend = function(n) 0.01 * (1:n) + rnorm(n),
#'   structural_break = function(n) c(rnorm(n/2), rnorm(n/2, mean = 2)),
#'   regime_switching = function(n) {
#'     states <- sample(1:2, n, replace = TRUE, prob = c(0.7, 0.3))
#'     rnorm(n, mean = ifelse(states == 1, -1, 1))
#'   }
#' )
#' 
#' # Compare models
#' comparison <- compare_models_on_different_series(
#'   sample_size = 500, 
#'   time_series_types = generators,
#'   models = c("GAS", "TVP", "Constant")
#' )
compare_models_on_different_series <- function(sample_size = 1000,
                                              time_series_types = NULL,
                                              K = 3,
                                              models = c("GAS", "TVP", "Exogenous", "Constant"),
                                              seed = 123) {
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Define default time series generators if none provided
  if (is.null(time_series_types)) {
    time_series_types <- list(
      stationary = function(n) rnorm(n),
      trend = function(n) 0.01 * (1:n) + rnorm(n),
      structural_break = function(n) c(rnorm(n/2), rnorm(n/2, mean = 2)),
      regime_switching = function(n) {
        states <- sample(1:3, n, replace = TRUE, prob = c(0.6, 0.3, 0.1))
        rnorm(n, mean = c(-1, 0, 1)[states], sd = c(0.5, 1, 1.5)[states])
      },
      volatility_clustering = function(n) {
        vol <- rep(0.1, n)
        for (i in 2:n) {
          vol[i] <- 0.1 + 0.8 * vol[i-1] + 0.1 * abs(rnorm(1))
        }
        rnorm(n, sd = vol)
      }
    )
  }
  
  # Initialize results dataframe
  results <- data.frame()
  
  # Generate an exogenous process for the Exogenous model
  if ("Exogenous" %in% models) {
    X_Exo <- rnorm(sample_size)
  } else {
    X_Exo <- NULL
  }
  
  # Loop through each time series type
  for (series_name in names(time_series_types)) {
    cat("Generating", series_name, "time series...\n")
    
    # Generate the time series
    generator <- time_series_types[[series_name]]
    y <- generator(sample_size)
    
    # Compare models on this time series
    cat("Comparing models on", series_name, "time series...\n")
    comparison <- compare_tvtp_models(
      y = y,
      X_Exo = X_Exo,
      K = K,
      models = models,
      B_burnin = 100,
      C = 50,
      verbose = FALSE
    )
    
    # Add time series type to the results
    comparison$TimeSeriesType <- series_name
    
    # Append to overall results
    results <- rbind(results, comparison)
  }
  
  # Create a summary table
  summary_table <- reshape(results,
                          idvar = c("TimeSeriesType", "Model"),
                          timevar = NULL,
                          direction = "wide",
                          drop = c("NumParams", "EstimationTime"))
  
  # Compute rankings for each time series type
  rankings <- data.frame()
  for (series_name in unique(results$TimeSeriesType)) {
    series_results <- results[results$TimeSeriesType == series_name, ]
    
    # Rank by AIC (lower is better)
    aic_ranks <- rank(series_results$AIC)
    
    # Create ranking dataframe
    series_ranks <- data.frame(
      TimeSeriesType = series_name,
      Model = series_results$Model,
      AIC_Rank = aic_ranks
    )
    
    # Append to overall rankings
    rankings <- rbind(rankings, series_ranks)
  }
  
  # Return both detailed results and rankings
  return(list(
    detailed_results = results,
    summary = summary_table,
    rankings = rankings
  ))
}

#' Generate and analyze a variety of test datasets
#'
#' @param N Length of each time series
#' @param seed Random seed for reproducibility
#' @return List of generated time series
#' @details
#' Creates a set of time series with various characteristics for testing models.
#' These include stationary series, series with trends, structural breaks,
#' regime switching, and volatility clustering.
generate_test_datasets <- function(N = 1000, seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Generate stationary time series
  stationary <- rnorm(N)
  
  # Generate time series with linear trend
  trend <- 0.01 * (1:N) + rnorm(N)
  
  # Generate time series with structural break
  structural_break <- c(rnorm(N/2), rnorm(N/2, mean = 2))
  
  # Generate regime-switching time series
  regime_switching <- numeric(N)
  states <- sample(1:3, N, replace = TRUE, prob = c(0.6, 0.3, 0.1))
  for (i in 1:N) {
    if (states[i] == 1) {
      regime_switching[i] <- rnorm(1, mean = -1, sd = 0.5)
    } else if (states[i] == 2) {
      regime_switching[i] <- rnorm(1, mean = 0, sd = 1)
    } else {
      regime_switching[i] <- rnorm(1, mean = 1, sd = 1.5)
    }
  }
  
  # Generate time series with volatility clustering
  volatility_clustering <- numeric(N)
  vol <- rep(0.1, N)
  for (i in 2:N) {
    vol[i] <- 0.1 + 0.8 * vol[i-1] + 0.1 * abs(rnorm(1))
  }
  for (i in 1:N) {
    volatility_clustering[i] <- rnorm(1, sd = vol[i])
  }
  
  # Generate cyclical time series
  cyclical <- sin(2 * pi * (1:N) / 100) + rnorm(N, sd = 0.5)
  
  # Generate time series with changing persistence
  changing_persistence <- numeric(N)
  for (i in 1:N) {
    if (i <= N/3) {
      # Low persistence
      changing_persistence[i] <- 0.3 * rnorm(1)
      if (i > 1) changing_persistence[i] <- changing_persistence[i] + 0.3 * changing_persistence[i-1]
    } else if (i <= 2*N/3) {
      # Medium persistence
      changing_persistence[i] <- 0.2 * rnorm(1)
      if (i > 1) changing_persistence[i] <- changing_persistence[i] + 0.6 * changing_persistence[i-1]
    } else {
      # High persistence
      changing_persistence[i] <- 0.1 * rnorm(1)
      if (i > 1) changing_persistence[i] <- changing_persistence[i] + 0.9 * changing_persistence[i-1]
    }
  }
  
  # Return all generated time series
  return(list(
    stationary = stationary,
    trend = trend,
    structural_break = structural_break,
    regime_switching = regime_switching,
    volatility_clustering = volatility_clustering,
    cyclical = cyclical,
    changing_persistence = changing_persistence
  ))
}

#' Compare parameter recovery ability of GAS vs TVP models
#'
#' @param num_repetitions Number of simulation runs
#' @param sample_size Length of each simulation
#' @param K Number of regimes
#' @param setup Parameter settings for data generation
#' @param seed Random seed
#' @return List with comparison results
#' @details
#' Compares how well the GAS and TVP models can recover the true parameter
#' values when data is generated from either model.
compare_gas_tvp_parameter_recovery <- function(num_repetitions = 20, 
                                              sample_size = 1000,
                                              K = 3,
                                              setup = NULL,
                                              seed = 123) {
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Configure default parameter setup if not provided
  if (is.null(setup)) {
    n_transition <- K*(K-1)
    setup <- list(
      # Define parameters for both models
      mu = seq(-2, 2, length.out = K),
      sigma2 = seq(0.1, 1, length.out = K),
      init_trans = rep(0.2, n_transition),
      # Parameters specific to each model
      GAS = list(
        A = rep(0.1, n_transition),
        B = rep(0.9, n_transition)
      ),
      TVP = list(
        A = rep(0.1, n_transition)
      )
    )
  }
  
  # Extract common parameters
  mu <- setup$mu
  sigma2 <- setup$sigma2
  init_trans <- setup$init_trans
  
  # Parameters specific to each model
  GAS_A <- setup$GAS$A
  GAS_B <- setup$GAS$B
  TVP_A <- setup$TVP$A
  
  # Initialize results storage
  results <- list(
    # Results for GAS data estimated with both models
    GAS_data = list(
      GAS_model = data.frame(),
      TVP_model = data.frame()
    ),
    # Results for TVP data estimated with both models
    TVP_data = list(
      GAS_model = data.frame(),
      TVP_model = data.frame()
    )
  )
  
  # Loop through repetitions
  for (rep in 1:num_repetitions) {
    cat("\nRunning repetition", rep, "of", num_repetitions, "\n")
    
    # Generate data from GAS model
    cat("  Generating data from GAS model...\n")
    GAS_data <- dataGASCD(1, sample_size, mu, sigma2, init_trans, GAS_A, GAS_B)[1,]
    
    # Generate data from TVP model
    cat("  Generating data from TVP model...\n")
    source("models/model_TVP.R")  # Ensure TVP model is loaded
    TVP_data <- dataTVPCD(1, sample_size, mu, sigma2, init_trans, TVP_A)[1,]
    
    # Estimate GAS model on GAS data
    cat("  Estimating GAS model on GAS data...\n")
    GAS_on_GAS <- tryCatch({
      model <- estimate_gas_model(GAS_data, K = K, verbose = FALSE)
      
      # Extract parameter estimates
      data.frame(
        rep = rep,
        datatype = "GAS",
        model = "GAS",
        loglik = model$diagnostics$loglik,
        aic = model$diagnostics$aic,
        mu_error = sqrt(mean((model$parameters$mu - mu)^2)),
        sigma2_error = sqrt(mean((model$parameters$sigma2 - sigma2)^2)),
        init_trans_error = sqrt(mean((model$parameters$init_trans - init_trans)^2)),
        A_error = sqrt(mean((model$parameters$A - GAS_A)^2)),
        B_error = sqrt(mean((model$parameters$B - GAS_B)^2))
      )
    }, error = function(e) {
      cat("    Error:", e$message, "\n")
      data.frame(
        rep = rep,
        datatype = "GAS",
        model = "GAS",
        loglik = NA,
        aic = NA,
        mu_error = NA,
        sigma2_error = NA,
        init_trans_error = NA,
        A_error = NA,
        B_error = NA
      )
    })
    
    # Estimate TVP model on GAS data
    cat("  Estimating TVP model on GAS data...\n")
    TVP_on_GAS <- tryCatch({
      model <- estimate_tvp_model(GAS_data, K = K, verbose = FALSE)
      
      # Extract parameter estimates
      data.frame(
        rep = rep,
        datatype = "GAS",
        model = "TVP",
        loglik = model$diagnostics$loglik,
        aic = model$diagnostics$aic,
        mu_error = sqrt(mean((model$parameters$mu - mu)^2)),
        sigma2_error = sqrt(mean((model$parameters$sigma2 - sigma2)^2)),
        init_trans_error = sqrt(mean((model$parameters$init_trans - init_trans)^2)),
        A_error = sqrt(mean((model$parameters$A - TVP_A)^2)),
        B_error = NA  # TVP model doesn't have B parameter
      )
    }, error = function(e) {
      cat("    Error:", e$message, "\n")
      data.frame(
        rep = rep,
        datatype = "GAS",
        model = "TVP",
        loglik = NA,
        aic = NA,
        mu_error = NA,
        sigma2_error = NA,
        init_trans_error = NA,
        A_error = NA,
        B_error = NA
      )
    })
    
    # Estimate GAS model on TVP data
    cat("  Estimating GAS model on TVP data...\n")
    GAS_on_TVP <- tryCatch({
      model <- estimate_gas_model(TVP_data, K = K, verbose = FALSE)
      
      # Extract parameter estimates
      data.frame(
        rep = rep,
        datatype = "TVP",
        model = "GAS",
        loglik = model$diagnostics$loglik,
        aic = model$diagnostics$aic,
        mu_error = sqrt(mean((model$parameters$mu - mu)^2)),
        sigma2_error = sqrt(mean((model$parameters$sigma2 - sigma2)^2)),
        init_trans_error = sqrt(mean((model$parameters$init_trans - init_trans)^2)),
        A_error = NA,  # No true A to compare with (different interpretation)
        B_error = NA   # No true B to compare with
      )
    }, error = function(e) {
      cat("    Error:", e$message, "\n")
      data.frame(
        rep = rep,
        datatype = "TVP",
        model = "GAS",
        loglik = NA,
        aic = NA,
        mu_error = NA,
        sigma2_error = NA,
        init_trans_error = NA,
        A_error = NA,
        B_error = NA
      )
    })
    
    # Estimate TVP model on TVP data
    cat("  Estimating TVP model on TVP data...\n")
    TVP_on_TVP <- tryCatch({
      model <- estimate_tvp_model(TVP_data, K = K, verbose = FALSE)
      
      # Extract parameter estimates
      data.frame(
        rep = rep,
        datatype = "TVP",
        model = "TVP",
        loglik = model$diagnostics$loglik,
        aic = model$diagnostics$aic,
        mu_error = sqrt(mean((model$parameters$mu - mu)^2)),
        sigma2_error = sqrt(mean((model$parameters$sigma2 - sigma2)^2)),
        init_trans_error = sqrt(mean((model$parameters$init_trans - init_trans)^2)),
        A_error = sqrt(mean((model$parameters$A - TVP_A)^2)),
        B_error = NA  # TVP model doesn't have B parameter
      )
    }, error = function(e) {
      cat("    Error:", e$message, "\n")
      data.frame(
        rep = rep,
        datatype = "TVP",
        model = "TVP",
        loglik = NA,
        aic = NA,
        mu_error = NA,
        sigma2_error = NA,
        init_trans_error = NA,
        A_error = NA,
        B_error = NA
      )
    })
    
    # Add results to the appropriate lists
    results$GAS_data$GAS_model <- rbind(results$GAS_data$GAS_model, GAS_on_GAS)
    results$GAS_data$TVP_model <- rbind(results$GAS_data$TVP_model, TVP_on_GAS)
    results$TVP_data$GAS_model <- rbind(results$TVP_data$GAS_model, GAS_on_TVP)
    results$TVP_data$TVP_model <- rbind(results$TVP_data$TVP_model, TVP_on_TVP)
  }
  
  # Summarize results
  summarize_results <- function(df) {
    if (nrow(df) == 0) return(NULL)
    data.frame(
      model = df$model[1],
      datatype = df$datatype[1],
      avg_loglik = mean(df$loglik, na.rm = TRUE),
      avg_aic = mean(df$aic, na.rm = TRUE),
      avg_mu_error = mean(df$mu_error, na.rm = TRUE),
      avg_sigma2_error = mean(df$sigma2_error, na.rm = TRUE),
      avg_init_trans_error = mean(df$init_trans_error, na.rm = TRUE),
      avg_A_error = mean(df$A_error, na.rm = TRUE),
      avg_B_error = mean(df$B_error, na.rm = TRUE)
    )
  }
  
  # Create summary
  summary <- rbind(
    summarize_results(results$GAS_data$GAS_model),
    summarize_results(results$GAS_data$TVP_model),
    summarize_results(results$TVP_data$GAS_model),
    summarize_results(results$TVP_data$TVP_model)
  )
  
  # Return results and summary
  return(list(
    results = results,
    summary = summary,
    setup = setup
  ))
}

#' Plot results from compare_gas_tvp_parameter_recovery
#'
#' @param results Results from compare_gas_tvp_parameter_recovery
#' @param output_file Path to save the output PDF (NULL for no saving)
#' @param width Width of the output PDF
#' @param height Height of the output PDF
#' @return NULL (creates plots as a side effect)
plot_parameter_recovery_comparison <- function(results, output_file = NULL, width = 10, height = 8) {
  # Extract data
  GAS_data_GAS_model <- results$results$GAS_data$GAS_model
  GAS_data_TVP_model <- results$results$GAS_data$TVP_model
  TVP_data_GAS_model <- results$results$TVP_data$GAS_model
  TVP_data_TVP_model <- results$results$TVP_data$TVP_model
  
  # Start plotting device if output_file is provided
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
    on.exit(dev.off())
  }
  
  # Set up plotting area with 2x2 panels
  par(mfrow = c(2, 2))
  
  # Plot 1: AIC comparison for GAS data
  boxplot(list(
    "GAS Model" = GAS_data_GAS_model$aic,
    "TVP Model" = GAS_data_TVP_model$aic
  ), main = "Model Fit on GAS-Generated Data",
  ylab = "AIC (lower is better)",
  col = c("lightblue", "lightgreen"))
  
  # Plot 2: AIC comparison for TVP data
  boxplot(list(
    "GAS Model" = TVP_data_GAS_model$aic,
    "TVP Model" = TVP_data_TVP_model$aic
  ), main = "Model Fit on TVP-Generated Data",
  ylab = "AIC (lower is better)",
  col = c("lightblue", "lightgreen"))
  
  # Plot 3: Parameter recovery for GAS data
  boxplot(list(
    "GAS μ Error" = GAS_data_GAS_model$mu_error,
    "TVP μ Error" = GAS_data_TVP_model$mu_error,
    "GAS σ² Error" = GAS_data_GAS_model$sigma2_error,
    "TVP σ² Error" = GAS_data_TVP_model$sigma2_error,
    "GAS Trans Error" = GAS_data_GAS_model$init_trans_error,
    "TVP Trans Error" = GAS_data_TVP_model$init_trans_error
  ), main = "Parameter Recovery on GAS-Generated Data",
  ylab = "RMSE",
  col = rep(c("lightblue", "lightgreen"), 3),
  las = 2)
  
  # Plot 4: Parameter recovery for TVP data
  boxplot(list(
    "GAS μ Error" = TVP_data_GAS_model$mu_error,
    "TVP μ Error" = TVP_data_TVP_model$mu_error,
    "GAS σ² Error" = TVP_data_GAS_model$sigma2_error,
    "TVP σ² Error" = TVP_data_TVP_model$sigma2_error,
    "GAS Trans Error" = TVP_data_GAS_model$init_trans_error,
    "TVP Trans Error" = TVP_data_TVP_model$init_trans_error
  ), main = "Parameter Recovery on TVP-Generated Data",
  ylab = "RMSE",
  col = rep(c("lightblue", "lightgreen"), 3),
  las = 2)
  
  # Reset plotting parameters if not saving to file
  if (is.null(output_file)) {
    par(mfrow = c(1, 1))
  }
}

#' Add GAS model to existing comparison results from Sendstad et al. (2025)
#'
#' @param empirical_data Time series to analyze
#' @param existing_results List with results from previous model comparisons
#' @param K Number of regimes
#' @param verbose Whether to print progress information
#' @return List with updated comparison results including GAS model
#' @details
#' Takes existing results from the Constant and TVP models and adds
#' results from the GAS model for a complete comparison.
analyze_empirical_data_with_gas <- function(empirical_data, existing_results = NULL, 
                                           K = 3, verbose = TRUE) {
  # Check if empirical data is provided
  if (is.null(empirical_data)) {
    stop("Empirical data must be provided")
  }
  
  if (verbose) {
    cat("Analyzing empirical data with GAS model (K =", K, ")...\n")
  }
  
  # Estimate GAS model
  gas_model <- estimate_gas_model(empirical_data, K = K, verbose = verbose)
  
  # Create model summary
  gas_summary <- data.frame(
    Model = "GAS",
    LogLik = gas_model$diagnostics$loglik,
    AIC = gas_model$diagnostics$aic,
    BIC = gas_model$diagnostics$bic,
    NumParams = gas_model$diagnostics$num_params
  )
  
  # If existing results are provided, combine with new results
  if (!is.null(existing_results)) {
    if (is.data.frame(existing_results)) {
      # If existing_results is a data frame, append our results
      combined_summary <- rbind(existing_results, gas_summary)
    } else if (is.list(existing_results) && "model_comparison" %in% names(existing_results)) {
      # If existing_results is a list with model_comparison, update it
      combined_summary <- rbind(existing_results$model_comparison, gas_summary)
    } else {
      warning("Existing results format not recognized. Returning GAS results only.")
      combined_summary <- gas_summary
    }
  } else {
    combined_summary <- gas_summary
  }
  
  # Sort by AIC
  combined_summary <- combined_summary[order(combined_summary$AIC), ]
  
  # Calculate regime persistence metrics
  persistence <- calculate_persistence(gas_model$filtered_probabilities)
  
  # Create a summary of GAS parameter estimates
  gas_params <- data.frame(
    Model = "GAS",
    Parameter = c(
      paste0("mu", 1:K),
      paste0("sigma2", 1:K),
      paste0("A_mean", 1),
      paste0("B_mean", 1)
    ),
    Value = c(
      gas_model$parameters$mu,
      gas_model$parameters$sigma2,
      mean(gas_model$parameters$A),
      mean(gas_model$parameters$B)
    )
  )
  
  # Return updated comparison results
  results <- list(
    model_comparison = combined_summary,
    gas_model = gas_model,
    gas_params = gas_params,
    gas_persistence = persistence
  )
  
  return(results)
}

#' Create plots comparing regime probabilities from different models
#'
#' @param empirical_data Original time series data
#' @param gas_model Estimated GAS model
#' @param tvp_model Estimated TVP model (optional)
#' @param constant_model Estimated Constant model (optional)
#' @param dates Vector of dates corresponding to the data (optional)
#' @param output_file Path to save the output PDF (NULL for no saving)
#' @param width Width of the output PDF
#' @param height Height of the output PDF
#' @return NULL (creates plots as a side effect)
#' @details
#' Creates plots comparing regime probabilities for GAS, TVP, and Constant models,
#' as well as the observed data and transition probabilities.
plot_model_comparison <- function(empirical_data, gas_model, tvp_model = NULL, 
                                 constant_model = NULL, dates = NULL,
                                 output_file = NULL, width = 12, height = 10) {
  # Start plotting device if output_file is provided
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
    on.exit(dev.off())
  }
  
  # Create x-axis values
  n <- length(empirical_data)
  if (is.null(dates)) {
    x_vals <- 1:n
    x_lab <- "Time"
  } else {
    x_vals <- dates
    x_lab <- "Date"
  }
  
  # Set up layout
  layout_matrix <- matrix(c(1, 1, 2, 3, 4, 5), ncol = 2, byrow = TRUE)
  layout(layout_matrix, heights = c(1, 2, 2))
  
  # Plot 1: Original data
  par(mar = c(4, 4, 2, 2))
  plot(x_vals, empirical_data, type = "l", xlab = x_lab, ylab = "Value",
       main = "Observed Data", col = "black")
  
  # Extract regime probabilities
  K <- ncol(gas_model$filtered_probabilities)
  
  # Define colors for regimes
  regime_colors <- rainbow(K)
  
  # Plot 2: GAS model regime probabilities
  par(mar = c(4, 4, 2, 2))
  matplot(x_vals, gas_model$filtered_probabilities, type = "l", 
          xlab = x_lab, ylab = "Probability",
          main = "GAS Model Regime Probabilities",
          col = regime_colors, lty = 1)
  legend("topright", legend = paste("Regime", 1:K), 
         col = regime_colors, lty = 1, cex = 0.8)
  
  # Plot 3: TVP model regime probabilities (if provided)
  par(mar = c(4, 4, 2, 2))
  if (!is.null(tvp_model)) {
    matplot(x_vals, tvp_model$filtered_probabilities, type = "l", 
            xlab = x_lab, ylab = "Probability",
            main = "TVP Model Regime Probabilities",
            col = regime_colors, lty = 1)
    legend("topright", legend = paste("Regime", 1:K), 
           col = regime_colors, lty = 1, cex = 0.8)
  } else {
    plot(0, 0, type = "n", xlab = "", ylab = "", main = "TVP Model Not Available")
  }
  
  # Plot 4: Constant model regime probabilities (if provided)
  par(mar = c(4, 4, 2, 2))
  if (!is.null(constant_model)) {
    matplot(x_vals, constant_model$filtered_probabilities, type = "l", 
            xlab = x_lab, ylab = "Probability",
            main = "Constant Model Regime Probabilities",
            col = regime_colors, lty = 1)
    legend("topright", legend = paste("Regime", 1:K), 
           col = regime_colors, lty = 1, cex = 0.8)
  } else {
    plot(0, 0, type = "n", xlab = "", ylab = "", main = "Constant Model Not Available")
  }
  
  # Plot 5: GAS model transition probabilities
  par(mar = c(4, 4, 2, 2))
  if (ncol(gas_model$transition_probabilities) == nrow(gas_model$filtered_probabilities)) {
    # Determine how many transition probabilities to plot (avoid clutter)
    K_states <- K
    n_transition <- K_states * (K_states - 1)
    
    # Only plot a subset of transition probabilities if there are many
    if (n_transition > 6) {
      # Select a few important transitions
      selected_transitions <- c(1, 2, n_transition-1, n_transition)
      trans_probs <- gas_model$transition_probabilities[selected_transitions, ]
      trans_names <- paste("Trans", selected_transitions)
    } else {
      trans_probs <- gas_model$transition_probabilities
      trans_names <- paste("Trans", 1:nrow(trans_probs))
    }
    
    matplot(x_vals, t(trans_probs), type = "l", 
            xlab = x_lab, ylab = "Probability",
            main = "GAS Model Transition Probabilities",
            col = 1:nrow(trans_probs), lty = 1)
    legend("topright", legend = trans_names, 
           col = 1:nrow(trans_probs), lty = 1, cex = 0.8)
  } else {
    plot(0, 0, type = "n", xlab = "", ylab = "", 
         main = "Transition Probabilities Not Available")
  }
  
  # Reset plotting parameters if not saving to file
  if (is.null(output_file)) {
    par(mfrow = c(1, 1))
  }
}
