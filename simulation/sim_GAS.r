#' Simulation Scripts for the GAS Model
#' 
#' This file provides simulation scripts and comparison tools for the score-driven
#' time-varying transition probability model (GAS model).

# Load required model implementations
source("models/model_GAS.R")
source("helpers/parallel_utils.R")

#' Run a comprehensive simulation study for the GAS model
#'
#' @param num_repetitions Number of simulation repetitions
#' @param sample_sizes Vector of sample sizes to test
#' @param K Number of regimes
#' @param param_settings List with parameter settings
#' @param seed Random seed for reproducibility
#' @param output_dir Directory to save detailed results (NULL for no saving)
#' @param n_starts_default Default number of starting points for estimation (default: 5)
#' @param equal_variances Whether to constrain all regime variances to be equal (default: FALSE)
#' @param sim_parallel Simulation-level parallelization: "auto", TRUE, FALSE (default: "auto")
#' @param max_cores Maximum cores to use (default: NULL, uses conservative detection)
#' @param reserve_cores Number of cores to reserve for system (default: 2)
#' @param use_fallback Whether to enable fallback to constant model on GAS failure (default: TRUE)
#' @param verbose Whether to print progress information (default: TRUE)
#' @return Data frame with simulation results
#' @details
#' Conducts a comprehensive simulation study for the GAS model with different
#' sample sizes and measures estimation accuracy and computational performance.
#' 
#' Uses conservative parallel processing that reserves cores for system use
#' and intelligently allocates between simulation-level and estimation-level
#' parallelization based on available resources.
#'
#' @examples
#' # Define parameter settings
#' param_settings <- list(
#'   mu = c(-2, 1, 2),
#'   sigma2 = c(0.02, 0.2, 0.6),
#'   init_trans = rep(0.2, 6),
#'   A = c(0.1, -0.1, 0.05, -0.05, 0.2, -0.2),
#'   B = c(0.9, 0.8, 0.85, 0.95, 0.75, 0.9)
#' )
#' 
#' # Run simulation study with auto parallelization
#' results <- run_gas_simulation_study(
#'   num_repetitions = 100,
#'   sample_sizes = c(500, 1000),
#'   K = 3,
#'   param_settings = param_settings
#' )
run_gas_simulation_study <- function(num_repetitions = 100, 
                                     sample_sizes = c(250, 500, 1000), 
                                     K = 3,
                                     param_settings = NULL,
                                     seed = 123,
                                     output_dir = NULL,
                                     n_starts_default = 5,
                                     equal_variances = FALSE,
                                     sim_parallel = "auto",
                                     max_cores = NULL,
                                     reserve_cores = 2,
                                     use_fallback = TRUE,
                                     verbose = TRUE) {
  
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Get conservative core count
  total_cores <- get_conservative_cores(max_cores, reserve_cores)
  
  # Get parallel allocation strategy
  allocation <- get_parallel_allocation(total_cores, num_repetitions, n_starts_default)
  
  # Determine if we should use simulation-level parallelization
  if (sim_parallel == "auto") {
    # Auto mode: use allocation strategy (no longer have coordinated mode)
    use_sim_parallel <- (allocation$sim_workers > 1)
  } else {
    use_sim_parallel <- as.logical(sim_parallel) && (allocation$sim_workers > 1)
  }
  
  # Override allocation if user explicitly disabled simulation parallelization
  if (!use_sim_parallel) {
    allocation <- list(
      sim_workers = 1L,
      cores_per_sim = total_cores
    )
  }
  
  if (verbose) {
    cat("=== GAS MODEL SIMULATION STUDY ===\n")
    cat("Repetitions:", num_repetitions, "\n")
    cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
    cat("Total cores available:", total_cores, "(", reserve_cores, "reserved)\n")
    cat("Parallel strategy:", describe_parallel_strategy(allocation), "\n")
    if (allocation$sim_workers > 1) {
      cat("Simulation workers:", allocation$sim_workers, "\n")
      cat("Cores per simulation:", allocation$cores_per_sim, "\n")
    }
    cat("Starting points per estimation:", n_starts_default, "\n")
    cat("Fallback to constant model:", ifelse(use_fallback, "enabled", "disabled"), "\n")
    cat("=================================\n\n")
  }
  
  # Set up parameter defaults
  if (is.null(param_settings)) {
    param_settings <- list(
      mu = seq(-3, 3, length.out = K),
      sigma2 = seq(0.1, 1, length.out = K),
      init_trans = rep(0.2, K*(K-1)),
      A = rep(0.1, K*(K-1)),
      B = rep(0.9, K*(K-1))
    )
  }
  
  # Parameter validation
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
  if (length(param_settings$B) != K*(K-1)) {
    stop("Parameter 'B' must have length K*(K-1)")
  }
  
  # Extract parameters
  mu <- param_settings$mu
  sigma2 <- param_settings$sigma2
  init_trans <- param_settings$init_trans
  A <- param_settings$A
  B <- param_settings$B
  
  # Set burn-in and cut-off
  B_param <- 100  # Rename to avoid confusion with B parameter
  C <- 50
  
  # Initialize results container
  all_results <- data.frame()
  
  # Define the single repetition function for potential parallelization
  run_single_repetition <- function(rep, N) {
    # Time the data generation
    start_time <- Sys.time()
    
    # Generate a single path
    data_sim <- dataGASCD(1, N + B_param + C, mu, sigma2, init_trans, A, B, burn_in = 0)
    
    # Extract the data
    y <- data_sim[1,]
    
    gen_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # Time the estimation
    start_time <- Sys.time()
    
    # Estimate the model with proper core allocation
    tryCatch({
      # Only enable estimation parallelization for estimation-only strategy
      use_estimation_parallel <- (allocation$sim_workers == 1 && allocation$cores_per_sim > 1)
      
      estimate <- estimate_gas_model(
        y = y,
        K = K,
        B = B_param,
        C = C,
        initial_params = NULL,
        bounds = NULL,
        n_starts = n_starts_default,
        equal_variances = equal_variances,
        parallel = use_estimation_parallel,
        cores = allocation$cores_per_sim,
        seed = seed + rep,  # Ensure reproducible but different seeds
        use_fallback = use_fallback,
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
      
      # Calculate regime persistence metrics if available
      transition_rate <- NA
      avg_duration <- NA
      if (!is.null(estimate$filtered_probabilities)) {
        # Simple transition rate calculation
        regime_sequence <- apply(estimate$filtered_probabilities, 1, which.max)
        transitions <- sum(diff(regime_sequence) != 0)
        transition_rate <- transitions / (length(regime_sequence) - 1)
        
        # Average duration calculation
        rle_result <- rle(regime_sequence)
        avg_duration <- mean(rle_result$lengths)
      }
      
      # Check if fallback was used
      fallback_used <- FALSE
      if (use_fallback && !is.null(estimate$model_info$fallback_used)) {
        fallback_used <- estimate$model_info$fallback_used
      }
      
      # Create result
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
        TransitionRate = transition_rate,
        AvgDuration = avg_duration,
        NStarts = n_starts_default,
        equal_variances = equal_variances,
        CoresPerSim = allocation$cores_per_sim,
        ParallelStrategy = describe_parallel_strategy(allocation),
        FallbackUsed = fallback_used
      )
      
      if (verbose && allocation$sim_workers == 1) {
        cat("Repetition", rep, "completed.\n")
      }
      
      return(result)
      
    }, error = function(e) {
      # Handle estimation failures gracefully
      if (verbose) {
        cat("Repetition", rep, "failed:", e$message, "\n")
      }
      
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
        NStarts = n_starts_default,
        equal_variances = equal_variances,
        CoresPerSim = allocation$cores_per_sim,
        ParallelStrategy = describe_parallel_strategy(allocation),
        FallbackUsed = NA
      )
      
      return(result)
    })
  }
  
  # Loop through sample sizes
  for (N in sample_sizes) {
    if (verbose) {
      cat("----- Sample size:", N, "-----\n")
    }
    
    # Run repetitions (either in parallel or sequential)
    if (allocation$sim_workers > 1 && requireNamespace("future.apply", quietly = TRUE)) {
      # Set up simulation-level parallelization
      future::plan(future::multisession, workers = allocation$sim_workers)
      on.exit(future::plan(future::sequential), add = TRUE)
      
      # Parallel execution across repetitions
      rep_results <- future.apply::future_lapply(1:num_repetitions, function(rep) {
        run_single_repetition(rep, N)
      }, future.seed = TRUE)
      
      # Combine results
      sample_results <- do.call(rbind, rep_results)
      
    } else {
      # Sequential execution
      sample_results <- data.frame()
      for (rep in 1:num_repetitions) {
        if (verbose && rep %% 10 == 0) {
          cat("Repetition", rep, "of", num_repetitions, "\n")
        }
        
        rep_result <- run_single_repetition(rep, N)
        sample_results <- rbind(sample_results, rep_result)
      }
    }
    
    # Add to overall results
    all_results <- rbind(all_results, sample_results)
    
    # Summarize results for current sample size
    valid_results <- sample_results[!is.na(sample_results$EstimationTime), ]
    
    if (verbose && nrow(valid_results) > 0) {
      cat("\nSummary for sample size", N, ":\n")
      cat("Successful estimations:", nrow(valid_results), "of", num_repetitions, "\n")
      cat("Mean estimation time:", 
          round(mean(valid_results$EstimationTime, na.rm = TRUE), 3), "seconds\n")
      cat("Mean parameter errors:\n")
      cat("  Mu:", round(mean(valid_results$MuError, na.rm = TRUE), 4), "\n")
      cat("  Sigma2:", round(mean(valid_results$Sigma2Error, na.rm = TRUE), 4), "\n")
      cat("  InitTrans:", round(mean(valid_results$InitTransError, na.rm = TRUE), 4), "\n")
      cat("  A:", round(mean(valid_results$AError, na.rm = TRUE), 4), "\n")
      cat("  B:", round(mean(valid_results$BError, na.rm = TRUE), 4), "\n")
      if (use_fallback) {
        fallback_rate <- mean(valid_results$FallbackUsed, na.rm = TRUE) * 100
        cat("Fallback usage:", round(fallback_rate, 1), "%\n")
      }
      cat("Mean LogLik:", round(mean(valid_results$LogLik, na.rm = TRUE), 2), "\n\n")
    }
  }
  
  if (verbose) {
    total_time <- sum(all_results$EstimationTime, na.rm = TRUE)
    cat("=== STUDY COMPLETED ===\n")
    cat("Total estimation time:", round(total_time, 2), "seconds\n")
    cat("Strategy used:", describe_parallel_strategy(allocation), "\n")
    if (!is.null(output_dir)) {
      cat("Results saved to:", output_dir, "\n")
    }
    cat("=======================\n")
  }
  
  # Save detailed results if requested
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Save main results
    results_file <- file.path(output_dir, paste0("gas_simulation_results_", 
                                                 format(Sys.Date(), "%Y%m%d"), ".csv"))
    write.csv(all_results, results_file, row.names = FALSE)
    
    # Save configuration
    config_file <- file.path(output_dir, paste0("gas_simulation_config_", 
                                                format(Sys.Date(), "%Y%m%d"), ".txt"))
    config_info <- list(
      timestamp = Sys.time(),
      num_repetitions = num_repetitions,
      sample_sizes = sample_sizes,
      K = K,
      param_settings = param_settings,
      seed = seed,
      total_cores = total_cores,
      reserve_cores = reserve_cores,
      allocation = allocation,
      parallel_strategy = describe_parallel_strategy(allocation),
      use_fallback = use_fallback
    )
    
    writeLines(capture.output(str(config_info)), config_file)
  }
  
  return(all_results)
}
