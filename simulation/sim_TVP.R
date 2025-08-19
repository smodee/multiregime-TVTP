#' Simulation Scripts for the TVP Model
#' 
#' This file provides simulation scripts and comparison tools for the time-varying
#' transition probability model with autoregressive dynamics (TVP).

# Load required model implementations
source("models/model_TVP.R")
source("helpers/parallel_utils.R")

#' Run a comprehensive simulation study for the TVP model
#'
#' @param num_repetitions Number of simulation repetitions
#' @param sample_sizes Vector of sample sizes to test
#' @param K Number of regimes
#' @param param_settings List with parameter settings
#' @param seed Random seed for reproducibility
#' @param output_dir Directory to save detailed results (NULL for no saving)
#' @param n_starts_default Default number of starting points for estimation (default: 3)
#' @param sim_parallel Simulation-level parallelization: "auto", TRUE, FALSE (default: "auto")
#' @param max_cores Maximum cores to use (default: NULL, uses conservative detection)
#' @param reserve_cores Number of cores to reserve for system (default: 2)
#' @param verbose Whether to print progress information (default: TRUE)
#' @return Data frame with simulation results
#' @details
#' Conducts a comprehensive simulation study for the TVP model with different
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
#'   A = c(0.1, -0.1, 0.05, -0.05, 0.2, -0.2)
#' )
#' 
#' # Run simulation study with auto parallelization
#' results <- run_tvp_simulation_study(
#'   num_repetitions = 100,
#'   sample_sizes = c(500, 1000),
#'   K = 3,
#'   param_settings = param_settings
#' )
run_tvp_simulation_study <- function(num_repetitions = 100, 
                                     sample_sizes = c(250, 500, 1000), 
                                     K = 3,
                                     param_settings = NULL,
                                     seed = 123,
                                     output_dir = NULL,
                                     n_starts_default = 3,
                                     sim_parallel = "auto",
                                     max_cores = NULL,
                                     reserve_cores = 2,
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
    cat("=== TVP SIMULATION STUDY CONFIGURATION ===\n")
    cat("Repetitions:", num_repetitions, "\n")
    cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
    cat("Total cores available:", total_cores, "(", reserve_cores, "reserved)\n")
    cat("Parallel strategy:", describe_parallel_strategy(allocation), "\n")
    if (allocation$sim_workers > 1) {
      cat("Simulation workers:", allocation$sim_workers, "\n")
      cat("Cores per simulation:", allocation$cores_per_sim, "\n")
    }
    cat("Starting points per estimation:", n_starts_default, "\n")
    cat("==========================================\n\n")
  }
  
  # Set up parameter defaults
  if (is.null(param_settings)) {
    param_settings <- list(
      mu = seq(-3, 3, length.out = K),
      sigma2 = seq(0.1, 1, length.out = K),
      init_trans = rep(0.2, K*(K-1)),
      A = rep(0.1, K*(K-1))
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
  
  # Extract parameters
  mu <- param_settings$mu
  sigma2 <- param_settings$sigma2
  init_trans <- param_settings$init_trans
  A <- param_settings$A
  
  # Set burn-in and cut-off
  B <- 100
  C <- 50
  
  # Initialize results container
  all_results <- data.frame()
  
  # Define the single repetition function for potential parallelization
  run_single_repetition <- function(rep, N) {
    # Time the data generation
    start_time <- Sys.time()
    
    # Generate a single path
    data_sim <- dataTVPCD(1, N + B + C, mu, sigma2, init_trans, A, burn_in = 0)
    
    # Extract the data
    y <- data_sim[1,]
    
    gen_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # Time the estimation
    start_time <- Sys.time()
    
    # Estimate the model with proper core allocation
    tryCatch({
      # Only enable estimation parallelization for estimation-only strategy
      use_estimation_parallel <- (allocation$sim_workers == 1 && allocation$cores_per_sim > 1)
      
      estimate <- estimate_tvp_model(
        y = y,
        K = K,
        B = B,
        C = C,
        initial_params = NULL,
        bounds = NULL,
        n_starts = n_starts_default,
        parallel = use_estimation_parallel,
        cores = allocation$cores_per_sim,
        seed = seed + rep,  # Ensure reproducible but different seeds
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
        NStarts = n_starts_default,
        CoresPerSim = allocation$cores_per_sim,
        ParallelStrategy = describe_parallel_strategy(allocation)
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
        NStarts = n_starts_default,
        CoresPerSim = allocation$cores_per_sim,
        ParallelStrategy = describe_parallel_strategy(allocation)
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
    results_file <- file.path(output_dir, paste0("tvp_simulation_results_", 
                                                 format(Sys.Date(), "%Y%m%d"), ".csv"))
    write.csv(all_results, results_file, row.names = FALSE)
    
    # Save configuration
    config_file <- file.path(output_dir, paste0("tvp_simulation_config_", 
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
      parallel_strategy = describe_parallel_strategy(allocation)
    )
    
    writeLines(capture.output(str(config_info)), config_file)
  }
  
  return(all_results)
}

#' Test various autoregressive coefficient specifications for TVP model
#'
#' @param y Observed data
#' @param A_values Vector of A coefficient values to test
#' @param K Number of regimes (default: 3)
#' @param B Burn-in period (default: 100)
#' @param C Cut-off period (default: 50)
#' @param n_starts Number of starting points for each estimation (default: 3)
#' @param parallel Whether to use parallel processing for estimation (default: TRUE)
#' @param cores Number of cores to use (default: NULL, auto-detect)
#' @param seed Random seed for reproducibility (default: 123)
#' @param verbose Whether to print progress information (default: TRUE)
#' @return Data frame comparing model fit for different A values
#' @details
#' Tests different autoregressive coefficient values for the TVP model and compares
#' their model fit using AIC, BIC and likelihood. Uses multiple starting points for
#' robust estimation and can leverage parallel processing.
#'
#' @examples
#' # Generate some sample data
#' y <- rnorm(1000)
#' 
#' # Test different A coefficient values
#' A_values <- c(0, 0.1, 0.2, 0.5)
#' comparison <- compare_tvp_coefficients(y, A_values, n_starts = 5)
#' 
#' # Quick testing mode for development
#' quick_comparison <- compare_tvp_coefficients(y, A_values, 
#'                                             n_starts = 1, parallel = FALSE)
compare_tvp_coefficients <- function(y, A_values = c(0, 0.1, 0.2, 0.5), 
                                     K = 3, B = 100, C = 50,
                                     n_starts = 3, parallel = TRUE, cores = NULL,
                                     seed = 123, verbose = TRUE) {
  
  # Input validation
  if (!is.numeric(y)) {
    stop("Input 'y' must be a numeric vector")
  }
  
  if (!is.numeric(A_values) || length(A_values) == 0) {
    stop("A_values must be a non-empty numeric vector")
  }
  
  if (length(y) < (B + C + 100)) {
    warning("Data length may be too short for reliable estimation with the specified B and C")
  }
  
  # Set up cores and parallelization
  if (is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 1)
  }
  cores <- min(cores, n_starts)
  
  # Determine parallelization strategy
  num_tests <- length(A_values)
  if (parallel && num_tests > 1 && cores >= 2) {
    # Use test-level parallelization if we have multiple tests and cores
    use_test_parallel <- TRUE
    cores_per_test <- max(1, cores %/% min(num_tests, cores))
    n_starts_per_test <- min(n_starts, cores_per_test)
    
    if (verbose) {
      cat("=== TVP Coefficient Comparison ===\n")
      cat("A values to test:", num_tests, "\n")
      cat("Parallelization: Test-level with", min(num_tests, cores), "workers\n")
      cat("Cores per test:", cores_per_test, "\n")
      cat("Starting points per test:", n_starts_per_test, "\n")
      cat("=================================\n\n")
    }
  } else {
    # Sequential test execution with full estimation parallelization
    use_test_parallel <- FALSE
    cores_per_test <- cores
    n_starts_per_test <- n_starts
    
    if (verbose) {
      cat("=== TVP Coefficient Comparison ===\n")
      cat("A values to test:", num_tests, "\n")
      cat("Parallelization: Sequential tests, parallel estimation\n")
      cat("Cores per estimation:", cores_per_test, "\n")
      cat("Starting points per estimation:", n_starts_per_test, "\n")
      cat("=================================\n\n")
    }
  }
  
  # Function to test a single A value
  test_single_A_value <- function(A_info) {
    A_val <- A_info$value
    A_idx <- A_info$index
    
    if (verbose && !use_test_parallel) {
      cat("Testing A coefficient:", A_val, "...")
    }
    
    start_time <- Sys.time()
    
    # Create A parameter vector (same value for all transitions)
    A_vec <- rep(A_val, K*(K-1))
    
    # Estimate the model
    tryCatch({
      estimate <- estimate_tvp_model(
        y = y,
        K = K,
        B = B,
        C = C,
        initial_params = list(A = A_vec),  # Use the specific A value as starting point
        bounds = NULL,
        n_starts = n_starts_per_test,
        parallel = (cores_per_test > 1),
        cores = cores_per_test,
        seed = seed + A_idx,  # Unique seed per test
        verbose = FALSE
      )
      
      estimation_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Extract convergence information if available
      n_converged <- NA
      n_failed <- NA
      if (n_starts_per_test > 1 && !is.null(estimate$multistart_info)) {
        n_converged <- estimate$multistart_info$n_converged
        n_failed <- estimate$multistart_info$n_failed
      }
      
      # Store results
      result <- data.frame(
        A_Value = A_val,
        LogLik = estimate$diagnostics$loglik,
        AIC = estimate$diagnostics$aic,
        BIC = estimate$diagnostics$bic,
        NumParams = estimate$diagnostics$num_params,
        EstimationTime = estimation_time,
        NStarts = n_starts_per_test,
        NConverged = n_converged,
        NFailed = n_failed,
        stringsAsFactors = FALSE
      )
      
      if (verbose && !use_test_parallel) {
        cat(" completed\n")
        cat("  LogLik:", round(estimate$diagnostics$loglik, 2), 
            ", AIC:", round(estimate$diagnostics$aic, 2), 
            ", BIC:", round(estimate$diagnostics$bic, 2), "\n")
        if (n_starts_per_test > 1 && !is.na(n_converged)) {
          cat("  Convergence:", n_converged, "of", n_starts_per_test, "starts succeeded\n")
        }
        cat("\n")
      }
      
      return(result)
      
    }, error = function(e) {
      estimation_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      if (verbose) {
        if (!use_test_parallel) cat(" failed\n")
        cat("  Error in estimation:", e$message, "\n\n")
      }
      
      # Store error result
      data.frame(
        A_Value = A_val,
        LogLik = NA,
        AIC = NA,
        BIC = NA,
        NumParams = NA,
        EstimationTime = estimation_time,
        NStarts = n_starts_per_test,
        NConverged = NA,
        NFailed = NA,
        ErrorMessage = e$message,
        stringsAsFactors = FALSE
      )
    })
  }
  
  # Prepare test information
  test_list <- lapply(seq_along(A_values), function(i) {
    list(value = A_values[i], index = i)
  })
  
  # Run the comparison (either in parallel or sequential)
  if (use_test_parallel && requireNamespace("future.apply", quietly = TRUE)) {
    # Set up parallel processing for tests
    original_plan <- future::plan()
    future::plan(future::multisession, workers = min(num_tests, cores))
    on.exit(future::plan(original_plan), add = TRUE)
    
    if (verbose) {
      cat("Running", num_tests, "A value tests in parallel...\n\n")
    }
    
    # Test A values in parallel
    results_list <- future.apply::future_lapply(test_list, test_single_A_value,
                                                future.seed = TRUE)
    
    # Combine results
    results <- do.call(rbind, results_list)
    
  } else {
    # Test A values sequentially
    results <- data.frame()
    
    for (test_info in test_list) {
      test_result <- test_single_A_value(test_info)
      results <- rbind(results, test_result)
    }
  }
  
  # Sort results by AIC (best first)
  valid_results <- results[!is.na(results$AIC), ]
  if (nrow(valid_results) > 0) {
    results <- results[order(results$AIC, na.last = TRUE), ]
  }
  
  # Summary output
  if (verbose) {
    cat("=== Comparison Results ===\n")
    successful_fits <- sum(!is.na(results$AIC))
    cat("Successful estimations:", successful_fits, "of", nrow(results), "\n")
    
    if (successful_fits > 0) {
      best_A <- results$A_Value[1]
      best_aic <- results$AIC[1]
      cat("Best A coefficient:", best_A, "(AIC =", round(best_aic, 2), ")\n")
      
      if (successful_fits > 1) {
        aic_diff <- results$AIC[2] - results$AIC[1]
        cat("AIC improvement over second best:", round(aic_diff, 2), "\n")
      }
      
      # Show timing summary
      total_time <- sum(results$EstimationTime, na.rm = TRUE)
      cat("Total estimation time:", round(total_time, 2), "seconds\n")
    }
    
    # Show failed estimations if any
    failed_tests <- results$A_Value[is.na(results$AIC)]
    if (length(failed_tests) > 0) {
      cat("Failed A values:", paste(failed_tests, collapse = ", "), "\n")
    }
    
    cat("=========================\n\n")
  }
  
  return(results)
}

#' Generate test datasets with known TVP characteristics
#'
#' @param N Length of the time series
#' @param K Number of regimes
#' @param A_strength Strength of autoregressive effects (default: 0.2)
#' @param seed Random seed for reproducibility
#' @return List of different TVP time series
#' @details
#' Generates a variety of time series with different TVP characteristics
#' that can be used for model testing and validation.
#'
#' @examples
#' # Generate TVP datasets of length 1000
#' datasets <- generate_tvp_test_datasets(1000)
#' 
#' # Plot the first few datasets
#' par(mfrow=c(2,2))
#' for(i in 1:4) {
#'   name <- names(datasets)[i]
#'   plot(datasets[[name]], type="l", main=name)
#' }
generate_tvp_test_datasets <- function(N = 1000, K = 3, A_strength = 0.2, seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Initialize list to store datasets
  datasets <- list()
  
  # Common parameters
  mu <- seq(-2, 2, length.out = K)
  sigma2 <- rep(0.5, K)
  init_trans <- rep(0.2, K*(K-1))
  
  # 1. Weak autoregressive effects
  A_weak <- rep(0.05, K*(K-1))
  datasets$weak_AR <- dataTVPCD(1, N, mu, sigma2, init_trans, A_weak)[1,]
  
  # 2. Medium autoregressive effects  
  A_medium <- rep(A_strength, K*(K-1))
  datasets$medium_AR <- dataTVPCD(1, N, mu, sigma2, init_trans, A_medium)[1,]
  
  # 3. Strong autoregressive effects
  A_strong <- rep(0.5, K*(K-1))
  datasets$strong_AR <- dataTVPCD(1, N, mu, sigma2, init_trans, A_strong)[1,]
  
  # 4. Asymmetric effects (different A for different transitions)
  A_asymmetric <- c(0.1, -0.1, 0.3, -0.2, 0.2, -0.3)[1:(K*(K-1))]
  datasets$asymmetric_AR <- dataTVPCD(1, N, mu, sigma2, init_trans, A_asymmetric)[1,]
  
  # 5. High persistence regimes
  mu_persistent <- c(-3, 0, 3)
  sigma2_persistent <- c(0.1, 0.1, 0.1)
  init_trans_persistent <- rep(0.05, K*(K-1))  # Low transition probabilities
  A_persistent <- rep(0.1, K*(K-1))
  datasets$high_persistence <- dataTVPCD(1, N, mu_persistent, sigma2_persistent, 
                                         init_trans_persistent, A_persistent)[1,]
  
  # 6. High switching frequency
  init_trans_frequent <- rep(0.4, K*(K-1))  # High transition probabilities
  A_frequent <- rep(0.2, K*(K-1))
  datasets$frequent_switching <- dataTVPCD(1, N, mu, sigma2, 
                                           init_trans_frequent, A_frequent)[1,]
  
  # 7. Heteroskedastic regimes
  mu_hetero <- c(-1, 0, 1)
  sigma2_hetero <- c(0.1, 0.5, 1.5)  # Different volatilities
  A_hetero <- rep(0.15, K*(K-1))
  datasets$heteroskedastic <- dataTVPCD(1, N, mu_hetero, sigma2_hetero, 
                                        init_trans, A_hetero)[1,]
  
  return(datasets)
}
