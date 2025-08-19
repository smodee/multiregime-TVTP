#' Simulation Scripts for the Constant Model
#' 
#' This file provides simulation scripts and comparison tools for the constant 
#' transition probability regime-switching model.

# Load required model implementations
source("models/model_constant.R")
source("helpers/parallel_utils.R")

#' Run a comprehensive simulation study for the constant model
#'
#' @param num_repetitions Number of simulation repetitions
#' @param sample_sizes Vector of sample sizes to test
#' @param K Number of regimes
#' @param param_settings List with parameter settings
#' @param seed Random seed for reproducibility
#' @param output_dir Directory to save detailed results (NULL for no saving)
#' @param n_starts_default Default number of starting points for estimation (default: 3)
#' @param equal_variances Whether to constrain all regime variances to be equal (default: FALSE)
#' @param sim_parallel Simulation-level parallelization: "auto", TRUE, FALSE (default: "auto")
#' @param max_cores Maximum cores to use (default: NULL, uses conservative detection)
#' @param reserve_cores Number of cores to reserve for system (default: 2)
#' @param verbose Whether to print progress information (default: TRUE)
#' @return Data frame with simulation results
#' @details
#' Conducts a comprehensive simulation study for the constant model with different
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
#'   trans_prob = c(0.85, 0.1, 0.05, 0.9, 0.05, 0.95)
#' )
#' 
#' # Run simulation study with auto parallelization
#' results <- run_const_simulation_study(
#'   num_repetitions = 100,
#'   sample_sizes = c(500, 1000),
#'   K = 3,
#'   param_settings = param_settings
#' )
run_const_simulation_study <- function(num_repetitions = 100, 
                                       sample_sizes = c(250, 500, 1000), 
                                       K = 3,
                                       param_settings = NULL,
                                       seed = 123,
                                       output_dir = NULL,
                                       n_starts_default = 3,
                                       equal_variances = FALSE,
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
    cat("=== CONSTANT MODEL SIMULATION STUDY ===\n")
    cat("Repetitions:", num_repetitions, "\n")
    cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
    cat("Total cores available:", total_cores, "(", reserve_cores, "reserved)\n")
    cat("Parallel strategy:", describe_parallel_strategy(allocation), "\n")
    if (allocation$sim_workers > 1) {
      cat("Simulation workers:", allocation$sim_workers, "\n")
      cat("Cores per simulation:", allocation$cores_per_sim, "\n")
    }
    cat("Starting points per estimation:", n_starts_default, "\n")
    cat("=======================================\n\n")
  }
  
  # Set up parameter defaults
  if (is.null(param_settings)) {
    # Create default transition probabilities for K regimes
    # High diagonal persistence (0.85) with equal off-diagonal probabilities
    n_transition <- K * (K - 1)
    trans_prob <- numeric(n_transition)
    
    # Fill transition probabilities
    idx <- 1
    for (i in 1:K) {
      off_diag_prob <- 0.15 / (K - 1)  # Equal off-diagonal probabilities
      for (j in 1:K) {
        if (i != j) {
          trans_prob[idx] <- off_diag_prob
          idx <- idx + 1
        }
      }
    }
    
    param_settings <- list(
      mu = seq(-3, 3, length.out = K),
      sigma2 = seq(0.1, 1, length.out = K),
      trans_prob = trans_prob
    )
  }
  
  # Parameter validation
  if (length(param_settings$mu) != K) {
    stop("Parameter 'mu' must have length K")
  }
  if (length(param_settings$sigma2) != K) {
    stop("Parameter 'sigma2' must have length K")
  }
  if (length(param_settings$trans_prob) != K*(K-1)) {
    stop("Parameter 'trans_prob' must have length K*(K-1)")
  }
  
  # Extract parameters
  mu <- param_settings$mu
  sigma2 <- param_settings$sigma2
  trans_prob <- param_settings$trans_prob
  
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
    data_sim <- dataConstCD(1, N + B + C, mu, sigma2, trans_prob, burn_in = 0)
    
    # Extract the data
    y <- data_sim[1,]
    
    gen_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # Time the estimation
    start_time <- Sys.time()
    
    # Estimate the model with proper core allocation
    tryCatch({
      # Only enable estimation parallelization for estimation-only strategy
      use_estimation_parallel <- (allocation$sim_workers == 1 && allocation$cores_per_sim > 1)
      
      estimate <- estimate_const_model(
        y = y,
        K = K,
        B = B,
        C = C,
        initial_params = NULL,
        bounds = NULL,
        n_starts = n_starts_default,
        equal_variances = equal_variances,
        parallel = use_estimation_parallel,
        cores = allocation$cores_per_sim,
        seed = seed + rep,  # Ensure reproducible but different seeds
        verbose = FALSE
      )
      
      est_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Extract estimated parameters
      mu_est <- estimate$parameters$mu
      sigma2_est <- estimate$parameters$sigma2
      trans_prob_est <- estimate$parameters$trans_prob
      
      # Calculate estimation errors
      mu_error <- sqrt(mean((mu_est - mu)^2))
      sigma2_error <- sqrt(mean((sigma2_est - sigma2)^2))
      trans_prob_error <- sqrt(mean((trans_prob_est - trans_prob)^2))
      
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
        TransProbError = trans_prob_error,
        NStarts = n_starts_default,
        equal_variances = equal_variances,
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
        TransProbError = NA,
        NStarts = n_starts_default,
        equal_variances = equal_variances,
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
      cat("  TransProb:", round(mean(valid_results$TransProbError, na.rm = TRUE), 4), "\n")
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
    results_file <- file.path(output_dir, paste0("const_simulation_results_", 
                                                 format(Sys.Date(), "%Y%m%d"), ".csv"))
    write.csv(all_results, results_file, row.names = FALSE)
    
    # Save configuration
    config_file <- file.path(output_dir, paste0("const_simulation_config_", 
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

#' Compare different starting point strategies for constant model estimation
#'
#' @param y Observed data
#' @param K Number of regimes (default: 3)
#' @param B Burn-in period (default: 100)
#' @param C Cut-off period (default: 50)
#' @param starting_strategies List of starting point strategies to test (default: common strategies)
#' @param parallel Whether to use parallel processing for estimation (default: TRUE)
#' @param cores Number of cores to use (default: NULL, auto-detect)
#' @param seed Random seed for reproducibility (default: 123)
#' @param verbose Whether to print progress information (default: TRUE)
#' @return Data frame comparing model fit for different starting strategies
#' @details
#' Tests different starting point strategies for the constant model estimation
#' and compares their effectiveness using model fit criteria and convergence rates.
#' This is useful for understanding the robustness of estimation to initialization.
#'
#' @examples
#' # Generate some sample data
#' y <- rnorm(1000)
#' 
#' # Compare starting strategies with default options
#' comparison <- compare_starting_strategies(y, K = 3)
#' 
#' # Compare with custom strategies
#' custom_strategies <- list(
#'   single_start = 1,
#'   few_starts = 3,
#'   many_starts = 10,
#'   extensive_starts = 20
#' )
#' result <- compare_starting_strategies(y, starting_strategies = custom_strategies)
compare_starting_strategies <- function(y, K = 3, B = 100, C = 50,
                                        starting_strategies = list(
                                          single_start = 1,
                                          few_starts = 3,
                                          moderate_starts = 5,
                                          many_starts = 10
                                        ),
                                        parallel = TRUE, cores = NULL,
                                        seed = 123, verbose = TRUE) {
  
  # Input validation
  if (!is.numeric(y)) {
    stop("Input 'y' must be a numeric vector")
  }
  
  if (!is.list(starting_strategies) || length(starting_strategies) == 0) {
    stop("starting_strategies must be a non-empty list")
  }
  
  if (length(y) < (B + C + 100)) {
    warning("Data length may be too short for reliable estimation with the specified B and C")
  }
  
  # Set up cores and parallelization
  if (is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Determine parallelization strategy
  num_strategies <- length(starting_strategies)
  max_starts <- max(unlist(starting_strategies))
  cores <- min(cores, max_starts)  # Don't use more cores than max starts
  
  if (parallel && num_strategies > 1 && cores >= 2) {
    # Use strategy-level parallelization if we have multiple strategies and cores
    use_strategy_parallel <- TRUE
    cores_per_strategy <- max(1, cores %/% min(num_strategies, cores))
    
    if (verbose) {
      cat("=== Starting Strategy Comparison ===\n")
      cat("Strategies to test:", num_strategies, "\n")
      cat("Parallelization: Strategy-level with", min(num_strategies, cores), "workers\n")
      cat("Cores per strategy:", cores_per_strategy, "\n")
      cat("====================================\n\n")
    }
  } else {
    # Sequential strategy testing with full estimation parallelization
    use_strategy_parallel <- FALSE
    cores_per_strategy <- cores
    
    if (verbose) {
      cat("=== Starting Strategy Comparison ===\n")
      cat("Strategies to test:", num_strategies, "\n")
      cat("Parallelization: Sequential strategies, parallel estimation\n")
      cat("Cores per estimation:", cores_per_strategy, "\n")
      cat("====================================\n\n")
    }
  }
  
  # Function to test a single starting strategy
  test_single_strategy <- function(strategy_info) {
    name <- strategy_info$name
    n_starts <- strategy_info$n_starts
    
    if (verbose && !use_strategy_parallel) {
      cat("Testing starting strategy:", name, "(", n_starts, "starts)...")
    }
    
    start_time <- Sys.time()
    
    # Estimate the model
    tryCatch({
      estimate <- estimate_const_model(
        y = y,
        K = K,
        B = B,
        C = C,
        initial_params = NULL,
        bounds = NULL,
        n_starts = n_starts,
        parallel = (cores_per_strategy > 1 && n_starts > 1),
        cores = min(cores_per_strategy, n_starts),
        seed = seed + which(names(starting_strategies) == name),  # Unique seed per strategy
        verbose = FALSE
      )
      
      estimation_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Extract convergence information if available
      n_converged <- NA
      n_failed <- NA
      if (n_starts > 1 && !is.null(estimate$multistart_info)) {
        n_converged <- estimate$multistart_info$n_converged
        n_failed <- estimate$multistart_info$n_failed
      }
      
      # Store results
      result <- data.frame(
        Strategy = name,
        NStarts = n_starts,
        LogLik = estimate$diagnostics$loglik,
        AIC = estimate$diagnostics$aic,
        BIC = estimate$diagnostics$bic,
        NumParams = estimate$diagnostics$num_params,
        EstimationTime = estimation_time,
        NConverged = n_converged,
        NFailed = n_failed,
        TimePerStart = estimation_time / n_starts,
        stringsAsFactors = FALSE
      )
      
      if (verbose && !use_strategy_parallel) {
        cat(" completed\n")
        cat("  LogLik:", round(estimate$diagnostics$loglik, 2), 
            ", AIC:", round(estimate$diagnostics$aic, 2), 
            ", Time:", round(estimation_time, 2), "s\n")
        if (n_starts > 1 && !is.na(n_converged)) {
          cat("  Convergence:", n_converged, "of", n_starts, "starts succeeded\n")
        }
        cat("\n")
      }
      
      return(result)
      
    }, error = function(e) {
      estimation_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      if (verbose) {
        if (!use_strategy_parallel) cat(" failed\n")
        cat("  Error in estimation:", e$message, "\n\n")
      }
      
      # Store error result
      data.frame(
        Strategy = name,
        NStarts = n_starts,
        LogLik = NA,
        AIC = NA,
        BIC = NA,
        NumParams = NA,
        EstimationTime = estimation_time,
        NConverged = NA,
        NFailed = NA,
        TimePerStart = NA,
        ErrorMessage = e$message,
        stringsAsFactors = FALSE
      )
    })
  }
  
  # Prepare strategy information for testing
  strategy_list <- lapply(names(starting_strategies), function(name) {
    list(name = name, n_starts = starting_strategies[[name]])
  })
  
  # Run the comparison (either in parallel or sequential)
  if (use_strategy_parallel && requireNamespace("future.apply", quietly = TRUE)) {
    # Set up parallel processing for strategies
    original_plan <- future::plan()
    future::plan(future::multisession, workers = min(num_strategies, cores))
    on.exit(future::plan(original_plan), add = TRUE)
    
    if (verbose) {
      cat("Running", num_strategies, "strategies in parallel...\n\n")
    }
    
    # Test strategies in parallel
    results_list <- future.apply::future_lapply(strategy_list, test_single_strategy,
                                                future.seed = TRUE)
    
    # Combine results
    results <- do.call(rbind, results_list)
    
  } else {
    # Test strategies sequentially
    results <- data.frame()
    
    for (strategy_info in strategy_list) {
      strategy_result <- test_single_strategy(strategy_info)
      results <- rbind(results, strategy_result)
    }
  }
  
  # Sort results by AIC (best first)
  valid_results <- results[!is.na(results$AIC), ]
  if (nrow(valid_results) > 0) {
    results <- results[order(results$AIC, na.last = TRUE), ]
  }
  
  # Summary output
  if (verbose) {
    cat("=== Strategy Comparison Results ===\n")
    successful_fits <- sum(!is.na(results$AIC))
    cat("Successful estimations:", successful_fits, "of", nrow(results), "\n")
    
    if (successful_fits > 0) {
      best_strategy <- results$Strategy[1]
      best_aic <- results$AIC[1]
      best_time <- results$EstimationTime[1]
      cat("Best strategy:", best_strategy, "(AIC =", round(best_aic, 2), 
          ", Time =", round(best_time, 2), "s)\n")
      
      if (successful_fits > 1) {
        aic_diff <- results$AIC[2] - results$AIC[1]
        cat("AIC improvement over second best:", round(aic_diff, 2), "\n")
      }
      
      # Show efficiency analysis
      if (nrow(valid_results) > 1) {
        cat("\nEfficiency analysis:\n")
        for (i in 1:nrow(valid_results)) {
          row <- valid_results[i, ]
          efficiency_score <- (best_aic - row$AIC) / row$EstimationTime
          cat("  ", row$Strategy, ": Efficiency score =", round(efficiency_score, 4), "\n")
        }
      }
      
      # Show timing summary
      total_time <- sum(results$EstimationTime, na.rm = TRUE)
      cat("Total estimation time:", round(total_time, 2), "seconds\n")
    }
    
    # Show failed estimations if any
    failed_strategies <- results$Strategy[is.na(results$AIC)]
    if (length(failed_strategies) > 0) {
      cat("Failed strategies:", paste(failed_strategies, collapse = ", "), "\n")
    }
    
    cat("===================================\n\n")
  }
  
  return(results)
}

#' Generate synthetic datasets with known constant transition probabilities
#'
#' @param N Length of each dataset
#' @param num_datasets Number of datasets to generate
#' @param K Number of regimes (default: 3)
#' @param regime_types List of regime characteristic types (default: common types)
#' @param seed Random seed for reproducibility (default: 123)
#' @return List of datasets with different regime characteristics
#' @details
#' Generates multiple datasets with different regime characteristics to test
#' the constant model's performance across various scenarios. Each dataset
#' has known parameters for validation purposes.
#'
#' @examples
#' # Generate datasets of length 1000
#' datasets <- generate_const_test_datasets(1000, 5)
#' 
#' # Generate datasets with custom regime types
#' custom_types <- list(
#'   similar_regimes = list(mu_spread = 0.5, sigma_spread = 0.1),
#'   distinct_regimes = list(mu_spread = 3.0, sigma_spread = 1.0)
#' )
#' datasets <- generate_const_test_datasets(1000, 2, regime_types = custom_types)
generate_const_test_datasets <- function(N, num_datasets = 5, K = 3,
                                         regime_types = list(
                                           balanced_persistent = list(
                                             mu_spread = 2, sigma_spread = 0.5, persistence = 0.9
                                           ),
                                           balanced_switching = list(
                                             mu_spread = 2, sigma_spread = 0.5, persistence = 0.7
                                           ),
                                           distinct_persistent = list(
                                             mu_spread = 4, sigma_spread = 1.0, persistence = 0.9
                                           ),
                                           similar_switching = list(
                                             mu_spread = 1, sigma_spread = 0.2, persistence = 0.6
                                           ),
                                           heteroskedastic = list(
                                             mu_spread = 2, sigma_spread = 2.0, persistence = 0.8
                                           )
                                         ),
                                         seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Initialize results list
  datasets <- list()
  
  # Generate datasets for each regime type
  for (type_name in names(regime_types)) {
    type_params <- regime_types[[type_name]]
    
    # Extract parameters with defaults
    mu_spread <- ifelse(is.null(type_params$mu_spread), 2, type_params$mu_spread)
    sigma_spread <- ifelse(is.null(type_params$sigma_spread), 0.5, type_params$sigma_spread)
    persistence <- ifelse(is.null(type_params$persistence), 0.8, type_params$persistence)
    
    # Generate multiple datasets of this type
    for (dataset_idx in 1:num_datasets) {
      # Create regime parameters
      mu <- seq(-mu_spread/2, mu_spread/2, length.out = K)
      sigma2 <- seq(0.1, 0.1 + sigma_spread, length.out = K)
      
      # Create transition probability matrix with specified persistence
      n_transition <- K * (K - 1)
      trans_prob <- numeric(n_transition)
      
      # Fill transition probabilities
      idx <- 1
      for (i in 1:K) {
        off_diag_prob <- (1 - persistence) / (K - 1)  # Equal off-diagonal probabilities
        for (j in 1:K) {
          if (i != j) {
            trans_prob[idx] <- off_diag_prob
            idx <- idx + 1
          }
        }
      }
      
      # Generate data
      data_sim <- dataConstCD(1, N, mu, sigma2, trans_prob, burn_in = 100)
      y <- data_sim[1,]
      
      # Store dataset with metadata
      dataset_name <- paste0(type_name, "_", dataset_idx)
      datasets[[dataset_name]] <- list(
        data = y,
        true_parameters = list(
          mu = mu,
          sigma2 = sigma2,
          trans_prob = trans_prob
        ),
        regime_type = type_name,
        dataset_index = dataset_idx,
        characteristics = list(
          mu_spread = mu_spread,
          sigma_spread = sigma_spread,
          persistence = persistence,
          K = K,
          N = N
        )
      )
    }
  }
  
  return(datasets)
}
