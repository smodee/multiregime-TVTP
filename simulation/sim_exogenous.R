#' Simulation Scripts for the Exogenous Model
#' 
#' This file provides simulation scripts and comparison tools for the exogenous 
#' time-varying transition probability model.

# Load required model implementations
source("models/model_exogenous.R")
source("helpers/parallel_utils.R")

#' Run a comprehensive simulation study for the exogenous model
#'
#' @param num_repetitions Number of simulation repetitions
#' @param sample_sizes Vector of sample sizes to test
#' @param K Number of regimes
#' @param param_settings List with parameter settings
#' @param exo_process Function to generate exogenous process or matrix of pre-generated processes
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
#' Conducts a comprehensive simulation study for the exogenous model with different
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
#' # Define exogenous process generator
#' exo_gen <- function(N) rnorm(N)
#' 
#' # Run simulation study with auto parallelization
#' results <- run_exo_simulation_study(
#'   num_repetitions = 100,
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
    cat("=== EXOGENOUS MODEL SIMULATION STUDY ===\n")
    cat("Repetitions:", num_repetitions, "\n")
    cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
    cat("Total cores available:", total_cores, "(", reserve_cores, "reserved)\n")
    cat("Parallel strategy:", describe_parallel_strategy(allocation), "\n")
    if (allocation$sim_workers > 1) {
      cat("Simulation workers:", allocation$sim_workers, "\n")
      cat("Cores per simulation:", allocation$cores_per_sim, "\n")
    }
    cat("Starting points per estimation:", n_starts_default, "\n")
    cat("========================================\n\n")
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
    
    # Generate exogenous process
    if (is.function(exo_process)) {
      X_exo <- exo_process(N + B + C)
    } else if (is.matrix(exo_process) && nrow(exo_process) >= rep) {
      X_exo <- exo_process[rep, 1:(N + B + C)]
    } else {
      stop("exo_process must be a function or matrix with sufficient rows")
    }
    
    # Generate a single path with exogenous process
    data_sim <- dataTVPXExoCD(1, N + B + C, mu, sigma2, init_trans, A, X_exo, burn_in = 0)
    
    # Extract the data
    y <- data_sim[1,]
    
    gen_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # Time the estimation
    start_time <- Sys.time()
    
    # Estimate the model with proper core allocation
    tryCatch({
      # Only enable estimation parallelization for estimation-only strategy
      use_estimation_parallel <- (allocation$sim_workers == 1 && allocation$cores_per_sim > 1)
      
      estimate <- estimate_exo_model(
        y = y,
        X_exo = X_exo,
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
        InitTransError = NA,
        AError = NA,
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
    results_file <- file.path(output_dir, paste0("exo_simulation_results_", 
                                                 format(Sys.Date(), "%Y%m%d"), ".csv"))
    write.csv(all_results, results_file, row.names = FALSE)
    
    # Save configuration
    config_file <- file.path(output_dir, paste0("exo_simulation_config_", 
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

#' Test various exogenous process specifications
#'
#' @param y Observed data
#' @param exo_processes List of exogenous processes to test
#' @param K Number of regimes (default: 3)
#' @param B Burn-in period (default: 100)
#' @param C Cut-off period (default: 50)
#' @param n_starts Number of starting points for each estimation (default: 3)
#' @param parallel Whether to use parallel processing for estimation (default: TRUE)
#' @param cores Number of cores to use (default: NULL, auto-detect)
#' @param seed Random seed for reproducibility (default: 123)
#' @param verbose Whether to print progress information (default: TRUE)
#' @return Data frame comparing model fit for different exogenous processes
#' @details
#' Tests different exogenous processes as drivers for the transition probabilities
#' and compares their model fit using AIC, BIC and likelihood. Uses multiple starting
#' points for robust estimation and can leverage parallel processing.
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
#' # Compare the processes with robust estimation
#' comparison <- compare_exo_processes(y, exo_processes, n_starts = 5)
#' 
#' # Compare with minimal resources for quick testing
#' quick_comparison <- compare_exo_processes(y, exo_processes, 
#'                                          n_starts = 1, parallel = FALSE)
compare_exo_processes <- function(y, exo_processes, K = 3, B = 100, C = 50,
                                  n_starts = 3, parallel = TRUE, cores = NULL,
                                  seed = 123, verbose = TRUE) {
  
  # Input validation
  if (!is.numeric(y)) {
    stop("Input 'y' must be a numeric vector")
  }
  
  if (!is.list(exo_processes)) {
    stop("Input 'exo_processes' must be a list of exogenous processes")
  }
  
  if (length(exo_processes) == 0) {
    stop("At least one exogenous process must be provided")
  }
  
  # Check that all exogenous processes have the correct length
  for (name in names(exo_processes)) {
    if (length(exo_processes[[name]]) != length(y)) {
      stop(paste("Exogenous process", name, "must have the same length as y"))
    }
  }
  
  # Set up cores and parallelization
  if (is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 1)
  }
  cores <- min(cores, n_starts)
  
  # Determine parallelization strategy
  num_processes <- length(exo_processes)
  if (parallel && num_processes > 1 && cores >= 2) {
    # Use process-level parallelization if we have multiple processes and cores
    use_process_parallel <- TRUE
    cores_per_process <- max(1, cores %/% min(num_processes, cores))
    n_starts_per_process <- min(n_starts, cores_per_process)
    
    if (verbose) {
      cat("=== Exogenous Process Comparison ===\n")
      cat("Processes to test:", num_processes, "\n")
      cat("Parallelization: Process-level with", min(num_processes, cores), "workers\n")
      cat("Cores per process:", cores_per_process, "\n")
      cat("Starting points per process:", n_starts_per_process, "\n")
      cat("===================================\n\n")
    }
  } else {
    # Sequential process testing with full estimation parallelization
    use_process_parallel <- FALSE
    cores_per_process <- cores
    n_starts_per_process <- n_starts
    
    if (verbose) {
      cat("=== Exogenous Process Comparison ===\n")
      cat("Processes to test:", num_processes, "\n")
      cat("Parallelization: Sequential processes, parallel estimation\n")
      cat("Cores per estimation:", cores_per_process, "\n")
      cat("Starting points per estimation:", n_starts_per_process, "\n")
      cat("===================================\n\n")
    }
  }
  
  # Function to test a single exogenous process
  test_single_process <- function(process_info) {
    name <- process_info$name
    X_Exo <- process_info$process
    
    if (verbose && !use_process_parallel) {
      cat("Testing exogenous process:", name, "...")
    }
    
    start_time <- Sys.time()
    
    # Estimate the model
    tryCatch({
      estimate <- estimate_exo_model(
        y = y,
        X_Exo = X_Exo,
        K = K,
        B = B,
        C = C,
        initial_params = NULL,
        bounds = NULL,
        n_starts = n_starts_per_process,
        parallel = (cores_per_process > 1),
        cores = cores_per_process,
        seed = seed + which(names(exo_processes) == name),  # Unique seed per process
        verbose = FALSE
      )
      
      estimation_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Extract convergence information if available
      n_converged <- NA
      n_failed <- NA
      if (n_starts_per_process > 1 && !is.null(estimate$multistart_info)) {
        n_converged <- estimate$multistart_info$n_converged
        n_failed <- estimate$multistart_info$n_failed
      }
      
      # Store results
      result <- data.frame(
        ExoProcess = name,
        LogLik = estimate$diagnostics$loglik,
        AIC = estimate$diagnostics$aic,
        BIC = estimate$diagnostics$bic,
        NumParams = estimate$diagnostics$num_params,
        EstimationTime = estimation_time,
        NStarts = n_starts_per_process,
        NConverged = n_converged,
        NFailed = n_failed,
        stringsAsFactors = FALSE
      )
      
      if (verbose && !use_process_parallel) {
        cat(" completed\n")
        cat("  LogLik:", round(estimate$diagnostics$loglik, 2), 
            ", AIC:", round(estimate$diagnostics$aic, 2), 
            ", BIC:", round(estimate$diagnostics$bic, 2), "\n")
        if (n_starts_per_process > 1 && !is.na(n_converged)) {
          cat("  Convergence:", n_converged, "of", n_starts_per_process, "starts succeeded\n")
        }
        cat("\n")
      }
      
      return(result)
      
    }, error = function(e) {
      estimation_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      if (verbose) {
        if (!use_process_parallel) cat(" failed\n")
        cat("  Error in estimation:", e$message, "\n\n")
      }
      
      # Store error result
      data.frame(
        ExoProcess = name,
        LogLik = NA,
        AIC = NA,
        BIC = NA,
        NumParams = NA,
        EstimationTime = estimation_time,
        NStarts = n_starts_per_process,
        NConverged = NA,
        NFailed = NA,
        ErrorMessage = e$message,
        stringsAsFactors = FALSE
      )
    })
  }
  
  # Prepare process information for testing
  process_list <- lapply(names(exo_processes), function(name) {
    list(name = name, process = exo_processes[[name]])
  })
  
  # Run the comparison (either in parallel or sequential)
  if (use_process_parallel && requireNamespace("future.apply", quietly = TRUE)) {
    # Set up parallel processing for processes
    original_plan <- future::plan()
    future::plan(future::multisession, workers = min(num_processes, cores))
    on.exit(future::plan(original_plan), add = TRUE)
    
    if (verbose) {
      cat("Running", num_processes, "processes in parallel...\n\n")
    }
    
    # Test processes in parallel
    results_list <- future.apply::future_lapply(process_list, test_single_process,
                                                future.seed = TRUE)
    
    # Combine results
    results <- do.call(rbind, results_list)
    
  } else {
    # Test processes sequentially
    results <- data.frame()
    
    for (process_info in process_list) {
      process_result <- test_single_process(process_info)
      results <- rbind(results, process_result)
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
      best_process <- results$ExoProcess[1]
      best_aic <- results$AIC[1]
      cat("Best process:", best_process, "(AIC =", round(best_aic, 2), ")\n")
      
      if (successful_fits > 1) {
        aic_diff <- results$AIC[2] - results$AIC[1]
        cat("AIC improvement over second best:", round(aic_diff, 2), "\n")
      }
      
      # Show timing summary
      total_time <- sum(results$EstimationTime, na.rm = TRUE)
      cat("Total estimation time:", round(total_time, 2), "seconds\n")
    }
    
    # Show failed estimations if any
    failed_processes <- results$ExoProcess[is.na(results$AIC)]
    if (length(failed_processes) > 0) {
      cat("Failed processes:", paste(failed_processes, collapse = ", "), "\n")
    }
    
    cat("=========================\n\n")
  }
  
  return(results)
}

#' Create and test various transformations of an exogenous process
#'
#' @param base_process Base exogenous process
#' @param y Observed data
#' @param K Number of regimes (default: 3)
#' @param transformations List of transformation functions (default: common transformations)
#' @param B Burn-in period (default: 100)
#' @param C Cut-off period (default: 50)
#' @param n_starts Number of starting points for each estimation (default: 3)
#' @param parallel Whether to use parallel processing (default: TRUE)
#' @param cores Number of cores to use (default: NULL, auto-detect)
#' @param seed Random seed for reproducibility (default: 123)
#' @param verbose Whether to print progress information (default: TRUE)
#' @return Data frame comparing model fit for different transformations
#' @details
#' Applies different transformations to a base exogenous process and compares
#' their effectiveness as drivers for transition probabilities using model fit criteria.
#' Uses multiple starting points for robust estimation and can leverage parallel processing.
#' 
#' The function tests common transformations like identity, absolute value, squared,
#' exponential (capped), and lagged versions of the base process. Custom transformations
#' can be provided via the transformations parameter.
#'
#' @examples
#' # Generate data and base process
#' y <- rnorm(1000)
#' base_process <- rnorm(1000)
#' 
#' # Test with default transformations and robust estimation
#' result <- test_exo_transformations(base_process, y, n_starts = 5)
#' 
#' # Test with custom transformations
#' custom_transformations <- list(
#'   identity = function(x) x,
#'   log_abs = function(x) log(abs(x) + 1e-6),
#'   tanh = function(x) tanh(x),
#'   differenced = function(x) c(0, diff(x))
#' )
#' result <- test_exo_transformations(base_process, y, 
#'                                   transformations = custom_transformations)
#' 
#' # Quick testing mode for development
#' quick_result <- test_exo_transformations(base_process, y, 
#'                                         n_starts = 1, parallel = FALSE)
test_exo_transformations <- function(base_process, y, K = 3, 
                                     transformations = list(
                                       identity = function(x) x,
                                       abs = function(x) abs(x),
                                       squared = function(x) x^2,
                                       exp = function(x) exp(pmin(x, 5)),
                                       lagged = function(x) c(NA, x[-length(x)])
                                     ),
                                     B = 100, C = 50,
                                     n_starts = 3, parallel = TRUE, cores = NULL,
                                     seed = 123, verbose = TRUE) {
  
  # Input validation
  if (length(base_process) != length(y)) {
    stop("base_process and y must have the same length")
  }
  
  if (!is.list(transformations) || length(transformations) == 0) {
    stop("transformations must be a non-empty list of functions")
  }
  
  if (length(base_process) < (B + C + 100)) {
    warning("Data length may be too short for reliable estimation with the specified B and C")
  }
  
  if (verbose) {
    cat("=== Testing Exogenous Process Transformations ===\n")
    cat("Base process length:", length(base_process), "\n")
    cat("Transformations to test:", length(transformations), "\n")
    cat("Transformation names:", paste(names(transformations), collapse = ", "), "\n")
    cat("===============================================\n\n")
  }
  
  # Apply transformations to create different exogenous processes
  exo_processes <- list()
  failed_transformations <- character(0)
  
  for (name in names(transformations)) {
    if (verbose) {
      cat("Applying transformation:", name, "...")
    }
    
    tryCatch({
      transformed_process <- transformations[[name]](base_process)
      
      # Validate transformation result
      if (length(transformed_process) != length(base_process)) {
        stop(paste("Transformation", name, "changed the length of the process"))
      }
      
      # Replace NAs and infinite values with appropriate substitutes
      if (any(is.na(transformed_process))) {
        n_na <- sum(is.na(transformed_process))
        if (verbose) cat(" (replacing", n_na, "NAs)")
        transformed_process[is.na(transformed_process)] <- 0
      }
      
      if (any(is.infinite(transformed_process))) {
        n_inf <- sum(is.infinite(transformed_process))
        if (verbose) cat(" (replacing", n_inf, "infinite values)")
        # Replace +Inf with 99th percentile, -Inf with 1st percentile
        finite_vals <- transformed_process[is.finite(transformed_process)]
        if (length(finite_vals) > 0) {
          p99 <- quantile(finite_vals, 0.99, na.rm = TRUE)
          p01 <- quantile(finite_vals, 0.01, na.rm = TRUE)
          transformed_process[transformed_process == Inf] <- p99
          transformed_process[transformed_process == -Inf] <- p01
        } else {
          transformed_process[is.infinite(transformed_process)] <- 0
        }
      }
      
      exo_processes[[name]] <- transformed_process
      
      if (verbose) {
        # Show basic statistics for the transformed process
        stats_summary <- sprintf("(range: [%.3f, %.3f], mean: %.3f)",
                                 min(transformed_process), max(transformed_process),
                                 mean(transformed_process))
        cat(" success", stats_summary, "\n")
      }
      
    }, error = function(e) {
      if (verbose) {
        cat(" failed -", e$message, "\n")
      }
      failed_transformations <- c(failed_transformations, name)
    })
  }
  
  # Check if we have any successful transformations
  if (length(exo_processes) == 0) {
    stop("All transformations failed. Please check your transformation functions.")
  }
  
  if (length(failed_transformations) > 0 && verbose) {
    cat("\nFailed transformations:", paste(failed_transformations, collapse = ", "), "\n")
  }
  
  if (verbose) {
    cat("\nSuccessfully created", length(exo_processes), "transformed processes\n")
    cat("Proceeding to model comparison...\n\n")
  }
  
  # Compare the transformed processes using the updated compare_exo_processes function
  tryCatch({
    results <- compare_exo_processes(
      y = y,
      exo_processes = exo_processes,
      K = K,
      B = B,
      C = C,
      n_starts = n_starts,
      parallel = parallel,
      cores = cores,
      seed = seed,
      verbose = verbose
    )
    
    # Add transformation information to results
    results$TransformationType <- results$ExoProcess
    
    # Add summary statistics for each transformation
    transform_stats <- data.frame(
      ExoProcess = names(exo_processes),
      TransformMean = sapply(exo_processes, mean, na.rm = TRUE),
      TransformSD = sapply(exo_processes, sd, na.rm = TRUE),
      TransformMin = sapply(exo_processes, min, na.rm = TRUE),
      TransformMax = sapply(exo_processes, max, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    # Merge with results
    results <- merge(results, transform_stats, by = "ExoProcess", all.x = TRUE)
    
    # Reorder columns for better readability
    col_order <- c("ExoProcess", "TransformationType", "LogLik", "AIC", "BIC", 
                   "NumParams", "EstimationTime", "NStarts", "NConverged", "NFailed",
                   "TransformMean", "TransformSD", "TransformMin", "TransformMax")
    
    # Only include columns that exist
    col_order <- col_order[col_order %in% names(results)]
    results <- results[, col_order]
    
    if (verbose) {
      cat("=== Transformation Comparison Summary ===\n")
      successful_fits <- sum(!is.na(results$AIC))
      cat("Successful transformations:", successful_fits, "of", nrow(results), "\n")
      
      if (successful_fits > 0) {
        best_transform <- results$ExoProcess[1]
        best_aic <- results$AIC[1]
        cat("Best transformation:", best_transform, "(AIC =", round(best_aic, 2), ")\n")
        
        # Show improvement over identity transformation if it exists
        identity_idx <- which(results$ExoProcess == "identity")
        if (length(identity_idx) > 0 && !is.na(results$AIC[identity_idx])) {
          identity_aic <- results$AIC[identity_idx]
          improvement <- identity_aic - best_aic
          if (improvement > 0) {
            cat("AIC improvement over identity:", round(improvement, 2), "\n")
          } else {
            cat("Identity transformation was optimal (or nearly so)\n")
          }
        }
      }
      
      # Show any failed estimations
      failed_estimations <- results$ExoProcess[is.na(results$AIC)]
      if (length(failed_estimations) > 0) {
        cat("Failed estimations:", paste(failed_estimations, collapse = ", "), "\n")
      }
      
      cat("========================================\n\n")
    }
    
    return(results)
    
  }, error = function(e) {
    stop("Error in model comparison: ", e$message)
  })
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
