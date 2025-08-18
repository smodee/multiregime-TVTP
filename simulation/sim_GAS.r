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
#' @param n_starts_default Default number of starting points for estimation (default: 5)
#' @param sim_parallel Simulation-level parallelization: "auto", TRUE, FALSE (default: "auto")
#' @param max_cores Maximum cores to use (default: NULL, uses detectCores() - 1)
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method (default: "simple")
#' @param use_fallback Whether to use constant model fallback for small A (default: TRUE)
#' @param verbose Whether to print progress information (default: TRUE)
#' @return Data frame with simulation results
#' @details
#' Conducts a comprehensive simulation study for the GAS model with different
#' sample sizes and measures estimation accuracy and computational performance.
#' 
#' Uses coordinated parallelization that intelligently allocates cores between
#' simulation-level and estimation-level parallelization based on study size.
#' 
#' The GAS model is more complex than other models due to score-driven dynamics,
#' so defaults to more starting points (5) for robust estimation.
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
#' # Run simulation study with auto parallelization
#' results <- run_gas_simulation_study(
#'   num_repetitions = 50,
#'   sample_sizes = c(500, 1000),
#'   K = 3,
#'   param_settings = param_settings,
#'   sim_parallel = "auto"
#' )
#' 
#' # Quick testing mode for development
#' quick_results <- run_gas_simulation_study(
#'   num_repetitions = 5,
#'   sample_sizes = c(500),
#'   K = 3,
#'   n_starts_default = 1,
#'   sim_parallel = FALSE
#' )
run_gas_simulation_study <- function(num_repetitions = 100, 
                                     sample_sizes = c(250, 500, 1000), 
                                     K = 3,
                                     param_settings = NULL,
                                     seed = 123,
                                     output_dir = NULL,
                                     n_starts_default = 5,
                                     sim_parallel = "auto",
                                     max_cores = NULL,
                                     n_nodes = 30,
                                     scaling_method = "simple",
                                     use_fallback = TRUE,
                                     verbose = TRUE) {
  
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Determine available cores
  if (is.null(max_cores)) {
    total_cores <- max(1, parallel::detectCores() - 1)
  } else {
    total_cores <- max(1, min(max_cores, parallel::detectCores() - 1))
  }
  
  # Intelligent parallelization strategy
  if (sim_parallel == "auto") {
    # Use coordinated parallelization for larger studies
    # GAS model is complex, so use parallel earlier than for simpler models
    use_sim_parallel <- (num_repetitions >= 15 && total_cores >= 4)
  } else {
    use_sim_parallel <- as.logical(sim_parallel)
  }
  
  # Calculate optimal core allocation
  if (use_sim_parallel && num_repetitions > 1) {
    # Parallelize both simulation and estimation levels
    sim_workers <- min(total_cores, num_repetitions, 
                       max(2, total_cores %/% 2))  # At least 2 cores per worker
    cores_per_sim <- max(1, total_cores %/% sim_workers)
    n_starts_per_sim <- min(n_starts_default, cores_per_sim)
    
    # Validation check: ensure we don't waste cores
    if (cores_per_sim > n_starts_per_sim && verbose) {
      cat("Note: Adjusting cores per simulation from", cores_per_sim, 
          "to", n_starts_per_sim, "to match starting points\n")
      cores_per_sim <- n_starts_per_sim
    }
    
    parallel_strategy <- "coordinated"
  } else {
    # Use only estimation-level parallelization
    sim_workers <- 1
    cores_per_sim <- total_cores
    n_starts_per_sim <- min(n_starts_default, cores_per_sim)
    parallel_strategy <- "estimation_only"
  }
  
  if (verbose) {
    cat("=== GAS Model Simulation Study Configuration ===\n")
    cat("Repetitions:", num_repetitions, "\n")
    cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
    cat("Total cores available:", total_cores, "\n")
    cat("Parallelization strategy:", parallel_strategy, "\n")
    if (use_sim_parallel) {
      cat("Simulation workers:", sim_workers, "\n")
      cat("Cores per simulation:", cores_per_sim, "\n")
    }
    cat("Starting points per estimation:", n_starts_per_sim, "\n")
    cat("Fallback to constant model:", ifelse(use_fallback, "enabled", "disabled"), "\n")
    cat("===============================================\n\n")
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
  B_burnin <- 100
  C <- 50
  
  # Define the single repetition function for potential parallelization
  run_single_repetition <- function(rep, N) {
    # Time the data generation
    start_time <- Sys.time()
    
    # Generate a single path using GAS data generation function
    data_sim <- dataGASCD(1, N + B_burnin + C, mu, sigma2, init_trans, A, B, burn_in = 0)
    
    # Extract the data
    y <- data_sim[1,]
    
    gen_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # Time the estimation
    start_time <- Sys.time()
    
    # Estimate the model with coordinated parallelization
    tryCatch({
      estimate <- estimate_gas_model(
        y = y,
        K = K,
        B_burnin = B_burnin,
        C = C,
        initial_params = NULL,
        bounds = NULL,
        n_nodes = n_nodes,
        scaling_method = scaling_method,
        use_fallback = use_fallback,
        n_starts = n_starts_per_sim,
        parallel = (cores_per_sim > 1),
        cores = cores_per_sim,
        seed = seed + rep,  # Ensure reproducible but different seeds
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
        NStarts = n_starts_per_sim,
        CoresUsed = cores_per_sim,
        ParallelStrategy = parallel_strategy,
        FallbackUsed = fallback_used
      )
      
      if (verbose && !use_sim_parallel) {
        cat("Repetition", rep, "completed. LogLik:", 
            round(estimate$diagnostics$loglik, 2))
        if (fallback_used) {
          cat(" (fallback used)")
        }
        cat("\n")
      }
      
      # Save detailed results if output directory is provided
      if (!is.null(output_dir)) {
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        
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
            est_time = est_time,
            parallel_strategy = parallel_strategy,
            n_starts = n_starts_per_sim,
            cores_used = cores_per_sim,
            fallback_used = fallback_used
          )
        )
        
        save_file <- file.path(output_dir, 
                               sprintf("gas_sim_N%d_rep%d_K%d.rds", N, rep, K))
        saveRDS(sim_details, save_file)
      }
      
      return(result)
      
    }, error = function(e) {
      if (verbose) {
        cat("Error in repetition", rep, ":", e$message, "\n")
      }
      
      # Return error result
      data.frame(
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
        NStarts = n_starts_per_sim,
        CoresUsed = cores_per_sim,
        ParallelStrategy = parallel_strategy,
        FallbackUsed = NA,
        ErrorMessage = e$message
      )
    })
  }
  
  # Prepare results storage
  all_results <- data.frame()
  
  # Set up parallel processing if using simulation-level parallelization
  if (use_sim_parallel && requireNamespace("future.apply", quietly = TRUE)) {
    original_plan <- future::plan()
    future::plan(future::multisession, workers = sim_workers)
    on.exit(future::plan(original_plan), add = TRUE)
    
    if (verbose) {
      cat("Simulation-level parallelization enabled with", sim_workers, "workers\n\n")
    }
  } else if (use_sim_parallel) {
    warning("future.apply package not available. Running simulations sequentially.")
    use_sim_parallel <- FALSE
  }
  
  # Loop through sample sizes
  for (N in sample_sizes) {
    if (verbose) {
      cat("----- Sample size:", N, "-----\n")
    }
    
    # Run repetitions (either in parallel or sequential)
    if (use_sim_parallel && requireNamespace("future.apply", quietly = TRUE)) {
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
      
      # Check fallback usage
      if ("FallbackUsed" %in% names(valid_results)) {
        fallback_count <- sum(valid_results$FallbackUsed, na.rm = TRUE)
        if (fallback_count > 0) {
          cat("Fallback to constant model used:", fallback_count, "times\n")
        }
      }
      
      cat("Mean parameter errors:\n")
      cat("  Mu:", round(mean(valid_results$MuError, na.rm = TRUE), 4), "\n")
      cat("  Sigma2:", round(mean(valid_results$Sigma2Error, na.rm = TRUE), 4), "\n")
      cat("  InitTrans:", round(mean(valid_results$InitTransError, na.rm = TRUE), 4), "\n")
      cat("  A:", round(mean(valid_results$AError, na.rm = TRUE), 4), "\n")
      cat("  B:", round(mean(valid_results$BError, na.rm = TRUE), 4), "\n")
      cat("Mean LogLik:", round(mean(valid_results$LogLik, na.rm = TRUE), 2), "\n")
      
      # Regime dynamics summary
      if ("TransitionRate" %in% names(valid_results)) {
        cat("Mean transition rate:", round(mean(valid_results$TransitionRate, na.rm = TRUE), 4), "\n")
        cat("Mean regime duration:", round(mean(valid_results$AvgDuration, na.rm = TRUE), 2), "periods\n")
      }
      cat("\n")
    }
  }
  
  if (verbose) {
    total_time <- sum(all_results$EstimationTime, na.rm = TRUE)
    cat("=== Study Completed ===\n")
    cat("Total estimation time:", round(total_time, 2), "seconds\n")
    cat("Strategy used:", parallel_strategy, "\n")
    if (use_sim_parallel) {
      cat("Theoretical speedup from coordination: ~", 
          round(sim_workers, 1), "x\n")
    }
    
    # Overall fallback summary
    if ("FallbackUsed" %in% names(all_results)) {
      total_fallbacks <- sum(all_results$FallbackUsed, na.rm = TRUE)
      if (total_fallbacks > 0) {
        cat("Total fallbacks to constant model:", total_fallbacks, 
            "of", nrow(all_results), "estimations\n")
      }
    }
    cat("======================\n")
  }
  
  return(all_results)
}
