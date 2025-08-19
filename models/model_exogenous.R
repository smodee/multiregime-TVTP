#' Time-Varying Transition Probability Models with Exogenous Variables
#' 
#' This file implements a regime-switching model where transition probabilities
#' depend on exogenous variables.
#'
#' Reference: Sendstad, Chronopoulos, & Li (2025) - The Value of Turning-Point 
#' Detection for Optimal Investment

# Load required helper functions
source("helpers/utility_functions.R")
source("helpers/transition_helpers.R")
source("helpers/parameter_transforms.R")

#' Generate data from an exogenous variable-driven regime-switching model
#'
#' @param M Number of simulation runs to be performed
#' @param N Length of the simulation runs (discretized time)
#' @param mu Vector of true means corresponding to each regime
#' @param sigma2 Vector of true variances corresponding to each regime
#' @param init_trans Initial transition probabilities for the latent process
#' @param A Exogenous factor weights, one for each transition probability
#' @param X_Exo Exogenous process that drives transition probability changes
#' @return Matrix of simulated data with M rows and N columns
#' @details 
#' Simulates data from a regime switching model where transition probabilities
#' depend on values of an exogenous process. The relationship is defined as:
#' f[t+1] = omega + A * X_Exo[t], followed by logistic transformation.
#'
#' @examples
#' # Generate data for a 3-regime model
#' mu_true <- c(-2, 1, 2)
#' sigma2_true <- c(0.02, 0.2, 0.6)
#' init_trans_true <- rep(0.2, 6)  # Off-diagonal elements for a 3x3 matrix
#' A_true <- c(0.1, -0.1, 0.05, -0.05, 0.2, -0.2)
#' X_Exo <- rnorm(1000)  # Generate exogenous process
#' data_sim <- dataexoCD(10, 1000, mu_true, sigma2_true, init_trans_true, A_true, X_Exo)
dataexoCD <- function(M, N, mu, sigma2, init_trans, A, X_Exo) {
  # Ensure that means and variances have been provided for all regimes
  if (length(mu) != length(sigma2)) {
    stop("Error: Unequal number of means and variances. Mean and variance have to be supplied for each regime.")
  }
  
  # Determine the number of regimes and transition parameters for the model
  K <- length(mu)
  n_transition <- K*(K-1)
  
  # Check that X_Exo is long enough
  if (length(X_Exo) < N) {
    stop("Error: Exogenous variable series is too short for the requested simulation length.")
  }
  
  # Ensure we have the right number of transition probabilities and A parameters
  if (length(init_trans) != n_transition) {
    stop(sprintf("Error: Expected %d transition probabilities for %d regimes, got %d.", 
                 n_transition, K, length(init_trans)))
  }
  if (length(A) != n_transition) {
    stop(sprintf("Error: Expected %d A parameters for %d regimes, got %d.", 
                 n_transition, K, length(A)))
  }
  
  # Set up a matrix to save the output
  data <- matrix(0, M, N)
  
  # Get baseline transition probabilities that are logit-transformed
  omega_LR <- logit(init_trans)
  omega <- omega_LR * (1 - A)  # This scaling ensures that when A=0, we get back to init_trans
  
  for (i in 1:M) {
    # Initialize all data structures
    eta <- matrix(0, nrow=K, ncol=N)     # Likelihood of each regime
    tot_lik <- numeric(N)                # Total likelihood
    X_t <- matrix(0, nrow=K, ncol=N)     # Filtered probabilities after observation
    X_tlag <- matrix(0, nrow=K, ncol=N)  # Predicted probabilities before observation
    S <- numeric(N)                      # Latent state
    y.sim <- numeric(N)                  # Random increments according to state
    
    # Initialize f values for transition probabilities and convert to valid probabilities
    f <- matrix(0, nrow=n_transition, ncol=N)
    p_trans <- matrix(0, nrow=n_transition, ncol=N)
    
    # Initial state probabilities (uniform distribution)
    X_t[,1] <- rep(1/K, K)
    
    # Set initial values
    f[,1] <- omega + A * X_Exo[1]
    p_trans_raw <- logistic(f[,1])
    p_trans[,1] <- convert_to_valid_probs(p_trans_raw, K)
    
    for (t in 1:(N-1)) {
      # Generate predicted probabilities
      Pmatrix <- transition_matrix(p_trans[,t], check_validity = FALSE)
      X_tlag[,t] <- Pmatrix %*% X_t[,t]
      
      # Sample a state based on the predicted probabilities and 
      # simulate data conditional on that state
      S[t] <- sample(1:K, 1, prob=X_tlag[,t])
      y.sim[t] <- rnorm(1, mu[S[t]], sqrt(sigma2[S[t]]))
      
      # Calculate likelihoods
      for (k in 1:K) {
        eta[k,t] <- dnorm(y.sim[t], mu[k], sqrt(sigma2[k]))
      }
      tot_lik[t] <- sum(eta[,t]*X_tlag[,t])
      
      # Calculate filtered probabilities
      X_t[,t+1] <- (eta[,t]*X_tlag[,t])/tot_lik[t]
      
      # Update transition probabilities based on the exogenous variable
      f[,t+1] <- omega + A * X_Exo[t+1]
      p_trans_raw <- logistic(f[,t+1])
      p_trans[,t+1] <- convert_to_valid_probs(p_trans_raw, K)
    }
    
    # For the last time point
    Pmatrix <- transition_matrix(p_trans[,N-1], check_validity = FALSE)
    X_tlag[,N] <- Pmatrix %*% X_t[,N-1]
    S[N] <- sample(1:K, 1, prob=X_tlag[,N])
    y.sim[N] <- rnorm(1, mu[S[N]], sqrt(sigma2[S[N]]))
    
    # Save the simulation run in the data matrix
    data[i,] <- y.sim
  }
  
  return(data)
}

#' Filter observed data through the exogenous variable-driven regime-switching model
#'
#' @param mu Vector of means corresponding to each regime
#' @param sigma2 Vector of variances corresponding to each regime
#' @param init_trans Initial transition probabilities for the latent process
#' @param A Exogenous factor weights, one for each transition probability
#' @param X_Exo Exogenous process that drives transition probability changes
#' @param y Observed time series increments
#' @param B Burn-in to be excluded at the beginning of the time series
#' @param C Cut-off to be excluded at the end of the time series
#' @return Negative log-likelihood of observed data under the model
#' @details
#' Filters observed data through the model to compute the likelihood.
#' Returns the negative log-likelihood for compatibility with optimization functions.
#'
#' @examples
#' # Calculate likelihood for a 3-regime model
#' mu <- c(-2, 1, 2)
#' sigma2 <- c(0.02, 0.2, 0.6)
#' init_trans <- rep(0.2, 6)
#' A <- c(0.1, -0.1, 0.05, -0.05, 0.2, -0.2)
#' X_Exo <- rnorm(1000)
#' y <- rnorm(1000) 
#' loglik <- Rfiltering_TVPXExo(mu, sigma2, init_trans, A, X_Exo, y, 100, 50)
Rfiltering_TVPXExo <- function(mu, sigma2, init_trans, A, X_Exo, y, B, C) {
  # Ensure that means and variances have been provided for all regimes
  if (length(mu) != length(sigma2)) {
    stop("Error: Unequal number of means and variances. Mean and variance have to be supplied for each regime.")
  }
  
  # Determine length of the time series
  M <- length(y)
  
  # Check that X_Exo is long enough
  if (length(X_Exo) < M) {
    stop("Error: Exogenous variable series is too short for the observed data.")
  }
  
  # Determine the number of regimes and transition parameters for the model
  K <- length(mu)
  n_transition <- K*(K-1)
  
  # Ensure we have the right number of transition probabilities and A parameters
  if (length(init_trans) != n_transition) {
    stop(sprintf("Error: Expected %d transition probabilities for %d regimes, got %d.", 
                 n_transition, K, length(init_trans)))
  }
  if (length(A) != n_transition) {
    stop(sprintf("Error: Expected %d A parameters for %d regimes, got %d.", 
                 n_transition, K, length(A)))
  }
  
  # Initialize variables
  eta <- matrix(0, nrow=K, ncol=M)     # Likelihood of each regime
  tot_lik <- numeric(M)                # Total likelihood
  X_t <- matrix(0, nrow=K, ncol=M)     # Filtered probabilities after observation
  X_tlag <- matrix(0, nrow=K, ncol=M)  # Predicted probabilities before observation
  
  # Get baseline transition probabilities that are logit-transformed
  omega_LR <- logit(init_trans)
  omega <- omega_LR * (1 - A)  # This scaling ensures that when A=0, we get back to init_trans
  
  # Initialize f values for transition probabilities
  f <- matrix(0, nrow=n_transition, ncol=M)
  p_trans <- matrix(0, nrow=n_transition, ncol=M)
  
  # Set initial values
  f[,1] <- omega + A * X_Exo[1]
  p_trans_raw <- logistic(f[,1])
  p_trans[,1] <- convert_to_valid_probs(p_trans_raw, K)
  
  # Initial state probabilities (using stationary distribution)
  X_t[,1] <- stat_dist(p_trans[,1], fallback_value = rep(1/K, K))
  
  for (t in 1:(M-1)) {
    # Generate predicted probabilities
    Pmatrix <- transition_matrix(p_trans[,t], check_validity = FALSE)
    X_tlag[,t] <- Pmatrix %*% X_t[,t]
    
    # Calculate likelihoods
    for (k in 1:K) {
      eta[k,t] <- dnorm(y[t], mu[k], sqrt(sigma2[k]))
    }
    tot_lik[t] <- sum(eta[,t]*X_tlag[,t])
    
    # Protect against numerical issues
    if (tot_lik[t] <= 0 || is.na(tot_lik[t])) {
      tot_lik[t] <- .Machine$double.eps
      X_t[,t+1] <- rep(1/K, K)  # Reset to uniform when we get an invalid likelihood
    } else {
      # Calculate filtered probabilities
      X_t[,t+1] <- (eta[,t]*X_tlag[,t])/tot_lik[t]
    }
    
    # Update transition probabilities based on the exogenous variable
    f[,t+1] <- omega + A * X_Exo[t+1]
    p_trans_raw <- logistic(f[,t+1])
    p_trans[,t+1] <- convert_to_valid_probs(p_trans_raw, K)
  }
  
  # Calculate likelihood for the last time point
  Pmatrix <- transition_matrix(p_trans[,M-1], check_validity = FALSE)
  X_tlag[,M] <- Pmatrix %*% X_t[,M-1]
  for (k in 1:K) {
    eta[k,M] <- dnorm(y[M], mu[k], sqrt(sigma2[k]))
  }
  tot_lik[M] <- sum(eta[,M]*X_tlag[,M])
  
  # Protect against numerical issues
  if (tot_lik[M] <= 0 || is.na(tot_lik[M])) {
    tot_lik[M] <- .Machine$double.eps
  }
  
  # Sum log-likelihoods, but exclude burn-in and cut-off
  valid_indices <- (B+1):(M-C)
  if (length(valid_indices) <= 0) {
    stop("Error: No valid data points after applying burn-in and cut-off.")
  }
  
  logLikSum <- sum(log(tot_lik[valid_indices]))
  
  # Return negative sum of log-likelihoods (for minimizing)
  res <- -logLikSum
  
  # Store additional information as attributes
  attr(res, "X.t") <- t(X_t)
  attr(res, "X.tlag") <- t(X_tlag)
  attr(res, "p_trans") <- p_trans
  
  return(res)
}

#' Estimate parameters for the exogenous variable-driven regime-switching model
#'
#' @param y Observed time series
#' @param X_Exo Exogenous variable that drives transition probabilities
#' @param K Number of regimes
#' @param B Burn-in period to exclude from likelihood calculation
#' @param C Cut-off period to exclude from likelihood calculation
#' @param initial_params Initial parameter guesses (optional, only used when n_starts=1)
#' @param bounds Parameter bounds for optimization in full parameter format (optional)
#' @param n_starts Number of starting points for optimization (default: 1)
#' @param parallel Whether to run multiple starts in parallel (default: FALSE)
#' @param cores Number of cores to use for parallel processing (default: detectCores() - 1)
#' @param seed Random seed for starting point generation (optional)
#' @param verbose Whether to print progress information (default: TRUE)
#' @param equal_variances Whether to constrain all regime variances to be equal (default: FALSE)
#' @return List with estimated parameters and model diagnostics
#' @details
#' Estimates model parameters using maximum likelihood estimation.
#' When n_starts > 1, generates multiple diverse starting points and runs
#' optimization in parallel (if parallel=TRUE) to improve robustness against
#' local optima. Returns the result with the highest likelihood.
#' 
#' Uses the future package for cross-platform parallel processing that works
#' on Windows, macOS, and Linux. The output format is identical regardless 
#' of single or multi-start estimation.
#' 
#' When equal_variances=TRUE, all regime variances are constrained to be equal
#' by compressing the variance parameters during optimization and expanding them
#' back afterward. Parameters and bounds should always be provided in full format
#' (with separate variance for each regime); compression is handled internally.
#'
#' @examples
#' # Single start (current behavior)
#' X_exo <- rnorm(1000)
#' results <- estimate_exogenous_model(data, X_exo, K=3)
#' 
#' # Multi-start with parallel processing
#' results <- estimate_exogenous_model(data, X_exo, K=3, n_starts=10, parallel=TRUE)
#' 
#' # Multi-start sequential (for debugging)
#' results <- estimate_exogenous_model(data, X_exo, K=3, n_starts=5, parallel=FALSE)
#' 
#' # Constrain variances to be equal (useful for comparing with benchmarks)
#' results <- estimate_exogenous_model(data, X_exo, K=3, equal_variances=TRUE)
estimate_exogenous_model <- function(y, X_Exo, K = 3, B = 100, C = 50, 
                                     initial_params = NULL, bounds = NULL,
                                     n_starts = 1, parallel = FALSE, cores = NULL,
                                     seed = NULL, verbose = TRUE, equal_variances = FALSE) {
  
  # Set up cores (don't use more cores than starts)
  if (is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 1)
  }
  cores <- min(cores, n_starts)
  
  if (verbose) {
    cat("Estimating Exogenous model with", K, "regimes")
    if (equal_variances) {
      cat(" (equal variances)")
    }
    cat("\n")
    if (n_starts > 1) {
      cat("Using", n_starts, "starting points")
      if (parallel) {
        cat(" with", cores, "cores in parallel\n")
      } else {
        cat(" sequentially\n")
      }
    } else {
      cat("Using single starting point\n")
    }
    start_time <- Sys.time()
  }
  
  # Set up parallel processing plan
  if (parallel && n_starts > 1 && requireNamespace("future.apply", quietly = TRUE)) {
    future::plan(future::multisession, workers = cores)
    on.exit(future::plan(future::sequential), add = TRUE)  # Cleanup
    if (verbose) {
      cat("Parallel processing enabled with future package\n")
    }
  } else {
    future::plan(future::sequential)
    if (parallel && n_starts > 1) {
      warning("future.apply package not available, running sequentially")
    }
  }
  
  # Generate starting points
  if (n_starts == 1 && !is.null(initial_params)) {
    # Use user-provided starting point
    starting_points <- list(initial_params)
    if (verbose) {
      cat("Using provided initial parameters\n")
    }
  } else {
    # Generate diverse starting points
    if (verbose && n_starts > 1) {
      cat("Generating", n_starts, "diverse starting points...\n")
    }
    starting_points <- generate_starting_points(y, K, "exogenous", n_starts, seed)
  }
  
  # Create default bounds if none provided
  if (is.null(bounds)) {
    n_transition <- K * (K - 1)
    
    lower_bounds <- c(rep(-Inf, K),               # No bounds on means
                      rep(-Inf, K),               # Variance >= 0 (log-transformed)
                      rep(-Inf, n_transition),    # Probabilities >= 0 (logit-transformed)
                      rep(-1, n_transition))      # A-coefficients bounded
    
    upper_bounds <- c(rep(Inf, K),                # No bounds on means
                      rep(Inf, K),                # Variance unbounded
                      rep(Inf, n_transition),     # Probabilities <= 1 (logit-transformed)
                      rep(1, n_transition))       # A-coefficients bounded
    
    bounds <- list(lower = lower_bounds, upper = upper_bounds)
  }
  
  # Compress bounds if using equal variances (applies to both user-provided and default bounds)
  if (equal_variances) {
    bounds$lower <- compress_variances(bounds$lower, K)
    bounds$upper <- compress_variances(bounds$upper, K)
  }
  
  # Define the optimization function for a single start
  optimize_single_start <- function(start_idx) {
    start_params <- starting_points[[start_idx]]
    
    if (verbose && !parallel) {
      cat("  Starting point", start_idx, "of", n_starts, "\n")
    }
    
    tryCatch({
      # 1. Transform parameters to unconstrained space for optimization
      transformed_params <- transform_parameters(start_params, "exogenous")
      
      # 2. Compress parameters (conditional)
      if (equal_variances) {
        transformed_params <- compress_variances(transformed_params, K)
      }
      
      # 3. Run optimization (same way regardless of equal_variances)
      trace_setting <- if (verbose > 1) 1 else 0
      optimization_result <- nlminb(
        start = transformed_params,
        objective = function(par_t, X_Exo, y, B, C) {
          # Handle compression inside the objective function
          working_par_t <- par_t
          if (equal_variances) {
            working_par_t <- expand_variances(par_t, K)
          }
          
          # Transform parameters back to original parameter space
          par <- untransform_parameters(working_par_t, "exogenous")
          
          mu <- mean_from_par(par, "exogenous")
          sigma2 <- sigma2_from_par(par, "exogenous")
          init_trans <- transp_from_par(par, "exogenous")
          A <- A_from_par(par, "exogenous")
          
          # Calculate likelihood and return it
          l <- Rfiltering_TVPXExo(mu, sigma2, init_trans, A, X_Exo, y, B, C)
          return(l[1])
        },
        lower = bounds$lower,
        upper = bounds$upper,
        X_Exo = X_Exo,
        y = y,
        B = B,
        C = C,
        control = list(eval.max = 1e6, iter.max = 1e6, trace = trace_setting)
      )
      
      # 4. Decompress parameters (conditional)
      result_params <- optimization_result$par
      if (equal_variances) {
        result_params <- expand_variances(result_params, K)
      }
      
      # 5. Untransform parameters (always happens)
      estimated_params <- untransform_parameters(result_params, "exogenous")
      
      # Store the final parameters in the optimization result for consistency
      optimization_result$final_par <- estimated_params
      
      # Return result with metadata
      list(
        result = optimization_result,
        start_idx = start_idx,
        convergence = optimization_result$convergence,
        objective = optimization_result$objective,
        start_params = start_params
      )
    }, error = function(e) {
      # Return error information
      list(
        result = NULL,
        start_idx = start_idx,
        convergence = 999,  # Error code
        objective = Inf,
        error = e$message,
        start_params = start_params
      )
    })
  }
  
  # Run optimization(s) using future package
  if (verbose && n_starts > 1) {
    cat("Running optimizations...\n")
    opt_start_time <- Sys.time()
  }
  
  results <- future.apply::future_lapply(
    X = 1:n_starts,
    FUN = optimize_single_start,
    future.seed = seed
  )
  
  if (verbose && n_starts > 1) {
    opt_end_time <- Sys.time()
    cat("Optimizations completed in", 
        format(difftime(opt_end_time, opt_start_time), digits = 4), "\n")
  }
  
  # Find the best result
  valid_results <- results[sapply(results, function(x) !is.null(x$result))]
  
  if (length(valid_results) == 0) {
    stop("All starting points failed to converge")
  }
  
  # Select result with lowest objective value (highest likelihood)
  objectives <- sapply(valid_results, function(x) x$objective)
  best_idx <- which.min(objectives)
  best_result <- valid_results[[best_idx]]
  
  # Count convergence issues
  convergence_codes <- sapply(results, function(x) x$convergence)
  n_converged <- sum(convergence_codes == 0)
  n_failed <- sum(convergence_codes == 999)
  
  if (verbose && n_starts > 1) {
    cat("Best result from starting point", best_result$start_idx, "\n")
    cat("Convergence summary:", n_converged, "converged,", 
        n_starts - n_converged - n_failed, "non-convergent,", n_failed, "failed\n")
    cat("Best negative log-likelihood:", best_result$objective, "\n")
  }
  
  # Extract the best optimization result
  optimization_result <- best_result$result
  
  # Get final parameters (already processed through our pipeline above)
  estimated_params <- optimization_result$final_par
  
  # Extract different parameter components
  mu_est <- mean_from_par(estimated_params, "exogenous")
  sigma2_est <- sigma2_from_par(estimated_params, "exogenous")
  init_trans_est <- transp_from_par(estimated_params, "exogenous")
  A_est <- A_from_par(estimated_params, "exogenous")
  
  # Calculate model diagnostics
  # Note: use the original parameter count for AIC/BIC calculation
  num_params_original <- length(estimated_params)  # Full parameter vector
  num_data_points <- length(y) - B - C
  
  aic <- 2 * optimization_result$objective + 2 * num_params_original
  bic <- 2 * optimization_result$objective + num_params_original * log(num_data_points)
  
  # Calculate filtered probabilities
  full_likelihood <- Rfiltering_TVPXExo(
    mu_est, sigma2_est, init_trans_est, A_est, X_Exo, y, B, C
  )
  
  filtered_probs <- attr(full_likelihood, "X.t")
  transition_probs <- attr(full_likelihood, "p_trans")
  
  # Prepare results (same format as original function)
  results_list <- list(
    parameters = list(
      mu = mu_est,
      sigma2 = sigma2_est,
      init_trans = init_trans_est,
      A = A_est
    ),
    diagnostics = list(
      loglik = -optimization_result$objective,
      aic = aic,
      bic = bic,
      num_params = num_params_original,
      num_data_points = num_data_points
    ),
    optimization = optimization_result,
    filtered_probabilities = filtered_probs,
    transition_probabilities = transition_probs,
    model_info = list(
      type = "exogenous",
      K = K,
      B = B,
      C = C,
      equal_variances = equal_variances
    )
  )
  
  # Add multi-start specific information (only if multiple starts)
  if (n_starts > 1) {
    results_list$multistart_info <- list(
      n_starts = n_starts,
      parallel_used = parallel && n_starts > 1 && requireNamespace("future.apply", quietly = TRUE),
      cores_used = cores,
      best_start_idx = best_result$start_idx,
      n_converged = n_converged,
      n_failed = n_failed,
      all_objectives = sapply(results, function(x) x$objective),
      convergence_codes = convergence_codes
    )
  }
  
  if (verbose) {
    end_time <- Sys.time()
    cat("Total estimation time:", 
        format(difftime(end_time, start_time), digits = 4), "\n")
    cat("AIC:", aic, "BIC:", bic, "\n")
    cat("Estimated means:", round(mu_est, 4), "\n")
    cat("Estimated variances:", round(sigma2_est, 4), "\n")
    if (equal_variances) {
      cat("Note: All variances constrained to be equal\n")
    }
  }
  
  return(results_list)
}
