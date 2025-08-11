#' Constant Transition Probability Regime-Switching Models
#' 
#' This file implements a regime-switching model with constant transition probabilities
#' between regimes.
#'
#' Reference: Sendstad, Chronopoulos, & Li (2025) - The Value of Turning-Point 
#' Detection for Optimal Investment

# Load required helper functions
source("helpers/utility_functions.R")
source("helpers/transition_helpers.R")
source("helpers/parameter_transforms.R")

#' Generate data from a constant transition probability regime-switching model
#'
#' @param M Number of simulation runs to be performed
#' @param N Length of the simulation runs (discretized time)
#' @param mu Vector of true means corresponding to each regime
#' @param sigma2 Vector of true variances corresponding to each regime
#' @param trans_prob Transition probabilities (off-diagonal elements) for the latent process
#' @param burn_in Number of burn-in observations to discard (default: 100)
#' @return Matrix of simulated data with M rows and N columns
#' @details 
#' Simulates data from a regime switching model where transition probabilities
#' remain constant over time.
#'
#' @examples
#' # Generate data for a 3-regime model
#' mu_true <- c(-2, 1, 2)
#' sigma2_true <- c(0.02, 0.2, 0.6)
#' trans_prob_true <- rep(0.2, 6)  # Off-diagonal elements for a 3x3 matrix
#' data_sim <- dataConstCD(10, 1000, mu_true, sigma2_true, trans_prob_true)
dataConstCD <- function(M, N, mu, sigma2, trans_prob, burn_in = 100) {
  # Parameters:
  # M           Number of simulation runs to be performed
  # N           Length of the simulation runs (discretized time)
  # mu          Vector of true means corresponding to each regime
  # sigma2      Vector of true variances corresponding to each regime
  # trans_prob  Transition probabilities (off-diagonal elements)
  # burn_in     Number of burn-in observations to discard
  
  # Ensure that means and variances have been provided for all regimes
  if (length(mu) != length(sigma2)) {
    stop("Error: Unequal number of means and variances. Mean and variance have to be supplied for each regime.")
  }
  
  # Determine the number of regimes and transition parameters for the model
  K <- length(mu)
  n_transition <- K*(K-1)
  
  # Ensure we have the right number of transition probabilities
  if (length(trans_prob) != n_transition) {
    stop(sprintf("Error: Expected %d transition probabilities for %d regimes, got %d.", 
                 n_transition, K, length(trans_prob)))
  }
  
  # Validate transition probabilities and convert to a valid stochastic matrix
  valid_trans_prob <- convert_to_valid_probs(trans_prob, K)
  transition_matrix <- transition_matrix(valid_trans_prob)
  
  # Set up a matrix to save the output
  data <- matrix(0, M, N)
  
  for (i in 1:M) {
    # Initialize all data structures including burn-in period
    total_length <- N + burn_in
    eta <- matrix(0, nrow=K, ncol=total_length)     # Likelihood of each regime
    tot_lik <- numeric(total_length)                # Total likelihood
    X_t <- matrix(0, nrow=K, ncol=total_length)     # Filtered probabilities after observation
    X_tlag <- matrix(0, nrow=K, ncol=total_length)  # Predicted probabilities before observation
    S <- numeric(total_length)                      # Latent state
    y.sim <- numeric(total_length)                  # Random increments according to state
    
    # Initial state probabilities (stationary distribution or uniform)
    start_dist <- try(stat_dist(transition_matrix), silent = TRUE)
    if (inherits(start_dist, "try-error")) {
      X_t[,1] <- rep(1/K, K)  # Fallback to uniform distribution
    } else {
      X_t[,1] <- start_dist
    }
    
    for (t in 1:(total_length-1)) {
      # Generate predicted probabilities using the constant transition matrix
      X_tlag[,t] <- transition_matrix %*% X_t[,t]
      
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
    }
    
    # For the last time point
    X_tlag[,total_length] <- transition_matrix %*% X_t[,total_length-1]
    S[total_length] <- sample(1:K, 1, prob=X_tlag[,total_length])
    y.sim[total_length] <- rnorm(1, mu[S[total_length]], sqrt(sigma2[S[total_length]]))
    
    # Remove burn-in and save the simulation run in the data matrix
    data[i,] <- y.sim[(burn_in+1):total_length]
  }
  
  return(data)
}

#' Filter observed data through the constant transition probability regime-switching model
#'
#' @param mu Vector of means corresponding to each regime
#' @param sigma2 Vector of variances corresponding to each regime
#' @param trans_prob Transition probabilities (off-diagonal elements)
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
#' trans_prob <- rep(0.2, 6)
#' y <- rnorm(1000) 
#' loglik <- Rfiltering_Const(mu, sigma2, trans_prob, y, 100, 50)
Rfiltering_Const <- function(mu, sigma2, trans_prob, y, B, C) {
  # Parameters passed to the function are:
  # mu          Vector of true means corresponding to each regime
  # sigma2      Vector of true variances corresponding to each regime
  # trans_prob  Transition probabilities (off-diagonal elements)
  # y           Observed time series increments
  # B           Burn-in to be excluded at the beginning of the time series
  # C           Cut-off to be excluded at the end of the time series
  
  # Ensure that means and variances have been provided for all regimes
  if (length(mu) != length(sigma2)) {
    stop("Error: Unequal number of means and variances. Mean and variance have to be supplied for each regime.")
  }
  
  # Determine length of the time series
  M <- length(y)
  
  # Determine the number of regimes and transition parameters for the model
  K <- length(mu)
  n_transition <- K*(K-1)
  
  # Ensure we have the right number of transition probabilities
  if (length(trans_prob) != n_transition) {
    stop(sprintf("Error: Expected %d transition probabilities for %d regimes, got %d.", 
                 n_transition, K, length(trans_prob)))
  }
  
  # Validate transition probabilities and convert to a valid stochastic matrix
  valid_trans_prob <- convert_to_valid_probs(trans_prob, K)
  transition_mat <- transition_matrix(valid_trans_prob)
  
  # Initialize variables
  eta <- matrix(0, nrow=K, ncol=M)     # Likelihood of each regime
  tot_lik <- numeric(M)                # Total likelihood
  X_t <- matrix(0, nrow=K, ncol=M)     # Filtered probabilities after observation
  X_tlag <- matrix(0, nrow=K, ncol=M)  # Predicted probabilities before observation
  
  # Initialize state probabilities (stationary distribution or uniform)
  start_dist <- try(stat_dist(transition_mat), silent = TRUE)
  if (inherits(start_dist, "try-error")) {
    X_t[,1] <- rep(1/K, K)  # Fallback to uniform distribution
  } else {
    X_t[,1] <- start_dist
  }
  
  for (t in 1:(M-1)) {
    # Generate predicted probabilities using the constant transition matrix
    X_tlag[,t] <- transition_mat %*% X_t[,t]
    
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
  }
  
  # Calculate likelihood for the last time point
  X_tlag[,M] <- transition_mat %*% X_t[,M-1]
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
  attr(res, "transition_matrix") <- transition_mat
  
  return(res)
}

#' Wrapper for the likelihood calculator to be used for max-likelihood estimation
#'
#' @param par_t Transformed parameters in unconstrained space
#' @param y Observed time series increments
#' @param B Burn-in to be excluded at the beginning of the time series
#' @param C Cut-off to be excluded at the end of the time series
#' @return Negative log-likelihood of observed data under the model
#' @details
#' Transforms parameters from unconstrained space back to natural space
#' and calls the filtering function to calculate the likelihood.
#'
#' @examples
#' # Optimize parameters for a 3-regime model
#' transformed_params <- transform_parameters(c(mu, sigma2, trans_prob), "constant")
#' result <- nlminb(transformed_params, Rfiltering.single.trasf_Const, 
#'                 y = y, B = 100, C = 50)
Rfiltering.single.trasf_Const <- function(par_t, y, B, C) {
  # Transform parameters back to original parameter space
  par <- untransform_parameters(par_t, "constant")
  
  K <- count_regime(par_t, "constant")
  
  mu <- mean_from_par(par, "constant")
  sigma2 <- sigma2_from_par(par, "constant")
  trans_prob <- transp_from_par(par, "constant")
  
  # Calculate likelihood and return it
  l <- Rfiltering_Const(mu, sigma2, trans_prob, y, B, C)
  return(l[1])
}

#' Estimate parameters for the constant transition probability regime-switching model
#'
#' @param y Observed time series
#' @param K Number of regimes
#' @param B Burn-in period to exclude from likelihood calculation
#' @param C Cut-off period to exclude from likelihood calculation
#' @param initial_params Initial parameter guesses (optional, only used when n_starts=1)
#' @param bounds Parameter bounds for optimization (optional)
#' @param n_starts Number of starting points for optimization (default: 1)
#' @param parallel Whether to run multiple starts in parallel (default: FALSE)
#' @param cores Number of cores to use for parallel processing (default: detectCores() - 1)
#' @param seed Random seed for starting point generation (optional)
#' @param verbose Whether to print progress information (default: TRUE)
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
#' @examples
#' # Single start (current behavior)
#' results <- estimate_const_model(data, K=3)
#' 
#' # Multi-start with parallel processing
#' results <- estimate_const_model(data, K=3, n_starts=10, parallel=TRUE)
#' 
#' # Multi-start sequential (for debugging)
#' results <- estimate_const_model(data, K=3, n_starts=5, parallel=FALSE)
estimate_const_model <- function(y, K = 3, B = 100, C = 50, 
                                 initial_params = NULL, bounds = NULL,
                                 n_starts = 1, parallel = FALSE, cores = NULL,
                                 seed = NULL, verbose = TRUE) {
  
  # Set up cores (don't use more cores than starts)
  if (is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 1)
  }
  cores <- min(cores, n_starts)
  
  if (verbose) {
    cat("Estimating constant transition probability model with", K, "regimes\n")
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
    starting_points <- generate_starting_points(y, K, "constant", n_starts, seed)
  }
  
  # Create default bounds if none provided
  if (is.null(bounds)) {
    lower_bounds <- c(rep(-Inf, K),           # No bounds on means
                      rep(-Inf, K),            # Variance >= 0 (log-transformed)
                      rep(-Inf, K*(K-1)))      # Probabilities >= 0 (logit-transformed)
    
    upper_bounds <- c(rep(Inf, K),            # No bounds on means
                      rep(Inf, K),             # Variance unbounded
                      rep(Inf, K*(K-1)))       # Probabilities <= 1 (logit-transformed)
    
    bounds <- list(lower = lower_bounds, upper = upper_bounds)
  }
  
  # Define the optimization function for a single start
  optimize_single_start <- function(start_idx) {
    start_params <- starting_points[[start_idx]]
    
    if (verbose && !parallel) {
      cat("  Starting point", start_idx, "of", n_starts, "\n")
    }
    
    tryCatch({
      # Transform parameters to unconstrained space for optimization
      transformed_params <- transform_parameters(start_params, "constant")
      
      # Run optimization
      trace_setting <- if (verbose > 1) 1 else 0
      optimization_result <- nlminb(
        start = transformed_params,
        objective = Rfiltering.single.trasf_Const,
        lower = bounds$lower,
        upper = bounds$upper,
        y = y,
        B = B,
        C = C,
        control = list(eval.max = 1e6, iter.max = 1e6, trace = trace_setting)
      )
      
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
  
  # Transform parameters back to natural space
  estimated_params <- untransform_parameters(optimization_result$par, "constant")
  
  # Extract different parameter components
  mu_est <- mean_from_par(estimated_params, "constant")
  sigma2_est <- sigma2_from_par(estimated_params, "constant")
  trans_prob_est <- transp_from_par(estimated_params, "constant")
  
  # Calculate model diagnostics
  num_params <- length(optimization_result$par)
  num_data_points <- length(y) - B - C
  
  aic <- 2 * optimization_result$objective + 2 * num_params
  bic <- 2 * optimization_result$objective + num_params * log(num_data_points)
  
  # Calculate filtered probabilities
  full_likelihood <- Rfiltering_Const(
    mu_est, sigma2_est, trans_prob_est, y, B, C
  )
  
  filtered_probs <- attr(full_likelihood, "X.t")
  transition_mat <- attr(full_likelihood, "transition_matrix")
  
  # Prepare results (same format as original function)
  results_list <- list(
    parameters = list(
      mu = mu_est,
      sigma2 = sigma2_est,
      trans_prob = trans_prob_est,
      transition_matrix = transition_mat
    ),
    diagnostics = list(
      loglik = -optimization_result$objective,
      aic = aic,
      bic = bic,
      num_params = num_params,
      num_data_points = num_data_points
    ),
    optimization = optimization_result,
    filtered_probabilities = filtered_probs,
    model_info = list(
      type = "Constant",
      K = K,
      B = B,
      C = C
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
  }
  
  return(results_list)
}
