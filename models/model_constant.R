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

#' Transform parameters from natural space to unconstrained space for constant model
#'
#' @param par Parameters in their natural space
#' @return Parameters in unconstrained space
#' @details
#' Applies log transformation to variances and logit to transition probabilities
#' to make parameters suitable for unconstrained optimization.
transform_Const <- function(par) {
  mu <- mean_from_par(par)
  sigma2 <- sigma2_from_par(par)
  trans_prob <- transp_from_par(par)
  
  # Check for invalid values
  if (any(sigma2 <= 0)) {
    stop("Variance parameters must be positive")
  }
  if (any(trans_prob <= 0) || any(trans_prob >= 1)) {
    stop("Transition probabilities must be between 0 and 1")
  }
  
  return(c(mu, log(sigma2), logit(trans_prob)))
}

#' Transform parameters from unconstrained space back to natural space for constant model
#'
#' @param par_t Parameters in unconstrained space
#' @return Parameters in their natural space
#' @details
#' Applies exp transformation to log-variances and logistic to logit-probabilities
#' to convert back to naturally bounded parameters.
untransform_Const <- function(par_t) {
  K <- count_regime_const(par_t)
  
  mu <- par_t[1:K]
  log_sigma2 <- par_t[(K+1):(2*K)]
  logit_trans_prob <- par_t[(2*K+1):length(par_t)]
  
  return(c(mu, exp(log_sigma2), logistic(logit_trans_prob)))
}

#' Determine the number of regimes from a parameter vector length for constant model
#'
#' @param par Parameter vector
#' @return Number of regimes (K)
count_regime_const <- function(par) {
  # For constant model, there are K means, K variances, and K(K-1) transition probabilities
  # So total length is 2K + K(K-1) = K(K+1)
  # We need to solve for K: K^2 + K - length(par) = 0
  
  # Quadratic formula: K = (-1 + sqrt(1 + 4*length(par)))/2
  K <- (-1 + sqrt(1 + 4*length(par)))/2
  
  if (abs(K - round(K)) > 1e-10) {
    stop("The parameter vector length is not compatible with any valid number of regimes.")
  }
  
  return(round(K))
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
#' transformed_params <- transform_Const(c(mu, sigma2, trans_prob))
#' result <- nlminb(transformed_params, Rfiltering.single.trasf_Const, 
#'                 y = y, B = 100, C = 50)
Rfiltering.single.trasf_Const <- function(par_t, y, B, C) {
  # Transform parameters back to original parameter space
  par <- untransform_Const(par_t)
  
  K <- count_regime_const(par_t)
  
  mu <- par[1:K]
  sigma2 <- par[(K+1):(2*K)]
  trans_prob <- par[(2*K+1):length(par)]
  
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
#' @param initial_params Initial parameter guesses (optional)
#' @param bounds Parameter bounds for optimization (optional)
#' @param verbose Whether to print progress information (default: TRUE)
#' @return List with estimated parameters and model diagnostics
#' @details
#' Estimates model parameters using maximum likelihood estimation.
#' Returns the estimated parameters and various diagnostics including
#' AIC, BIC, filtered probabilities, and optimization details.
#'
#' @examples
#' # Estimate a 3-regime model
#' data <- rnorm(1000)
#' results <- estimate_const_model(data, K=3)
estimate_const_model <- function(y, K = 3, B = 100, C = 50, 
                                 initial_params = NULL, bounds = NULL,
                                 verbose = TRUE) {
  if (verbose) {
    cat("Estimating constant transition probability model with", K, "regimes\n")
    start_time <- Sys.time()
  }
  
  # Create default initial parameters if none provided
  if (is.null(initial_params)) {
    # Create sensible initial guesses based on data characteristics
    y_mean <- mean(y, na.rm = TRUE)
    y_sd <- sd(y, na.rm = TRUE)
    
    # Create regime means that span the range of the data
    spread <- 2 * y_sd
    mu_guess <- seq(y_mean - spread, y_mean + spread, length.out = K)
    
    # Create regime variances that increase with regime number
    sigma2_guess <- seq(y_sd^2 * 0.5, y_sd^2 * 1.5, length.out = K)
    
    # Initialize transition probabilities
    trans_prob_guess <- rep(0.2, K*(K-1))
    
    initial_params <- c(mu_guess, sigma2_guess, trans_prob_guess)
    
    if (verbose) {
      cat("Generated initial parameter guesses:\n")
      cat("Means:", round(mu_guess, 4), "\n")
      cat("Variances:", round(sigma2_guess, 4), "\n")
    }
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
  
  # Transform parameters to unconstrained space for optimization
  transformed_params <- transform_Const(initial_params)
  
  # Optimize parameters
  if (verbose) {
    cat("Starting optimization...\n")
    opt_start_time <- Sys.time()
  }
  
  optimization_result <- nlminb(
    start = transformed_params,
    objective = Rfiltering.single.trasf_Const,
    lower = bounds$lower,
    upper = bounds$upper,
    y = y,
    B = B,
    C = C,
    control = list(eval.max = 1e6, iter.max = 1e6, trace = ifelse(verbose, 1, 0))
  )
  
  if (verbose) {
    opt_end_time <- Sys.time()
    cat("Optimization completed in", 
        format(difftime(opt_end_time, opt_start_time), digits = 4), "\n")
    cat("Optimization convergence code:", optimization_result$convergence, "\n")
    cat("Final negative log-likelihood:", optimization_result$objective, "\n")
  }
  
  # Transform parameters back to natural space
  estimated_params <- untransform_Const(optimization_result$par)
  
  # Extract different parameter components
  mu_est <- estimated_params[1:K]
  sigma2_est <- estimated_params[(K+1):(2*K)]
  trans_prob_est <- estimated_params[(2*K+1):length(estimated_params)]
  
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
  
  # Prepare results
  results <- list(
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
  
  if (verbose) {
    end_time <- Sys.time()
    cat("Total estimation time:", 
        format(difftime(end_time, start_time), digits = 4), "\n")
    cat("AIC:", aic, "BIC:", bic, "\n")
    cat("Estimated means:", round(mu_est, 4), "\n")
    cat("Estimated variances:", round(sigma2_est, 4), "\n")
  }
  
  return(results)
}

#' Simulate and estimate a constant transition probability model with specified parameters
#'
#' @param M Number of simulation paths
#' @param N Length of each simulation path
#' @param mu True regime means
#' @param sigma2 True regime variances
#' @param trans_prob True transition probabilities
#' @param B Burn-in period
#' @param C Cut-off period
#' @param verbose Whether to print progress information
#' @return List with simulation results and estimation results
#' @details
#' This function simulates data from a constant transition probability model and then estimates 
#' the model parameters, providing a way to test the estimation procedure with known parameters.
#'
#' @examples
#' # Test estimation with 10 simulated paths
#' mu_true <- c(-2, 1, 2)
#' sigma2_true <- c(0.02, 0.2, 0.6)
#' trans_prob_true <- rep(0.2, 6)
#' results <- example_Const_simulation(10, 1000, mu_true, sigma2_true, trans_prob_true)
example_Const_simulation <- function(M = 10, N = 1000, 
                                     mu = c(-2, 1, 2), 
                                     sigma2 = c(0.02, 0.2, 0.6), 
                                     trans_prob = rep(0.2, 6),
                                     B = 100, C = 50, verbose = TRUE) {
  # Parameters
  K <- length(mu)    # Number of regimes
  
  if (verbose) {
    cat("Setting up parameters for constant transition probability simulation...\n")
    cat("Regimes:", K, "\n")
    cat("Means:", mu, "\n")
    cat("Variances:", sigma2, "\n")
  }
  
  # Generate simulation data
  if (verbose) {
    cat("Generating simulation data for", M, "paths, each with length", N, "...\n")
    start_time <- Sys.time()
  }
  
  data_Const_sim <- dataConstCD(M, N, mu, sigma2, trans_prob)
  
  if (verbose) {
    end_time <- Sys.time()
    cat("Data generation completed in", 
        format(difftime(end_time, start_time), digits = 4), "\n")
  }
  
  # Storage for estimation results
  param_estimates <- matrix(0, M, K*(K+1))
  colnames(param_estimates) <- c(paste0("mu", 1:K), 
                                 paste0("sigma2", 1:K), 
                                 paste0("trans", 1:(K*(K-1))))
  diagnostics <- matrix(0, M, 3)
  colnames(diagnostics) <- c("loglik", "aic", "bic")
  
  # Estimate parameters for each simulation path
  if (verbose) {
    cat("Starting parameter estimation for", M, "simulation paths...\n")
    total_start_time <- Sys.time()
  }
  
  for (i in 1:M) {
    if (verbose) {
      path_start_time <- Sys.time()
      cat("Processing path", i, "of", M, "... ")
    }
    
    # Extract the current simulation path
    current_data <- data_Const_sim[i,]
    
    # Initial parameter guesses for optimization
    mu_guess <- mu * 0.9 + rnorm(K, 0, 0.1*abs(mu))  # Add noise to true values
    sigma2_guess <- sigma2 * 0.9 + rnorm(K, 0, 0.1*sigma2)  # Add noise to true values
    trans_prob_guess <- pmax(pmin(trans_prob + rnorm(K*(K-1), 0, 0.05), 0.95), 0.05)  # Add noise, keep in (0,1)
    
    par_guess <- c(mu_guess, sigma2_guess, trans_prob_guess)
    
    # Estimate the model
    Const_est <- try(nlminb(
      start = transform_Const(par_guess),
      objective = Rfiltering.single.trasf_Const, 
      lower = bounds$lower,
      upper = bounds$upper,
      y = current_data,
      B = B, 
      C = C,
      control = list(eval.max = 1e6, iter.max = 1e6, trace = 0)
    ), silent = TRUE)
    
    if (!inherits(Const_est, "try-error")) {
      # Transform parameters back to natural space
      Const_par_fin <- untransform_Const(Const_est$par)
      
      # Store the parameter estimates
      param_estimates[i,] <- Const_par_fin
      
      # Store diagnostics
      num_params <- length(Const_est$par)
      diagnostics[i, 1] <- -Const_est$objective
      diagnostics[i, 2] <- 2 * Const_est$objective + 2 * num_params
      diagnostics[i, 3] <- 2 * Const_est$objective + num_params * log(N - B - C)
      
      if (verbose) {
        path_end_time <- Sys.time()
        path_time <- difftime(path_end_time, path_start_time, units = "secs")
        cat("completed in", round(path_time, 2), "seconds (", 
            formatC(Const_est$objective, digits = 8, format = "f"), ")\n")
      }
    } else {
      # Handle estimation errors
      param_estimates[i,] <- NA
      diagnostics[i,] <- NA
      
      if (verbose) {
        cat("FAILED - Optimization error\n")
      }
    }
    
    # Estimate remaining time if verbose
    if (verbose && i > 1) {
      elapsed <- difftime(Sys.time(), total_start_time, units = "secs")
      time_per_path <- as.numeric(elapsed) / i
      remaining_paths <- M - i
      est_remaining <- time_per_path * remaining_paths
      
      cat("Overall progress:", i, "/", M, "paths -", 
          round(i/M*100), "% complete. Est. remaining time:", 
          round(est_remaining/60, 1), "minutes\n")
    }
  }
  
  if (verbose) {
    total_end_time <- Sys.time()
    total_time <- difftime(total_end_time, total_start_time, units = "mins")
    cat("\nParameter estimation completed in", 
        format(total_time, digits = 4), "\n")
  }
  
  # Calculate summary statistics for parameter estimates
  estimates_summary <- data.frame(
    Parameter = colnames(param_estimates),
    True_Value = c(mu, sigma2, trans_prob),
    Mean_Estimate = colMeans(param_estimates, na.rm = TRUE),
    SD_Estimate = apply(param_estimates, 2, sd, na.rm = TRUE),
    Bias = colMeans(param_estimates, na.rm = TRUE) - c(mu, sigma2, trans_prob),
    Rel_Bias_Pct = 100 * (colMeans(param_estimates, na.rm = TRUE) - c(mu, sigma2, trans_prob)) / 
      ifelse(abs(c(mu, sigma2, trans_prob)) > 1e-10, 
             c(mu, sigma2, trans_prob), 
             1e-10)
  )
  
  # Compare estimated parameters with true parameters
  if (verbose) {
    cat("\nParameter Estimation Summary:\n")
    print(estimates_summary)
    
    cat("\nMean Log-Likelihood:", mean(diagnostics[,1], na.rm = TRUE), "\n")
    cat("Mean AIC:", mean(diagnostics[,2], na.rm = TRUE), "\n")
    cat("Mean BIC:", mean(diagnostics[,3], na.rm = TRUE), "\n")
  }
  
  # Prepare return values
  results <- list(
    simulation = list(
      data = data_Const_sim,
      parameters = list(
        mu = mu,
        sigma2 = sigma2,
        trans_prob = trans_prob
      )
    ),
    estimation = list(
      parameters = param_estimates,
      diagnostics = diagnostics,
      summary = estimates_summary
    ),
    settings = list(
      M = M,
      N = N,
      K = K,
      B = B,
      C = C
    )
  )
  
  return(results)
}
