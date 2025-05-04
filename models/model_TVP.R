#' Time-Varying Transition Probability Models with Autoregressive Dynamics
#' 
#' This file implements a regime-switching model where transition probabilities
#' depend on the process itself (autoregressive dynamics).
#'
#' Reference: Sendstad, Chronopoulos, & Li (2025) - The Value of Turning-Point 
#' Detection for Optimal Investment

# Load required helper functions
source("helpers/utility_functions.R")
source("helpers/transition_helpers.R")
source("helpers/parameter_transforms.R")

#' Generate data from an autoregressive regime-switching model
#'
#' @param M Number of simulation runs to be performed
#' @param N Length of the simulation runs (discretized time)
#' @param mu Vector of true means corresponding to each regime
#' @param sigma2 Vector of true variances corresponding to each regime
#' @param init_trans Initial transition probabilities for the latent process
#' @param A Autoregressive factor weights, one for each transition probability
#' @param burn_in Number of burn-in observations to discard (default: 100)
#' @return Matrix of simulated data with M rows and N columns
#' @details 
#' Simulates data from a regime switching model where transition probabilities
#' depend on previous values of the process itself. The relationship is defined as:
#' f[t+1] = omega + A * y[t], followed by logistic transformation.
#'
#' @examples
#' # Generate data for a 3-regime model
#' mu_true <- c(-2, 1, 2)
#' sigma2_true <- c(0.02, 0.2, 0.6)
#' init_trans_true <- rep(0.2, 6)  # Off-diagonal elements for a 3x3 matrix
#' A_true <- c(0.1, -0.1, 0.05, -0.05, 0.2, -0.2)
#' data_sim <- dataTVPCD(10, 1000, mu_true, sigma2_true, init_trans_true, A_true)
dataTVPCD <- function(M, N, mu, sigma2, init_trans, A, burn_in = 100) {
  # Parameters:
  # M           Number of simulation runs to be performed
  # N           Length of the simulation runs (discretized time)
  # mu          Vector of true means corresponding to each regime
  # sigma2      Vector of true variances corresponding to each regime
  # init_trans  Initial transition probabilities for the latent process
  # A           Autoregressive factor weights, one for each transition probability
  # burn_in     Number of burn-in observations to discard
  
  # Ensure that means and variances have been provided for all regimes
  if (length(mu) != length(sigma2)) {
    stop("Error: Unequal number of means and variances. Mean and variance have to be supplied for each regime.")
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
  
  # Set up a matrix to save the output
  data <- matrix(0, M, N)
  
  # Get baseline transition probabilities that are logit-transformed
  omega_LR <- logit(init_trans)
  omega <- omega_LR * (1 - A)  # This scaling ensures that when A=0, we get back to init_trans
  
  for (i in 1:M) {
    # Initialize all data structures including burn-in period
    total_length <- N + burn_in
    eta <- matrix(0, nrow=K, ncol=total_length)     # Likelihood of each regime
    tot_lik <- numeric(total_length)                # Total likelihood
    X_t <- matrix(0, nrow=K, ncol=total_length)     # Filtered probabilities after observation
    X_tlag <- matrix(0, nrow=K, ncol=total_length)  # Predicted probabilities before observation
    S <- numeric(total_length)                      # Latent state
    y.sim <- numeric(total_length)                  # Random increments according to state
    
    # Initialize f values for transition probabilities and convert to valid probabilities
    f <- matrix(0, nrow=n_transition, ncol=total_length)
    p_trans <- matrix(0, nrow=n_transition, ncol=total_length)
    
    # Initial state probabilities (uniform distribution)
    X_t[,1] <- rep(1/K, K)
    
    # Set initial values
    f[,1] <- omega
    p_trans_raw <- logistic(f[,1])
    p_trans[,1] <- convert_to_valid_probs(p_trans_raw, K)
    
    for (t in 1:(total_length-1)) {
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
      
      # Update transition probabilities based on the observed y.sim[t]
      # This is the key difference from the exogenous model
      f[,t+1] <- omega + A * y.sim[t]
      p_trans_raw <- logistic(f[,t+1])
      p_trans[,t+1] <- convert_to_valid_probs(p_trans_raw, K)
    }
    
    # For the last time point
    Pmatrix <- transition_matrix(p_trans[,total_length-1], check_validity = FALSE)
    X_tlag[,total_length] <- Pmatrix %*% X_t[,total_length-1]
    S[total_length] <- sample(1:K, 1, prob=X_tlag[,total_length])
    y.sim[total_length] <- rnorm(1, mu[S[total_length]], sqrt(sigma2[S[total_length]]))
    
    # Remove burn-in and save the simulation run in the data matrix
    data[i,] <- y.sim[(burn_in+1):total_length]
  }
  
  return(data)
}

#' Filter observed data through the autoregressive regime-switching model
#'
#' @param mu Vector of means corresponding to each regime
#' @param sigma2 Vector of variances corresponding to each regime
#' @param init_trans Initial transition probabilities for the latent process
#' @param A Autoregressive factor weights, one for each transition probability
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
#' y <- rnorm(1000) 
#' loglik <- Rfiltering_TVP(mu, sigma2, init_trans, A, y, 100, 50)
Rfiltering_TVP <- function(mu, sigma2, init_trans, A, y, B, C) {
  # Parameters passed to the function are:
  # mu          Vector of true means corresponding to each regime
  # sigma2      Vector of true variances corresponding to each regime
  # init_trans  Initial transition probabilities for the latent process
  # A           Autoregressive factor weights, one for each transition probability
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
  f[,1] <- omega
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
    
    # Update transition probabilities based on the observed y[t]
    # Note: we're using the actual observations now, not simulated values
    f[,t+1] <- omega + A * y[t]
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
#' transformed_params <- transform_TVP(c(mu, sigma2, init_trans, A))
#' result <- nlminb(transformed_params, Rfiltering.single.trasf_TVP, 
#'                 y = y, B = 100, C = 50)
Rfiltering.single.trasf_TVP <- function(par_t, y, B, C) {
  # Transform parameters back to original parameter space
  par <- untransform_TVP(par_t)
  
  mu <- mean_from_par(par)
  sigma2 <- sigma2_from_par(par)
  init_trans <- transp_from_par(par)
  A <- A_from_par(par)
  
  # Calculate likelihood and return it
  l <- Rfiltering_TVP(mu, sigma2, init_trans, A, y, B, C)
  return(l[1])
}

#' Estimate parameters for the autoregressive regime-switching model
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
#' results <- estimate_tvp_model(data, K=3)
estimate_tvp_model <- function(y, K = 3, B = 100, C = 50, 
                               initial_params = NULL, bounds = NULL,
                               verbose = TRUE) {
  if (verbose) {
    cat("Estimating TVP model with", K, "regimes\n")
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
    
    # Initialize transition probabilities and A coefficients
    init_trans_guess <- rep(0.2, K*(K-1))
    A_guess <- rep(0, K*(K-1))
    
    initial_params <- c(mu_guess, sigma2_guess, init_trans_guess, A_guess)
    
    if (verbose) {
      cat("Generated initial parameter guesses:\n")
      cat("Means:", round(mu_guess, 4), "\n")
      cat("Variances:", round(sigma2_guess, 4), "\n")
    }
  }
  
  # Create default bounds if none provided
  if (is.null(bounds)) {
    D <- 1e-3  # Small value to avoid boundary issues
    lower_bounds <- c(rep(-Inf, K),           # No bounds on means
                      rep(-Inf, K),            # Variance >= 0 (log-transformed)
                      rep(-Inf, K*(K-1)),      # Probabilities >= 0 (logit-transformed)
                      rep(-1, K*(K-1)))        # A-coefficients bounded
    
    upper_bounds <- c(rep(Inf, K),            # No bounds on means
                      rep(Inf, K),             # Variance unbounded
                      rep(Inf, K*(K-1)),       # Probabilities <= 1 (logit-transformed)
                      rep(1, K*(K-1)))         # A-coefficients bounded
    
    bounds <- list(lower = lower_bounds, upper = upper_bounds)
  }
  
  # Transform parameters to unconstrained space for optimization
  transformed_params <- transform_TVP(initial_params)
  
  # Optimize parameters
  if (verbose) {
    cat("Starting optimization...\n")
    opt_start_time <- Sys.time()
  }
  
  optimization_result <- nlminb(
    start = transformed_params,
    objective = Rfiltering.single.trasf_TVP,
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
  estimated_params <- untransform_TVP(optimization_result$par)
  
  # Extract different parameter components
  mu_est <- mean_from_par(estimated_params)
  sigma2_est <- sigma2_from_par(estimated_params)
  init_trans_est <- transp_from_par(estimated_params)
  A_est <- A_from_par(estimated_params)
  
  # Calculate model diagnostics
  num_params <- length(optimization_result$par)
  num_data_points <- length(y) - B - C
  
  aic <- 2 * optimization_result$objective + 2 * num_params
  bic <- 2 * optimization_result$objective + num_params * log(num_data_points)
  
  # Calculate filtered probabilities
  full_likelihood <- Rfiltering_TVP(
    mu_est, sigma2_est, init_trans_est, A_est, y, B, C
  )
  
  filtered_probs <- attr(full_likelihood, "X.t")
  transition_probs <- attr(full_likelihood, "p_trans")
  
  # Prepare results
  results <- list(
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
      num_params = num_params,
      num_data_points = num_data_points
    ),
    optimization = optimization_result,
    filtered_probabilities = filtered_probs,
    transition_probabilities = transition_probs,
    model_info = list(
      type = "TVP",
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

#' Simulate and estimate a TVP model with specified parameters
#'
#' @param M Number of simulation paths
#' @param N Length of each simulation path
#' @param mu True regime means
#' @param sigma2 True regime variances
#' @param init_trans True initial transition probabilities
#' @param A True autoregressive coefficients
#' @param B Burn-in period
#' @param C Cut-off period
#' @param verbose Whether to print progress information
#' @return List with simulation results and estimation results
#' @details
#' This function simulates data from a TVP model and then estimates 
#' the model parameters, providing a way to test the estimation procedure
#' with known parameters.
#'
#' @examples
#' # Test estimation with 10 simulated paths
#' mu_true <- c(-2, 1, 2)
#' sigma2_true <- c(0.02, 0.2, 0.6)
#' init_trans_true <- rep(0.2, 6)
#' A_true <- c(0.1, -0.1, 0.05, -0.05, 0.2, -0.2)
#' results <- example_TVP_simulation(10, 1000, mu_true, sigma2_true, init_trans_true, A_true)
example_TVP_simulation <- function(M = 10, N = 1000, 
                                   mu = c(-2, 1, 2), 
                                   sigma2 = c(0.02, 0.2, 0.6), 
                                   init_trans = rep(0.2, 6), 
                                   A = c(0.1, -0.1, 0.05, -0.05, 0.2, -0.2),
                                   B = 100, C = 50, verbose = TRUE) {
  # Parameters
  K <- length(mu)    # Number of regimes
  
  if (verbose) {
    cat("Setting up parameters for TVP simulation...\n")
    cat("Regimes:", K, "\n")
    cat("Means:", mu, "\n")
    cat("Variances:", sigma2, "\n")
  }
  
  # Generate simulation data
  if (verbose) {
    cat("Generating simulation data for", M, "paths, each with length", N, "...\n")
    start_time <- Sys.time()
  }
  
  data_TVP_sim <- dataTVPCD(M, N, mu, sigma2, init_trans, A)
  
  if (verbose) {
    end_time <- Sys.time()
    cat("Data generation completed in", 
        format(difftime(end_time, start_time), digits = 4), "\n")
  }
  
  # Storage for estimation results
  param_estimates <- matrix(0, M, 2*K^2)
  colnames(param_estimates) <- c(paste0("mu", 1:K), 
                                 paste0("sigma2", 1:K), 
                                 paste0("trans", 1:(K*(K-1))), 
                                 paste0("A", 1:(K*(K-1))))
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
    current_data <- data_TVP_sim[i,]
    
    # Initial parameter guesses for optimization
    mu_guess <- mu * 0.9 + rnorm(K, 0, 0.1*abs(mu))  # Add noise to true values
    sigma2_guess <- sigma2 * 0.9 + rnorm(K, 0, 0.1*sigma2)  # Add noise to true values
    init_trans_guess <- pmax(pmin(init_trans + rnorm(K*(K-1), 0, 0.05), 0.95), 0.05)  # Add noise, keep in (0,1)
    A_guess <- rep(0, K*(K-1))  # Start with no effect
    
    par_guess <- c(mu_guess, sigma2_guess, init_trans_guess, A_guess)
    
    # Set parameter bounds
    D <- 1e-3  # Small value to avoid boundary issues
    lower_bounds <- c(rep(-Inf, K),           # No bounds on means
                      rep(-Inf, K),           # Variance >= 0 (log-transformed)
                      rep(-Inf, K*(K-1)),     # Probabilities >= 0 (logit-transformed)
                      rep(-1, K*(K-1)))       # A-coefficients bounded
    
    upper_bounds <- c(rep(Inf, K),            # No bounds on means
                      rep(Inf, K),            # Variance unbounded
                      rep(Inf, K*(K-1)),      # Probabilities <= 1 (logit-transformed)
                      rep(1, K*(K-1)))        # A-coefficients bounded
    
    # Estimate the model
    TVP_est <- try(nlminb(
      start = transform_TVP(par_guess),
      objective = Rfiltering.single.trasf_TVP, 
      lower = lower_bounds,
      upper = upper_bounds,
      y = current_data,
      B = B, 
      C = C,
      control = list(eval.max = 1e6, iter.max = 1e6, trace = 0)
    ), silent = TRUE)
    
    if (!inherits(TVP_est, "try-error")) {
      # Transform parameters back to natural space
      TVP_par_fin <- untransform_TVP(TVP_est$par)
      
      # Store the parameter estimates
      param_estimates[i,] <- TVP_par_fin
      
      # Store diagnostics
      num_params <- length(TVP_est$par)
      diagnostics[i, 1] <- -TVP_est$objective
      diagnostics[i, 2] <- 2 * TVP_est$objective + 2 * num_params
      diagnostics[i, 3] <- 2 * TVP_est$objective + num_params * log(N - B - C)
      
      if (verbose) {
        path_end_time <- Sys.time()
        path_time <- difftime(path_end_time, path_start_time, units = "secs")
        cat("completed in", round(path_time, 2), "seconds (", 
            formatC(TVP_est$objective, digits = 8, format = "f"), ")\n")
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
    True_Value = c(mu, sigma2, init_trans, A),
    Mean_Estimate = colMeans(param_estimates, na.rm = TRUE),
    SD_Estimate = apply(param_estimates, 2, sd, na.rm = TRUE),
    Bias = colMeans(param_estimates, na.rm = TRUE) - c(mu, sigma2, init_trans, A),
    Rel_Bias_Pct = 100 * (colMeans(param_estimates, na.rm = TRUE) - c(mu, sigma2, init_trans, A)) / 
      ifelse(abs(c(mu, sigma2, init_trans, A)) > 1e-10, 
             c(mu, sigma2, init_trans, A), 
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
      data = data_TVP_sim,
      parameters = list(
        mu = mu,
        sigma2 = sigma2,
        init_trans = init_trans,
        A = A
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
