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
#' @param burn_in Number of burn-in observations to discard (default: 100)
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
dataexoCD <- function(M, N, mu, sigma2, init_trans, A, X_Exo, burn_in = 100) {
  # Ensure that means and variances have been provided for all regimes
  if (length(mu) != length(sigma2)) {
    stop("Error: Unequal number of means and variances. Mean and variance have to be supplied for each regime.")
  }
  
  # Determine the number of regimes and transition parameters for the model
  K <- length(mu)
  n_transition <- K*(K-1)
  
  # Check that X_Exo is long enough
  if (length(X_Exo) < N + burn_in) {
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
    f[,1] <- omega + A * X_Exo[1]
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
      
      # Update transition probabilities based on the exogenous variable
      f[,t+1] <- omega + A * X_Exo[t+1]
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

#' Wrapper for the likelihood calculator to be used for max-likelihood estimation
#'
#' @param par_t Transformed parameters in unconstrained space
#' @param X_Exo Exogenous process that drives transition probability changes
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
#' transformed_params <- transform_TVPXExo(c(mu, sigma2, init_trans, A))
#' result <- nlminb(transformed_params, Rfiltering.single.trasf_TVPXExo, 
#'                 X_Exo = X_Exo, y = y, B = 100, C = 50)
Rfiltering.single.trasf_TVPXExo <- function(par_t, X_Exo, y, B, C) {
  # Transform parameters back to original parameter space
  par <- untransform_TVPXExo(par_t)
  
  mu <- mean_from_par(par)
  sigma2 <- sigma2_from_par(par)
  init_trans <- transp_from_par(par)
  A <- A_from_par(par)
  
  # Calculate likelihood and return it
  l <- Rfiltering_TVPXExo(mu, sigma2, init_trans, A, X_Exo, y, B, C)
  return(l[1])
}

#' Estimate parameters for the exogenous variable-driven regime-switching model
#'
#' @param y Observed time series
#' @param X_Exo Exogenous variable that drives transition probabilities
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
#' exog <- rnorm(1000) 
#' results <- estimate_exo_model(data, exog, K=3)
estimate_exo_model <- function(y, X_Exo, K = 3, B = 100, C = 50, 
                               initial_params = NULL, bounds = NULL,
                               verbose = TRUE) {
  if (verbose) {
    cat("Estimating exogenous model with", K, "regimes\n")
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
    init_trans_guess <- rep(1/K, K*(K-1))
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
    lower_bounds <- c(rep(-Inf, K),       # No bounds on means
                      rep(-Inf, K),       # Variance >= 0 (log-transformed)
                      rep(-Inf, K*(K-1)), # Probabilities >= 0 (logit-transformed)
                      rep(-1, K*(K-1)))   # A-coefficients bounded
    
    upper_bounds <- c(rep(Inf, K),        # No bounds on means
                      rep(Inf, K),        # Variance unbounded
                      rep(Inf, K*(K-1)),  # Probabilities <= 1 (logit-transformed)
                      rep(1, K*(K-1)))    # A-coefficients bounded
    
    bounds <- list(lower = lower_bounds, upper = upper_bounds)
  }
  
  # Transform parameters to unconstrained space for optimization
  transformed_params <- transform_TVPXExo(initial_params)
  
  # Optimize parameters
  if (verbose) {
    cat("Starting optimization...\n")
    opt_start_time <- Sys.time()
  }
  
  optimization_result <- nlminb(
    start = transformed_params,
    objective = Rfiltering.single.trasf_TVPXExo,
    lower = bounds$lower,
    upper = bounds$upper,
    X_Exo = X_Exo,
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
  estimated_params <- untransform_TVPXExo(optimization_result$par)
  
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
  full_likelihood <- Rfiltering_TVPXExo(
    mu_est, sigma2_est, init_trans_est, A_est, X_Exo, y, B, C
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
      type = "exogenous",
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
