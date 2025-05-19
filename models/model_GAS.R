#' Score-Driven (GAS) Time-Varying Transition Probability Models
#' 
#' This file implements a regime-switching model where transition probabilities
#' vary over time according to the score of the predictive likelihood function.
#'
#' Reference: Sendstad, Chronopoulos, & Li (2025) - The Value of Turning-Point 
#' Detection for Optimal Investment
#' Based on: Bazzi et al. (2017) - Time-Varying Transition Probabilities for Markov Regime
#' Switching Models

# Load required helper functions
source("helpers/utility_functions.R")
source("helpers/transition_helpers.R")
source("helpers/parameter_transforms.R")

#' Generate data from a score-driven regime-switching model (GAS)
#'
#' @param M Number of simulation runs to be performed
#' @param N Length of the simulation runs (discretized time)
#' @param mu Vector of true means corresponding to each regime
#' @param sigma2 Vector of true variances corresponding to each regime
#' @param init_trans Initial transition probabilities for the latent process
#' @param A Scale parameter for score updates (sensitivity)
#' @param B Persistence parameter for score updates (memory)
#' @param burn_in Number of burn-in observations to discard (default: 100)
#' @return Matrix of simulated data with M rows and N columns
#' @details 
#' Simulates data from a regime switching model where transition probabilities
#' are updated using score-driven dynamics based on the predictive likelihood.
#' The update equation is f[t+1] = omega + A*s[t] + B*(f[t] - omega).
#'
#' @examples
#' # Generate data for a 3-regime model
#' mu_true <- c(-2, 1, 2)
#' sigma2_true <- c(0.02, 0.2, 0.6)
#' init_trans_true <- rep(0.2, 6)  # Off-diagonal elements for a 3x3 matrix
#' A_true <- rep(0.1, 6)  # Score scaling parameters
#' B_true <- rep(0.9, 6)  # Persistence parameters
#' data_sim <- dataGASCD(10, 1000, mu_true, sigma2_true, init_trans_true, A_true, B_true)
dataGASCD <- function(M, N, mu, sigma2, init_trans, A, B, burn_in = 100) {
  # Parameters:
  # M           Number of simulation runs to be performed
  # N           Length of the simulation runs (discretized time)
  # mu          Vector of true means corresponding to each regime
  # sigma2      Vector of true variances corresponding to each regime
  # init_trans  Initial transition probabilities for the latent process
  # A           Scaling parameters for score updates (sensitivity)
  # B           Persistence parameters for score updates (memory)
  # burn_in     Number of burn-in observations to discard
  
  # Ensure that means and variances have been provided for all regimes
  if (length(mu) != length(sigma2)) {
    stop("Error: Unequal number of means and variances. Mean and variance have to be supplied for each regime.")
  }
  
  # Determine the number of regimes and transition parameters for the model
  K <- length(mu)
  n_transition <- K*(K-1)
  
  # Ensure we have the right number of parameters
  if (length(init_trans) != n_transition) {
    stop(sprintf("Error: Expected %d transition probabilities for %d regimes, got %d.", 
                 n_transition, K, length(init_trans)))
  }
  if (length(A) != n_transition) {
    stop(sprintf("Error: Expected %d A parameters for %d regimes, got %d.", 
                 n_transition, K, length(A)))
  }
  if (length(B) != n_transition) {
    stop(sprintf("Error: Expected %d B parameters for %d regimes, got %d.", 
                 n_transition, K, length(B)))
  }
  
  # Set up a matrix to save the output
  data <- matrix(0, M, N)
  
  # Get baseline transition probabilities that are logit-transformed
  omega <- logit(init_trans)
  
  # For Gauss-Hermite quadrature (used in scaling score updates)
  GH_points <- 30
  
  for (i in 1:M) {
    # Initialize all data structures including burn-in period
    total_length <- N + burn_in
    
    # Arrays for storing various quantities throughout the simulation
    eta <- matrix(0, nrow=K, ncol=total_length)      # Likelihood of each regime
    tot_lik <- numeric(total_length)                 # Total likelihood
    X_t <- matrix(0, nrow=K, ncol=total_length)      # Filtered probabilities after observation
    X_tlag <- matrix(0, nrow=K, ncol=total_length)   # Predicted probabilities before observation
    S <- numeric(total_length)                       # Latent state
    y.sim <- numeric(total_length)                   # Random increments according to state
    
    # Arrays for score-driven dynamics
    f <- matrix(0, nrow=n_transition, ncol=total_length)       # Time-varying parameters (logit scale)
    p_trans <- matrix(0, nrow=n_transition, ncol=total_length) # Transition probabilities
    score <- matrix(0, nrow=n_transition, ncol=total_length)   # Score vectors
    score_scaled <- matrix(0, nrow=n_transition, ncol=total_length) # Scaled score vectors
    info <- numeric(total_length)                              # Fisher information
    
    # Initial state probabilities (uniform distribution)
    X_t[,1] <- rep(1/K, K)
    
    # Set initial values for transition parameters
    f[,1] <- omega
    p_trans_raw <- logistic(f[,1])
    p_trans[,1] <- convert_to_valid_probs(p_trans_raw, K)
    
    # Initialize setup for Gauss-Hermite quadrature (for scaling scores)
    # This is used to compute expectations in the Fisher information matrix
    gh_setup <- NULL
    
    for (t in 1:(total_length-1)) {
      # Generate predicted probabilities using current transition probabilities
      Pmatrix <- transition_matrix(p_trans[,t], check_validity = FALSE)
      X_tlag[,t] <- Pmatrix %*% X_t[,t]
      
      # Sample a state based on the predicted probabilities and 
      # simulate data conditional on that state
      S[t] <- sample(1:K, 1, prob=X_tlag[,t])
      y.sim[t] <- rnorm(1, mu[S[t]], sqrt(sigma2[S[t]]))
      
      # Calculate likelihoods for each regime
      for (k in 1:K) {
        eta[k,t] <- dnorm(y.sim[t], mu[k], sqrt(sigma2[k]))
      }
      tot_lik[t] <- sum(eta[,t]*X_tlag[,t])
      
      # Calculate filtered probabilities
      X_t[,t+1] <- (eta[,t]*X_tlag[,t])/tot_lik[t]
      
      # Calculate the score vector (differences in likelihoods between regimes)
      # Formula from the paper: r_t = (p(y_t|θ_0) - p(y_t|θ_1)) / p(y_t|θ) * g(f_t, θ, I_{t-1})
      
      # Compute likelihood differences for regimes
      delta_like <- numeric(n_transition)
      idx <- 1
      for (r in 1:K) {
        for (c in 1:K) {
          if (r != c) {
            # This calculates density difference impact for each transition probability
            # Following Bazzi et al. (2017), equation (13)
            delta_like[idx] <- (eta[r,t] - eta[c,t]) / tot_lik[t]
            idx <- idx + 1
          }
        }
      }
      
      # Compute the g vector - how filtered probabilities affect transitions
      # Following Bazzi et al. (2017), equation (14)
      g_vector <- numeric(n_transition)
      idx <- 1
      for (r in 1:K) {
        for (c in 1:K) {
          if (r != c) {
            # Impact of filtered probability on transition probability
            g_vector[idx] <- X_t[r,t] * p_trans[idx,t] * (1 - p_trans[idx,t])
            idx <- idx + 1
          }
        }
      }
      
      # Unscaled score is product of likelihood differences and probability impacts
      score[,t] <- delta_like * g_vector
      
      # Compute Fisher information scaling for better numerical stability
      # This is an approximation using Gauss-Hermite quadrature
      # The scaling is important for stability as described in Bazzi et al. (2017)
      
      # In a real implementation, we would compute Fisher information matrix
      # and use the Moore-Penrose pseudo-inverse as in equation (15)
      # Here we use a simplification where we approximate the scaling
      
      # Normalize g_vector for numerical stability
      g_norm <- sqrt(sum(g_vector^2))
      if (g_norm > 0) {
        g_normalized <- g_vector / g_norm
      } else {
        g_normalized <- g_vector
      }
      
      # For simplicity in simulation, we're using a simplified scale factor
      # In a real implementation, this would involve integrating over the
      # distribution of y_t to get the expected Fisher information
      
      # Use a reasonable scaling factor based on the likelihood
      scale_factor <- 1.0 / sqrt(1.0 + abs(delta_like))
      
      # Scaled score vector
      score_scaled[,t] <- scale_factor * score[,t]
      
      # Update transition parameters using the score-driven dynamics
      # f_{t+1} = ω + A*s_t + B*(f_t - ω)
      f[,t+1] <- omega + A * score_scaled[,t] + B * (f[,t] - omega)
      
      # Transform back to probability space
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

#' Filter observed data through the score-driven (GAS) regime-switching model
#'
#' @param mu Vector of means corresponding to each regime
#' @param sigma2 Vector of variances corresponding to each regime
#' @param init_trans Initial transition probabilities for the latent process
#' @param A Scale parameter for score updates (sensitivity)
#' @param B Persistence parameter for score updates (memory)
#' @param y Observed time series increments
#' @param B_burnin Burn-in to be excluded at the beginning of the time series
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
#' A <- rep(0.1, 6)
#' B <- rep(0.9, 6)
#' y <- rnorm(1000) 
#' loglik <- Rfiltering_GAS(mu, sigma2, init_trans, A, B, y, 100, 50)
Rfiltering_GAS <- function(mu, sigma2, init_trans, A, B, y, B_burnin, C) {
  # Parameters passed to the function are:
  # mu          Vector of true means corresponding to each regime
  # sigma2      Vector of true variances corresponding to each regime
  # init_trans  Initial transition probabilities for the latent process
  # A           Scale parameter for score updates (sensitivity)
  # B           Persistence parameter for score updates (memory)
  # y           Observed time series increments
  # B_burnin    Burn-in to be excluded at the beginning of the time series
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
  
  # Ensure we have the right number of parameters
  if (length(init_trans) != n_transition) {
    stop(sprintf("Error: Expected %d transition probabilities for %d regimes, got %d.", 
                 n_transition, K, length(init_trans)))
  }
  if (length(A) != n_transition) {
    stop(sprintf("Error: Expected %d A parameters for %d regimes, got %d.", 
                 n_transition, K, length(A)))
  }
  if (length(B) != n_transition) {
    stop(sprintf("Error: Expected %d B parameters for %d regimes, got %d.", 
                 n_transition, K, length(B)))
  }
  
  # Initialize variables
  eta <- matrix(0, nrow=K, ncol=M)      # Likelihood of each regime
  tot_lik <- numeric(M)                 # Total likelihood
  X_t <- matrix(0, nrow=K, ncol=M)      # Filtered probabilities after observation
  X_tlag <- matrix(0, nrow=K, ncol=M)   # Predicted probabilities before observation
  
  # Arrays for score-driven dynamics
  f <- matrix(0, nrow=n_transition, ncol=M)       # Time-varying parameters (logit scale)
  p_trans <- matrix(0, nrow=n_transition, ncol=M) # Transition probabilities
  score <- matrix(0, nrow=n_transition, ncol=M)   # Score vectors
  score_scaled <- matrix(0, nrow=n_transition, ncol=M) # Scaled score vectors
  
  # Get baseline transition probabilities that are logit-transformed
  omega <- logit(init_trans)
  
  # Set initial values for transition parameters
  f[,1] <- omega
  p_trans_raw <- logistic(f[,1])
  p_trans[,1] <- convert_to_valid_probs(p_trans_raw, K)
  
  # Initial state probabilities (stationary distribution or uniform)
  X_t[,1] <- stat_dist(p_trans[,1], fallback_value = rep(1/K, K))
  
  for (t in 1:(M-1)) {
    # Generate predicted probabilities using current transition probabilities
    Pmatrix <- transition_matrix(p_trans[,t], check_validity = FALSE)
    X_tlag[,t] <- Pmatrix %*% X_t[,t]
    
    # Calculate likelihoods for each regime
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
    
    # Calculate the score vector (following Bazzi et al. 2017)
    # Compute likelihood differences for regimes
    delta_like <- numeric(n_transition)
    idx <- 1
    for (r in 1:K) {
      for (c in 1:K) {
        if (r != c) {
          delta_like[idx] <- (eta[r,t] - eta[c,t]) / tot_lik[t]
          idx <- idx + 1
        }
      }
    }
    
    # Compute the g vector
    g_vector <- numeric(n_transition)
    idx <- 1
    for (r in 1:K) {
      for (c in 1:K) {
        if (r != c) {
          g_vector[idx] <- X_t[r,t] * p_trans[idx,t] * (1 - p_trans[idx,t])
          idx <- idx + 1
        }
      }
    }
    
    # Unscaled score
    score[,t] <- delta_like * g_vector
    
    # Normalize g_vector for numerical stability
    g_norm <- sqrt(sum(g_vector^2))
    if (g_norm > 0) {
      g_normalized <- g_vector / g_norm
    } else {
      g_normalized <- g_vector
    }
    
    # Use a reasonable scaling factor
    scale_factor <- 1.0 / sqrt(1.0 + abs(delta_like))
    
    # Scaled score vector
    score_scaled[,t] <- scale_factor * score[,t]
    
    # Update transition parameters using the score-driven dynamics
    # f_{t+1} = ω + A*s_t + B*(f_t - ω)
    f[,t+1] <- omega + A * score_scaled[,t] + B * (f[,t] - omega)
    
    # Transform back to probability space
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
  valid_indices <- (B_burnin+1):(M-C)
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
  attr(res, "score") <- score
  attr(res, "score_scaled") <- score_scaled
  
  return(res)
}

#' Transform parameters from natural space to unconstrained space for GAS model
#'
#' @param par Parameters in their natural space
#' @return Parameters in unconstrained space
#' @details
#' Applies log transformation to variances, logit to transition probabilities,
#' and logit to A and B coefficient to make parameters suitable for unconstrained optimization.
transform_GAS <- function(par) {
  K <- sqrt(length(par)/3)  # For GAS model, we have K means, K variances, K*(K-1) transitions, K*(K-1) A's, and K*(K-1) B's
  
  if (K != as.integer(K)) {
    stop("The parameter vector length is not compatible with any valid number of regimes.")
  }
  
  K <- as.integer(K)
  n_transition <- K*(K-1)
  
  mu <- par[1:K]
  sigma2 <- par[(K+1):(2*K)]
  init_trans <- par[(2*K+1):(2*K+n_transition)]
  A <- par[(2*K+n_transition+1):(2*K+2*n_transition)]
  B <- par[(2*K+2*n_transition+1):(2*K+3*n_transition)]
  
  # Check for invalid values
  if (any(sigma2 <= 0)) {
    stop("Variance parameters must be positive")
  }
  if (any(init_trans <= 0) || any(init_trans >= 1)) {
    stop("Transition probabilities must be between 0 and 1")
  }
  if (any(A < 0) || any(A > 1)) {
    stop("A parameters must be between 0 and 1")
  }
  if (any(B < 0) || any(B > 1)) {
    stop("B parameters must be between 0 and 1")
  }
  
  # Transform parameters to unconstrained space
  mu_t <- mu  # No transformation for means
  sigma2_t <- log(sigma2)  # Log transformation for variances
  init_trans_t <- logit(init_trans)  # Logit transformation for probabilities
  
  # For A and B coefficients, we use logit transformation to map from (0,1) to (-Inf,Inf)
  A_t <- logit(A)
  B_t <- logit(B)
  
  return(c(mu_t, sigma2_t, init_trans_t, A_t, B_t))
}

#' Transform parameters from unconstrained space back to natural space for GAS model
#'
#' @param par_t Parameters in unconstrained space
#' @return Parameters in their natural space
#' @details
#' Applies exp transformation to log-variances, logistic to logit-probabilities,
#' and logistic to logit-coefficients to convert back to naturally bounded parameters.
untransform_GAS <- function(par_t) {
  K <- sqrt(length(par_t)/3)  # For GAS model, we have K means, K variances, K*(K-1) transitions, K*(K-1) A's, and K*(K-1) B's
  
  if (K != as.integer(K)) {
    stop("The parameter vector length is not compatible with any valid number of regimes.")
  }
  
  K <- as.integer(K)
  n_transition <- K*(K-1)
  
  mu_t <- par_t[1:K]
  sigma2_t <- par_t[(K+1):(2*K)]
  init_trans_t <- par_t[(2*K+1):(2*K+n_transition)]
  A_t <- par_t[(2*K+n_transition+1):(2*K+2*n_transition)]
  B_t <- par_t[(2*K+2*n_transition+1):(2*K+3*n_transition)]
  
  # Transform parameters back to natural space
  mu <- mu_t  # No transformation for means
  sigma2 <- exp(sigma2_t)  # Exp transformation for variances
  init_trans <- logistic(init_trans_t)  # Logistic transformation for probabilities
  
  # For A and B coefficients, we use logistic transformation to map from (-Inf,Inf) to (0,1)
  A <- logistic(A_t)
  B <- logistic(B_t)
  
  return(c(mu, sigma2, init_trans, A, B))
}

#' Count the number of regimes from a parameter vector for GAS model
#'
#' @param par Parameter vector
#' @return Number of regimes
count_regime_GAS <- function(par) {
  # For GAS model, we have K means, K variances, K*(K-1) transitions, K*(K-1) A's, and K*(K-1) B's
  # So total length is K + K + K*(K-1) + K*(K-1) + K*(K-1) = 2K + 3K*(K-1) = 2K + 3K^2 - 3K = 3K^2 - K
  # We need to solve for K: 3K^2 - K - length(par) = 0
  
  # Quadratic formula: K = (1 + sqrt(1 + 12*length(par)))/6
  K <- (1 + sqrt(1 + 12*length(par)))/6
  
  if (abs(K - round(K)) > 1e-10) {
    stop("The parameter vector length is not compatible with any valid number of regimes.")
  }
  
  return(round(K))
}

#' Wrapper for the likelihood calculator to be used for max-likelihood estimation
#'
#' @param par_t Transformed parameters in unconstrained space
#' @param y Observed time series increments
#' @param B_burnin Burn-in to be excluded at the beginning of the time series
#' @param C Cut-off to be excluded at the end of the time series
#' @return Negative log-likelihood of observed data under the model
#' @details
#' Transforms parameters from unconstrained space back to natural space
#' and calls the filtering function to calculate the likelihood.
#'
#' @examples
#' # Optimize parameters for a 3-regime model
#' transformed_params <- transform_GAS(c(mu, sigma2, init_trans, A, B))
#' result <- nlminb(transformed_params, Rfiltering.single.trasf_GAS, 
#'                 y = y, B_burnin = 100, C = 50)
Rfiltering.single.trasf_GAS <- function(par_t, y, B_burnin, C) {
  # Transform parameters back to original parameter space
  par <- untransform_GAS(par_t)
  
  # Determine the number of regimes
  K <- count_regime_GAS(par)
  n_transition <- K*(K-1)
  
  # Extract the components
  mu <- par[1:K]
  sigma2 <- par[(K+1):(2*K)]
  init_trans <- par[(2*K+1):(2*K+n_transition)]
  A <- par[(2*K+n_transition+1):(2*K+2*n_transition)]
  B <- par[(2*K+2*n_transition+1):(2*K+3*n_transition)]
  
  # Calculate likelihood and return it
  l <- Rfiltering_GAS(mu, sigma2, init_trans, A, B, y, B_burnin, C)
  return(l[1])
}

#' Estimate parameters for the score-driven (GAS) regime-switching model
#'
#' @param y Observed time series
#' @param K Number of regimes
#' @param B_burnin Burn-in period to exclude from likelihood calculation
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
#' results <- estimate_gas_model(data, K=3)
estimate_gas_model <- function(y, K = 3, B_burnin = 100, C = 50, 
                               initial_params = NULL, bounds = NULL,
                               verbose = TRUE) {
  if (verbose) {
    cat("Estimating GAS model with", K, "regimes\n")
    start_time <- Sys.time()
  }
  
  # Number of transition probabilities
  n_transition <- K*(K-1)
  
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
    
    # Initialize transition probabilities, A and B coefficients
    init_trans_guess <- rep(0.2, n_transition)
    A_guess <- rep(0.1, n_transition)
    B_guess <- rep(0.9, n_transition)
    
    initial_params <- c(mu_guess, sigma2_guess, init_trans_guess, A_guess, B_guess)
    
    if (verbose) {
      cat("Generated initial parameter guesses:\n")
      cat("Means:", round(mu_guess, 4), "\n")
      cat("Variances:", round(sigma2_guess, 4), "\n")
    }
  }
  
  # Create default bounds if none provided
  if (is.null(bounds)) {
    lower_bounds <- c(rep(-Inf, K),               # No bounds on means
                      rep(-Inf, K),               # Variance >= 0 (log-transformed)
                      rep(-Inf, n_transition),    # Probabilities >= 0 (logit-transformed)
                      rep(-Inf, n_transition),    # A-coefficients >= 0 (logit-transformed)
                      rep(-Inf, n_transition))    # B-coefficients >= 0 (logit-transformed)
    
    upper_bounds <- c(rep(Inf, K),                # No bounds on means
                      rep(Inf, K),                # Variance unbounded
                      rep(Inf, n_transition),     # Probabilities <= 1 (logit-transformed)
                      rep(Inf, n_transition),     # A-coefficients <= 1 (logit-transformed)
                      rep(Inf, n_transition))     # B-coefficients <= 1 (logit-transformed)
    
    bounds <- list(lower = lower_bounds, upper = upper_bounds)
  }
  
  # Transform parameters to unconstrained space for optimization
  transformed_params <- transform_GAS(initial_params)
  
  # Optimize parameters
  if (verbose) {
    cat("Starting optimization...\n")
    opt_start_time <- Sys.time()
  }
  
  optimization_result <- nlminb(
    start = transformed_params,
    objective = Rfiltering.single.trasf_GAS,
    lower = bounds$lower,
    upper = bounds$upper,
    y = y,
    B_burnin = B_burnin,
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
  estimated_params <- untransform_GAS(optimization_result$par)
  
  # Extract different parameter components
  mu_est <- estimated_params[1:K]
  sigma2_est <- estimated_params[(K+1):(2*K)]
  init_trans_est <- estimated_params[(2*K+1):(2*K+n_transition)]
  A_est <- estimated_params[(2*K+n_transition+1):(2*K+2*n_transition)]
  B_est <- estimated_params[(2*K+2*n_transition+1):(2*K+3*n_transition)]
  
  # Calculate model diagnostics
  num_params <- length(optimization_result$par)
  num_data_points <- length(y) - B_burnin - C
  
  aic <- 2 * optimization_result$objective + 2 * num_params
  bic <- 2 * optimization_result$objective + num_params * log(num_data_points)
  
  # Calculate filtered probabilities and other model outputs
  full_likelihood <- Rfiltering_GAS(
    mu_est, sigma2_est, init_trans_est, A_est, B_est, y, B_burnin, C
  )
  
  filtered_probs <- attr(full_likelihood, "X.t")
  transition_probs <- attr(full_likelihood, "p_trans")
  scores <- attr(full_likelihood, "score")
  scaled_scores <- attr(full_likelihood, "score_scaled")
  
  # Prepare results
  results <- list(
    parameters = list(
      mu = mu_est,
      sigma2 = sigma2_est,
      init_trans = init_trans_est,
      A = A_est,
      B = B_est
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
    scores = scores,
    scaled_scores = scaled_scores,
    model_info = list(
      type = "GAS",
      K = K,
      B_burnin = B_burnin,
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
    cat("Estimated A coefficients (mean):", round(mean(A_est), 4), "\n")
    cat("Estimated B coefficients (mean):", round(mean(B_est), 4), "\n")
  }
  
  return(results)
}

#' Simulate and estimate a GAS model with specified parameters
#'
#' @param M Number of simulation paths
#' @param N Length of each simulation path
#' @param mu True regime means
#' @param sigma2 True regime variances
#' @param init_trans True initial transition probabilities
#' @param A True score scaling coefficients
#' @param B True persistence coefficients
#' @param B_burnin Burn-in period
#' @param C Cut-off period
#' @param verbose Whether to print progress information
#' @return List with simulation results and estimation results
#' @details
#' This function simulates data from a GAS model and then estimates 
#' the model parameters, providing a way to test the estimation procedure
#' with known parameters.
#'
#' @examples
#' # Test estimation with 10 simulated paths
#' mu_true <- c(-2, 1, 2)
#' sigma2_true <- c(0.02, 0.2, 0.6)
#' init_trans_true <- rep(0.2, 6)
#' A_true <- rep(0.1, 6)
#' B_true <- rep(0.9, 6)
#' results <- example_GAS_simulation(10, 1000, mu_true, sigma2_true, init_trans_true, A_true, B_true)
example_GAS_simulation <- function(M = 10, N = 1000, 
                                   mu = c(-2, 1, 2), 
                                   sigma2 = c(0.02, 0.2, 0.6), 
                                   init_trans = rep(0.2, 6), 
                                   A = rep(0.1, 6),
                                   B = rep(0.9, 6),
                                   B_burnin = 100, C = 50, verbose = TRUE) {
  # Parameters
  K <- length(mu)    # Number of regimes
  n_transition <- K*(K-1)
  
  if (verbose) {
    cat("Setting up parameters for GAS simulation...\n")
    cat("Regimes:", K, "\n")
    cat("Means:", mu, "\n")
    cat("Variances:", sigma2, "\n")
    cat("Score scaling (A) mean:", mean(A), "\n")
    cat("Persistence (B) mean:", mean(B), "\n")
  }
  
  # Generate simulation data
  if (verbose) {
    cat("Generating simulation data for", M, "paths, each with length", N, "...\n")
    start_time <- Sys.time()
  }
  
  data_GAS_sim <- dataGASCD(M, N, mu, sigma2, init_trans, A, B)
  
  if (verbose) {
    end_time <- Sys.time()
    cat("Data generation completed in", 
        format(difftime(end_time, start_time), digits = 4), "\n")
  }
  
  # Storage for estimation results
  param_estimates <- matrix(0, M, 2*K + 3*n_transition)
  colnames(param_estimates) <- c(paste0("mu", 1:K), 
                                paste0("sigma2", 1:K), 
                                paste0("trans", 1:n_transition), 
                                paste0("A", 1:n_transition),
                                paste0("B", 1:n_transition))
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
    current_data <- data_GAS_sim[i,]
    
    # Initial parameter guesses for optimization
    mu_guess <- mu * 0.9 + rnorm(K, 0, 0.1*abs(mu))  # Add noise to true values
    sigma2_guess <- sigma2 * 0.9 + rnorm(K, 0, 0.1*sigma2)  # Add noise to true values
    init_trans_guess <- pmax(pmin(init_trans + rnorm(n_transition, 0, 0.05), 0.95), 0.05)  # Add noise, keep in (0,1)
    A_guess <- rep(0.1, n_transition)  # Start with moderate values
    B_guess <- rep(0.8, n_transition)  # Start with high persistence
    
    par_guess <- c(mu_guess, sigma2_guess, init_trans_guess, A_guess, B_guess)
    
    # Estimate the model
    GAS_est <- try(nlminb(
      start = transform_GAS(par_guess),
      objective = Rfiltering.single.trasf_GAS, 
      y = current_data,
      B_burnin = B_burnin, 
      C = C,
      control = list(eval.max = 1e6, iter.max = 1e6, trace = 0)
    ), silent = TRUE)
    
    if (!inherits(GAS_est, "try-error")) {
      # Transform parameters back to natural space
      GAS_par_fin <- untransform_GAS(GAS_est$par)
      
      # Store the parameter estimates
      param_estimates[i,] <- GAS_par_fin
      
      # Store diagnostics
      num_params <- length(GAS_est$par)
      diagnostics[i, 1] <- -GAS_est$objective
      diagnostics[i, 2] <- 2 * GAS_est$objective + 2 * num_params
      diagnostics[i, 3] <- 2 * GAS_est$objective + num_params * log(N - B_burnin - C)
      
      if (verbose) {
        path_end_time <- Sys.time()
        path_time <- difftime(path_end_time, path_start_time, units = "secs")
        cat("completed in", round(path_time, 2), "seconds (", 
            formatC(GAS_est$objective, digits = 8, format = "f"), ")\n")
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
  true_params <- c(mu, sigma2, init_trans, A, B)
  estimates_summary <- data.frame(
    Parameter = colnames(param_estimates),
    True_Value = true_params,
    Mean_Estimate = colMeans(param_estimates, na.rm = TRUE),
    SD_Estimate = apply(param_estimates, 2, sd, na.rm = TRUE),
    Bias = colMeans(param_estimates, na.rm = TRUE) - true_params,
    Rel_Bias_Pct = 100 * (colMeans(param_estimates, na.rm = TRUE) - true_params) / 
      ifelse(abs(true_params) > 1e-10, true_params, 1e-10)
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
      data = data_GAS_sim,
      parameters = list(
        mu = mu,
        sigma2 = sigma2,
        init_trans = init_trans,
        A = A,
        B = B
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
      B_burnin = B_burnin,
      C = C
    )
  )
  
  return(results)
}

#' Compare TVP, Exogenous, and GAS Models
#'
#' @param y Observed time series
#' @param X_Exo Optional exogenous variable for the Exogenous model
#' @param K Number of regimes
#' @param models Character vector specifying which models to compare
#' @param B_burnin Burn-in period
#' @param C Cut-off period
#' @param verbose Whether to print progress information
#' @return Data frame comparing model performance
#' @details
#' Estimates multiple time-varying transition probability models and 
#' compares their performance using likelihood-based criteria.
#'
#' @examples
#' # Compare all three models
#' data <- rnorm(1000)
#' exo <- rnorm(1000)
#' comparison <- compare_tvtp_models(data, exo, K=3, 
#'                                  models=c("TVP", "Exogenous", "GAS"))
compare_tvtp_models <- function(y, X_Exo = NULL, K = 3, 
                               models = c("Constant", "TVP", "Exogenous", "GAS"),
                               B_burnin = 100, C = 50, verbose = TRUE) {
  # Check if required models are available
  available_models <- c("Constant", "TVP", "Exogenous", "GAS")
  models <- match.arg(models, available_models, several.ok = TRUE)
  
  # Check if exogenous variable is provided when needed
  if ("Exogenous" %in% models && is.null(X_Exo)) {
    stop("Exogenous model requires X_Exo to be specified")
  }
  
  # Make sure exogenous variable has the right length
  if (!is.null(X_Exo) && length(X_Exo) != length(y)) {
    stop("X_Exo must have the same length as y")
  }
  
  # Prepare results
  results <- data.frame()
  
  # Estimate each model
  for (model_type in models) {
    if (verbose) {
      cat("\n--- Estimating", model_type, "model ---\n")
    }
    
    start_time <- Sys.time()
    
    # Estimate the specified model
    if (model_type == "Constant") {
      # Source the constant model file if not already loaded
      if (!exists("estimate_const_model")) {
        source("models/model_constant.R")
      }
      
      model_result <- estimate_const_model(y, K, B_burnin, C, verbose = verbose)
    } else if (model_type == "TVP") {
      # Source the TVP model file if not already loaded
      if (!exists("estimate_tvp_model")) {
        source("models/model_TVP.R")
      }
      
      model_result <- estimate_tvp_model(y, K, B_burnin, C, verbose = verbose)
    } else if (model_type == "Exogenous") {
      # Source the exogenous model file if not already loaded
      if (!exists("estimate_exo_model")) {
        source("models/model_exogenous.R")
      }
      
      model_result <- estimate_exo_model(y, X_Exo, K, B_burnin, C, verbose = verbose)
    } else if (model_type == "GAS") {
      model_result <- estimate_gas_model(y, K, B_burnin, C, verbose = verbose)
    }
    
    end_time <- Sys.time()
    elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Record results
    model_info <- data.frame(
      Model = model_type,
      LogLik = model_result$diagnostics$loglik,
      AIC = model_result$diagnostics$aic,
      BIC = model_result$diagnostics$bic,
      NumParams = model_result$diagnostics$num_params,
      EstimationTime = elapsed
    )
    
    results <- rbind(results, model_info)
    
    if (verbose) {
      cat("Estimation completed in", format_time(elapsed), "\n")
      cat("LogLik:", round(model_result$diagnostics$loglik, 2), 
          "AIC:", round(model_result$diagnostics$aic, 2),
          "BIC:", round(model_result$diagnostics$bic, 2), "\n")
    }
  }
  
  # Sort by AIC (best model first)
  results <- results[order(results$AIC), ]
  
  return(results)
}

#' Calculate regime persistence metrics
#'
#' @param filtered_probs Matrix of filtered probabilities (rows are time, columns are regimes)
#' @param threshold Probability threshold for regime identification (default: 0.5)
#' @return List with persistence metrics
#' @details
#' Calculates various persistence metrics for the estimated regimes
#' to better understand the model dynamics.
calculate_persistence <- function(filtered_probs, threshold = 0.5) {
  # Determine the most likely regime at each point in time
  most_likely_regime <- apply(filtered_probs, 1, which.max)
  
  # Calculate the number of transitions
  num_transitions <- sum(diff(most_likely_regime) != 0)
  transition_rate <- num_transitions / (length(most_likely_regime) - 1)
  
  # Calculate average duration in each regime
  regimes <- unique(most_likely_regime)
  durations <- list()
  
  for (r in regimes) {
    # Find consecutive stretches of this regime
    in_regime <- most_likely_regime == r
    run_lengths <- rle(in_regime)
    lengths <- run_lengths$lengths[run_lengths$values]
    
    if (length(lengths) > 0) {
      durations[[as.character(r)]] <- lengths
    } else {
      durations[[as.character(r)]] <- 0
    }
  }
  
  # Calculate average durations
  avg_durations <- sapply(durations, mean)
  
  # Calculate regime occupancy percentages
  regime_counts <- table(most_likely_regime)
  regime_percentages <- 100 * regime_counts / length(most_likely_regime)
  
  return(list(
    most_likely_regime = most_likely_regime,
    num_transitions = num_transitions,
    transition_rate = transition_rate,
    durations = durations,
    avg_durations = avg_durations,
    regime_counts = regime_counts,
    regime_percentages = regime_percentages
  ))
}

#' Helper function to format time in a human-readable way
#'
#' @param seconds Time in seconds
#' @return String with formatted time
#' @examples
#' format_time(65)  # Returns "1m 5s"
format_time <- function(seconds) {
  if (seconds < 60) {
    return(sprintf("%.1fs", seconds))
  } else if (seconds < 3600) {
    minutes <- floor(seconds / 60)
    secs <- round(seconds - minutes * 60)
    return(sprintf("%dm %ds", minutes, secs))
  } else {
    hours <- floor(seconds / 3600)
    minutes <- floor((seconds - hours * 3600) / 60)
    return(sprintf("%dh %dm", hours, minutes))
  }
}
