#' Generalized Autoregressive Score Models with Time-Varying Transition Probabilities
#' 
#' This file implements a regime-switching model where transition probabilities
#' are driven by the score of the predictive likelihood function (GAS dynamics).
#'
#' Reference: Bazzi, M., Blasques, F., Koopman, S.J., & Lucas, A. (2017).
#' Time-Varying Transition Probabilities for Markov Regime Switching Models.
#' Journal of Time Series Analysis, 38(3), 458-478.

# Load required helper functions
source("helpers/utility_functions.R")
source("helpers/transition_helpers.R")
source("helpers/parameter_transforms.R")
source("helpers/score_functions.R")
source("models/model_constant.R")

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
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method ("moore_penrose", "simple", "normalized", or "original")
#' @param quad_sample_size Sample size for creating representative data for quadrature setup (default: 1000)
#' @return Matrix of simulated data with M rows and N columns
#' @details 
#' Simulates data from a regime switching model where transition probabilities
#' are updated using score-driven dynamics based on the predictive likelihood.
#' The update equation is f[t+1] = omega + A*s[t] + B*(f[t] - omega).
#' 
#' This implementation uses proper GAS score scaling following Bazzi et al. (2017)
#' with Gauss-Hermite quadrature and Moore-Penrose pseudo-inverse scaling.
#' 
#' The quadrature setup is based on representative data generated from the regime
#' parameters, ensuring the integration domain matches the data-generating process.
#'
#' @examples
#' # Generate data for a 3-regime model
#' mu_true <- c(-2, 1, 2)
#' sigma2_true <- c(0.02, 0.2, 0.6)
#' init_trans_true <- rep(0.2, 6)  # Off-diagonal elements for a 3x3 matrix
#' A_true <- rep(0.1, 6)  # Score scaling parameters
#' B_true <- rep(0.9, 6)  # Persistence parameters
#' data_sim <- dataGASCD(10, 1000, mu_true, sigma2_true, init_trans_true, A_true, B_true)
dataGASCD <- function(M, N, mu, sigma2, init_trans, A, B, burn_in = 100,
                      n_nodes = 30, scaling_method = "moore_penrose", quad_sample_size = 1000) {
  # Parameters:
  # M                 Number of simulation runs to be performed
  # N                 Length of the simulation runs (discretized time)
  # mu                Vector of true means corresponding to each regime
  # sigma2            Vector of true variances corresponding to each regime
  # init_trans        Initial transition probabilities for the latent process
  # A                 Scaling parameters for score updates (sensitivity)
  # B                 Persistence parameters for score updates (memory)
  # burn_in           Number of burn-in observations to discard
  # n_nodes           Number of Gauss-Hermite quadrature nodes
  # scaling_method    Score scaling method to use
  # quad_sample_size  Sample size for quadrature setup
  
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
  
  # Validate scaling method
  scaling_method <- match.arg(scaling_method, c("moore_penrose", "simple", "normalized", "original"))
  
  # Setup Gauss-Hermite quadrature based on regime parameters
  # Create representative data from regime parameters for quadrature setup
  tryCatch({
    # Generate sample sizes for each regime (equal weighting)
    regime_weights <- rep(1/K, K)
    sample_sizes <- round(quad_sample_size * regime_weights)
    
    # Ensure we have at least some samples for each regime
    sample_sizes <- pmax(sample_sizes, 50)
    
    # Generate representative data from each regime
    typical_data <- c()
    for (k in 1:K) {
      regime_data <- rnorm(sample_sizes[k], mu[k], sqrt(sigma2[k]))
      typical_data <- c(typical_data, regime_data)
    }
    
    # Setup quadrature based on this representative data
    gh_setup <- setup_gauss_hermite_quadrature(typical_data, n_nodes = n_nodes, method = "data_based")
    
  }, error = function(e) {
    warning("Failed to setup quadrature with regime parameters. Using standardized setup: ", e$message)
    # Fallback to standardized setup
    gh_setup <<- setup_gauss_hermite_quadrature(rnorm(quad_sample_size), n_nodes = n_nodes, method = "standardized")
  })
  
  # Set up a matrix to save the output
  data <- matrix(0, M, N)
  
  # Get baseline transition probabilities that are logit-transformed
  omega <- logit(init_trans)
  
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
    score_scaled <- matrix(0, nrow=n_transition, ncol=total_length) # Scaled score vectors
    
    # For debugging and diagnostics
    fisher_info_series <- numeric(total_length)
    score_norms <- numeric(total_length)
    
    # Initial state probabilities (uniform distribution)
    X_t[,1] <- rep(1/K, K)
    
    # Set initial values for transition parameters
    f[,1] <- omega
    p_trans_raw <- logistic(f[,1])
    p_trans[,1] <- convert_to_valid_probs(p_trans_raw, K)
    
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
      if (tot_lik[t] <= 0 || is.na(tot_lik[t])) {
        # Handle numerical issues with warning
        warning(paste("Invalid total likelihood", tot_lik[t], "at time", t, 
                      "in simulation", i, "- resetting to uniform probabilities and zero scores"))
        tot_lik[t] <- .Machine$double.eps
        X_t[,t+1] <- rep(1/K, K)
        score_scaled[,t] <- rep(0, n_transition)
        fisher_info_series[t] <- NA
        score_norms[t] <- 0
      } else {
        X_t[,t+1] <- (eta[,t]*X_tlag[,t])/tot_lik[t]
        
        # Calculate properly scaled score using our new helper functions
        tryCatch({
          # Normalize probabilities to ensure they sum to 1 (fix numerical precision issues)
          X_tlag_sum <- sum(X_tlag[,t])
          if (X_tlag_sum <= 0 || !is.finite(X_tlag_sum)) {
            # Fallback to uniform distribution if X_tlag is invalid
            warning(paste("Invalid X_tlag probabilities at time", t, "in simulation", i,
                          "with sum =", X_tlag_sum, ". Using uniform fallback."))
            X_tlag_normalized <- rep(1/K, K)
          } else {
            X_tlag_normalized <- X_tlag[,t] / X_tlag_sum
          }
          
          # Also ensure X_t probabilities are valid
          X_t_sum <- sum(X_t[,t])
          if (X_t_sum <= 0 || !is.finite(X_t_sum)) {
            warning(paste("Invalid X_t probabilities at time", t, "in simulation", i,
                          "with sum =", X_t_sum, ". Using uniform fallback."))
            X_t_normalized <- rep(1/K, K)
          } else {
            X_t_normalized <- X_t[,t] / X_t_sum
          }
          
          # Use proper GAS score calculation from helpers
          scaled_score_t <- calculate_gas_score(
            y_t = y.sim[t],
            mu = mu,
            sigma2 = sigma2,
            X_t_lag = X_tlag_normalized,  # Use normalized version
            X_t_prev = X_t_normalized,   # Use normalized version
            p_trans = p_trans[,t],
            gh_setup = gh_setup,
            K = K,
            scaling_method = scaling_method
          )
          
          score_scaled[,t] <- scaled_score_t
          
          # Extract diagnostics if available
          calc_info <- attr(scaled_score_t, "calculation_info")
          if (!is.null(calc_info)) {
            fisher_info_series[t] <- calc_info$fisher_info
            score_norms[t] <- calc_info$scaled_score_norm
          } else {
            fisher_info_series[t] <- NA
            score_norms[t] <- sqrt(sum(scaled_score_t^2))
          }
          
        }, error = function(e) {
          # Fallback to zero score if calculation fails - with warning
          warning(paste("GAS score calculation failed at time", t, "in simulation", i, 
                        "- using zero scores. Error:", e$message))
          score_scaled[,t] <<- rep(0, n_transition)
          fisher_info_series[t] <<- NA
          score_norms[t] <<- 0
        })
      }
      
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
    
    # Check for excessive zero scores and warn if problematic
    zero_score_count <- sum(apply(score_scaled[, 1:(total_length-1)], 2, function(x) all(x == 0)))
    if (zero_score_count > (total_length-1) * 0.1) {  # Warn if >10% zero scores
      warning(paste("Simulation", i, "used zero scores for", zero_score_count, 
                    "out of", total_length-1, "time points (", 
                    round(100*zero_score_count/(total_length-1), 1), "%)"))
    }
    
    # Optional: Store additional simulation diagnostics for debugging
    if (getOption("store_gas_simulation_details", FALSE)) {
      attr(data, paste0("sim_", i, "_diagnostics")) <- list(
        S = S[(burn_in+1):total_length],
        X_t = X_t[, (burn_in+1):total_length],
        p_trans = p_trans[, (burn_in+1):total_length],
        fisher_info = fisher_info_series[(burn_in+1):total_length],
        score_norms = score_norms[(burn_in+1):total_length]
      )
    }
  }
  
  # Store simulation metadata
  attr(data, "simulation_info") <- list(
    M = M,
    N = N,
    K = K,
    burn_in = burn_in,
    parameters = list(
      mu = mu,
      sigma2 = sigma2,
      init_trans = init_trans,
      A = A,
      B = B
    ),
    gas_settings = list(
      n_nodes = n_nodes,
      scaling_method = scaling_method,
      quad_sample_size = quad_sample_size
    ),
    gh_setup = gh_setup
  )
  
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
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method ("moore_penrose", "simple", "normalized", or "original")
#' @return Negative log-likelihood of observed data under the model
#' @details
#' Filters observed data through the model to compute the likelihood using proper
#' GAS score scaling following Bazzi et al. (2017). This implementation uses
#' the sophisticated score calculation helpers with Gauss-Hermite quadrature 
#' and Moore-Penrose pseudo-inverse scaling.
#'
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
Rfiltering_GAS <- function(mu, sigma2, init_trans, A, B, y, B_burnin, C, 
                           n_nodes = 30, scaling_method = "moore_penrose") {
  # Parameters passed to the function are:
  # mu          Vector of true means corresponding to each regime
  # sigma2      Vector of true variances corresponding to each regime
  # init_trans  Initial transition probabilities for the latent process
  # A           Scale parameter for score updates (sensitivity)
  # B           Persistence parameter for score updates (memory)
  # y           Observed time series increments
  # B_burnin    Burn-in to be excluded at the beginning of the time series
  # C           Cut-off to be excluded at the end of the time series
  # n_nodes     Number of Gauss-Hermite quadrature nodes
  # scaling_method Score scaling method to use
  
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
  
  # Validate scaling method
  scaling_method <- match.arg(scaling_method, c("moore_penrose", "simple", "normalized", "original"))
  
  # Setup Gauss-Hermite quadrature based on observed data
  # This is the key improvement - using actual data characteristics
  tryCatch({
    gh_setup <- setup_gauss_hermite_quadrature(y, n_nodes = n_nodes, method = "data_based")
  }, error = function(e) {
    warning("Failed to setup quadrature with data. Using standardized setup: ", e$message)
    # Fallback to standardized setup
    gh_setup <<- setup_gauss_hermite_quadrature(rnorm(1000), n_nodes = n_nodes, method = "standardized")
  })
  
  # Initialize variables
  eta <- matrix(0, nrow=K, ncol=M)      # Likelihood of each regime
  tot_lik <- numeric(M)                 # Total likelihood
  X_t <- matrix(0, nrow=K, ncol=M)      # Filtered probabilities after observation
  X_tlag <- matrix(0, nrow=K, ncol=M)   # Predicted probabilities before observation
  
  # Arrays for score-driven dynamics
  f <- matrix(0, nrow=n_transition, ncol=M)       # Time-varying parameters (logit scale)
  p_trans <- matrix(0, nrow=n_transition, ncol=M) # Transition probabilities
  score_scaled <- matrix(0, nrow=n_transition, ncol=M) # Scaled score vectors
  
  # For debugging and analysis
  fisher_info_series <- numeric(M)
  score_norms <- numeric(M)
  
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
      warning(paste("Invalid total likelihood", tot_lik[t], "at time", t, 
                    "- resetting to uniform probabilities and zero scores"))
      tot_lik[t] <- .Machine$double.eps
      X_t[,t+1] <- rep(1/K, K)  # Reset to uniform when we get an invalid likelihood
      
      # Use zero score for problematic observations
      score_scaled[,t] <- rep(0, n_transition)
      fisher_info_series[t] <- NA
      score_norms[t] <- 0
      
    } else {
      # Calculate filtered probabilities
      X_t[,t+1] <- (eta[,t]*X_tlag[,t])/tot_lik[t]
      
      # Calculate properly scaled score using our helper functions
      tryCatch({
        # Normalize X_tlag to ensure it sums to 1 (fix numerical precision issues)
        X_tlag_sum <- sum(X_tlag[,t])
        if (X_tlag_sum <= 0 || !is.finite(X_tlag_sum)) {
          # Fallback to uniform distribution if X_tlag is invalid
          warning(paste("Invalid X_tlag probabilities at time", t, 
                        "with sum =", X_tlag_sum, ". Using uniform fallback."))
          X_tlag_normalized <- rep(1/K, K)
        } else {
          X_tlag_normalized <- X_tlag[,t] / X_tlag_sum
        }
        
        # Also ensure X_t probabilities are valid
        X_t_sum <- sum(X_t[,t])
        if (X_t_sum <= 0 || !is.finite(X_t_sum)) {
          warning(paste("Invalid X_t probabilities at time", t, 
                        "with sum =", X_t_sum, ". Using uniform fallback."))
          X_t_normalized <- rep(1/K, K)
        } else {
          X_t_normalized <- X_t[,t] / X_t_sum
        }
        
        # Use proper GAS score calculation from helpers
        scaled_score_t <- calculate_gas_score(
          y_t = y[t],
          mu = mu,
          sigma2 = sigma2,
          X_t_lag = X_tlag_normalized,  # Use normalized version
          X_t_prev = X_t_normalized,   # Use normalized version
          p_trans = p_trans[,t],
          gh_setup = gh_setup,
          K = K,
          scaling_method = scaling_method
        )
        
        score_scaled[,t] <- scaled_score_t
        
        # Extract diagnostics if available
        calc_info <- attr(scaled_score_t, "calculation_info")
        if (!is.null(calc_info)) {
          fisher_info_series[t] <- calc_info$fisher_info
          score_norms[t] <- calc_info$scaled_score_norm
        } else {
          fisher_info_series[t] <- NA
          score_norms[t] <- sqrt(sum(scaled_score_t^2))
        }
        
      }, error = function(e) {
        # Fallback to zero score if calculation fails - with warning
        warning(paste("GAS score calculation failed at time", t, 
                      "- using zero scores. Error:", e$message))
        score_scaled[,t] <<- rep(0, n_transition)
        fisher_info_series[t] <<- NA
        score_norms[t] <<- 0
      })
    }
    
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
  
  # Check for excessive zero scores and warn if problematic
  valid_score_indices <- B_burnin:(M-1-C)  # Exclude burn-in and cut-off from score check
  if (length(valid_score_indices) > 0) {
    zero_score_count <- sum(apply(score_scaled[, valid_score_indices, drop=FALSE], 2, function(x) all(x == 0)))
    if (zero_score_count > length(valid_score_indices) * 0.1) {  # Warn if >10% zero scores
      warning(paste("GAS filtering used zero scores for", zero_score_count, 
                    "out of", length(valid_score_indices), "time points (", 
                    round(100*zero_score_count/length(valid_score_indices), 1), "%)"))
    }
  }
  
  # Return negative sum of log-likelihoods (for minimizing)
  res <- -logLikSum
  
  # Store additional information as attributes
  attr(res, "X.t") <- t(X_t)
  attr(res, "X.tlag") <- t(X_tlag)
  attr(res, "p_trans") <- p_trans
  attr(res, "score_scaled") <- score_scaled
  attr(res, "f") <- f  # Time-varying parameters in logit space
  
  # Store GAS-specific diagnostics
  attr(res, "gas_diagnostics") <- list(
    fisher_info_series = fisher_info_series,
    score_norms = score_norms,
    gh_setup = gh_setup,
    scaling_method = scaling_method,
    n_nodes = n_nodes,
    mean_fisher_info = mean(fisher_info_series, na.rm = TRUE),
    mean_score_norm = mean(score_norms, na.rm = TRUE),
    zero_score_percentage = if(length(valid_score_indices) > 0) {
      100 * sum(apply(score_scaled[, valid_score_indices, drop=FALSE], 2, function(x) all(x == 0))) / length(valid_score_indices)
    } else { 0 }
  )
  
  return(res)
}

#' Transform parameters from natural space to unconstrained space for GAS model
#'
#' @param par Parameters in their natural space
#' @return Parameters in unconstrained space
#' @details
#' Applies log transformation to variances, logit to transition probabilities,
#' and logit to A and B coefficients to make parameters suitable for unconstrained optimization.
transform_GAS <- function(par) {
  K <- count_regime_GAS(par)
  n_transition <- K*(K-1)
  
  mu <- par[1:K]
  sigma2 <- par[(K+1):(2*K)]
  init_trans <- par[(2*K+1):(2*K+n_transition)]
  A <- par[(2*K+n_transition+1):(2*K+2*n_transition)]
  B <- par[(2*K+2*n_transition+1):(2*K+3*n_transition)]
  
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
  
  return(c(mu, log(sigma2), logit(init_trans), logit(A), logit(B)))
}

#' Transform parameters from unconstrained space back to natural space for GAS model
#'
#' @param par_t Parameters in unconstrained space
#' @return Parameters in their natural space
#' @details
#' Applies exp transformation to log-variances, logistic to logit-probabilities,
#' and logistic to logit-coefficients to convert back to naturally bounded parameters.
untransform_GAS <- function(par_t) {
  K <- count_regime_GAS(par_t)
  n_transition <- K*(K-1)
  
  mu_t <- par_t[1:K]
  sigma2_t <- par_t[(K+1):(2*K)]
  init_trans_t <- par_t[(2*K+1):(2*K+n_transition)]
  A_t <- par_t[(2*K+n_transition+1):(2*K+2*n_transition)]
  B_t <- par_t[(2*K+2*n_transition+1):(2*K+3*n_transition)]
  
  return(c(mu_t, exp(sigma2_t), logistic(init_trans_t), logistic(A_t), logistic(B_t)))
}

#' Count the number of regimes from a parameter vector for GAS model
#'
#' @param par Parameter vector
#' @return Number of regimes
count_regime_GAS <- function(par) {
  # For GAS model: K + K + K*(K-1) + K*(K-1) + K*(K-1) = 3K^2 - K
  # So: 3K^2 - K - length(par) = 0
  # Using quadratic formula: K = (1 + sqrt(1 + 12*length(par)))/6
  discriminant <- 1 + 12*length(par)
  if (discriminant < 0) {
    stop("Invalid parameter vector length for GAS model")
  }
  
  K <- (1 + sqrt(discriminant))/6
  
  if (abs(K - round(K)) > 1e-10) {
    stop(paste("The parameter vector length", length(par), 
               "is not compatible with any valid number of regimes."))
  }
  
  return(round(K))
}

#' Wrapper for the likelihood calculator to be used for max-likelihood estimation
#'
#' @param par_t Transformed parameters in unconstrained space
#' @param y Observed time series increments
#' @param B_burnin Burn-in to be excluded at the beginning of the time series
#' @param C Cut-off to be excluded at the end of the time series
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method (default: "moore_penrose")
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
Rfiltering.single.trasf_GAS <- function(par_t, y, B_burnin, C, 
                                        n_nodes = 30, scaling_method = "moore_penrose") {
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
  l <- Rfiltering_GAS(mu, sigma2, init_trans, A, B, y, B_burnin, C, n_nodes, scaling_method)
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
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method (default: "moore_penrose")
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
                               n_nodes = 30, scaling_method = "moore_penrose",
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
      cat("A coefficients (sensitivity):", round(A_guess, 4), "\n")
      cat("B coefficients (persistence):", round(B_guess, 4), "\n")
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
    n_nodes = n_nodes,
    scaling_method = scaling_method,
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
    mu_est, sigma2_est, init_trans_est, A_est, B_est, y, B_burnin, C, n_nodes, scaling_method
  )
  
  filtered_probs <- attr(full_likelihood, "X.t")
  transition_probs <- attr(full_likelihood, "p_trans")
  scaled_scores <- attr(full_likelihood, "score_scaled")
  gas_diagnostics_raw <- attr(full_likelihood, "gas_diagnostics")
  
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
    scaled_scores = scaled_scores,
    gas_diagnostics = gas_diagnostics_raw,
    model_info = list(
      type = "GAS",
      K = K,
      B_burnin = B_burnin,
      C = C,
      n_nodes = n_nodes,
      scaling_method = scaling_method
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
    
    # GAS-specific diagnostics
    if (!is.null(gas_diagnostics_raw)) {
      cat("GAS diagnostics - Mean Fisher Info:", round(gas_diagnostics_raw$mean_fisher_info, 4), "\n")
      cat("GAS diagnostics - Mean Score Norm:", round(gas_diagnostics_raw$mean_score_norm, 4), "\n")
    }
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
#' @param n_nodes Number of Gauss-Hermite quadrature nodes
#' @param scaling_method Score scaling method
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
                                   B_burnin = 100, C = 50, 
                                   n_nodes = 30, scaling_method = "moore_penrose",
                                   verbose = TRUE) {
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
    cat("GAS settings - Nodes:", n_nodes, "Scaling:", scaling_method, "\n")
  }
  
  # Generate simulation data
  if (verbose) {
    cat("Generating simulation data for", M, "paths, each with length", N, "...\n")
    start_time <- Sys.time()
  }
  
  data_GAS_sim <- dataGASCD(M, N, mu, sigma2, init_trans, A, B, 
                            burn_in = 0, n_nodes = n_nodes, 
                            scaling_method = scaling_method)
  
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
  
  # Storage for GAS-specific diagnostics
  gas_diagnostics <- list()
  
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
      n_nodes = n_nodes,
      scaling_method = scaling_method,
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
      
      # Store GAS-specific diagnostics (simplified for batch processing)
      gas_diagnostics[[i]] <- list(
        converged = GAS_est$convergence == 0,
        iterations = GAS_est$iterations,
        scaling_method = scaling_method
      )
      
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
      gas_diagnostics[[i]] <- NULL
      
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
    
    # GAS-specific summary
    successful_runs <- !is.na(diagnostics[,1])
    if (sum(successful_runs) > 0) {
      convergence_rate <- mean(sapply(gas_diagnostics[successful_runs], function(x) x$converged), na.rm = TRUE)
      cat("Convergence rate:", round(convergence_rate * 100, 1), "%\n")
    }
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
    gas_diagnostics = gas_diagnostics,
    settings = list(
      M = M,
      N = N,
      K = K,
      B_burnin = B_burnin,
      C = C,
      n_nodes = n_nodes,
      scaling_method = scaling_method
    )
  )
  
  return(results)
}

#' Helper function to calculate regime persistence metrics
#'
#' @param filtered_probs Matrix of filtered probabilities (T x K)
#' @return List with persistence metrics
#' @details
#' Calculates various metrics about regime persistence and transitions
#' from the filtered probability sequences.
calculate_persistence <- function(filtered_probs) {
  if (!is.matrix(filtered_probs)) {
    stop("filtered_probs must be a matrix")
  }
  
  T_obs <- nrow(filtered_probs)
  K <- ncol(filtered_probs)
  
  # Get most likely regime at each time point
  regime_sequence <- apply(filtered_probs, 1, which.max)
  
  # Count transitions
  transitions <- sum(diff(regime_sequence) != 0)
  transition_rate <- transitions / (T_obs - 1)
  
  # Calculate average durations in each regime
  rle_result <- rle(regime_sequence)
  regime_durations <- split(rle_result$lengths, rle_result$values)
  
  avg_durations <- numeric(K)
  for (k in 1:K) {
    if (k %in% names(regime_durations)) {
      avg_durations[k] <- mean(regime_durations[[as.character(k)]])
    } else {
      avg_durations[k] <- 0
    }
  }
  
  # Calculate regime occupancy percentages
  regime_counts <- table(factor(regime_sequence, levels = 1:K))
  regime_percentages <- as.numeric(regime_counts) / T_obs * 100
  
  return(list(
    num_transitions = transitions,
    transition_rate = transition_rate,
    avg_durations = avg_durations,
    regime_percentages = regime_percentages,
    regime_sequence = regime_sequence
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



#' Create parameter bounds for L-BFGS-B optimization in natural space
#'
#' @param K Number of regimes
#' @param y_mean Mean of observed data (for setting sensible bounds on mu)
#' @param y_sd Standard deviation of observed data (for setting sensible bounds)
#' @return List with lower and upper bounds for natural parameters
create_gas_bounds_natural <- function(K, y_mean = 0, y_sd = 1) {
  n_transition <- K * (K - 1)
  
  # Parameter order: mu, sigma2, init_trans, A, B
  
  # Bounds for means (allow wide range around data characteristics)
  mu_lower <- rep(y_mean - 10 * y_sd, K)
  mu_upper <- rep(y_mean + 10 * y_sd, K)
  
  # Bounds for variances (must be positive, reasonable upper bound)
  sigma2_lower <- rep(1e-6, K)
  sigma2_upper <- rep(100 * y_sd^2, K)
  
  # Bounds for transition probabilities (must be in (0,1))
  trans_lower <- rep(1e-6, n_transition)
  trans_upper <- rep(1 - 1e-6, n_transition)
  
  # Bounds for A parameters (sensitivity, should be non-negative and bounded)
  A_lower <- rep(0, n_transition)
  A_upper <- rep(1, n_transition)
  
  # Bounds for B parameters (persistence, should be in [0,1) for stationarity)
  B_lower <- rep(0, n_transition)
  B_upper <- rep(1 - 1e-6, n_transition)
  
  return(list(
    lower = c(mu_lower, sigma2_lower, trans_lower, A_lower, B_lower),
    upper = c(mu_upper, sigma2_upper, trans_upper, A_upper, B_upper)
  ))
}

#' Validate parameters in natural space
#'
#' @param par Parameter vector in natural space
#' @param K Number of regimes
#' @return TRUE if valid, throws error otherwise
validate_gas_params_natural <- function(par, K) {
  n_transition <- K * (K - 1)
  expected_length <- K + K + n_transition + n_transition + n_transition  # mu + sigma2 + trans + A + B
  
  if (length(par) != expected_length) {
    stop(sprintf("Parameter vector has length %d, expected %d for K=%d regimes", 
                 length(par), expected_length, K))
  }
  
  # Extract components
  mu <- par[1:K]
  sigma2 <- par[(K+1):(2*K)]
  init_trans <- par[(2*K+1):(2*K+n_transition)]
  A <- par[(2*K+n_transition+1):(2*K+2*n_transition)]
  B <- par[(2*K+2*n_transition+1):(2*K+3*n_transition)]
  
  # Check for invalid values
  if (any(is.na(par)) || any(is.infinite(par))) {
    stop("Parameter vector contains NA or infinite values")
  }
  
  # Check variance bounds
  if (any(sigma2 <= 0)) {
    stop("All variance parameters must be positive")
  }
  
  # Check transition probability bounds
  if (any(init_trans <= 0) || any(init_trans >= 1)) {
    stop("All transition probabilities must be in (0,1)")
  }
  
  # Check A parameter bounds
  if (any(A < 0) || any(A > 1)) {
    stop("All A parameters must be in [0,1]")
  }
  
  # Check B parameter bounds
  if (any(B < 0) || any(B >= 1)) {
    stop("All B parameters must be in [0,1)")
  }
  
  return(TRUE)
}



#' Wrapper for the likelihood calculator using natural parameters with L-BFGS-B
#'
#' @param par Natural parameters (not transformed)
#' @param y Observed time series increments
#' @param B_burnin Burn-in to be excluded at the beginning of the time series
#' @param C Cut-off to be excluded at the end of the time series
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method (default: "moore_penrose")
#' @return Negative log-likelihood of observed data under the model
#' @details
#' This wrapper works directly with natural parameters for L-BFGS-B optimization.
#' It provides informative error messages for debugging rather than silent penalties.
Rfiltering.single.natural_GAS <- function(par, y, B_burnin, C, 
                                          n_nodes = 30, scaling_method = "moore_penrose") {
  
  # Determine the number of regimes
  K <- count_regime_GAS(par)
  n_transition <- K*(K-1)
  
  # Validate parameters with informative error messages
  validate_gas_params_natural(par, K)
  
  # Extract the components (already in natural space)
  mu <- par[1:K]
  sigma2 <- par[(K+1):(2*K)]
  init_trans <- par[(2*K+1):(2*K+n_transition)]
  A <- par[(2*K+n_transition+1):(2*K+2*n_transition)]
  B <- par[(2*K+2*n_transition+1):(2*K+3*n_transition)]
  
  # Check for transition probabilities that would create invalid transition matrices
  p_trans_test <- convert_to_valid_probs(init_trans, K)
  test_matrix <- transition_matrix(p_trans_test, check_validity = TRUE)
  
  # Calculate likelihood
  l <- Rfiltering_GAS(mu, sigma2, init_trans, A, B, y, B_burnin, C, n_nodes, scaling_method)
  
  # Check for invalid likelihood
  if (!is.finite(l) || is.na(l)) {
    stop(sprintf("Likelihood calculation returned invalid value: %s. 
                 Parameters: mu=%s, sigma2=%s, A_mean=%.4f, B_mean=%.4f", 
                 l, 
                 paste(round(mu, 4), collapse=","),
                 paste(round(sigma2, 4), collapse=","),
                 mean(A), mean(B)))
  }
  
  return(l[1])
}



#' Estimate parameters for the score-driven (GAS) regime-switching model using L-BFGS-B
#'
#' @param y Observed time series
#' @param K Number of regimes
#' @param B_burnin Burn-in period to exclude from likelihood calculation
#' @param C Cut-off period to exclude from likelihood calculation
#' @param initial_params Initial parameter guesses (optional)
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method (default: "moore_penrose")
#' @param max_iterations Maximum number of optimization iterations (default: 1000)
#' @param factr Controls precision of L-BFGS-B (default: 1e7)
#' @param verbose Whether to print progress information (default: TRUE)
#' @return List with estimated parameters and model diagnostics
estimate_gas_model_lbfgsb <- function(y, K = 3, B_burnin = 100, C = 50, 
                                      initial_params = NULL,
                                      n_nodes = 30, scaling_method = "moore_penrose",
                                      max_iterations = 1000, factr = 1e7,
                                      verbose = TRUE) {
  if (verbose) {
    cat("Estimating GAS model with", K, "regimes using L-BFGS-B\n")
    start_time <- Sys.time()
  }
  
  # Number of transition probabilities
  n_transition <- K*(K-1)
  
  # Get data characteristics for bounds
  y_mean <- mean(y, na.rm = TRUE)
  y_sd <- sd(y, na.rm = TRUE)
  
  # Create parameter bounds
  bounds <- create_gas_bounds_natural(K, y_mean, y_sd)
  
  # Create default initial parameters if none provided
  if (is.null(initial_params)) {
    # Create MUCH better initial guesses using data characteristics
    y_mean <- mean(y, na.rm = TRUE)
    y_sd <- sd(y, na.rm = TRUE)
    
    # Use more conservative initial regime means (closer to data mean)
    spread <- 1.0 * y_sd  # Reduced from 2.0 * y_sd
    mu_guess <- seq(y_mean - spread, y_mean + spread, length.out = K)
    
    # Use variance closer to data variance
    sigma2_guess <- rep(y_sd^2, K)  # Start all regimes with same variance
    
    # Conservative transition probabilities and GAS parameters
    init_trans_guess <- rep(0.2, n_transition)  # More conservative
    A_guess <- rep(0.01, n_transition)  # Much smaller for stability
    B_guess <- rep(0.98, n_transition)  # Higher persistence
    
    initial_params <- c(mu_guess, sigma2_guess, init_trans_guess, A_guess, B_guess)
    
    if (verbose) {
      cat("Generated conservative initial parameter guesses:\n")
      cat("Means:", round(mu_guess, 4), "\n")
      cat("Variances:", round(sigma2_guess, 4), "\n")
      cat("A coefficients (sensitivity):", round(A_guess, 4), "\n")
      cat("B coefficients (persistence):", round(B_guess, 4), "\n")
    }
  }
  
  # Validate initial parameters
  validate_gas_params_natural(initial_params, K)
  
  # Ensure initial parameters are within bounds
  initial_params <- pmax(initial_params, bounds$lower)
  initial_params <- pmin(initial_params, bounds$upper)
  
  # Test the likelihood function with initial parameters to catch early errors
  if (verbose) {
    cat("Testing initial parameters...\n")
  }
  
  test_result <- tryCatch({
    Rfiltering.single.natural_GAS(
      par = initial_params,
      y = y,
      B_burnin = B_burnin,
      C = C,
      n_nodes = n_nodes,
      scaling_method = scaling_method
    )
  }, error = function(e) {
    stop(paste("Initial parameters failed likelihood calculation:", e$message))
  })
  
  if (verbose) {
    cat("Initial negative log-likelihood:", round(test_result, 4), "\n")
  }
  
  # Optimize parameters using L-BFGS-B
  if (verbose) {
    cat("Starting L-BFGS-B optimization...\n")
    opt_start_time <- Sys.time()
  }
  
  # Try multiple starting points for robustness
  n_starts <- 3  # Try 3 different starting points
  best_result <- NULL
  best_likelihood <- Inf
  
  for (start_i in 1:n_starts) {
    if (verbose && n_starts > 1) {
      cat("Trying starting point", start_i, "of", n_starts, "...\n")
    }
    
    # Create slightly different starting points
    if (start_i == 1) {
      current_start <- initial_params
    } else {
      # Add small random perturbations for additional starts
      set.seed(start_i * 100)
      perturbation <- c(
        rnorm(K, 0, 0.1 * y_sd),           # Small mu perturbations
        rnorm(K, 0, 0.1 * y_sd^2),         # Small sigma2 perturbations  
        rnorm(n_transition, 0, 0.05),      # Small transition prob perturbations
        rnorm(n_transition, 0, 0.005),     # Tiny A perturbations
        rnorm(n_transition, 0, 0.01)       # Small B perturbations
      )
      current_start <- initial_params + perturbation
      
      # Ensure bounds are respected
      current_start <- pmax(current_start, bounds$lower)
      current_start <- pmin(current_start, bounds$upper)
    }
    
    # Test this starting point
    test_result <- tryCatch({
      Rfiltering.single.natural_GAS(
        par = current_start,
        y = y,
        B_burnin = B_burnin,
        C = C,
        n_nodes = n_nodes,
        scaling_method = scaling_method
      )
    }, error = function(e) {
      if (verbose) cat("Starting point", start_i, "failed initial test:", e$message, "\n")
      return(Inf)
    })
    
    if (verbose && n_starts > 1) {
      cat("Starting point", start_i, "initial likelihood:", round(test_result, 4), "\n")
    }
    
    # Only optimize if initial test passed
    if (is.finite(test_result)) {
      current_result <- tryCatch({
        optim(
          par = current_start,
          fn = Rfiltering.single.natural_GAS,
          method = "L-BFGS-B",
          lower = bounds$lower,
          upper = bounds$upper,
          y = y,
          B_burnin = B_burnin,
          C = C,
          n_nodes = n_nodes,
          scaling_method = scaling_method,
          control = list(
            maxit = max_iterations,
            factr = factr,
            trace = 0,  # Silent for multiple starts
            REPORT = 100
          )
        )
      }, error = function(e) {
        if (verbose) cat("Optimization failed for starting point", start_i, ":", e$message, "\n")
        return(NULL)
      })
      
      # Check if this is the best result so far
      if (!is.null(current_result) && is.finite(current_result$value) && current_result$value < best_likelihood) {
        best_result <- current_result
        best_likelihood <- current_result$value
        if (verbose && n_starts > 1) {
          cat("New best result from starting point", start_i, "- likelihood:", round(best_likelihood, 4), "\n")
        }
      }
    }
  }
  
  # Use the best result
  if (is.null(best_result)) {
    stop("All starting points failed. Try different initial values or check your data.")
  }
  
  optimization_result <- best_result
  
  if (verbose) {
    opt_end_time <- Sys.time()
    cat("Optimization completed in", 
        format(difftime(opt_end_time, opt_start_time), digits = 4), "\n")
    cat("Optimization convergence code:", optimization_result$convergence, "\n")
    cat("Final negative log-likelihood:", optimization_result$value, "\n")
    cat("Function evaluations:", optimization_result$counts[1], "\n")
    cat("Gradient evaluations:", optimization_result$counts[2], "\n")
  }
  
  # Extract estimated parameters (already in natural space)
  estimated_params <- optimization_result$par
  
  # Validate final parameters
  validate_gas_params_natural(estimated_params, K)
  
  # Extract different parameter components
  mu_est <- estimated_params[1:K]
  sigma2_est <- estimated_params[(K+1):(2*K)]
  init_trans_est <- estimated_params[(2*K+1):(2*K+n_transition)]
  A_est <- estimated_params[(2*K+n_transition+1):(2*K+2*n_transition)]
  B_est <- estimated_params[(2*K+2*n_transition+1):(2*K+3*n_transition)]
  
  # Calculate model diagnostics
  num_params <- length(estimated_params)
  num_data_points <- length(y) - B_burnin - C
  
  aic <- 2 * optimization_result$value + 2 * num_params
  bic <- 2 * optimization_result$value + num_params * log(num_data_points)
  
  # Calculate filtered probabilities and other model outputs
  full_likelihood <- tryCatch({
    Rfiltering_GAS(
      mu_est, sigma2_est, init_trans_est, A_est, B_est, y, B_burnin, C, n_nodes, scaling_method
    )
  }, error = function(e) {
    stop(paste("Failed to calculate final model outputs with estimated parameters:", e$message))
  })
  
  filtered_probs <- attr(full_likelihood, "X.t")
  transition_probs <- attr(full_likelihood, "p_trans")
  scaled_scores <- attr(full_likelihood, "score_scaled")
  gas_diagnostics_raw <- attr(full_likelihood, "gas_diagnostics")
  
  # Calculate parameter standard errors using numerical Hessian
  if (verbose) {
    cat("Calculating standard errors...\n")
  }
  
  # Approximate Hessian using finite differences
  hessian_result <- tryCatch({
    numDeriv::hessian(
      func = Rfiltering.single.natural_GAS,
      x = estimated_params,
      y = y,
      B_burnin = B_burnin,
      C = C,
      n_nodes = n_nodes,
      scaling_method = scaling_method
    )
  }, error = function(e) {
    warning("Failed to compute Hessian: ", e$message)
    return(NULL)
  })
  
  # Calculate standard errors from Hessian
  standard_errors <- if (!is.null(hessian_result)) {
    tryCatch({
      sqrt(diag(solve(hessian_result)))
    }, error = function(e) {
      warning("Failed to invert Hessian for standard errors: ", e$message)
      rep(NA, num_params)
    })
  } else {
    rep(NA, num_params)
  }
  
  # Prepare results
  results <- list(
    parameters = list(
      mu = mu_est,
      sigma2 = sigma2_est,
      init_trans = init_trans_est,
      A = A_est,
      B = B_est
    ),
    standard_errors = list(
      mu = standard_errors[1:K],
      sigma2 = standard_errors[(K+1):(2*K)],
      init_trans = standard_errors[(2*K+1):(2*K+n_transition)],
      A = standard_errors[(2*K+n_transition+1):(2*K+2*n_transition)],
      B = standard_errors[(2*K+2*n_transition+1):(2*K+3*n_transition)]
    ),
    diagnostics = list(
      loglik = -optimization_result$value,
      aic = aic,
      bic = bic,
      num_params = num_params,
      num_data_points = num_data_points
    ),
    optimization = optimization_result,
    hessian = hessian_result,
    filtered_probabilities = filtered_probs,
    transition_probabilities = transition_probs,
    scaled_scores = scaled_scores,
    gas_diagnostics = gas_diagnostics_raw,
    model_info = list(
      type = "GAS",
      K = K,
      B_burnin = B_burnin,
      C = C,
      n_nodes = n_nodes,
      scaling_method = scaling_method,
      optimization_method = "L-BFGS-B"
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
    
    # GAS-specific diagnostics
    if (!is.null(gas_diagnostics_raw)) {
      cat("GAS diagnostics - Mean Fisher Info:", round(gas_diagnostics_raw$mean_fisher_info, 4), "\n")
      cat("GAS diagnostics - Mean Score Norm:", round(gas_diagnostics_raw$mean_score_norm, 4), "\n")
      cat("GAS diagnostics - Zero Score %:", round(gas_diagnostics_raw$zero_score_percentage, 1), "%\n")
    }
    
    # Convergence diagnostics
    convergence_msg <- switch(as.character(optimization_result$convergence),
                              "0" = "successful convergence",
                              "1" = "iteration limit reached",
                              "51" = "warning from L-BFGS-B",
                              "52" = "error from L-BFGS-B",
                              paste("convergence code", optimization_result$convergence))
    cat("Convergence:", convergence_msg, "\n")
  }
  
  return(results)
}



#' Filter observed data through GAS model with automatic constant model fallback
#'
#' @param mu Vector of means corresponding to each regime
#' @param sigma2 Vector of variances corresponding to each regime
#' @param init_trans Initial transition probabilities for the latent process
#' @param A Scale parameter for score updates (sensitivity)
#' @param B Persistence parameter for score updates (memory)
#' @param y Observed time series increments
#' @param B_burnin Burn-in to be excluded at the beginning of the time series
#' @param C Cut-off to be excluded at the end of the time series
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method (default: "moore_penrose")
#' @param A_threshold Threshold below which to use constant model (default: 1e-4)
#' @param verbose Whether to report fallback usage (default: FALSE)
#' @return Negative log-likelihood of observed data under the model
#' @details
#' Automatically falls back to constant transition probability model when A parameters
#' are effectively zero, providing both speed improvements and numerical stability.
Rfiltering_GAS_with_fallback <- function(mu, sigma2, init_trans, A, B, y, B_burnin, C, 
                                         n_nodes = 30, scaling_method = "moore_penrose",
                                         A_threshold = 1e-4, verbose = FALSE) {
  
  # Check if A parameters are effectively zero
  max_A <- max(abs(A))
  use_constant_model <- max_A < A_threshold
  
  if (use_constant_model) {
    if (verbose) {
      cat("Using constant model fallback (max|A| =", round(max_A, 6), "< threshold =", A_threshold, ")\n")
    }
    
    # Use fast constant model filtering
    result <- Rfiltering_Const(mu, sigma2, init_trans, y, B_burnin, C)
    
    # Add attributes to match GAS output format
    attr(result, "model_type") <- "constant_fallback"
    attr(result, "max_A") <- max_A
    attr(result, "A_threshold") <- A_threshold
    
    # Create dummy GAS-specific attributes for compatibility
    attr(result, "score_scaled") <- matrix(0, nrow = length(A), ncol = length(y))
    attr(result, "f") <- matrix(logit(init_trans), nrow = length(A), ncol = length(y))
    attr(result, "gas_diagnostics") <- list(
      model_type = "constant_fallback",
      max_A = max_A,
      zero_score_percentage = 100,
      mean_fisher_info = NA,
      mean_score_norm = 0
    )
    
  } else {
    if (verbose) {
      cat("Using full GAS model (max|A| =", round(max_A, 6), ">= threshold =", A_threshold, ")\n")
    }
    
    # Use full GAS filtering
    result <- Rfiltering_GAS(mu, sigma2, init_trans, A, B, y, B_burnin, C, n_nodes, scaling_method)
    
    # Add model type information
    attr(result, "model_type") <- "full_gas"
    attr(result, "max_A") <- max_A
    attr(result, "A_threshold") <- A_threshold
  }
  
  return(result)
}

#' Wrapper for likelihood calculator with automatic fallback
#'
#' @param par Natural parameters (not transformed)
#' @param y Observed time series increments
#' @param B_burnin Burn-in to be excluded at the beginning of the time series
#' @param C Cut-off to be excluded at the end of the time series
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method (default: "moore_penrose")
#' @param A_threshold Threshold below which to use constant model (default: 1e-4)
#' @return Negative log-likelihood of observed data under the model
Rfiltering.single.natural_GAS_with_fallback <- function(par, y, B_burnin, C, 
                                                        n_nodes = 30, scaling_method = "moore_penrose",
                                                        A_threshold = 1e-4) {
  
  # Determine the number of regimes
  K <- count_regime_GAS(par)
  n_transition <- K*(K-1)
  
  # Validate parameters with informative error messages
  validate_gas_params_natural(par, K)
  
  # Extract the components (already in natural space)
  mu <- par[1:K]
  sigma2 <- par[(K+1):(2*K)]
  init_trans <- par[(2*K+1):(2*K+n_transition)]
  A <- par[(2*K+n_transition+1):(2*K+2*n_transition)]
  B <- par[(2*K+2*n_transition+1):(2*K+3*n_transition)]
  
  # Check for transition probabilities that would create invalid transition matrices
  p_trans_test <- convert_to_valid_probs(init_trans, K)
  test_matrix <- transition_matrix(p_trans_test, check_validity = TRUE)
  
  # Calculate likelihood using fallback-enabled filtering
  l <- Rfiltering_GAS_with_fallback(mu, sigma2, init_trans, A, B, y, B_burnin, C, 
                                    n_nodes, scaling_method, A_threshold, verbose = FALSE)
  
  # Check for invalid likelihood
  if (!is.finite(l) || is.na(l)) {
    stop(sprintf("Likelihood calculation returned invalid value: %s. 
                 Parameters: mu=%s, sigma2=%s, A_mean=%.4f, B_mean=%.4f, model_type=%s", 
                 l, 
                 paste(round(mu, 4), collapse=","),
                 paste(round(sigma2, 4), collapse=","),
                 mean(A), mean(B),
                 attr(l, "model_type")))
  }
  
  return(l[1])
}




#' Estimate GAS model with automatic constant model fallback
#'
#' @param y Observed time series
#' @param K Number of regimes
#' @param B_burnin Burn-in period to exclude from likelihood calculation
#' @param C Cut-off period to exclude from likelihood calculation
#' @param initial_params Initial parameter guesses (optional)
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method (default: "moore_penrose")
#' @param max_iterations Maximum number of optimization iterations (default: 1000)
#' @param factr Controls precision of L-BFGS-B (default: 1e7)
#' @param A_threshold Threshold below which to use constant model (default: 1e-4)
#' @param verbose Whether to print progress information (default: TRUE)
#' @return List with estimated parameters and model diagnostics
estimate_gas_model_with_fallback <- function(y, K = 3, B_burnin = 100, C = 50, 
                                             initial_params = NULL,
                                             n_nodes = 30, scaling_method = "moore_penrose",
                                             max_iterations = 1000, factr = 1e7,
                                             A_threshold = 1e-4, verbose = TRUE) {
  if (verbose) {
    cat("Estimating GAS model with automatic fallback (A threshold =", A_threshold, ")\n")
    start_time <- Sys.time()
  }
  
  # Number of transition probabilities
  n_transition <- K*(K-1)
  
  # Get data characteristics for bounds
  y_mean <- mean(y, na.rm = TRUE)
  y_sd <- sd(y, na.rm = TRUE)
  
  # Create parameter bounds
  bounds <- create_gas_bounds_natural(K, y_mean, y_sd)
  
  # Create default initial parameters if none provided
  if (is.null(initial_params)) {
    # Create conservative initial guesses based on data characteristics
    y_mean <- mean(y, na.rm = TRUE)
    y_sd <- sd(y, na.rm = TRUE)
    
    # Use conservative initial regime means (closer to data mean)
    spread <- 1.0 * y_sd
    mu_guess <- seq(y_mean - spread, y_mean + spread, length.out = K)
    
    # Use variance closer to data variance
    sigma2_guess <- rep(y_sd^2, K)
    
    # Conservative transition probabilities and GAS parameters
    init_trans_guess <- rep(0.2, n_transition)
    A_guess <- rep(0.01, n_transition)  # Start small
    B_guess <- rep(0.98, n_transition)  # High persistence
    
    initial_params <- c(mu_guess, sigma2_guess, init_trans_guess, A_guess, B_guess)
    
    if (verbose) {
      cat("Generated conservative initial parameter guesses:\n")
      cat("Means:", round(mu_guess, 4), "\n")
      cat("Variances:", round(sigma2_guess, 4), "\n")
      cat("A coefficients (sensitivity):", round(A_guess, 4), "\n")
      cat("B coefficients (persistence):", round(B_guess, 4), "\n")
    }
  }
  
  # Validate initial parameters
  validate_gas_params_natural(initial_params, K)
  
  # Ensure initial parameters are within bounds
  initial_params <- pmax(initial_params, bounds$lower)
  initial_params <- pmin(initial_params, bounds$upper)
  
  # Test the likelihood function with initial parameters
  if (verbose) {
    cat("Testing initial parameters...\n")
  }
  
  test_result <- tryCatch({
    Rfiltering.single.natural_GAS_with_fallback(
      par = initial_params,
      y = y,
      B_burnin = B_burnin,
      C = C,
      n_nodes = n_nodes,
      scaling_method = scaling_method,
      A_threshold = A_threshold
    )
  }, error = function(e) {
    stop(paste("Initial parameters failed likelihood calculation:", e$message))
  })
  
  if (verbose) {
    cat("Initial negative log-likelihood:", round(test_result, 4), "\n")
  }
  
  # Try multiple starting points for robustness
  n_starts <- 3
  best_result <- NULL
  best_likelihood <- Inf
  
  for (start_i in 1:n_starts) {
    if (verbose && n_starts > 1) {
      cat("Trying starting point", start_i, "of", n_starts, "...\n")
    }
    
    # Create slightly different starting points
    if (start_i == 1) {
      current_start <- initial_params
    } else {
      # Add small random perturbations for additional starts
      set.seed(start_i * 100)
      perturbation <- c(
        rnorm(K, 0, 0.1 * y_sd),           # Small mu perturbations
        rnorm(K, 0, 0.1 * y_sd^2),         # Small sigma2 perturbations  
        rnorm(n_transition, 0, 0.05),      # Small transition prob perturbations
        rnorm(n_transition, 0, 0.005),     # Tiny A perturbations
        rnorm(n_transition, 0, 0.01)       # Small B perturbations
      )
      current_start <- initial_params + perturbation
      
      # Ensure bounds are respected
      current_start <- pmax(current_start, bounds$lower)
      current_start <- pmin(current_start, bounds$upper)
    }
    
    # Test this starting point
    test_result <- tryCatch({
      Rfiltering.single.natural_GAS_with_fallback(
        par = current_start,
        y = y,
        B_burnin = B_burnin,
        C = C,
        n_nodes = n_nodes,
        scaling_method = scaling_method,
        A_threshold = A_threshold
      )
    }, error = function(e) {
      if (verbose) cat("Starting point", start_i, "failed initial test:", e$message, "\n")
      return(Inf)
    })
    
    if (verbose && n_starts > 1) {
      cat("Starting point", start_i, "initial likelihood:", round(test_result, 4), "\n")
    }
    
    # Only optimize if initial test passed
    if (is.finite(test_result)) {
      current_result <- tryCatch({
        optim(
          par = current_start,
          fn = Rfiltering.single.natural_GAS_with_fallback,
          method = "L-BFGS-B",
          lower = bounds$lower,
          upper = bounds$upper,
          y = y,
          B_burnin = B_burnin,
          C = C,
          n_nodes = n_nodes,
          scaling_method = scaling_method,
          A_threshold = A_threshold,
          control = list(
            maxit = max_iterations,
            factr = factr,
            trace = 0,  # Silent for multiple starts
            REPORT = 100
          )
        )
      }, error = function(e) {
        if (verbose) cat("Optimization failed for starting point", start_i, ":", e$message, "\n")
        return(NULL)
      })
      
      # Check if this is the best result so far
      if (!is.null(current_result) && is.finite(current_result$value) && current_result$value < best_likelihood) {
        best_result <- current_result
        best_likelihood <- current_result$value
        if (verbose && n_starts > 1) {
          cat("New best result from starting point", start_i, "- likelihood:", round(best_likelihood, 4), "\n")
        }
      }
    }
  }
  
  # Use the best result
  if (is.null(best_result)) {
    stop("All starting points failed. Try different initial values or check your data.")
  }
  
  optimization_result <- best_result
  
  if (verbose) {
    cat("Optimization completed\n")
    cat("Optimization convergence code:", optimization_result$convergence, "\n")
    cat("Final negative log-likelihood:", optimization_result$value, "\n")
    cat("Function evaluations:", optimization_result$counts[1], "\n")
    cat("Gradient evaluations:", optimization_result$counts[2], "\n")
  }
  
  # Extract estimated parameters (already in natural space)
  estimated_params <- optimization_result$par
  
  # Validate final parameters
  validate_gas_params_natural(estimated_params, K)
  
  # Extract different parameter components
  mu_est <- estimated_params[1:K]
  sigma2_est <- estimated_params[(K+1):(2*K)]
  init_trans_est <- estimated_params[(2*K+1):(2*K+n_transition)]
  A_est <- estimated_params[(2*K+n_transition+1):(2*K+2*n_transition)]
  B_est <- estimated_params[(2*K+2*n_transition+1):(2*K+3*n_transition)]
  
  # Calculate model diagnostics
  num_params <- length(estimated_params)
  num_data_points <- length(y) - B_burnin - C
  
  aic <- 2 * optimization_result$value + 2 * num_params
  bic <- 2 * optimization_result$value + num_params * log(num_data_points)
  
  # Calculate filtered probabilities and other model outputs
  full_likelihood <- tryCatch({
    Rfiltering_GAS_with_fallback(
      mu_est, sigma2_est, init_trans_est, A_est, B_est, y, B_burnin, C, 
      n_nodes, scaling_method, A_threshold, verbose = TRUE
    )
  }, error = function(e) {
    stop(paste("Failed to calculate final model outputs with estimated parameters:", e$message))
  })
  
  filtered_probs <- attr(full_likelihood, "X.t")
  transition_probs <- attr(full_likelihood, "p_trans")
  scaled_scores <- attr(full_likelihood, "score_scaled")
  gas_diagnostics_raw <- attr(full_likelihood, "gas_diagnostics")
  model_type <- attr(full_likelihood, "model_type")
  
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
      loglik = -optimization_result$value,
      aic = aic,
      bic = bic,
      num_params = num_params,
      num_data_points = num_data_points
    ),
    optimization = optimization_result,
    filtered_probabilities = filtered_probs,
    transition_probabilities = transition_probs,
    scaled_scores = scaled_scores,
    gas_diagnostics = gas_diagnostics_raw,
    model_info = list(
      type = "GAS_with_fallback",
      actual_model_used = model_type,
      K = K,
      B_burnin = B_burnin,
      C = C,
      n_nodes = n_nodes,
      scaling_method = scaling_method,
      A_threshold = A_threshold,
      optimization_method = "L-BFGS-B",
      max_A_estimated = max(abs(A_est))
    )
  )
  
  if (verbose) {
    end_time <- Sys.time()
    cat("Total estimation time:", 
        format(difftime(end_time, start_time), digits = 4), "\n")
    cat("Final model type used:", model_type, "\n")
    cat("Max |A| estimated:", round(max(abs(A_est)), 6), "\n")
    cat("AIC:", aic, "BIC:", bic, "\n")
    cat("Estimated means:", round(mu_est, 4), "\n")
    cat("Estimated variances:", round(sigma2_est, 4), "\n")
    cat("Estimated A coefficients (mean):", round(mean(A_est), 6), "\n")
    cat("Estimated B coefficients (mean):", round(mean(B_est), 4), "\n")
    
    # GAS-specific diagnostics
    if (!is.null(gas_diagnostics_raw)) {
      if (model_type == "constant_fallback") {
        cat("Used constant model fallback - no GAS score calculations performed\n")
      } else {
        cat("GAS diagnostics - Mean Fisher Info:", round(gas_diagnostics_raw$mean_fisher_info, 4), "\n")
        cat("GAS diagnostics - Mean Score Norm:", round(gas_diagnostics_raw$mean_score_norm, 4), "\n")
        cat("GAS diagnostics - Zero Score %:", round(gas_diagnostics_raw$zero_score_percentage, 1), "%\n")
      }
    }
  }
  
  return(results)
}
