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
#' @param use_fallback Whether to automatically use constant model fallback when A is small (default: TRUE)
#' @param A_threshold Threshold below which to use constant model fallback (default: 1e-4)
#' @param verbose Whether to report fallback usage (default: FALSE)
#' @return Negative log-likelihood of observed data under the model
#' @details
#' Filters observed data through the model to compute the likelihood using proper
#' GAS score scaling following Bazzi et al. (2017). When use_fallback=TRUE,
#' automatically falls back to constant transition probability model when A parameters
#' are effectively zero, providing both speed improvements and numerical stability.
#'
#' Returns the negative log-likelihood for compatibility with optimization functions.
#'
#' @examples
#' # Calculate likelihood for a 3-regime model
#' mu <- c(-2, 1, 2)
#' sigma2 <- c(0.02, 0.2, 0.6)
#' init_trans <- rep(0.2, 6)
#' A <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
#' B <- c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9)
#' y <- rnorm(1000)
#' loglik <- Rfiltering_GAS(mu, sigma2, init_trans, A, B, y, 100, 50)
Rfiltering_GAS <- function(mu, sigma2, init_trans, A, B, y, B_burnin, C, 
                           n_nodes = 30, scaling_method = "moore_penrose",
                           use_fallback = TRUE, A_threshold = 1e-4, verbose = FALSE) {
  
  # Check fallback condition if enabled
  if (use_fallback) {
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
      
      return(result)
    }
    
    if (verbose) {
      cat("Using full GAS model (max|A| =", round(max_A, 6), ">= threshold =", A_threshold, ")\n")
    }
  }

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
  
  # Convert to valid probabilities and create transition matrix
  p_trans <- convert_to_valid_probs(init_trans, K)
  
  # Set up Gauss-Hermite quadrature for score scaling
  gh_setup <- setup_gauss_hermite_quadrature(y, n_nodes)
  
  # Initialize storage for filtered probabilities and other outputs
  X_t <- matrix(0, K, M)
  X_tlag <- matrix(0, K, M)
  p_trans_series <- matrix(0, n_transition, M)
  f_series <- matrix(0, n_transition, M)
  score_scaled_series <- matrix(0, n_transition, M)
  
  # Initialize transition parameters in logit space
  f_t <- logit(init_trans)
  omega <- logit(init_trans)  # Long-run mean
  
  # Storage for GAS diagnostics
  fisher_info_series <- numeric(M)
  score_norms <- numeric(M)
  zero_score_count <- 0
  
  # Calculate log-likelihood
  logLik <- numeric(M)
  
  for (t in 1:M) {
    # Current transition probabilities
    current_trans <- logistic(f_t)
    p_trans_series[, t] <- current_trans
    f_series[, t] <- f_t
    
    # Create transition matrix and calculate predicted probabilities
    if (t == 1) {
      # Initial probabilities
      Pmatrix <- transition_matrix(current_trans, check_validity = FALSE)
      X_tlag[, t] <- solve(diag(K) - Pmatrix + matrix(1, K, K)) %*% rep(1, K)
    } else {
      Pmatrix <- transition_matrix(current_trans, check_validity = FALSE)
      X_tlag[, t] <- Pmatrix %*% X_t[, t-1]
    }
    
    # Always normalize, but warn if error is suspiciously large
    X_tlag_sum <- sum(X_tlag[, t])
    if (abs(X_tlag_sum - 1) > 0.01) {  # Warn if error > 1%
      warning(paste("Large probability normalization at time", t, 
                    "- sum was", round(X_tlag_sum, 4), "instead of 1.0"))
    }
    # Always normalize regardless
    X_tlag[, t] <- X_tlag[, t] / X_tlag_sum
    
    # Calculate regime likelihoods
    eta <- numeric(K)
    for (k in 1:K) {
      eta[k] <- dnorm(y[t], mu[k], sqrt(sigma2[k]))
    }
    
    # Calculate total likelihood and filtered probabilities
    tot_lik <- sum(eta * X_tlag[, t])
    logLik[t] <- log(tot_lik)
    
    X_t[, t] <- (eta * X_tlag[, t]) / tot_lik
    
    # Always normalize, but warn if error is suspiciously large
    X_t_sum <- sum(X_t[, t])
    if (abs(X_t_sum - 1) > 0.01) {  # Warn if error > 1%
      warning(paste("Large X_t probability normalization at time", t, 
                    "- sum was", round(X_t_sum, 4), "instead of 1.0"))
    }
    # Always normalize
    X_t[, t] <- X_t[, t] / X_t_sum
    
    # Calculate scaled score for updating f
    if (t < M) {  # Don't update after last observation
      score_scaled <- calculate_gas_score(
        y[t], mu, sigma2, X_tlag[, t], 
        if (t == 1) X_tlag[, t] else X_t[, t-1],
        current_trans, gh_setup, K, scaling_method
      )
      
      score_scaled_series[, t] <- score_scaled
      
      # Store diagnostics - safely handle missing attributes
      fisher_attr <- attr(score_scaled, "fisher_info")
      if (!is.null(fisher_attr) && length(fisher_attr) > 0) {
        fisher_info_series[t] <- fisher_attr
      } else {
        fisher_info_series[t] <- 1.0  # Default fallback value
      }
      
      if (all(abs(score_scaled) < .Machine$double.eps)) {
        zero_score_count <- zero_score_count + 1
      }
      
      # Update f using GAS dynamics: f[t+1] = omega + A*s[t] + B*(f[t] - omega)
      f_t <- omega + A * score_scaled + B * (f_t - omega)
    }
  }
  
  # Calculate final likelihood excluding burn-in and cut-off
  valid_indices <- (B_burnin + 1):(M - C)
  if (length(valid_indices) <= 0) {
    stop("Error: No valid time points remain after burn-in and cut-off.")
  }
  
  logLikSum <- sum(logLik[valid_indices])
  
  # Return negative sum of log-likelihoods (for minimizing)
  res <- -logLikSum
  
  # Store additional information as attributes
  attr(res, "X.t") <- t(X_t)
  attr(res, "X.tlag") <- t(X_tlag)
  attr(res, "p_trans") <- t(p_trans_series)
  attr(res, "f") <- t(f_series)
  attr(res, "score_scaled") <- t(score_scaled_series)
  
  # Add model type information
  attr(res, "model_type") <- "full_gas"
  if (use_fallback) {
    attr(res, "max_A") <- max(abs(A))
    attr(res, "A_threshold") <- A_threshold
  }
  
  # GAS-specific diagnostics
  valid_diag_indices <- intersect(1:(M-1), valid_indices)
  attr(res, "gas_diagnostics") <- list(
    model_type = "full_gas",
    zero_score_percentage = 100 * zero_score_count / length(valid_diag_indices),
    mean_fisher_info = mean(fisher_info_series[valid_diag_indices], na.rm = TRUE),
    mean_score_norm = mean(score_norms[valid_diag_indices], na.rm = TRUE)
  )
  
  return(res)
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
#' @param use_fallback Whether to automatically use constant model fallback when A is small (default: TRUE)
#' @param A_threshold Threshold below which to use constant model fallback (default: 1e-4)
#' @param verbose Whether to print progress information (default: TRUE)
#' @return List with estimated parameters and model diagnostics
#' @details
#' Estimates model parameters using maximum likelihood estimation.
#' When use_fallback=TRUE, automatically falls back to constant transition 
#' probability model when A parameters are effectively zero, providing both 
#' speed improvements and numerical stability.
#' Returns the estimated parameters and various diagnostics including
#' AIC, BIC, filtered probabilities, and optimization details.
#'
#' @examples
#' # Estimate a 3-regime model
#' data <- rnorm(1000)
#' results <- estimate_gas_model(data, K=3)
#' 
#' # Estimate without fallback
#' results <- estimate_gas_model(data, K=3, use_fallback=FALSE)
estimate_gas_model <- function(y, K = 3, B_burnin = 100, C = 50, 
                               initial_params = NULL, bounds = NULL,
                               n_nodes = 30, scaling_method = "moore_penrose",
                               use_fallback = TRUE, A_threshold = 1e-4,
                               verbose = TRUE) {
  if (verbose) {
    cat("Estimating GAS model with", K, "regimes")
    if (use_fallback) {
      cat(" (fallback enabled, A threshold =", A_threshold, ")")
    }
    cat("\n")
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
    objective = function(par_t, y, B_burnin, C, n_nodes, scaling_method, use_fallback, A_threshold) {
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
      l <- Rfiltering_GAS(mu, sigma2, init_trans, A, B, y, B_burnin, C, 
                          n_nodes, scaling_method, use_fallback, A_threshold, verbose = FALSE)
      return(l[1])
    },
    lower = bounds$lower,
    upper = bounds$upper,
    y = y,
    B_burnin = B_burnin,
    C = C,
    n_nodes = n_nodes,
    scaling_method = scaling_method,
    use_fallback = use_fallback,
    A_threshold = A_threshold,
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
    mu_est, sigma2_est, init_trans_est, A_est, B_est, y, B_burnin, C, 
    n_nodes, scaling_method, use_fallback, A_threshold, verbose = TRUE
  )
  
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
      actual_model_used = model_type,
      K = K,
      B_burnin = B_burnin,
      C = C,
      n_nodes = n_nodes,
      scaling_method = scaling_method,
      use_fallback = use_fallback,
      A_threshold = A_threshold
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
