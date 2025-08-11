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
                      n_nodes = 30, scaling_method = "simple", quad_sample_size = 1000) {
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
                           n_nodes = 30, scaling_method = "simple",
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
    if (abs(X_tlag_sum - 1) > 0.1) {  # Warn if error > 10%
      warning(paste("Large probability normalization at time", t, 
                    "- sum was", round(X_tlag_sum, 4), "instead of 1.0"))
    }
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
#' @param scaling_method Score scaling method (default: "simple")
#' @return Negative log-likelihood of observed data under the model
#' @details
#' Transforms parameters from unconstrained space back to natural space
#' and calls the filtering function to calculate the likelihood.
#'
#' @examples
#' # Optimize parameters for a 3-regime model
#' transformed_params <- transform_parameters(c(mu, sigma2, init_trans, A, B), "gas")
#' result <- nlminb(transformed_params, Rfiltering.single.trasf_GAS, 
#'                 y = y, B_burnin = 100, C = 50)
Rfiltering.single.trasf_GAS <- function(par_t, y, B_burnin, C, 
                                        n_nodes = 30, scaling_method = "simple") {
  # Transform parameters back to original parameter space
  par <- untransform_parameters(par_t, "gas")
  
  # Determine the number of regimes
  K <- count_regime(par, "gas")
  n_transition <- K*(K-1)
  
  # Extract the components
  mu <- mean_from_par(par, "gas")
  sigma2 <- sigma2_from_par(par, "gas")
  init_trans <- transp_from_par(par, "gas")
  A <- A_from_par(par, "gas")
  B <- B_from_par(par, "gas")
  
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
#' @param initial_params Initial parameter guesses (optional, only used when n_starts=1)
#' @param bounds Parameter bounds for optimization (optional)
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method (default: "simple")
#' @param use_fallback Whether to automatically use constant model fallback when A is small (default: TRUE)
#' @param A_threshold Threshold below which to use constant model fallback (default: 1e-4)
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
#' When use_fallback=TRUE, automatically falls back to constant transition 
#' probability model when A parameters are effectively zero, providing both 
#' speed improvements and numerical stability.
#' 
#' Uses the future package for cross-platform parallel processing that works
#' on Windows, macOS, and Linux. The output format is identical regardless 
#' of single or multi-start estimation.
#'
#' @examples
#' # Single start (current behavior)
#' results <- estimate_gas_model(data, K=3)
#' 
#' # Multi-start with parallel processing
#' results <- estimate_gas_model(data, K=3, n_starts=10, parallel=TRUE)
#' 
#' # Multi-start sequential (for debugging)
#' results <- estimate_gas_model(data, K=3, n_starts=5, parallel=FALSE)
#' 
#' # Multi-start without fallback
#' results <- estimate_gas_model(data, K=3, n_starts=10, parallel=TRUE, use_fallback=FALSE)
estimate_gas_model <- function(y, K = 3, B_burnin = 100, C = 50, 
                               initial_params = NULL, bounds = NULL,
                               n_nodes = 30, scaling_method = "simple",
                               use_fallback = TRUE, A_threshold = 1e-4,
                               n_starts = 1, parallel = FALSE, cores = NULL,
                               seed = NULL, verbose = TRUE) {
  
  # Set up cores (don't use more cores than starts)
  if (is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 1)
  }
  cores <- min(cores, n_starts)
  
  if (verbose) {
    cat("Estimating GAS model with", K, "regimes")
    if (use_fallback) {
      cat(" (fallback enabled, A threshold =", A_threshold, ")")
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
    starting_points <- generate_starting_points(y, K, "gas", n_starts, seed)
  }
  
  # Create default bounds if none provided
  if (is.null(bounds)) {
    n_transition <- K * (K - 1)
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
  
  # Define the optimization function for a single start
  optimize_single_start <- function(start_idx) {
    start_params <- starting_points[[start_idx]]
    
    if (verbose && !parallel) {
      cat("  Starting point", start_idx, "of", n_starts, "\n")
    }
    
    tryCatch({
      # Transform parameters to unconstrained space for optimization
      transformed_params <- transform_parameters(start_params, "gas")
      
      # Run optimization using the wrapper function that handles GAS-specific parameters
      trace_setting <- if (verbose > 1) 1 else 0
      optimization_result <- nlminb(
        start = transformed_params,
        objective = function(par_t, y, B_burnin, C, n_nodes, scaling_method, use_fallback, A_threshold) {
          # Transform parameters back to original parameter space
          par <- untransform_parameters(par_t, "gas")
          
          # Extract the components
          mu <- mean_from_par(par, "gas")
          sigma2 <- sigma2_from_par(par, "gas")
          init_trans <- transp_from_par(par, "gas")
          A <- A_from_par(par, "gas")
          B <- B_from_par(par, "gas")
          
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
  estimated_params <- untransform_parameters(optimization_result$par, "gas")
  
  # Extract different parameter components
  mu_est <- mean_from_par(estimated_params, "gas")
  sigma2_est <- sigma2_from_par(estimated_params, "gas")
  init_trans_est <- transp_from_par(estimated_params, "gas")
  A_est <- A_from_par(estimated_params, "gas")
  B_est <- B_from_par(estimated_params, "gas")
  
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
  
  # Prepare results (same format as original function)
  results_list <- list(
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
  
  return(results_list)
}
