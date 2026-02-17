#' Generalized Autoregressive Score Models with Time-Varying Transition Probabilities
#' 
#' This file implements a regime-switching model where transition probabilities
#' are driven by the score of the predictive likelihood function (GAS dynamics).
#'
#' Reference: Bazzi, M., Blasques, F., Koopman, S.J., & Lucas, A. (2017).
#' Time-Varying Transition Probabilities for Markov Regime Switching Models.
#' Journal of Time Series Analysis, 38(3), 458-478.

#' Generate data from a score-driven regime-switching model (GAS)
#'
#' @param M Number of simulation runs to be performed
#' @param N Length of the simulation runs (discretized time)
#' @param par Parameter vector with attributes (mu, sigma2, init_trans, A, B)
#' @param burn_in Number of burn-in observations to discard (default: 100)
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method ("moore_penrose", "simple", "normalized", or "factored")
#' @param quad_sample_size Sample size for creating representative data for quadrature setup (default: 1000)
#' @return Matrix of simulated data with M rows and N columns
#' @details 
#' Automatically detects configuration (diagonal vs off-diagonal, equal variances)
#' from parameter attributes.
#' 
#' Simulates data from a regime switching model where transition probabilities
#' are updated using score-driven dynamics based on the predictive likelihood.
#' The update equation is f\[t+1\] = omega + A*s\[t\] + B*(f\[t\] - omega).
#' 
#' This implementation uses proper GAS score scaling following Bazzi et al. (2017)
#' with Gauss-Hermite quadrature and Moore-Penrose pseudo-inverse scaling.
#'
#' @examples
#' # Generate data using diagonal parameterization
#' par_diag <- c(-1, 1, 0.5, 0.6, 0.8, 0.9, 0.1, 0.2, 0.8, 0.9)  # mu, sigma2, p11, p22, A, B
#' par_diag <- set_parameter_attributes(par_diag, K=2, model_type="gas",
#'                                      diag_probs=TRUE, equal_variances=FALSE)
#' data_sim <- dataGASCD(3, 200, par_diag)
#' @export
dataGASCD <- function(M, N, par, burn_in = 100, n_nodes = 30,
                      scaling_method = NULL, quad_sample_size = 1000) {
  
  # Validate parameter vector and extract configuration
  validate_parameter_attributes(par)
  
  K <- attr(par, "K")
  diag_probs <- attr(par, "diag_probs")
  equal_variances <- attr(par, "equal_variances")
  
  # Extract parameter components using attribute-based extraction
  mu <- extract_parameter_component(par, "mu")
  sigma2 <- extract_parameter_component(par, "sigma2")
  init_trans <- extract_parameter_component(par, "trans_prob")
  A <- extract_parameter_component(par, "A")
  B <- extract_parameter_component(par, "B")
  
  # Handle equal variances: expand single variance to K variances for simulation
  if (equal_variances && length(sigma2) == 1) {
    sigma2 <- rep(sigma2, K)
  }
  
  # Determine the number of transition parameters based on parameterization
  n_transition <- length(init_trans)

  # Setup Gauss-Hermite quadrature for score calculation
  # Use statmod::gauss.quad.prob to match the original C code exactly

  # Create representative sample from regime parameters (matching original approach)
  regime_sample <- c(rnorm(quad_sample_size/2, mu[1], sqrt(sigma2[1])),
                     rnorm(quad_sample_size/2, mu[K], sqrt(sigma2[K])))
  sample_median <- median(regime_sample)
  sample_sd <- sd(regime_sample)

  GQ <- statmod::gauss.quad.prob(n_nodes, "normal", mu = sample_median, sigma = sample_sd)

  # Create gh_setup object compatible with calculate_gas_score
  gh_setup <- list(
    nodes = GQ$nodes,
    weights = GQ$weights,
    n_nodes = n_nodes,
    mu = sample_median,
    sigma = sample_sd,
    mu_quad = sample_median,
    sigma_quad = sample_sd,
    method = "statmod_exact"
  )
  class(gh_setup) <- c("gauss_hermite_setup", "list")
  
  # Set up a matrix to save the output
  data <- matrix(0, M, N)
  
  # Get baseline transition parameters in f-space
  if (diag_probs) {
    params_to_f <- logit
    f_to_params <- logistic
    f_to_params_loop <- logistic_clamped
  } else {
    params_to_f <- logit
    f_to_params <- logistic
    f_to_params_loop <- logistic_clamped
  }
  omega <- params_to_f(init_trans)

  for (i in 1:M) {
    # Initialize all data structures including burn-in period
    total_length <- N + burn_in
    eta <- matrix(0, nrow=K, ncol=total_length)     # Likelihood of each regime
    tot_lik <- numeric(total_length)                # Total likelihood
    X_t <- matrix(0, nrow=K, ncol=total_length)     # Filtered probabilities after observation
    X_tlag <- matrix(0, nrow=K, ncol=total_length)  # Predicted probabilities before observation
    S <- numeric(total_length)                      # Latent state
    y.sim <- numeric(total_length)                  # Random increments according to state
    
    # Initialize GAS-specific variables
    f <- matrix(0, nrow=n_transition, ncol=total_length)
    p_trans <- matrix(0, nrow=n_transition, ncol=total_length)
    score_scaled <- matrix(0, nrow=n_transition, ncol=total_length)
    fisher_info_series <- numeric(total_length)
    score_norms <- numeric(total_length)
    
    # Initial state probabilities (uniform distribution)
    X_t[,1] <- rep(1/K, K)
    
    # Set initial f values
    f[,1] <- omega
    # Use clamped logistic to match original HMMGAS C implementation
    # This maps f directly to [1e-10, 1-1e-10] for numerical stability
    p_trans[,1] <- f_to_params_loop(f[,1])
    
    # Initialize score as zero
    score_scaled[,1] <- 0
    
    # Track zero scores for diagnostics
    zero_score_count <- 0
    
    for (t in 1:(total_length-1)) {
      # Generate predicted probabilities using generalized formula (works for any K>=2)
      X_tlag[, t] <- compute_predicted_probs(X_t[, t], p_trans[, t], diag_probs = diag_probs)

      # Sample a state based on the predicted probabilities and
      # simulate data conditional on that state
      S[t] <- sample(1:K, 1, prob=X_tlag[, t])
      y.sim[t] <- rnorm(1, mu[S[t]], sqrt(sigma2[S[t]]))
      
      # Calculate likelihoods
      for (k in 1:K) {
        eta[k,t] <- dnorm(y.sim[t], mu[k], sqrt(sigma2[k]))
      }
      tot_lik[t] <- sum(eta[,t]*X_tlag[,t])
      
      # Calculate filtered probabilities
      X_t[,t+1] <- (eta[,t]*X_tlag[,t])/tot_lik[t]
      # Warn if drift is pathological, then re-normalize
      if (abs(sum(X_t[,t+1]) - 1) > 1e-4) {
        warning(sprintf("Filtered probabilities at t=%d sum to %.6f, re-normalizing", t+1, sum(X_t[,t+1])))
      }
      X_t[,t+1] <- pmax(X_t[,t+1], .Machine$double.eps)
      X_t[,t+1] <- X_t[,t+1] / sum(X_t[,t+1])

      # Calculate GAS score for time t
      score_result <- calculate_gas_score(
        y_obs = y.sim[t],
        mu = mu,
        sigma2 = sigma2,
        trans_prob = p_trans[,t],
        diag_probs = diag_probs,
        X_pred = X_tlag[,t],
        gh_setup = gh_setup,
        scaling_method = scaling_method
      )
      
      # Store score and diagnostics
      score_scaled[,t+1] <- score_result$scaled_score
      fisher_info_series[t] <- score_result$fisher_info
      score_norms[t] <- sqrt(sum(score_result$raw_score^2))
      
      # Count zero scores for diagnostics
      if (all(abs(score_result$scaled_score) < 1e-10)) {
        zero_score_count <- zero_score_count + 1
      }
      
      # Update f values using GAS dynamics: f[t+1] = omega + A*s[t] + B*(f[t] - omega)
      f[,t+1] <- omega + A * score_scaled[,t+1] + B * (f[,t] - omega)
      
      # Convert f to transition parameters
      p_trans[,t+1] <- f_to_params_loop(f[,t+1])
    }

    # For the last time point - use generalized formula
    X_tlag[, total_length] <- compute_predicted_probs(X_t[, total_length-1], p_trans[, total_length-1], diag_probs = diag_probs)
    S[total_length] <- sample(1:K, 1, prob=X_tlag[, total_length])
    y.sim[total_length] <- rnorm(1, mu[S[total_length]], sqrt(sigma2[S[total_length]]))
    
    # Remove burn-in and save the simulation run in the data matrix
    data[i,] <- y.sim[(burn_in+1):total_length]
    
    # Optional warning about excessive zero scores
    if (zero_score_count > 0.5 * (total_length - 1)) {
      warning(paste("Simulation", i, "used zero scores for", zero_score_count, 
                    "out of", total_length-1, "time points (", 
                    round(100*zero_score_count/(total_length-1), 1), "%)"))
    }
  }
  
  # Store simulation metadata
  attr(data, "simulation_info") <- list(
    M = M,
    N = N,
    K = K,
    burn_in = burn_in,
    diag_probs = diag_probs,
    equal_variances = equal_variances,
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
#' @param par Parameter vector with attributes (mu, sigma2, init_trans, A, B)
#' @param y Observed time series increments
#' @param B_burnin Burn-in to be excluded at the beginning of the time series
#' @param C Cut-off to be excluded at the end of the time series
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method ("moore_penrose", "simple", "normalized", or "factored")
#' @param use_fallback Whether to automatically use constant model fallback when A is small (default: TRUE)
#' @param A_threshold Threshold below which to use constant model fallback (default: 1e-4)
#' @param diagnostics If TRUE, include detailed diagnostic information (default: FALSE)
#' @param verbose Whether to report fallback usage (default: FALSE)
#' @return Negative log-likelihood of observed data under the model
#' @details
#' Automatically reads configuration (diagonal vs off-diagonal, equal variances)
#' from parameter attributes.
#' 
#' Filters observed data through the model to compute the likelihood using proper
#' GAS score scaling following Bazzi et al. (2017). When use_fallback=TRUE,
#' automatically falls back to constant transition probability model when A parameters
#' are effectively zero, providing both speed improvements and numerical stability.
#' 
#' The diagnostics parameter controls whether detailed information is attached
#' as attributes (expensive during optimization, useful for final results).
#'
#' @examples
#' # Filter data using diagonal parameterization
#' par_diag <- c(-1, 1, 0.5, 0.6, 0.8, 0.9, 0.1, 0.2, 0.8, 0.9)
#' par_diag <- set_parameter_attributes(par_diag, K=2, model_type="gas",
#'                                      diag_probs=TRUE, equal_variances=FALSE)
#' y <- rnorm(200)
#' loglik <- Rfiltering_GAS(par_diag, y, 20, 10)
#' @export
Rfiltering_GAS <- function(par, y, B_burnin, C, n_nodes = 30, scaling_method = NULL,
                           use_fallback = TRUE, A_threshold = 1e-4, diagnostics = FALSE, verbose = FALSE) {
  
  # Only validate if diagnostics are enabled (performance optimization)
  if (diagnostics) {
    validate_parameter_attributes(par)
  }
  
  K <- attr(par, "K")
  diag_probs <- attr(par, "diag_probs")
  equal_variances <- attr(par, "equal_variances")
  
  # Extract parameter components using attribute-based extraction
  mu <- extract_parameter_component(par, "mu")
  sigma2 <- extract_parameter_component(par, "sigma2")
  init_trans <- extract_parameter_component(par, "trans_prob")
  A <- extract_parameter_component(par, "A")
  B <- extract_parameter_component(par, "B")
  
  # Handle equal variances: expand single variance to K variances for filtering
  if (equal_variances && length(sigma2) == 1) {
    sigma2 <- rep(sigma2, K)
  }
  
  # Check fallback condition if enabled
  if (use_fallback) {
    max_A <- max(abs(A))
    use_constant_model <- max_A < A_threshold
    
    if (use_constant_model) {
      if (verbose) {
        cat("Using constant model fallback (max|A| =", round(max_A, 6), "< threshold =", A_threshold, ")\n")
      }
      
      # Create constant model parameter vector with proper attributes
      # Use original (unexpanded) sigma2 to match equal_variances attribute
      original_sigma2 <- extract_parameter_component(par, "sigma2")
      const_par <- c(mu, original_sigma2, init_trans)
      const_par <- set_parameter_attributes(
        par = const_par,
        K = K,
        model_type = "constant",
        diag_probs = diag_probs,
        equal_variances = equal_variances
      )
      
      # Use fast constant model filtering
      result <- Rfiltering_Const(const_par, y, B_burnin, C, diagnostics = diagnostics)
      
      # Add GAS-specific attributes if diagnostics enabled
      if (diagnostics) {
        attr(result, "model_type") <- "constant_fallback"
        attr(result, "max_A") <- max_A
        attr(result, "A_threshold") <- A_threshold
        attr(result, "score_scaled") <- matrix(0, nrow = length(A), ncol = length(y))
        f_init <- if (diag_probs) logit(init_trans) else init_trans
        attr(result, "f") <- matrix(f_init, nrow = length(A), ncol = length(y))
        attr(result, "gas_diagnostics") <- list(
          model_type = "constant_fallback",
          max_A = max_A,
          zero_score_percentage = 100,
          mean_fisher_info = NA,
          mean_score_norm = 0
        )
      }
      
      return(result)
    }
    
    if (verbose) {
      cat("Using full GAS model (max|A| =", round(max_A, 6), ">= threshold =", A_threshold, ")\n")
    }
  }
  
  # Full GAS model implementation
  M <- length(y)
  n_transition <- length(init_trans)

  # Setup Gauss-Hermite quadrature for score calculation
  # CRITICAL: Use statmod::gauss.quad.prob to match the original C code exactly
  # The C code receives nodes/weights from R via gauss.quad.prob(30, "normal", mu, sigma)
  # where mu = median(representative_sample), sigma = sd(representative_sample)

  # Create representative sample from regime parameters (matching original approach)
  regime_sample <- c(rnorm(500, mu[1], sqrt(sigma2[1])),
                     rnorm(500, mu[K], sqrt(sigma2[K])))
  sample_median <- median(regime_sample)
  sample_sd <- sd(regime_sample)

  # Use statmod::gauss.quad.prob for exact match with C code
  GQ <- statmod::gauss.quad.prob(n_nodes, "normal", mu = sample_median, sigma = sample_sd)

  # Create gh_setup object compatible with calculate_gas_score
  gh_setup <- list(
    nodes = GQ$nodes,
    weights = GQ$weights,
    n_nodes = n_nodes,
    mu = sample_median,
    sigma = sample_sd,
    mu_quad = sample_median,
    sigma_quad = sample_sd,
    method = "statmod_exact"
  )
  class(gh_setup) <- c("gauss_hermite_setup", "list")
  
  # Initialize variables
  tot_lik <- numeric(M)                # Total likelihood
  X_t <- matrix(0, nrow=K, ncol=M)     # Filtered probabilities after observation
  X_tlag <- matrix(0, nrow=K, ncol=M)  # Predicted probabilities before observation

  # =========================================================================
  # PRE-COMPUTE EMISSION LIKELIHOODS (vectorized for performance)
  # =========================================================================
  eta <- matrix(0, nrow=K, ncol=M)
  sqrt_sigma2 <- sqrt(sigma2)
  for (k in 1:K) {
    eta[k, ] <- dnorm(y, mu[k], sqrt_sigma2[k])
  }

  # Initialize GAS-specific variables
  f <- matrix(0, nrow=n_transition, ncol=M)
  p_trans <- matrix(0, nrow=n_transition, ncol=M)
  score_scaled <- matrix(0, nrow=n_transition, ncol=M)
  fisher_info_series <- numeric(M)
  score_norms <- numeric(M)
  
  # Get baseline transition parameters in f-space
  # NOTE: For GAS model, omega = omega_LR (no (1-A) factor like TVP/Exogenous)
  if (diag_probs) {
    params_to_f <- logit
    f_to_params <- logistic
    f_to_params_loop <- logistic_clamped
  } else {
    params_to_f <- logit
    f_to_params <- logistic
    f_to_params_loop <- logistic_clamped
  }
  omega <- params_to_f(init_trans)
  omega_LR <- omega  # They are the same for GAS

  # =============================================================================
  # INITIALIZATION (t=1): Generalized for K>=2
  # =============================================================================
  # Set f[,1] = omega_LR
  f[, 1] <- omega_LR

  # Compute transition probabilities at t=1
  p_trans[, 1] <- f_to_params(f[, 1])

  # Compute stationary distribution using eigenvalue method (works for any K>=2)
  stationary_probs <- stat_dist(p_trans[, 1], diag_probs = diag_probs)

  # Compute initial predicted probabilities: X_tlag[,1] = P^T * stationary
  X_tlag[, 1] <- compute_initial_predicted_probs(p_trans[, 1], diag_probs = diag_probs)

  # Note: eta[,1] already computed in vectorized pre-computation above
  tot_lik[1] <- sum(eta[, 1] * X_tlag[, 1])

  # Protect against numerical issues at t=1
  if (tot_lik[1] <= 0 || is.na(tot_lik[1])) {
    tot_lik[1] <- .Machine$double.eps
    X_t[, 1] <- rep(1/K, K)
  } else {
    # Calculate filtered probabilities X_t[,1]
    X_t[, 1] <- (eta[, 1] * X_tlag[, 1]) / tot_lik[1]
    # Warn if drift is pathological, then re-normalize
    if (abs(sum(X_t[, 1]) - 1) > 1e-4) {
      warning(sprintf("Filtered probabilities at t=1 sum to %.6f, re-normalizing", sum(X_t[, 1])))
    }
    X_t[, 1] <- pmax(X_t[, 1], .Machine$double.eps)
    X_t[, 1] <- X_t[, 1] / sum(X_t[, 1])
  }

  # Initialize score as zero
  score_scaled[, 1] <- 0

  # Track zero scores for diagnostics
  zero_score_count <- 0

  # =============================================================================
  # Calculate score for t=1 and set f[,2]
  # =============================================================================
  # For the GAS model, at t=1 we use stationary probabilities for the score calculation
  # (the stationary_probs were already computed above for any K>=2)

  # Calculate GAS score for time t=1 using STATIONARY probs for g-vector
  # CRITICAL: The C code uses:
  #   - logLik[0] = log(eta[0]*X_tlag[0]+eta[1]*X_tlag[1]) for S denominator
  #   - g[0] = p1[0]*p11*(1-p11), g[1] = -p2[0]*p22*(1-p22) for g-vector (uses stationary probs)
  # So we pass:
  #   - X_pred = stationary_probs (for g-vector calculation)
  #   - tot_lik_external = tot_lik[1] (computed using X_tlag, for S denominator)
  score_result_1 <- calculate_gas_score(
    y_obs = y[1],
    mu = mu,
    sigma2 = sigma2,
    trans_prob = p_trans[, 1],
    diag_probs = diag_probs,
    X_pred = stationary_probs,  # Use STATIONARY probs for g-vector
    gh_setup = gh_setup,
    scaling_method = scaling_method,
    tot_lik_external = tot_lik[1]  # Use tot_lik based on X_tlag for S denominator
  )

  score_scaled[, 2] <- score_result_1$scaled_score
  fisher_info_series[1] <- score_result_1$fisher_info
  score_norms[1] <- sqrt(sum(score_result_1$raw_score^2))

  if (all(abs(score_result_1$scaled_score) < 1e-10)) {
    zero_score_count <- zero_score_count + 1
  }

  # Set f[,2] using GAS dynamics
  if (M >= 2) {
    f[, 2] <- omega + A * score_scaled[, 2] + B * (f[, 1] - omega)
    # Use clamped logistic for p_trans[,2] to match C code (lines 125-126)
    p_trans[, 2] <- f_to_params_loop(f[, 2])
  }

  # =============================================================================
  # Main loop: t = 2 to M - Generalized for K>=2
  # =============================================================================
  for (t in 2:M) {
    # Generate predicted probabilities using generalized formula
    X_tlag[, t] <- compute_predicted_probs(X_t[, t-1], p_trans[, t], diag_probs = diag_probs)

    # Note: eta[,t] already computed in vectorized pre-computation above
    tot_lik[t] <- sum(eta[, t] * X_tlag[, t])

    # Protect against numerical issues
    if (tot_lik[t] <= 0 || is.na(tot_lik[t])) {
      tot_lik[t] <- .Machine$double.eps
      X_t[, t] <- rep(1/K, K)
    } else {
      # Calculate filtered probabilities
      X_t[, t] <- (eta[, t] * X_tlag[, t]) / tot_lik[t]
      # Warn if drift is pathological, then re-normalize
      if (abs(sum(X_t[, t]) - 1) > 1e-4) {
        warning(sprintf("Filtered probabilities at t=%d sum to %.6f, re-normalizing", t, sum(X_t[, t])))
      }
      X_t[, t] <- pmax(X_t[, t], .Machine$double.eps)
      X_t[, t] <- X_t[, t] / sum(X_t[, t])
    }

    # Only compute score and update f/p_trans if NOT the last time point
    # C code: if((t<(T[0]-1))) { ... }
    if (t < M) {
      # Calculate GAS score for time t
      # CRITICAL: The C code uses two different probability sets:
      #   - S denominator uses exp(logLik[t]) = sum(eta * X_tlag) (predicted probs)
      #   - g-vector uses X_t[t-1] (filtered probs from t-1)
      # Original C code (line 169-171):
      #   S = (eta[t*2]-eta[t*2+1])/exp(logLik[t])   <- Uses logLik based on X_tlag
      #   g[t*2]   = X_t[(t-1)*2]*p11[t]*(1-p11[t])   <- Uses X_t from t-1
      #   g[t*2+1] = -X_t[(t-1)*2+1]*p22[t]*(1-p22[t])
      score_result <- calculate_gas_score(
        y_obs = y[t],
        mu = mu,
        sigma2 = sigma2,
        trans_prob = p_trans[, t],
        diag_probs = diag_probs,
        X_pred = X_t[, t-1],  # Use FILTERED probs from t-1 for g-vector
        gh_setup = gh_setup,
        scaling_method = scaling_method,
        tot_lik_external = tot_lik[t]  # Use tot_lik based on X_tlag for S denominator
      )

      # Store score and diagnostics
      score_scaled[, t+1] <- score_result$scaled_score
      fisher_info_series[t] <- score_result$fisher_info
      score_norms[t] <- sqrt(sum(score_result$raw_score^2))

      # Count zero scores for diagnostics
      if (all(abs(score_result$scaled_score) < 1e-10)) {
        zero_score_count <- zero_score_count + 1
      }

      # Update f values using GAS dynamics: f[t+1] = omega + A*s[t] + B*(f[t] - omega)
      f[, t+1] <- omega + A * score_scaled[, t+1] + B * (f[, t] - omega)

      # Convert f to transition probabilities using clamped logistic
      p_trans[, t+1] <- f_to_params_loop(f[, t+1])
    }
  }

  # Protect against numerical issues at the last time point
  if (tot_lik[M] <= 0 || is.na(tot_lik[M])) {
    tot_lik[M] <- .Machine$double.eps
  }
  
  # Sum log-likelihoods, but exclude burn-in and cut-off
  valid_indices <- (B_burnin+1):(M-C)
  if (length(valid_indices) <= 0) {
    stop("Error: No valid data points after applying burn-in and cut-off.")
  }
  
  log_lik_values <- log(tot_lik[valid_indices])
  
  # Handle any remaining invalid values
  log_lik_values[!is.finite(log_lik_values)] <- log(.Machine$double.eps)
  
  # Calculate total negative log-likelihood
  neg_log_lik <- -sum(log_lik_values)
  
  # Store additional information as attributes for diagnostics (if requested)
  if (diagnostics) {
    attr(neg_log_lik, "X.t") <- X_t
    attr(neg_log_lik, "X.tlag") <- X_tlag
    attr(neg_log_lik, "eta") <- eta
    attr(neg_log_lik, "tot.lik") <- tot_lik
    attr(neg_log_lik, "f") <- f
    attr(neg_log_lik, "p_trans") <- p_trans
    attr(neg_log_lik, "score_scaled") <- score_scaled
    attr(neg_log_lik, "log_lik_values") <- log_lik_values
    attr(neg_log_lik, "valid_indices") <- valid_indices
    attr(neg_log_lik, "model_info") <- list(
      K = K,
      diag_probs = diag_probs,
      equal_variances = equal_variances,
      model_type = "gas",
      n_transition = n_transition
    )
    
    # GAS-specific diagnostics
    valid_diag_indices <- intersect(1:(M-1), valid_indices)
    attr(neg_log_lik, "gas_diagnostics") <- list(
      model_type = "full_gas",
      zero_score_percentage = 100 * zero_score_count / length(valid_diag_indices),
      mean_fisher_info = mean(fisher_info_series[valid_diag_indices], na.rm = TRUE),
      mean_score_norm = mean(score_norms[valid_diag_indices], na.rm = TRUE)
    )
  }
  
  return(neg_log_lik)
}

#' Estimate GAS regime-switching model
#'
#' @param y Observed time series data
#' @param K Number of regimes
#' @param diag_probs If TRUE, use diagonal transition probability parameterization
#' @param equal_variances If TRUE, constrain all regimes to have equal variances
#' @param n_starts Number of random starting points for optimization (default: 10)
#' @param B_burnin Burn-in observations to exclude (default: 100)
#' @param C Cut-off observations to exclude (default: 50)
#' @param bounds Optional list with lower and upper parameter bounds
#' @param n_nodes Number of Gauss-Hermite quadrature nodes (default: 30)
#' @param scaling_method Score scaling method (default: "simple")
#' @param use_fallback Whether to automatically use constant model fallback when A is small (default: TRUE)
#' @param A_threshold Threshold below which to use constant model fallback (default: 1e-4)
#' @param early_stopping Enable early stopping for diverging starts (default: FALSE)
#' @param early_stop_patience Evaluations without improvement before stopping (default: 500)
#' @param early_stop_max_evals Maximum evaluations per start (default: 50000)
#' @param parallel Enable parallel processing for multiple starts (default: TRUE)
#' @param cores Number of cores for parallel processing (default: future::availableCores()-1)
#' @param seed Random seed for reproducibility (optional)
#' @param verbose Verbosity level (0=silent, 1=basic, 2=detailed)
#' @return List with estimation results including parameters, diagnostics, and metadata
#' @details
#' UPDATED to support both diagonal and off-diagonal transition probability parameterizations
#' using the new attribute-based parameter system. This enables exact validation against
#' the original simulation.R implementation when diag_probs=TRUE.
#'
#' @examples
#' \donttest{
#' # Estimate model with diagonal probabilities
#' y <- rnorm(200)
#' result_diag <- estimate_gas_model(y, K=2, diag_probs=TRUE, n_starts=3,
#'                                   B_burnin=20, C=10)
#'
#' # Estimate model with off-diagonal probabilities
#' result_offdiag <- estimate_gas_model(y, K=2, diag_probs=FALSE, n_starts=3,
#'                                      B_burnin=20, C=10)
#' }
#' @export
estimate_gas_model <- function(y, K, diag_probs = TRUE, equal_variances = FALSE,
                               n_starts = 10, B_burnin = 100, C = 50, bounds = NULL,
                               n_nodes = 30, scaling_method = NULL,
                               use_fallback = TRUE, A_threshold = 1e-4,
                               early_stopping = FALSE,
                               early_stop_patience = 500L,
                               early_stop_max_evals = 50000L,
                               parallel = TRUE, cores = NULL, seed = NULL, verbose = 1) {
  
  # Input validation
  if (!is.numeric(y) || length(y) == 0) {
    stop("y must be a non-empty numeric vector")
  }
  if (!is.numeric(K) || K < 2 || K != as.integer(K)) {
    stop("K must be an integer >= 2")
  }
  if (length(y) <= B_burnin + C + K) {
    stop("Time series too short for specified burn-in and cut-off")
  }
  
  # Setup parallel processing
  if (parallel) {
    if (is.null(cores)) {
      cores <- min(n_starts, future::availableCores() - 1, 4)  # Reasonable default
    }
    future::plan(future::multisession, workers = cores)
  } else {
    future::plan(future::sequential)
  }
  
  # Build early stopping configuration
  early_stop_config <- create_early_stop_config(
    enabled = early_stopping,
    patience = early_stop_patience,
    max_evals = early_stop_max_evals
  )

  if (verbose >= 1) {
    cat("Estimating GAS (score-driven) regime-switching model\n")
    cat("==================================================\n")
    cat("K:", K, "regimes\n")
    cat("Data points:", length(y), "(using", length(y) - B_burnin - C, "after burn-in/cut-off)\n")
    cat("Parameterization:", ifelse(diag_probs, "diagonal", "off-diagonal"), "transition probabilities\n")
    cat("Variances:", ifelse(equal_variances, "equal (shared)", "separate"), "\n")
    cat("Starting points:", n_starts, "\n")
    effective_scaling <- if (is.null(scaling_method)) {
      if (diag_probs) "factored (auto)" else "simple (auto)"
    } else {
      scaling_method
    }
    cat("GAS settings: scaling_method =", effective_scaling, ", n_nodes =", n_nodes, "\n")
    cat("Fallback enabled:", use_fallback, "(threshold =", A_threshold, ")\n")
    if (early_stopping) {
      cat("Early stopping: enabled (patience=", early_stop_patience,
          ", max_evals=", early_stop_max_evals, ")\n")
    }
    if (parallel) cat("Parallel processing:", cores, "cores\n")
    cat("\n")
  }
  
  # Generate diverse starting points using updated function
  if (verbose >= 1) cat("Generating starting points...\n")
  starting_points <- generate_starting_points(
    y = y,
    K = K,
    model_type = "gas",
    n_starts = n_starts,
    diag_probs = diag_probs,
    equal_variances = equal_variances,
    seed = seed
  )
  
  # Create parameter bounds if none provided
  if (is.null(bounds)) {
    # Calculate expected parameter count
    expected_length <- calculate_expected_length(K, "gas", diag_probs, equal_variances)
    
    # Create default bounds for transformed parameters
    n_mu <- K
    n_sigma2 <- ifelse(equal_variances, 1, K)
    n_trans <- ifelse(diag_probs, K, K*(K-1))
    n_A <- n_trans  # A coefficients match transition structure
    n_B <- n_trans  # B coefficients match transition structure
    
    lower_bounds <- c(rep(-Inf, n_mu),           # No bounds on means
                      rep(log(1e-10), n_sigma2), # Log-variances: floor at sigma2=1e-10
                      rep(-Inf, n_trans),        # Logit-probabilities
                      rep(-Inf, n_A),            # Logit-A coefficients
                      rep(-Inf, n_B))            # Logit-B coefficients
    
    upper_bounds <- c(rep(Inf, n_mu),       # No bounds on means
                      rep(Inf, n_sigma2),   # Log-variances
                      rep(Inf, n_trans),    # Logit-probabilities
                      rep(Inf, n_A),        # Logit-A coefficients
                      rep(Inf, n_B))        # Logit-B coefficients
    
    bounds <- list(lower = lower_bounds, upper = upper_bounds)
  }
  
  # Define the optimization function for a single start
  optimize_single_start <- function(start_idx) {
    start_params <- starting_points[[start_idx]]

    if (verbose >= 2 && !parallel) {
      cat("  Starting point", start_idx, "of", n_starts, "\n")
    }

    tryCatch({
      # Transform parameters to unconstrained space for optimization
      transformed_params <- transform_parameters(start_params)

      # Define the base objective function
      base_objective <- function(par_t) {
        par_t_with_attrs <- transformed_params
        par_t_with_attrs[] <- par_t
        attr(par_t_with_attrs, "parameterization") <- "transformed"
        par_natural <- untransform_parameters(par_t_with_attrs)
        neg_log_lik <- Rfiltering_GAS(par_natural, y, B_burnin, C,
                                      n_nodes = n_nodes, scaling_method = scaling_method,
                                      use_fallback = use_fallback, A_threshold = A_threshold,
                                      diagnostics = FALSE, verbose = FALSE)
        return(neg_log_lik)
      }

      # Wrap with early stopping if enabled
      if (early_stop_config$enabled) {
        es <- create_early_stop_objective(base_objective, early_stop_config)
        objective_fn <- es$objective
        es_tracker <- es$tracker
      } else {
        objective_fn <- base_objective
        es_tracker <- NULL
      }

      # Run optimization
      trace_setting <- if (verbose >= 2) 1 else 0
      optimization_result <- suppressWarnings(nlminb(
        start = transformed_params,
        objective = objective_fn,
        lower = bounds$lower,
        upper = bounds$upper,
        control = list(eval.max = 1e6, iter.max = 1e6, trace = trace_setting)
      ))

      # Transform final parameters back to natural space
      final_par_t <- transformed_params  # Copy attributes
      final_par_t[] <- optimization_result$par  # Update values
      attr(final_par_t, "parameterization") <- "transformed"

      estimated_params <- untransform_parameters(final_par_t)

      # Store the final parameters in the optimization result
      optimization_result$final_par <- estimated_params

      # Check early stopping status
      early_stopped <- if (!is.null(es_tracker)) es_tracker$early_stopped else FALSE
      stop_reason <- if (!is.null(es_tracker)) es_tracker$stop_reason else ""
      eval_count <- if (!is.null(es_tracker)) es_tracker$eval_count else NA_integer_

      # Return result with metadata
      list(
        result = optimization_result,
        start_idx = start_idx,
        convergence = optimization_result$convergence,
        objective = if (early_stopped) Inf else optimization_result$objective,
        start_params = start_params,
        early_stopped = early_stopped,
        stop_reason = stop_reason,
        eval_count = eval_count
      )
    }, error = function(e) {
      # Return error information
      list(
        result = NULL,
        start_idx = start_idx,
        convergence = 999,  # Error code
        objective = Inf,
        error = e$message,
        start_params = start_params,
        early_stopped = FALSE,
        stop_reason = "error",
        eval_count = NA_integer_
      )
    })
  }
  
  # Run optimization(s) using future package
  if (verbose >= 1 && n_starts > 1) {
    cat("Running optimizations...\n")
    opt_start_time <- Sys.time()
  }
  
  results <- future.apply::future_lapply(
    X = 1:n_starts,
    FUN = optimize_single_start,
    future.seed = seed
  )
  
  if (verbose >= 1 && n_starts > 1) {
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
  n_early_stopped <- sum(sapply(results, function(x) isTRUE(x$early_stopped)))

  if (verbose >= 1 && n_starts > 1) {
    cat("Best result from starting point", best_result$start_idx, "\n")
    cat("Convergence summary:", n_converged, "converged,",
        n_starts - n_converged - n_failed - n_early_stopped, "non-convergent,",
        n_early_stopped, "early-stopped,",
        n_failed, "failed\n")
    cat("Best negative log-likelihood:", sprintf("%.6f", best_result$objective), "\n")

    if (verbose >= 2 && n_early_stopped > 0) {
      es_results <- results[sapply(results, function(x) isTRUE(x$early_stopped))]
      for (r in es_results) {
        cat("  Start", r$start_idx, "early-stopped:", r$stop_reason,
            "after", r$eval_count, "evaluations\n")
      }
    }
  }

  # Extract the best optimization result
  optimization_result <- best_result$result
  estimated_params <- optimization_result$final_par

  # Calculate model diagnostics
  num_params <- length(estimated_params)
  num_data_points <- length(y) - B_burnin - C

  aic <- 2 * optimization_result$objective + 2 * num_params
  bic <- 2 * optimization_result$objective + num_params * log(num_data_points)
  
  # Calculate filtered probabilities and additional diagnostics (with full diagnostics)
  full_likelihood_result <- Rfiltering_GAS(estimated_params, y, B_burnin, C,
                                           n_nodes = n_nodes, scaling_method = scaling_method,
                                           use_fallback = use_fallback, A_threshold = A_threshold,
                                           diagnostics = TRUE, verbose = (verbose >= 2))
  
  # Extract GAS-specific diagnostics
  gas_diagnostics <- attr(full_likelihood_result, "gas_diagnostics")
  model_type_used <- if (!is.null(gas_diagnostics)) gas_diagnostics$model_type else "full_gas"
  
  # Clean up parallel processing
  if (parallel) {
    future::plan(future::sequential)
  }
  
  # Compile final results
  results_list <- list(
    parameters = estimated_params,
    optimization = optimization_result,
    diagnostics = list(
      neg_log_likelihood = optimization_result$objective,
      log_likelihood = -optimization_result$objective,
      aic = aic,
      bic = bic,
      num_parameters = num_params,
      num_data_points = num_data_points,
      convergence_code = optimization_result$convergence,
      n_starts = n_starts,
      n_converged = n_converged,
      n_early_stopped = n_early_stopped,
      n_failed = n_failed,
      early_stopping_enabled = early_stopping
    ),
    model_info = list(
      K = K,
      model_type = "gas",
      model_type_used = model_type_used,
      diag_probs = diag_probs,
      equal_variances = equal_variances,
      parameterization = "natural"
    ),
    gas_settings = list(
      n_nodes = n_nodes,
      scaling_method = scaling_method,
      use_fallback = use_fallback,
      A_threshold = A_threshold
    ),
    data_info = list(
      n_obs = length(y),
      burn_in = B_burnin,
      cut_off = C,
      n_used = num_data_points
    ),
    filtered_probabilities = attr(full_likelihood_result, "X.t"),
    time_varying_probs = attr(full_likelihood_result, "p_trans"),
    f_values = attr(full_likelihood_result, "f"),
    score_values = attr(full_likelihood_result, "score_scaled"),
    gas_diagnostics = gas_diagnostics,
    likelihood_components = list(
      log_lik_values = attr(full_likelihood_result, "log_lik_values"),
      total_likelihood = attr(full_likelihood_result, "tot.lik")
    ),
    all_results = results  # For detailed diagnostics if needed
  )
  
  # Add convenient parameter extraction
  results_list$mu_est <- extract_parameter_component(estimated_params, "mu")
  results_list$sigma2_est <- extract_parameter_component(estimated_params, "sigma2")
  results_list$init_trans_est <- extract_parameter_component(estimated_params, "trans_prob")
  results_list$A_est <- extract_parameter_component(estimated_params, "A")
  results_list$B_est <- extract_parameter_component(estimated_params, "B")
  
  if (verbose >= 1) {
    cat("\nEstimation completed successfully!\n")
    cat("Final model type used:", model_type_used, "\n")
    cat("Final log-likelihood:", sprintf("%.6f", -optimization_result$objective), "\n")
    cat("AIC:", sprintf("%.2f", aic), "BIC:", sprintf("%.2f", bic), "\n")
    
    # GAS-specific summary
    if (model_type_used == "constant_fallback") {
      cat("Used constant model fallback - A coefficients below threshold\n")
    } else if (!is.null(gas_diagnostics)) {
      cat("Max |A| estimated:", sprintf("%.6f", max(abs(results_list$A_est))), "\n")
      cat("Mean |B| estimated:", sprintf("%.4f", mean(abs(results_list$B_est))), "\n")
      cat("GAS score diagnostics:\n")
      cat("  Mean Fisher Info:", sprintf("%.4f", gas_diagnostics$mean_fisher_info), "\n")
      cat("  Mean Score Norm:", sprintf("%.4f", gas_diagnostics$mean_score_norm), "\n")
      cat("  Zero Score %:", sprintf("%.1f%%", gas_diagnostics$zero_score_percentage), "\n")
    }
  }
  
  return(results_list)
}
