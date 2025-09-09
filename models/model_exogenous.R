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
#' @param par Parameter vector with attributes (mu, sigma2, init_trans, A)
#' @param X_Exo Exogenous process that drives transition probability changes
#' @param burn_in Number of burn-in observations to discard (default: 100)
#' @return Matrix of simulated data with M rows and N columns
#' @details 
#' Automatically detects configuration (diagonal vs off-diagonal, equal variances)
#' from parameter attributes.
#' 
#' Simulates data from a regime switching model where transition probabilities
#' depend on values of an exogenous process. The relationship is defined as:
#' f[t+1] = omega + A * X_Exo[t], followed by logistic transformation.
#'
#' @examples
#' # Generate data using diagonal parameterization
#' par_diag <- c(-1, 1, 0.5, 0.6, 0.8, 0.9, 0.1, 0.2)  # mu, sigma2, p11, p22, A1, A2
#' par_diag <- set_parameter_attributes(par_diag, K=2, model_type="exogenous", 
#'                                      diag_probs=TRUE, equal_variances=FALSE)
#' X_Exo <- rnorm(1000)
#' data_sim <- dataTVPXExoCD(10, 1000, par_diag, X_Exo)
dataTVPXExoCD <- function(M, N, par, X_Exo, burn_in = 100) {
  
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
  
  # Handle equal variances: expand single variance to K variances for simulation
  if (equal_variances && length(sigma2) == 1) {
    sigma2 <- rep(sigma2, K)
  }
  
  # Check that X_Exo is long enough
  total_length <- N + burn_in
  if (length(X_Exo) < total_length) {
    stop("Error: Exogenous variable series is too short for the requested simulation length.")
  }
  
  # Determine the number of transition parameters based on parameterization
  n_transition <- length(init_trans)
  
  # Set up a matrix to save the output
  data <- matrix(0, M, N)
  
  # Get baseline transition probabilities that are logit-transformed
  omega_LR <- logit(init_trans)
  omega <- omega_LR * (1 - A)  # This scaling ensures that when A=0, we get back to init_trans
  
  for (i in 1:M) {
    # Initialize all data structures including burn-in period
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
    
    # Set initial values using the first exogenous variable value
    f[,1] <- omega + A * X_Exo[1]
    p_trans_raw <- logistic(f[,1])
    p_trans[,1] <- convert_to_valid_probs(p_trans_raw, diag_probs = diag_probs)
    
    for (t in 1:(total_length-1)) {
      # Generate predicted probabilities
      Pmatrix <- transition_matrix(p_trans[,t], diag_probs = diag_probs, check_validity = FALSE)
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
      # This is the key feature of the exogenous model
      f[,t+1] <- omega + A * X_Exo[t+1]
      p_trans_raw <- logistic(f[,t+1])
      p_trans[,t+1] <- convert_to_valid_probs(p_trans_raw, diag_probs = diag_probs)
    }
    
    # For the last time point
    Pmatrix <- transition_matrix(p_trans[,total_length-1], diag_probs = diag_probs, check_validity = FALSE)
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
#' @param par Parameter vector with attributes (mu, sigma2, init_trans, A)
#' @param X_Exo Exogenous process that drives transition probability changes
#' @param y Observed time series increments
#' @param B Burn-in to be excluded at the beginning of the time series
#' @param C Cut-off to be excluded at the end of the time series
#' @param diagnostics If TRUE, include detailed diagnostic information (default: TRUE)
#' @return Negative log-likelihood of observed data under the model
#' @details
#' Automatically reads configuration (diagonal vs off-diagonal, equal variances)
#' from parameter attributes.
#' 
#' Filters observed data through the model to compute the likelihood.
#' Returns the negative log-likelihood for compatibility with optimization functions.
#' 
#' The diagnostics parameter controls whether detailed information is attached
#' as attributes (expensive during optimization, useful for final results).
#'
#' @examples
#' # Filter data using diagonal parameterization
#' par_diag <- c(-1, 1, 0.5, 0.6, 0.8, 0.9, 0.1, 0.2)
#' par_diag <- set_parameter_attributes(par_diag, K=2, model_type="exogenous",
#'                                      diag_probs=TRUE, equal_variances=FALSE)
#' X_Exo <- rnorm(1000)
#' y <- rnorm(1000)
#' loglik <- Rfiltering_TVPXExo(par_diag, X_Exo, y, 100, 50)
Rfiltering_TVPXExo <- function(par, X_Exo, y, B, C, diagnostics = FALSE) {
  
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
  
  # Handle equal variances: expand single variance to K variances for filtering
  if (equal_variances && length(sigma2) == 1) {
    sigma2 <- rep(sigma2, K)
  }
  
  # Determine length of the time series
  M <- length(y)
  
  # Check that X_Exo is long enough
  if (length(X_Exo) < M) {
    stop("Error: Exogenous variable series is too short for the observed data.")
  }
  
  # Determine the number of transition parameters based on parameterization
  n_transition <- length(init_trans)
  
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
  
  # Set initial values using the first exogenous variable value
  f[,1] <- omega + A * X_Exo[1]
  p_trans_raw <- logistic(f[,1])
  p_trans[,1] <- convert_to_valid_probs(p_trans_raw, diag_probs = diag_probs)
  
  # Initial state probabilities (using stationary distribution)
  initial_probs <- tryCatch({
    stat_dist(p_trans[,1], diag_probs = diag_probs)
  }, error = function(e) {
    rep(1/K, K)  # Fallback to uniform
  })
  X_t[,1] <- initial_probs
  
  for (t in 1:(M-1)) {
    # Generate predicted probabilities
    Pmatrix <- transition_matrix(p_trans[,t], diag_probs = diag_probs, check_validity = FALSE)
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
    # This is the key feature of the exogenous model
    f[,t+1] <- omega + A * X_Exo[t+1]
    p_trans_raw <- logistic(f[,t+1])
    p_trans[,t+1] <- convert_to_valid_probs(p_trans_raw, diag_probs = diag_probs)
  }
  
  # Calculate likelihood for the last time point
  Pmatrix <- transition_matrix(p_trans[,M-1], diag_probs = diag_probs, check_validity = FALSE)
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
    attr(neg_log_lik, "log_lik_values") <- log_lik_values
    attr(neg_log_lik, "valid_indices") <- valid_indices
    attr(neg_log_lik, "model_info") <- list(
      K = K,
      diag_probs = diag_probs,
      equal_variances = equal_variances,
      model_type = "exogenous",
      n_transition = n_transition
    )
  }
  
  return(neg_log_lik)
}

#' Estimate exogenous regime-switching model
#'
#' @param y Observed time series data
#' @param X_Exo Exogenous process that drives transition probability changes
#' @param K Number of regimes
#' @param diag_probs If TRUE, use diagonal transition probability parameterization
#' @param equal_variances If TRUE, constrain all regimes to have equal variances
#' @param n_starts Number of random starting points for optimization (default: 10)
#' @param B Burn-in observations to exclude (default: 100)
#' @param C Cut-off observations to exclude (default: 50)
#' @param bounds Optional list with lower and upper parameter bounds
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
#' # Estimate model with diagonal probabilities
#' y <- rnorm(1000)
#' X_Exo <- rnorm(1000)
#' result_diag <- estimate_exo_model(y, X_Exo, K=2, diag_probs=TRUE, n_starts=5)
#' 
#' # Estimate model with off-diagonal probabilities
#' result_offdiag <- estimate_exo_model(y, X_Exo, K=2, diag_probs=FALSE, n_starts=5)
estimate_exo_model <- function(y, X_Exo, K, diag_probs = FALSE, equal_variances = FALSE,
                               n_starts = 10, B = 100, C = 50, bounds = NULL,
                               parallel = TRUE, cores = NULL, seed = NULL, verbose = 1) {
  
  # Input validation
  if (!is.numeric(y) || length(y) == 0) {
    stop("y must be a non-empty numeric vector")
  }
  if (!is.numeric(X_Exo) || length(X_Exo) == 0) {
    stop("X_Exo must be a non-empty numeric vector")
  }
  if (length(X_Exo) < length(y)) {
    stop("X_Exo must be at least as long as y")
  }
  if (!is.numeric(K) || K < 2 || K != as.integer(K)) {
    stop("K must be an integer >= 2")
  }
  if (length(y) <= B + C + K) {
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
  
  if (verbose >= 1) {
    cat("Estimating exogenous regime-switching model\n")
    cat("==========================================\n")
    cat("K:", K, "regimes\n")
    cat("Data points:", length(y), "(using", length(y) - B - C, "after burn-in/cut-off)\n")
    cat("Exogenous variable length:", length(X_Exo), "\n")
    cat("Parameterization:", ifelse(diag_probs, "diagonal", "off-diagonal"), "transition probabilities\n")
    cat("Variances:", ifelse(equal_variances, "equal (shared)", "separate"), "\n")
    cat("Starting points:", n_starts, "\n")
    if (parallel) cat("Parallel processing:", cores, "cores\n")
    cat("\n")
  }
  
  # Generate diverse starting points using updated function
  if (verbose >= 1) cat("Generating starting points...\n")
  starting_points <- generate_starting_points(
    y = y,
    K = K,
    model_type = "exogenous",
    n_starts = n_starts,
    diag_probs = diag_probs,
    equal_variances = equal_variances,
    seed = seed
  )
  
  # Create parameter bounds if none provided
  if (is.null(bounds)) {
    # Calculate expected parameter count
    expected_length <- calculate_expected_length(K, "exogenous", diag_probs, equal_variances)
    
    # Create default bounds for transformed parameters
    n_mu <- K
    n_sigma2 <- ifelse(equal_variances, 1, K)
    n_trans <- ifelse(diag_probs, K, K*(K-1))
    n_A <- n_trans  # A coefficients match transition structure
    
    lower_bounds <- c(rep(-Inf, n_mu),      # No bounds on means
                      rep(-Inf, n_sigma2),  # Log-variances
                      rep(-Inf, n_trans),   # Logit-probabilities  
                      rep(-Inf, n_A))       # Logit-A coefficients
    
    upper_bounds <- c(rep(Inf, n_mu),       # No bounds on means
                      rep(Inf, n_sigma2),   # Log-variances
                      rep(Inf, n_trans),    # Logit-probabilities
                      rep(Inf, n_A))        # Logit-A coefficients
    
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
      
      # Run optimization
      trace_setting <- if (verbose >= 2) 1 else 0
      optimization_result <- nlminb(
        start = transformed_params,
        objective = function(par_t) {
          # Transform parameters back to natural space with proper attributes
          par_t_with_attrs <- transformed_params  # Copy attributes from start
          par_t_with_attrs[] <- par_t  # Update values
          attr(par_t_with_attrs, "parameterization") <- "transformed"
          
          par_natural <- untransform_parameters(par_t_with_attrs)
          
          # Calculate negative log-likelihood (no diagnostics for speed)
          neg_log_lik <- Rfiltering_TVPXExo(par_natural, X_Exo, y, B, C, diagnostics = FALSE)
          return(neg_log_lik)
        },
        lower = bounds$lower,
        upper = bounds$upper,
        control = list(eval.max = 1e6, iter.max = 1e6, trace = trace_setting)
      )
      
      # Transform final parameters back to natural space
      final_par_t <- transformed_params  # Copy attributes
      final_par_t[] <- optimization_result$par  # Update values
      attr(final_par_t, "parameterization") <- "transformed"
      
      estimated_params <- untransform_parameters(final_par_t)
      
      # Store the final parameters in the optimization result
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
  
  if (verbose >= 1 && n_starts > 1) {
    cat("Best result from starting point", best_result$start_idx, "\n")
    cat("Convergence summary:", n_converged, "converged,", 
        n_starts - n_converged - n_failed, "non-convergent,", n_failed, "failed\n")
    cat("Best negative log-likelihood:", sprintf("%.6f", best_result$objective), "\n")
  }
  
  # Extract the best optimization result
  optimization_result <- best_result$result
  estimated_params <- optimization_result$final_par
  
  # Calculate model diagnostics
  num_params <- length(estimated_params)
  num_data_points <- length(y) - B - C
  
  aic <- 2 * optimization_result$objective + 2 * num_params
  bic <- 2 * optimization_result$objective + num_params * log(num_data_points)
  
  # Calculate filtered probabilities and additional diagnostics (with full diagnostics)
  full_likelihood_result <- Rfiltering_TVPXExo(estimated_params, X_Exo, y, B, C, diagnostics = TRUE)
  
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
      n_failed = n_failed
    ),
    model_info = list(
      K = K,
      model_type = "exogenous",
      diag_probs = diag_probs,
      equal_variances = equal_variances,
      parameterization = "natural"
    ),
    data_info = list(
      n_obs = length(y),
      burn_in = B,
      cut_off = C,
      n_used = num_data_points,
      exo_length = length(X_Exo)
    ),
    filtered_probabilities = attr(full_likelihood_result, "X.t"),
    time_varying_probs = attr(full_likelihood_result, "p_trans"),
    f_values = attr(full_likelihood_result, "f"),
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
  
  if (verbose >= 1) {
    cat("\nEstimation completed successfully!\n")
    cat("Final log-likelihood:", sprintf("%.6f", -optimization_result$objective), "\n")
    cat("AIC:", sprintf("%.2f", aic), "BIC:", sprintf("%.2f", bic), "\n")
  }
  
  return(results_list)
}
