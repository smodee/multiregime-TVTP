#' Constant Transition Probability Regime-Switching Models
#' 
#' This file implements a regime-switching model with constant transition probabilities
#' between regimes. UPDATED to support both diagonal and off-diagonal transition
#' probability parameterizations using the new attribute-based parameter system.
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
#' @param par Parameter vector with attributes (mu, sigma2, trans_prob)
#' @param burn_in Number of burn-in observations to discard (default: 100)
#' @return Matrix of simulated data with M rows and N columns
#' @details 
#' Automatically detects configuration (diagonal vs off-diagonal, equal variances)
#' from parameter attributes.
#' 
#' Simulates data from a regime switching model where transition probabilities
#' remain constant over time.
#'
#' @examples
#' # Generate data using diagonal parameterization (original style)
#' par_diag <- c(-1, 1, 0.5, 0.5, 0.8, 0.9)  # mu1, mu2, sig1, sig2, p11, p22
#' par_diag <- set_parameter_attributes(par_diag, K=2, model_type="constant", 
#'                                      diag_probs=TRUE, equal_variances=FALSE)
#' data_sim <- dataConstCD(10, 1000, par_diag)
#' 
#' # Generate data using off-diagonal parameterization (new style)
#' par_offdiag <- c(-1, 1, 0.5, 0.5, 0.2, 0.1)  # mu1, mu2, sig1, sig2, p12, p21
#' par_offdiag <- set_parameter_attributes(par_offdiag, K=2, model_type="constant",
#'                                         diag_probs=FALSE, equal_variances=FALSE)
#' data_sim <- dataConstCD(10, 1000, par_offdiag)
dataConstCD <- function(M, N, par, burn_in = 100) {
  
  # Validate parameter vector and extract configuration
  validate_parameter_attributes(par)
  
  K <- attr(par, "K")
  diag_probs <- attr(par, "diag_probs")
  equal_variances <- attr(par, "equal_variances")
  
  # Extract parameter components using attribute-based extraction
  mu <- extract_parameter_component(par, "mu")
  sigma2 <- extract_parameter_component(par, "sigma2")
  trans_prob <- extract_parameter_component(par, "trans_prob")
  
  # Handle equal variances: expand single variance to K variances for simulation
  if (equal_variances && length(sigma2) == 1) {
    sigma2 <- rep(sigma2, K)
  }
  
  # Create transition matrix using the appropriate parameterization
  transition_mat <- transition_matrix(trans_prob, diag_probs = diag_probs, check_validity = TRUE)
  
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
    start_dist <- tryCatch({
      stat_dist(transition_mat)
    }, error = function(e) {
      warning("Could not calculate stationary distribution, using uniform")
      rep(1/K, K)
    })
    X_t[,1] <- start_dist

    # Initialize X_tlag[,1] using stationary distribution and transition matrix
    # For row-stochastic P, prediction is: X_pred = P^T %*% X_prev
    X_tlag[,1] <- as.vector(t(transition_mat) %*% start_dist)

    for (t in 1:(total_length-1)) {
      # Generate predicted probabilities using the constant transition matrix
      # Use transpose because P is row-stochastic: P[i,j] = P(to j | from i)
      X_tlag[,t] <- as.vector(t(transition_mat) %*% X_t[,t])
      
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
    X_tlag[,total_length] <- as.vector(t(transition_mat) %*% X_t[,total_length-1])
    S[total_length] <- sample(1:K, 1, prob=X_tlag[,total_length])
    y.sim[total_length] <- rnorm(1, mu[S[total_length]], sqrt(sigma2[S[total_length]]))
    
    # Remove burn-in and save the simulation run in the data matrix
    data[i,] <- y.sim[(burn_in+1):total_length]
  }
  
  return(data)
}

#' Filter observed data through the constant transition probability regime-switching model
#'
#' @param par Parameter vector with attributes (mu, sigma2, trans_prob)
#' @param y Observed time series increments
#' @param B Burn-in to be excluded at the beginning of the time series
#' @param C Cut-off to be excluded at the end of the time series
#' @param diagnostics Choose whether to store calculation diagnostics or not
#' @return Negative log-likelihood of observed data under the model
#' @details
#' Automatically reads configuration (diagonal vs off-diagonal, equal variances)
#' from parameter attributes.
#' 
#' Filters observed data through the model to compute the likelihood.
#' Returns the negative log-likelihood for compatibility with optimization functions.
#'
#' @examples
#' # Filter data using diagonal parameterization
#' par_diag <- c(-1, 1, 0.5, 0.5, 0.8, 0.9)
#' par_diag <- set_parameter_attributes(par_diag, K=2, model_type="constant",
#'                                      diag_probs=TRUE, equal_variances=FALSE)
#' y <- rnorm(1000)
#' loglik <- Rfiltering_Const(par_diag, y, 100, 50)
Rfiltering_Const <- function(par, y, B, C, diagnostics = FALSE) {
  
  # Validate parameter vector and extract configuration
  if (diagnostics) {
    validate_parameter_attributes(par)
  }
  
  K <- attr(par, "K")
  diag_probs <- attr(par, "diag_probs")
  equal_variances <- attr(par, "equal_variances")
  
  # Extract parameter components using attribute-based extraction
  mu <- extract_parameter_component(par, "mu")
  sigma2 <- extract_parameter_component(par, "sigma2")
  trans_prob <- extract_parameter_component(par, "trans_prob")
  
  # Handle equal variances: expand single variance to K variances for filtering
  if (equal_variances && length(sigma2) == 1) {
    sigma2 <- rep(sigma2, K)
  }
  
  # Determine length of the time series
  M <- length(y)
  
  # Create transition matrix using the appropriate parameterization
  transition_mat <- transition_matrix(trans_prob, diag_probs = diag_probs, check_validity = TRUE)

  # Pre-compute transpose for prediction step (avoids repeated transpose in loop)
  transition_mat_T <- t(transition_mat)

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

  # Initial state probabilities (stationary distribution or uniform)
  start_dist <- tryCatch({
    stat_dist(transition_mat)
  }, error = function(e) {
    warning("Could not calculate stationary distribution, using uniform")
    rep(1/K, K)
  })
  X_t[,1] <- start_dist

  # Initialize X_tlag[,1] using stationary distribution and transition matrix
  # For row-stochastic P, prediction is: X_pred = P^T %*% X_prev
  X_tlag[,1] <- as.vector(transition_mat_T %*% start_dist)

  for (t in 1:(M-1)) {
    # Generate predicted probabilities using the constant transition matrix
    # Use pre-computed transpose for efficiency
    X_tlag[,t] <- as.vector(transition_mat_T %*% X_t[,t])

    # Note: eta[,t] already computed in vectorized pre-computation above
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
  X_tlag[,M] <- as.vector(transition_mat_T %*% X_t[,M-1])
  # Note: eta[,M] already computed in vectorized pre-computation above
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
  
  if (diagnostics) {
    # Store additional information as attributes for diagnostics
    attr(neg_log_lik, "X.t") <- X_t
    attr(neg_log_lik, "X.tlag") <- X_tlag
    attr(neg_log_lik, "eta") <- eta
    attr(neg_log_lik, "tot.lik") <- tot_lik
    attr(neg_log_lik, "transition_matrix") <- transition_mat
    attr(neg_log_lik, "log_lik_values") <- log_lik_values
    attr(neg_log_lik, "valid_indices") <- valid_indices
    attr(neg_log_lik, "model_info") <- list(
      K = K,
      diag_probs = diag_probs,
      equal_variances = equal_variances,
      model_type = "constant"
    )
  }
  
  return(neg_log_lik)
}

#' Estimate constant transition probability regime-switching model
#'
#' @param y Observed time series data
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
#' # Estimate model with diagonal probabilities (original style)
#' y <- rnorm(1000)
#' result_diag <- estimate_constant_model(y, K=2, diag_probs=TRUE, n_starts=5)
#' 
#' # Estimate model with off-diagonal probabilities (new style)  
#' result_offdiag <- estimate_constant_model(y, K=2, diag_probs=FALSE, n_starts=5)
estimate_constant_model <- function(y, K, diag_probs = TRUE, equal_variances = FALSE,
                                    n_starts = 10, B = 100, C = 50, bounds = NULL,
                                    parallel = TRUE, cores = NULL, seed = NULL, verbose = 1) {
  
  # Input validation
  if (!is.numeric(y) || length(y) == 0) {
    stop("y must be a non-empty numeric vector")
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
    cat("Estimating constant transition probability model\n")
    cat("==============================================\n")
    cat("K:", K, "regimes\n")
    cat("Data points:", length(y), "(using", length(y) - B - C, "after burn-in/cut-off)\n")
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
    model_type = "constant",
    n_starts = n_starts,
    diag_probs = diag_probs,
    equal_variances = equal_variances,
    seed = seed
  )
  
  # Create parameter bounds if none provided
  if (is.null(bounds)) {
    # Calculate expected parameter count
    expected_length <- calculate_expected_length(K, "constant", diag_probs, equal_variances)
    
    # Create default bounds for transformed parameters
    lower_bounds <- c(rep(-Inf, K),                 # No bounds on means
                      rep(-Inf, ifelse(equal_variances, 1, K)),  # Log-variances
                      rep(-Inf, ifelse(diag_probs, K, K*(K-1))))  # Logit-probabilities
    
    upper_bounds <- c(rep(Inf, K),                  # No bounds on means  
                      rep(Inf, ifelse(equal_variances, 1, K)),   # Log-variances
                      rep(Inf, ifelse(diag_probs, K, K*(K-1))))  # Logit-probabilities
    
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
          
          # Calculate negative log-likelihood
          neg_log_lik <- Rfiltering_Const(par_natural, y, B, C)
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
  
  # Calculate filtered probabilities and additional diagnostics
  full_likelihood_result <- Rfiltering_Const(estimated_params, y, B, C)
  
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
      model_type = "constant",
      diag_probs = diag_probs,
      equal_variances = equal_variances,
      parameterization = "natural"
    ),
    data_info = list(
      n_obs = length(y),
      burn_in = B,
      cut_off = C,
      n_used = num_data_points
    ),
    filtered_probabilities = attr(full_likelihood_result, "X.t"),
    transition_matrix = attr(full_likelihood_result, "transition_matrix"),
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
  
  if (verbose >= 1) {
    cat("\nEstimation completed successfully!\n")
    cat("Final log-likelihood:", sprintf("%.6f", -optimization_result$objective), "\n")
    cat("AIC:", sprintf("%.2f", aic), "BIC:", sprintf("%.2f", bic), "\n")
  }
  
  return(results_list)
}
