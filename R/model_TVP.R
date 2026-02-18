#' Time-Varying Transition Probability Models with Autoregressive Dynamics
#' 
#' This file implements a regime-switching model where transition probabilities
#' depend on the process itself (autoregressive dynamics).
#'
#' Extended to support arbitrary K-regime models with both diagonal and
#' off-diagonal transition probability parameterizations.

#' Generate data from an autoregressive regime-switching model
#'
#' @param M Number of simulation runs to be performed
#' @param N Length of the simulation runs (discretized time)
#' @param par Parameter vector with attributes (mu, sigma2, init_trans, A)
#' @param burn_in Number of burn-in observations to discard (default: 100)
#' @return Matrix of simulated data with M rows and N columns
#' @details 
#' Automatically detects configuration (diagonal vs off-diagonal, equal variances)
#' from parameter attributes.
#' 
#' Simulates data from a regime switching model where transition probabilities
#' depend on previous values of the process itself. The relationship is defined as:
#' f\[t+1\] = omega + A * y\[t\], followed by logistic transformation.
#'
#' @examples
#' # Generate data using diagonal parameterization (original style)
#' par_diag <- c(-1, 1, 0.5, 0.5, 0.8, 0.9, 0.1, 0.1)  # mu, sigma2, p11, p22, A1, A2
#' par_diag <- set_parameter_attributes(par_diag, K=2, model_type="tvp",
#'                                      diag_probs=TRUE, equal_variances=FALSE)
#' data_sim <- dataTVPCD(3, 200, par_diag)
#'
#' # Generate data using off-diagonal parameterization (new style)
#' par_offdiag <- c(-1, 1, 0.5, 0.5, 0.2, 0.1, 0.1, 0.1)  # mu, sigma2, p12, p21, A12, A21
#' par_offdiag <- set_parameter_attributes(par_offdiag, K=2, model_type="tvp",
#'                                         diag_probs=FALSE, equal_variances=FALSE)
#' data_sim <- dataTVPCD(3, 200, par_offdiag)
#' @export
dataTVPCD <- function(M, N, par, burn_in = 100) {

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

  # Determine the number of transition parameters based on parameterization
  n_transition <- length(init_trans)

  # Set up a matrix to save the output
  data <- matrix(0, M, N)

  # Get baseline transition parameters in f-space
  # omega_LR is the long-run (unconditional) value of f
  # omega is the intercept in the AR(1) process: f[t+1] = omega + A * y[t]
  if (diag_probs) {
    params_to_f <- logit
    f_to_params <- logistic
  } else {
    params_to_f <- logit
    f_to_params <- logistic
  }
  omega_LR <- params_to_f(init_trans)
  omega <- omega_LR * (1 - A)

  for (i in 1:M) {
    # Total length including burn-in
    total_length <- N + burn_in

    # Initialize data structures
    f <- matrix(0, nrow = n_transition, ncol = total_length)
    p_trans <- matrix(0, nrow = n_transition, ncol = total_length)
    X_tlag <- matrix(0, nrow = K, ncol = total_length)
    X_t <- matrix(0, nrow = K, ncol = total_length)
    eta <- matrix(0, nrow = K, ncol = total_length)
    S <- numeric(total_length)
    y.sim <- numeric(total_length)

    # =========================================================================
    # INITIALIZATION (t=1): Generalized for K>=2
    # =========================================================================
    # Set f[1] = omega_LR (long-run value)
    f[, 1] <- omega_LR
    p_trans[, 1] <- f_to_params(f[, 1])
    if (!diag_probs) p_trans[, 1] <- clamp_offdiag_rowsums(p_trans[, 1], K)

    # Compute initial predicted probabilities using stationary distribution
    # This works for any K>=2 using eigenvalue-based stationary distribution
    X_tlag[, 1] <- compute_initial_predicted_probs(p_trans[, 1], diag_probs = diag_probs)

    # Sample first state based on X_tlag (predicted probabilities)
    S[1] <- sample(1:K, 1, prob = X_tlag[, 1])
    y.sim[1] <- rnorm(1, mu[S[1]], sqrt(sigma2[S[1]]))

    # Compute eta and filtered probabilities for t=1
    for (k in 1:K) {
      eta[k, 1] <- dnorm(y.sim[1], mu[k], sqrt(sigma2[k]))
    }
    tot_lik_1 <- sum(eta[, 1] * X_tlag[, 1])
    X_t[, 1] <- (eta[, 1] * X_tlag[, 1]) / tot_lik_1

    # =========================================================================
    # MAIN LOOP (t=1 to total_length-1): Matches original exactly
    # Original loop: for(t in 1:(N-1)) with f[t+1] = omega + A*y[t] at START
    # =========================================================================
    for (t in 1:(total_length - 1)) {
      # Update f for next time step: f[t+1] = omega + A * y[t]
      # This is done at the BEGINNING of each iteration (matching original)
      f[, t + 1] <- omega + A * y.sim[t]
      p_trans[, t + 1] <- f_to_params(f[, t + 1])
      if (!diag_probs) p_trans[, t + 1] <- clamp_offdiag_rowsums(p_trans[, t + 1], K)

      # Compute predicted probabilities X_tlag[t+1] using generalized formula
      # Works for any K>=2
      if (t == 1) {
        X_tlag[, t + 1] <- compute_predicted_probs(X_t[, 1], p_trans[, t + 1], diag_probs = diag_probs)
      } else {
        X_tlag[, t + 1] <- compute_predicted_probs(X_t[, t], p_trans[, t + 1], diag_probs = diag_probs)
      }

      # Sample state and generate observation
      S[t + 1] <- sample(1:K, 1, prob = X_tlag[, t + 1])
      y.sim[t + 1] <- rnorm(1, mu[S[t + 1]], sqrt(sigma2[S[t + 1]]))

      # Compute eta and filtered probabilities
      for (k in 1:K) {
        eta[k, t + 1] <- dnorm(y.sim[t + 1], mu[k], sqrt(sigma2[k]))
      }
      tot_lik_t <- sum(eta[, t + 1] * X_tlag[, t + 1])
      X_t[, t + 1] <- (eta[, t + 1] * X_tlag[, t + 1]) / tot_lik_t
    }

    # Remove burn-in and save the simulation run
    data[i, ] <- y.sim[(burn_in + 1):total_length]
  }

  return(data)
}

#' Filter observed data through the autoregressive regime-switching model
#'
#' @param par Parameter vector with attributes (mu, sigma2, init_trans, A)
#' @param y Observed time series increments
#' @param n_burnin Burn-in to be excluded at the beginning of the time series
#' @param n_cutoff Cut-off to be excluded at the end of the time series
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
#' par_diag <- c(-1, 1, 0.5, 0.5, 0.8, 0.9, 0.1, 0.1)
#' par_diag <- set_parameter_attributes(par_diag, K=2, model_type="tvp",
#'                                      diag_probs=TRUE, equal_variances=FALSE)
#' y <- rnorm(200)
#' loglik <- Rfiltering_TVP(par_diag, y, n_burnin=20, n_cutoff=10)
#' @export
Rfiltering_TVP <- function(par, y, n_burnin, n_cutoff, diagnostics = FALSE) {

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

  # Determine the number of transition parameters based on parameterization
  n_transition <- length(init_trans)

  # Initialize variables
  tot_lik <- numeric(M)                # Total likelihood
  X_t <- matrix(0, nrow=K, ncol=M)     # Filtered probabilities after observation
  X_tlag <- matrix(0, nrow=K, ncol=M)  # Predicted probabilities before observation

  # =========================================================================
  # PRE-COMPUTE EMISSION LIKELIHOODS (vectorized for performance)
  # =========================================================================
  # Instead of computing dnorm inside the loop, compute all at once
  eta <- matrix(0, nrow=K, ncol=M)
  sqrt_sigma2 <- sqrt(sigma2)  # Pre-compute sqrt once
  for (k in 1:K) {
    eta[k, ] <- dnorm(y, mu[k], sqrt_sigma2[k])
  }

  # Get baseline transition parameters in f-space
  # omega_LR is the long-run (unconditional) value of f
  # omega is the intercept in the AR(1) process: f[t+1] = omega + A * y[t]
  if (diag_probs) {
    params_to_f <- logit
    f_to_params <- logistic
  } else {
    params_to_f <- logit
    f_to_params <- logistic
  }
  omega_LR <- params_to_f(init_trans)
  omega <- omega_LR * (1 - A)  # This scaling ensures that when A=0, we get back to init_trans

  # Initialize f values for transition probabilities
  f <- matrix(0, nrow=n_transition, ncol=M)
  p_trans <- matrix(0, nrow=n_transition, ncol=M)

  # =========================================================================
  # INITIALIZATION (t=1): Generalized for K>=2
  # =========================================================================
  # At t=1, use omega_LR (long-run value) for initial transition probabilities
  f[, 1] <- omega_LR
  p_trans[, 1] <- f_to_params(f[, 1])
  if (!diag_probs) p_trans[, 1] <- clamp_offdiag_rowsums(p_trans[, 1], K)

  # Compute initial predicted probabilities using stationary distribution
  # This works for any K>=2 using eigenvalue-based stationary distribution
  X_tlag[, 1] <- compute_initial_predicted_probs(p_trans[, 1], diag_probs = diag_probs)

  # Note: eta[,1] already computed in vectorized pre-computation above

  # Compute total likelihood and filtered probabilities for t=1
  tot_lik[1] <- sum(eta[, 1] * X_tlag[, 1])
  if (tot_lik[1] <= 0 || is.na(tot_lik[1])) {
    tot_lik[1] <- .Machine$double.eps
    X_t[, 1] <- c(0.5, 0.5)
  } else {
    X_t[, 1] <- (eta[, 1] * X_tlag[, 1]) / tot_lik[1]
  }

  # =========================================================================
  # CRITICAL: Set f[,2] = omega (NOT omega + A*y[1])
  # This matches the original C code lines 60-65 which sets f1[1]=omega1, f2[1]=omega2
  # BEFORE entering the main loop. This is a special initialization step.
  # =========================================================================
  if (M >= 2) {
    f[, 2] <- omega
    p_trans[, 2] <- f_to_params(f[, 2])
    if (!diag_probs) p_trans[, 2] <- clamp_offdiag_rowsums(p_trans[, 2], K)
  }

  # =========================================================================
  # MAIN LOOP (t=2 to M-1): Standard TVP filtering - Generalized for K>=2
  # =========================================================================
  if (M >= 2) {
    for (t in 2:(M-1)) {
      # Generate predicted probabilities using generalized formula
      X_tlag[, t] <- compute_predicted_probs(X_t[, t-1], p_trans[, t], diag_probs = diag_probs)

      # Note: eta[,t] already computed in vectorized pre-computation above
      tot_lik[t] <- sum(eta[, t] * X_tlag[, t])

      # Compute filtered probabilities
      if (tot_lik[t] <= 0 || is.na(tot_lik[t])) {
        tot_lik[t] <- .Machine$double.eps
        X_t[, t] <- rep(1/K, K)  # Uniform fallback for any K
      } else {
        X_t[, t] <- (eta[, t] * X_tlag[, t]) / tot_lik[t]
      }

      # Update transition probabilities for next time step
      # f[t+1] = omega + A * y[t]
      if (t < M - 1) {
        f[, t+1] <- omega + A * y[t]
        p_trans[, t+1] <- f_to_params(f[, t+1])
        if (!diag_probs) p_trans[, t+1] <- clamp_offdiag_rowsums(p_trans[, t+1], K)
      }
    }
  }

  # =========================================================================
  # LAST TIME POINT (t=M): Process final observation - Generalized for K>=2
  # =========================================================================
  if (M >= 2) {
    # Need to set f[,M] if not already set
    if (M > 2) {
      f[, M] <- omega + A * y[M-1]
      p_trans[, M] <- f_to_params(f[, M])
      if (!diag_probs) p_trans[, M] <- clamp_offdiag_rowsums(p_trans[, M], K)
    }

    # Generate predicted probabilities using generalized formula
    X_tlag[, M] <- compute_predicted_probs(X_t[, M-1], p_trans[, M], diag_probs = diag_probs)

    # Note: eta[,M] already computed in vectorized pre-computation above
    tot_lik[M] <- sum(eta[, M] * X_tlag[, M])

    if (tot_lik[M] <= 0 || is.na(tot_lik[M])) {
      tot_lik[M] <- .Machine$double.eps
    }
  }

  # Sum log-likelihoods, but exclude burn-in and cut-off
  valid_indices <- (n_burnin+1):(M-n_cutoff)
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
      model_type = "tvp",
      n_transition = n_transition
    )
  }
  
  return(neg_log_lik)
}

#' Estimate TVP regime-switching model
#'
#' @param y Observed time series data
#' @param K Number of regimes
#' @param diag_probs If TRUE, use diagonal transition probability parameterization
#' @param equal_variances If TRUE, constrain all regimes to have equal variances
#' @param n_starts Number of random starting points for optimization (default: 10)
#' @param n_burnin Burn-in observations to exclude (default: 100)
#' @param n_cutoff Cut-off observations to exclude (default: 50)
#' @param bounds Optional list with lower and upper parameter bounds
#' @param early_stopping Enable early stopping for diverging starts (default: FALSE)
#' @param early_stop_patience Evaluations without improvement before stopping
#'   (default: 10000, calibrated from K=3 estimation on real data)
#' @param early_stop_max_evals Maximum evaluations per start (default: 60000)
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
#' # Estimate model with diagonal probabilities (original style)
#' y <- rnorm(200)
#' result_diag <- estimate_tvp_model(y, K=2, diag_probs=TRUE, n_starts=3,
#'                                   n_burnin=20, n_cutoff=10)
#'
#' # Estimate model with off-diagonal probabilities (new style)
#' result_offdiag <- estimate_tvp_model(y, K=2, diag_probs=FALSE, n_starts=3,
#'                                      n_burnin=20, n_cutoff=10)
#' }
#' @export
estimate_tvp_model <- function(y, K, diag_probs = TRUE, equal_variances = FALSE,
                               n_starts = 10, n_burnin = 100, n_cutoff = 50, bounds = NULL,
                               early_stopping = FALSE,
                               early_stop_patience = 10000L,
                               early_stop_max_evals = 60000L,
                               parallel = TRUE, cores = NULL, seed = NULL, verbose = 1) {

  # Input validation
  if (!is.numeric(y) || length(y) == 0) {
    stop("y must be a non-empty numeric vector")
  }
  if (!is.numeric(K) || K < 2 || K != as.integer(K)) {
    stop("K must be an integer >= 2")
  }
  if (length(y) <= n_burnin + n_cutoff + K) {
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
    cat("Estimating TVP (autoregressive) regime-switching model\n")
    cat("===================================================\n")
    cat("K:", K, "regimes\n")
    cat("Data points:", length(y), "(using", length(y) - n_burnin - n_cutoff, "after burn-in/cut-off)\n")
    cat("Parameterization:", ifelse(diag_probs, "diagonal", "off-diagonal"), "transition probabilities\n")
    cat("Variances:", ifelse(equal_variances, "equal (shared)", "separate"), "\n")
    cat("Starting points:", n_starts, "\n")
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
    model_type = "tvp",
    n_starts = n_starts,
    diag_probs = diag_probs,
    equal_variances = equal_variances,
    seed = seed
  )
  
  # Create parameter bounds if none provided
  if (is.null(bounds)) {
    # Calculate expected parameter count
    expected_length <- calculate_expected_length(K, "tvp", diag_probs, equal_variances)
    
    # Create default bounds for transformed parameters
    n_mu <- K
    n_sigma2 <- ifelse(equal_variances, 1, K)
    n_trans <- ifelse(diag_probs, K, K*(K-1))
    n_A <- n_trans  # A coefficients match transition structure
    
    lower_bounds <- c(rep(-Inf, n_mu),           # No bounds on means
                      rep(log(1e-10), n_sigma2), # Log-variances: floor at sigma2=1e-10
                      rep(-Inf, n_trans),        # Logit-probabilities
                      rep(-Inf, n_A))            # Logit-A coefficients
    
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

      # Define the base objective function
      base_objective <- function(par_t) {
        par_t_with_attrs <- transformed_params
        par_t_with_attrs[] <- par_t
        attr(par_t_with_attrs, "parameterization") <- "transformed"
        par_natural <- untransform_parameters(par_t_with_attrs)
        neg_log_lik <- Rfiltering_TVP(par_natural, y, n_burnin, n_cutoff, diagnostics = FALSE)
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
      optimization_result <- nlminb(
        start = transformed_params,
        objective = objective_fn,
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
  num_data_points <- length(y) - n_burnin - n_cutoff

  aic <- 2 * optimization_result$objective + 2 * num_params
  bic <- 2 * optimization_result$objective + num_params * log(num_data_points)
  
  # Calculate filtered probabilities and additional diagnostics (with full diagnostics)
  full_likelihood_result <- Rfiltering_TVP(estimated_params, y, n_burnin, n_cutoff, diagnostics = TRUE)
  
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
      model_type = "tvp",
      diag_probs = diag_probs,
      equal_variances = equal_variances,
      parameterization = "natural"
    ),
    data_info = list(
      n_obs = length(y),
      burn_in = n_burnin,
      cut_off = n_cutoff,
      n_used = num_data_points
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
