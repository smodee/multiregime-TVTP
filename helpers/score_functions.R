#' GAS Score-Driven Helper Functions
#' 
#' IMPORTANT: This file contains both PRODUCTION-READY and EXPERIMENTAL functions.
#' 
#' PRODUCTION-READY (recommended for parameter estimation):
#' - calculate_gas_score() with scaling_method = "simple" (default)
#' - apply_moore_penrose_scaling() with method = "simple" or "normalized"  
#' - setup_gauss_hermite_quadrature()
#' 
#' EXPERIMENTAL (research only, not for estimation):
#' - calculate_gas_score_robust() - causes optimization instability
#' - calculate_fisher_information_robust() - computationally intensive
#' - Any functions using Moore-Penrose pseudo-inverse scaling
#' 
#' The experimental functions implement theoretically correct Moore-Penrose 
#' scaling from Bazzi et al. (2017) but are numerically unstable for optimization.
#' Use simple scaling methods for reliable parameter estimation.

# Load required helper functions
source("helpers/utility_functions.R")
source("helpers/transition_helpers.R")

# Check for required package
if (!require(fastGHQuad, quietly = TRUE)) {
  stop("Package 'fastGHQuad' is required for GAS score calculations. Please install it with: install.packages('fastGHQuad')")
}

#' Setup Gauss-Hermite quadrature for GAS score scaling
#'
#' @param y Observed time series data
#' @param n_nodes Number of quadrature nodes (default: 30)
#' @param method Method for determining quadrature parameters ("data_based" or "standardized")
#' @return List containing quadrature setup information
#' @details 
#' Sets up Gauss-Hermite quadrature nodes and weights for numerical integration
#' in the Fisher Information calculation. Following the original implementation,
#' the quadrature is centered around the median of the data with scale based on
#' the standard deviation.
#'
#' IMPORTANT: The integration domain is limited to ±4 standard deviations to avoid
#' numerical issues where regime densities become negligible (< 1e-10). This prevents
#' quadrature node failures that can bias parameter estimation toward zero.
#'
#' The "data_based" method uses median(y) and sd(y) as in the original implementation.
#' The "standardized" method uses mean=0, sd=1 for standardized integration.
#'
#' @examples
#' # Setup quadrature based on data characteristics
#' y <- rnorm(1000)
#' gh_setup <- setup_gauss_hermite_quadrature(y)
#' 
#' # Setup with more nodes for higher precision
#' gh_setup <- setup_gauss_hermite_quadrature(y, n_nodes = 50)
#' 
#' # Use standardized approach
#' gh_setup <- setup_gauss_hermite_quadrature(y, method = "standardized")
setup_gauss_hermite_quadrature <- function(y, n_nodes = 30, method = "data_based") {
  # Validate inputs
  if (!is.numeric(y) || length(y) == 0) {
    stop("Input 'y' must be a non-empty numeric vector")
  }
  
  if (!is.numeric(n_nodes) || n_nodes < 2 || n_nodes != as.integer(n_nodes)) {
    stop("Input 'n_nodes' must be an integer >= 2")
  }
  
  method <- match.arg(method, c("data_based", "standardized"))
  
  # Check for invalid values in y and clean data
  if (has_invalid_values(y)) {
    warning("Input data contains NA, NaN, or infinite values. These will be removed.")
    y_clean <- y[!is.na(y) & !is.nan(y) & is.finite(y)]
    
    if (length(y_clean) == 0) {
      stop("No valid data points remaining after removing invalid values")
    }
  } else {
    y_clean <- y
  }
  
  # Generate standard Gauss-Hermite quadrature nodes and weights
  # This gives us nodes and weights for standard normal N(0,1)
  gh_data <- gaussHermiteData(n_nodes)
  
  # Determine quadrature parameters based on method
  if (method == "data_based") {
    # Use data characteristics for centering and scaling
    mu_quad <- median(y_clean, na.rm = TRUE)
    sigma_quad <- sd(y_clean, na.rm = TRUE)
    
    # Protect against zero variance
    if (sigma_quad <= 0 || is.na(sigma_quad)) {
      warning("Data has zero or invalid variance. Using unit variance for quadrature.")
      sigma_quad <- 1.0
    }
    
    # CRITICAL FIX: Limit the integration domain to ±4 standard deviations
    # This prevents numerical issues where regime densities become negligible
    # Beyond ±4σ, normal densities drop below 1e-6, causing Fisher Information
    # calculation failures that bias A parameters toward zero
    node_range <- 4.0  # Covers 99.99% of probability mass
    
    # Transform nodes from standard to limited data-appropriate range
    # Instead of letting nodes go to ±∞, we constrain them to reasonable values
    raw_nodes <- gh_data$x  # These are from standard normal
    node_scale <- node_range / max(abs(raw_nodes))  # Scale factor to limit range
    
    # Apply the transformation: center at mu_quad, scale by sigma_quad, but limit range
    nodes <- mu_quad + sigma_quad * node_scale * raw_nodes
    weights <- gh_data$w / sqrt(pi)  # Standard Gauss-Hermite normalization
    
  } else {  # standardized method
    # Use standardized parameters without range limitation for comparison
    mu_quad <- 0.0
    sigma_quad <- 1.0
    
    # Standard transformation for comparison/reference
    nodes <- sqrt(2) * gh_data$x
    weights <- gh_data$w / sqrt(pi)
  }
  
  # Store additional information for diagnostics
  mu <- if (method == "data_based") mu_quad else 0.0
  sigma <- if (method == "data_based") sigma_quad else 1.0
  
  # Create the setup object with comprehensive information
  setup <- list(
    nodes = nodes,
    weights = weights,
    n_nodes = n_nodes,
    mu = mu,              # Keep as 'mu' for backward compatibility
    sigma = sigma,        # Keep as 'sigma' for backward compatibility  
    mu_quad = mu_quad,    # Explicit quadrature parameters
    sigma_quad = sigma_quad,
    method = method,
    node_range_limit = if (method == "data_based") node_range else NULL,
    data_characteristics = list(
      n_obs = length(y_clean),
      mean_y = mean(y_clean, na.rm = TRUE),
      median_y = median(y_clean, na.rm = TRUE),
      sd_y = sd(y_clean, na.rm = TRUE),
      min_y = min(y_clean, na.rm = TRUE),
      max_y = max(y_clean, na.rm = TRUE)
    )
  )
  
  # Add class for potential method dispatch
  class(setup) <- c("gauss_hermite_setup", "list")
  
  return(setup)
}

#' Print method for Gauss-Hermite setup objects
#'
#' @param x A gauss_hermite_setup object
#' @param ... Additional arguments (ignored)
#' @return NULL (prints information as side effect)
#' @details
#' Provides a summary of the quadrature setup including number of nodes,
#' integration parameters, and data characteristics.
#'
#' @examples
#' y <- rnorm(1000)
#' gh_setup <- setup_gauss_hermite_quadrature(y)
#' print(gh_setup)
print.gauss_hermite_setup <- function(x, ...) {
  cat("Gauss-Hermite Quadrature Setup\n")
  cat("==============================\n")
  cat("Number of nodes:", x$n_nodes, "\n")
  cat("Method:", x$method, "\n")
  cat("Integration center (μ):", round(x$mu_quad, 4), "\n")
  cat("Integration scale (σ):", round(x$sigma_quad, 4), "\n")
  cat("Node range: [", round(min(x$nodes), 4), ", ", round(max(x$nodes), 4), "]\n", sep = "")
  
  cat("\nData Characteristics:\n")
  cat("  Observations:", x$data_characteristics$n_obs, "\n")
  cat("  Mean:", round(x$data_characteristics$mean_y, 4), "\n")
  cat("  Median:", round(x$data_characteristics$median_y, 4), "\n")
  cat("  Std Dev:", round(x$data_characteristics$sd_y, 4), "\n")
  cat("  Range: [", round(x$data_characteristics$min_y, 4), ", ", 
      round(x$data_characteristics$max_y, 4), "]\n", sep = "")
}

#' Validate Gauss-Hermite setup object
#'
#' @param gh_setup A gauss_hermite_setup object to validate
#' @return TRUE if valid, otherwise throws an error with explanation
#' @details
#' Checks that the Gauss-Hermite setup object has all required components
#' and that they are valid for use in score calculations.
#'
#' @examples
#' y <- rnorm(1000)
#' gh_setup <- setup_gauss_hermite_quadrature(y)
#' validate_gauss_hermite_setup(gh_setup)  # Returns TRUE
validate_gauss_hermite_setup <- function(gh_setup) {
  # Check if object has the right class
  if (!inherits(gh_setup, "gauss_hermite_setup")) {
    stop("Object is not a valid gauss_hermite_setup object")
  }
  
  # Check required components
  required_components <- c("nodes", "weights", "n_nodes", "mu_quad", "sigma_quad", "method")
  missing_components <- setdiff(required_components, names(gh_setup))
  
  if (length(missing_components) > 0) {
    stop("Missing required components: ", paste(missing_components, collapse = ", "))
  }
  
  # Validate nodes and weights
  if (!is.numeric(gh_setup$nodes) || !is.numeric(gh_setup$weights)) {
    stop("Nodes and weights must be numeric")
  }
  
  if (length(gh_setup$nodes) != length(gh_setup$weights)) {
    stop("Nodes and weights must have the same length")
  }
  
  if (length(gh_setup$nodes) != gh_setup$n_nodes) {
    stop("Length of nodes does not match n_nodes")
  }
  
  # Check for invalid values
  if (has_invalid_values(gh_setup$nodes) || has_invalid_values(gh_setup$weights)) {
    stop("Nodes or weights contain invalid values (NA, NaN, Inf)")
  }
  
  # Check that weights are positive
  if (any(gh_setup$weights <= 0)) {
    stop("All weights must be positive")
  }
  
  # Check parameters
  if (!is.finite(gh_setup$mu_quad) || !is.finite(gh_setup$sigma_quad)) {
    stop("Quadrature parameters (mu_quad, sigma_quad) must be finite")
  }
  
  if (gh_setup$sigma_quad <= 0) {
    stop("Quadrature scale parameter (sigma_quad) must be positive")
  }
  
  # Check method
  if (!gh_setup$method %in% c("data_based", "standardized")) {
    stop("Method must be either 'data_based' or 'standardized'")
  }
  
  return(TRUE)
}

#' Calculate properly scaled GAS score vector
#'
#' @param y_t Current observation at time t
#' @param mu Vector of regime means (length K)
#' @param sigma2 Vector of regime variances (length K)
#' @param X_t_lag Predicted probabilities before observation (length K)
#' @param X_t_prev Filtered probabilities from previous time step (length K)
#' @param p_trans Current transition probabilities (length K*(K-1))
#' @param gh_setup Gauss-Hermite quadrature setup from setup_gauss_hermite_quadrature()
#' @param K Number of regimes
#' @param scaling_method Method for scaling ("moore_penrose", "simple", "normalized", or "original")
#' @return Properly scaled score vector (length K*(K-1))
#' @details 
#' Main interface function that combines all components of the GAS score calculation
#' following Bazzi et al. (2017). This function:
#' 
#' 1. Calculates regime likelihoods for the current observation
#' 2. Computes the raw score vector (equation 13 in Bazzi et al.)
#' 3. Calculates Fisher Information using Gauss-Hermite quadrature
#' 4. Applies Moore-Penrose pseudo-inverse scaling (equation 15 in Bazzi et al.)
#' 
#' The resulting scaled score vector is used to update the time-varying parameters
#' in the GAS model: f_{t+1} = ω + A*s_t + B*(f_t - ω)
#' 
#' This implementation generalizes the original 2-regime methodology to K regimes
#' while maintaining the mathematical foundation and numerical stability.
#'
#' @examples
#' # Setup for a 3-regime model
#' y <- rnorm(1000)
#' gh_setup <- setup_gauss_hermite_quadrature(y)
#' 
#' # Model parameters
#' K <- 3
#' mu <- c(-1, 0, 1)
#' sigma2 <- c(0.5, 1, 1.5)
#' X_t_lag <- c(0.3, 0.4, 0.3)
#' X_t_prev <- c(0.25, 0.45, 0.3)
#' p_trans <- rep(0.2, 6)
#' 
#' # Current observation
#' y_t <- 0.5
#' 
#' # Calculate scaled score
#' scaled_score <- calculate_gas_score(y_t, mu, sigma2, X_t_lag, X_t_prev, 
#'                                    p_trans, gh_setup, K)
#' 
#' # Use alternative scaling method
#' scaled_score <- calculate_gas_score(y_t, mu, sigma2, X_t_lag, X_t_prev, 
#'                                    p_trans, gh_setup, K, 
#'                                    scaling_method = "original")
calculate_gas_score <- function(y_t, mu, sigma2, X_t_lag, X_t_prev, p_trans, 
                                gh_setup, K, scaling_method = "simple") {
  if (scaling_method == "moore_penrose") {
    warning(
      paste(
        "WARNING: scaling_method='moore_penrose' is EXPERIMENTAL.",
        "For parameter estimation, use 'simple' (default) instead.",
        "Moore-Penrose scaling can cause optimization failures."
      ),
      call. = FALSE
    )
  }
  
  # Validate all inputs
  if (!is.numeric(y_t) || length(y_t) != 1) {
    stop("y_t must be a numeric scalar")
  }
  
  if (has_invalid_values(y_t)) {
    stop("y_t contains invalid values (NA, NaN, or Inf)")
  }
  
  validate_fisher_inputs(mu, sigma2, X_t_prev, p_trans, gh_setup, K)
  
  if (!is.numeric(X_t_lag) || length(X_t_lag) != K) {
    stop("X_t_lag must be a numeric vector of length K")
  }
  
  if (any(X_t_lag < 0) || any(X_t_lag > 1) || abs(sum(X_t_lag) - 1) > 1e-10) {
    stop("X_t_lag must be valid probabilities that sum to 1")
  }
  
  scaling_method <- match.arg(scaling_method, c("simple", "normalized", "original", "moore_penrose"))
  
  # Step 1: Calculate regime likelihoods for the current observation
  eta <- numeric(K)
  for (k in 1:K) {
    eta[k] <- dnorm(y_t, mu[k], sqrt(sigma2[k]))
  }
  
  # Calculate total likelihood
  tot_lik <- sum(eta * X_t_lag)
  
  # Protect against numerical issues
  if (tot_lik <= .Machine$double.eps || !is.finite(tot_lik)) {
    warning("Total likelihood is too small or invalid. Using fallback calculation.")
    tot_lik <- .Machine$double.eps
    # In case of numerical issues, return zero score
    return(rep(0, K*(K-1)))
  }
  
  # Step 2: Calculate raw score vector
  raw_score <- calculate_raw_score_vector(eta, tot_lik, X_t_prev, p_trans, K)
  
  # If raw score is essentially zero, return it as-is
  if (all(abs(raw_score) < .Machine$double.eps)) {
    attr(raw_score, "calculation_info") <- list(
      method = "zero_score",
      y_t = y_t,
      tot_lik = tot_lik,
      fisher_info = NA,
      scaling_method = scaling_method
    )
    return(raw_score)
  }
  
  # Step 3: Calculate Fisher Information
  fisher_info <- calculate_fisher_information(mu, sigma2, X_t_prev, p_trans, gh_setup, K)
  
  # Step 4: Apply appropriate scaling method
  if (scaling_method == "original") {
    # Use the original 2-regime methodology
    g_vector <- attr(raw_score, "g_components")
    if (is.null(g_vector)) {
      # Recalculate g_vector if not available
      g_vector <- numeric(K*(K-1))
      idx <- 1
      for (i in 1:K) {
        for (j in 1:K) {
          if (i != j) {
            g_vector[idx] <- X_t_prev[i] * p_trans[idx] * (1 - p_trans[idx])
            idx <- idx + 1
          }
        }
      }
    }
    scaled_score <- calculate_original_scaling(raw_score, fisher_info, g_vector)
  } else {
    # Use Moore-Penrose or other modern scaling methods
    scaled_score <- apply_moore_penrose_scaling(raw_score, fisher_info, scaling_method)
  }
  
  # Step 5: Final validation and cleanup
  if (has_invalid_values(scaled_score)) {
    warning("Scaled score contains invalid values. Using conservative fallback.")
    # Conservative fallback: small scaled version of raw score
    score_norm <- sqrt(sum(raw_score^2))
    if (score_norm > 0) {
      scaled_score <- raw_score / score_norm * 0.01  # Very conservative scaling
    } else {
      scaled_score <- rep(0, K*(K-1))
    }
  }
  
  # Add comprehensive attributes for debugging and analysis
  attr(scaled_score, "calculation_info") <- list(
    y_t = y_t,
    tot_lik = tot_lik,
    fisher_info = fisher_info,
    scaling_method = scaling_method,
    raw_score_norm = sqrt(sum(raw_score^2)),
    scaled_score_norm = sqrt(sum(scaled_score^2)),
    regime_likelihoods = eta,
    most_likely_regime = which.max(eta),
    likelihood_spread = max(eta) - min(eta)
  )
  
  attr(scaled_score, "regime_diagnostics") <- list(
    K = K,
    filtered_prob_entropy = -sum(X_t_prev * log(pmax(X_t_prev, .Machine$double.eps))),
    predicted_prob_entropy = -sum(X_t_lag * log(pmax(X_t_lag, .Machine$double.eps))),
    max_transition_prob = max(p_trans),
    min_transition_prob = min(p_trans)
  )
  
  return(scaled_score)
}

#' Batch calculate GAS scores for multiple time points
#'
#' @param y_series Time series of observations (length T)
#' @param mu Vector of regime means (length K)
#' @param sigma2 Vector of regime variances (length K)
#' @param X_t_lag_series Matrix of predicted probabilities (T x K)
#' @param X_t_prev_series Matrix of previous filtered probabilities (T x K)
#' @param p_trans_series Matrix of transition probabilities (T x K*(K-1))
#' @param gh_setup Gauss-Hermite quadrature setup
#' @param K Number of regimes
#' @param scaling_method Scaling method to use
#' @param verbose Whether to show progress
#' @return Matrix of scaled scores (T x K*(K-1))
#' @details
#' Efficiently calculates GAS scores for a series of observations. Useful for
#' batch processing and when working with long time series. Includes progress
#' reporting and error handling for individual time points.
#'
#' @examples
#' # Calculate scores for an entire time series
#' y_series <- rnorm(100)
#' # ... setup matrices for X_t_lag_series, X_t_prev_series, p_trans_series ...
#' # scores <- calculate_gas_scores_batch(y_series, mu, sigma2, ...)
calculate_gas_scores_batch <- function(y_series, mu, sigma2, X_t_lag_series, 
                                       X_t_prev_series, p_trans_series, gh_setup, K,
                                       scaling_method = "moore_penrose", verbose = FALSE) {
  # Validate inputs
  T_obs <- length(y_series)
  n_transition <- K * (K - 1)
  
  if (!is.matrix(X_t_lag_series) || nrow(X_t_lag_series) != T_obs || ncol(X_t_lag_series) != K) {
    stop("X_t_lag_series must be a T x K matrix")
  }
  
  if (!is.matrix(X_t_prev_series) || nrow(X_t_prev_series) != T_obs || ncol(X_t_prev_series) != K) {
    stop("X_t_prev_series must be a T x K matrix")
  }
  
  if (!is.matrix(p_trans_series) || nrow(p_trans_series) != T_obs || ncol(p_trans_series) != n_transition) {
    stop("p_trans_series must be a T x K*(K-1) matrix")
  }
  
  # Initialize output
  scaled_scores <- matrix(0, nrow = T_obs, ncol = n_transition)
  
  # Calculate scores for each time point
  for (t in 1:T_obs) {
    if (verbose && t %% 100 == 0) {
      cat("Processing time point", t, "of", T_obs, "\n")
    }
    
    tryCatch({
      scaled_scores[t, ] <- calculate_gas_score(
        y_t = y_series[t],
        mu = mu,
        sigma2 = sigma2,
        X_t_lag = X_t_lag_series[t, ],
        X_t_prev = X_t_prev_series[t, ],
        p_trans = p_trans_series[t, ],
        gh_setup = gh_setup,
        K = K,
        scaling_method = scaling_method
      )
    }, error = function(e) {
      if (verbose) {
        cat("Error at time point", t, ":", e$message, "\n")
      }
      # Use zero score for failed calculations
      scaled_scores[t, ] <<- rep(0, n_transition)
    })
  }
  
  return(scaled_scores)
}

#' Print comprehensive summary of GAS score calculation
#'
#' @param scaled_score Scaled score vector from calculate_gas_score()
#' @param ... Additional arguments (ignored)
#' @return NULL (prints summary as side effect)
#' @details
#' Provides detailed information about the GAS score calculation including
#' input parameters, intermediate results, and diagnostics.
print_gas_score_summary <- function(scaled_score, ...) {
  if (!is.numeric(scaled_score)) {
    stop("Input must be a numeric vector from calculate_gas_score()")
  }
  
  calc_info <- attr(scaled_score, "calculation_info")
  regime_diag <- attr(scaled_score, "regime_diagnostics")
  
  cat("GAS Score Calculation Summary\n")
  cat("=============================\n")
  
  if (!is.null(calc_info)) {
    cat("Observation (y_t):", round(calc_info$y_t, 4), "\n")
    cat("Total likelihood:", format(calc_info$tot_lik, scientific = TRUE, digits = 4), "\n")
    cat("Fisher Information:", format(calc_info$fisher_info, scientific = TRUE, digits = 4), "\n")
    cat("Scaling method:", calc_info$scaling_method, "\n")
    cat("Raw score norm:", format(calc_info$raw_score_norm, digits = 4), "\n")
    cat("Scaled score norm:", format(calc_info$scaled_score_norm, digits = 4), "\n")
    cat("Most likely regime:", calc_info$most_likely_regime, "\n")
    cat("Likelihood spread:", format(calc_info$likelihood_spread, digits = 4), "\n")
  }
  
  if (!is.null(regime_diag)) {
    cat("\nRegime Diagnostics:\n")
    cat("Number of regimes (K):", regime_diag$K, "\n")
    cat("Filtered prob entropy:", format(regime_diag$filtered_prob_entropy, digits = 4), "\n")
    cat("Predicted prob entropy:", format(regime_diag$predicted_prob_entropy, digits = 4), "\n")
    cat("Transition prob range: [", format(regime_diag$min_transition_prob, digits = 4), 
        ", ", format(regime_diag$max_transition_prob, digits = 4), "]\n", sep = "")
  }
  
  cat("\nScaled Score Vector:\n")
  cat(format(scaled_score, digits = 4), "\n")
}

#' Apply Moore-Penrose pseudo-inverse scaling to raw score vector
#'
#' @param raw_score Raw score vector from calculate_raw_score_vector() (length K*(K-1))
#' @param fisher_info Fisher Information scalar from calculate_fisher_information()
#' @param method Scaling method ("moore_penrose", "simple", or "normalized")
#' @return Scaled score vector (length K*(K-1))
#' @details 
#' Applies proper scaling to the raw score vector using the Moore-Penrose 
#' pseudo-inverse approach described in Bazzi et al. (2017), equation (15).
#' 
#' Since the Fisher Information matrix is singular by design (as noted in the paper),
#' we use the Moore-Penrose pseudo-inverse for scaling. The original 2-regime
#' implementation uses: s_t = S_t * r_t, where S_t is related to the square root
#' of the pseudo-inverse.
#' 
#' Following the original implementation: S_star = S1 / sqrt(I), where:
#' - S1 is the raw score component
#' - I is the Fisher Information
#' 
#' For the K-regime case, we generalize this by:
#' 1. Normalizing the raw score vector
#' 2. Scaling by the inverse square root of Fisher Information
#' 3. Applying additional normalization for numerical stability
#'
#' @examples
#' # Apply scaling to a raw score vector
#' raw_score <- c(0.1, -0.2, 0.05, -0.1, 0.15, -0.05)  # 6 elements for 3-regime model
#' fisher_info <- 2.5
#' 
#' scaled_score <- apply_moore_penrose_scaling(raw_score, fisher_info)
#' 
#' # Use simple scaling method
#' scaled_score <- apply_moore_penrose_scaling(raw_score, fisher_info, method = "simple")
apply_moore_penrose_scaling <- function(raw_score, fisher_info, method = "simple") {
  if (method == "moore_penrose") {
    warning(
      paste(
        "WARNING: method='moore_penrose' is EXPERIMENTAL and NOT recommended for estimation.",
        "It can cause optimization algorithms to fail due to numerical instability.",
        "For parameter estimation, use method='simple' (default) or 'normalized' instead.",
        "Moore-Penrose scaling is kept for research/comparison purposes only."
      ),
      call. = FALSE
    )
  }
  
  # Validate inputs
  if (!is.numeric(raw_score) || length(raw_score) == 0) {
    stop("raw_score must be a non-empty numeric vector")
  }
  
  if (!is.numeric(fisher_info) || length(fisher_info) != 1 || fisher_info <= 0) {
    stop("fisher_info must be a positive scalar")
  }
  
  if (has_invalid_values(raw_score)) {
    stop("raw_score contains invalid values (NA, NaN, or Inf)")
  }
  
  if (has_invalid_values(fisher_info)) {
    stop("fisher_info contains invalid values (NA, NaN, or Inf)")
  }
  
  method <- match.arg(method, c("moore_penrose", "simple", "normalized"))
  
  # Get dimensions
  n_transition <- length(raw_score)
  
  # Handle zero score vector
  if (all(abs(raw_score) < .Machine$double.eps)) {
    return(raw_score)  # Return zero vector as-is
  }
  
  # Apply scaling based on method
  if (method == "moore_penrose") {
    # Following Bazzi et al. (2017) and original implementation
    # This implements the Moore-Penrose pseudo-inverse scaling
    
    # Step 1: Calculate the norm of the raw score vector
    # This represents the "magnitude" of the score
    score_norm <- sqrt(sum(raw_score^2))
    
    if (score_norm < .Machine$double.eps) {
      # Handle near-zero score vector
      scaled_score <- raw_score
    } else {
      # Step 2: Normalize the score vector (unit vector direction)
      # This is equivalent to g/||g|| in the original formulation
      normalized_score <- raw_score / score_norm
      
      # Step 3: Apply Fisher Information scaling
      # Following original: S_star = S1 / sqrt(I)
      # where S1 is the score component and I is Fisher Information
      fisher_scale <- 1.0 / sqrt(fisher_info)
      
      # Step 4: Combine the scaling factors
      # This gives us the properly scaled score vector
      scaled_score <- fisher_scale * normalized_score * score_norm
      
      # Step 5: Additional numerical stability
      # Apply a reasonable bound to prevent extreme values
      max_scale <- 10.0  # Prevent excessive scaling
      scale_factor <- min(1.0, max_scale / max(abs(scaled_score)))
      scaled_score <- scaled_score * scale_factor
    }
    
  } else if (method == "simple") {
    # Simple inverse square root scaling (less sophisticated)
    fisher_scale <- 1.0 / sqrt(fisher_info)
    scaled_score <- raw_score * fisher_scale
    
  } else if (method == "normalized") {
    # Normalized scaling with Fisher Information
    score_norm <- sqrt(sum(raw_score^2))
    
    if (score_norm < .Machine$double.eps) {
      scaled_score <- raw_score
    } else {
      # Normalize to unit vector and scale by Fisher Information
      fisher_scale <- 1.0 / sqrt(fisher_info)
      scaled_score <- (raw_score / score_norm) * fisher_scale
    }
  }
  
  # Final validation and cleanup
  # Ensure no invalid values in output
  if (has_invalid_values(scaled_score)) {
    warning("Scaling produced invalid values. Using fallback scaling.")
    # Fallback: simple normalization
    score_norm <- sqrt(sum(raw_score^2))
    if (score_norm > 0) {
      scaled_score <- raw_score / score_norm * 0.1  # Conservative scaling
    } else {
      scaled_score <- raw_score
    }
  }
  
  # Add attributes for debugging and analysis
  attr(scaled_score, "scaling_info") <- list(
    method = method,
    fisher_info = fisher_info,
    raw_score_norm = sqrt(sum(raw_score^2)),
    scaled_score_norm = sqrt(sum(scaled_score^2)),
    scaling_factor = if (sqrt(sum(raw_score^2)) > 0) {
      sqrt(sum(scaled_score^2)) / sqrt(sum(raw_score^2))
    } else {
      1.0
    }
  )
  
  return(scaled_score)
}

#' Calculate score vector using original 2-regime methodology
#'
#' @param raw_score Raw score vector (length K*(K-1))
#' @param fisher_info Fisher Information scalar
#' @param g_vector G-vector from score calculation (length K*(K-1))
#' @return Scaled score vector following original methodology
#' @details
#' Alternative scaling method that more closely follows the original 2-regime
#' implementation. This function implements the exact methodology from the
#' original simulation.R code:
#' 
#' 1. Calculate g_mod = sqrt(sum(g^2))
#' 2. Normalize g if g_mod > 0: g = g / g_mod
#' 3. Calculate S_star = S1 / sqrt(I)
#' 4. Final score = S_star * g
#' 
#' This provides an alternative to the Moore-Penrose method for comparison.
calculate_original_scaling <- function(raw_score, fisher_info, g_vector) {
  # Validate inputs
  if (!is.numeric(raw_score) || !is.numeric(g_vector)) {
    stop("raw_score and g_vector must be numeric")
  }
  
  if (length(raw_score) != length(g_vector)) {
    stop("raw_score and g_vector must have the same length")
  }
  
  if (!is.numeric(fisher_info) || fisher_info <= 0) {
    stop("fisher_info must be a positive scalar")
  }
  
  # Following the original implementation exactly
  
  # Step 1: Calculate g_mod (norm of g-vector)
  g_mod <- sqrt(sum(g_vector^2))
  
  # Step 2: Normalize g-vector if possible
  if (g_mod > .Machine$double.eps) {
    g_normalized <- g_vector / g_mod
  } else {
    g_normalized <- g_vector
  }
  
  # Step 3: Calculate S_star following original formula
  # In the original: S_star = S1 / sqrt(I)
  # where S1 is (eta[1] - eta[2]) / tot_lik
  
  # For K-regime case, we use the magnitude of the score differences
  S1 <- sqrt(sum(raw_score^2))  # Magnitude of raw score
  
  if (fisher_info > .Machine$double.eps) {
    S_star <- S1 / sqrt(fisher_info)
  } else {
    S_star <- 0.0
  }
  
  # Step 4: Final scaled score
  scaled_score <- S_star * g_normalized
  
  return(scaled_score)
}

#' Validate scaling inputs and parameters
#'
#' @param raw_score Raw score vector to validate
#' @param fisher_info Fisher Information value to validate
#' @param method Scaling method to validate
#' @return TRUE if valid, otherwise throws an error
#' @details
#' Comprehensive validation for score scaling inputs. Ensures all parameters
#' are mathematically valid and within reasonable ranges.
validate_scaling_inputs <- function(raw_score, fisher_info, method = "moore_penrose") {
  # Check raw score
  if (!is.numeric(raw_score)) {
    stop("raw_score must be numeric")
  }
  
  if (length(raw_score) == 0) {
    stop("raw_score must not be empty")
  }
  
  if (has_invalid_values(raw_score)) {
    stop("raw_score contains invalid values (NA, NaN, or Inf)")
  }
  
  # Check Fisher Information
  if (!is.numeric(fisher_info)) {
    stop("fisher_info must be numeric")
  }
  
  if (length(fisher_info) != 1) {
    stop("fisher_info must be a scalar")
  }
  
  if (has_invalid_values(fisher_info)) {
    stop("fisher_info contains invalid values (NA, NaN, or Inf)")
  }
  
  if (fisher_info <= 0) {
    stop("fisher_info must be positive")
  }
  
  # Check method
  valid_methods <- c("moore_penrose", "simple", "normalized")
  if (!method %in% valid_methods) {
    stop("method must be one of: ", paste(valid_methods, collapse = ", "))
  }
  
  # Check for extreme values that might cause numerical issues
  if (fisher_info > 1e10) {
    warning("fisher_info is very large, scaling may be unstable")
  }
  
  if (fisher_info < 1e-10) {
    warning("fisher_info is very small, scaling may be unstable")
  }
  
  max_score <- max(abs(raw_score))
  if (max_score > 1e5) {
    warning("raw_score contains very large values, scaling may be unstable")
  }
  
  return(TRUE)
}

#' Calculate Fisher Information matrix using Gauss-Hermite quadrature
#'
#' @param mu Vector of regime means (length K)
#' @param sigma2 Vector of regime variances (length K)
#' @param X_t_prev Previous filtered probabilities (length K)
#' @param p_trans Current transition probabilities (length K*(K-1))
#' @param gh_setup Gauss-Hermite quadrature setup from setup_gauss_hermite_quadrature()
#' @param K Number of regimes
#' @return Fisher Information value (scalar)
#' @details 
#' Calculates the Fisher Information matrix for the GAS model using numerical
#' integration via Gauss-Hermite quadrature. Following the methodology in 
#' Bazzi et al. (2017), this involves integrating the squared score function
#' over all possible values of y_t.
#' 
#' The Fisher Information is used for scaling the score vector to ensure
#' proper statistical properties. Since the matrix is singular by design
#' (as noted in Bazzi et al.), we compute a scalar measure that captures
#' the overall information content.
#' 
#' The integration follows the original 2-regime implementation:
#' I = ∫ [(d1-d2)²/den] * φ(y; μ_quad, σ_quad²) dy
#' 
#' where d1, d2 are regime densities and den is the total density.
#'
#' @examples
#' # Setup quadrature and calculate Fisher Information
#' y <- rnorm(1000)
#' gh_setup <- setup_gauss_hermite_quadrature(y)
#' 
#' K <- 3
#' mu <- c(-1, 0, 1)
#' sigma2 <- c(0.5, 1, 1.5)
#' X_t_prev <- c(0.3, 0.4, 0.3)
#' p_trans <- rep(0.2, 6)
#' 
#' fisher_info <- calculate_fisher_information(mu, sigma2, X_t_prev, p_trans, gh_setup, K)
calculate_fisher_information <- function(mu, sigma2, X_t_prev, p_trans, gh_setup, K) {
  # Validate inputs
  validate_score_inputs(c(rep(1, K)), 1, X_t_prev, p_trans, K)  # Use dummy eta and tot_lik
  validate_gauss_hermite_setup(gh_setup)
  
  if (!is.numeric(mu) || length(mu) != K) {
    stop("mu must be a numeric vector of length K")
  }
  
  if (!is.numeric(sigma2) || length(sigma2) != K || any(sigma2 <= 0)) {
    stop("sigma2 must be a numeric vector of length K with positive values")
  }
  
  # Extract quadrature nodes and weights
  nodes <- gh_setup$nodes
  weights <- gh_setup$weights
  n_nodes <- gh_setup$n_nodes
  
  # Initialize Fisher Information accumulator
  fisher_info <- 0.0
  
  # Parameters for the quadrature weight function
  # This represents the distribution we're integrating over
  mu_weights <- gh_setup$mu_quad
  sigma_weights <- gh_setup$sigma_quad
  
  # Numerical integration using Gauss-Hermite quadrature
  # Following the original implementation in simulation.R
  for (j in 1:n_nodes) {
    # Current quadrature node (y-value)
    y_node <- nodes[j]
    current_weight <- weights[j]
    
    # Calculate regime densities at this node
    regime_densities <- numeric(K)
    for (k in 1:K) {
      regime_densities[k] <- dnorm(y_node, mu[k], sqrt(sigma2[k]))
    }
    
    # Calculate predicted probabilities for this hypothetical observation
    # This requires reconstructing the transition matrix from p_trans
    p_trans_valid <- convert_to_valid_probs(p_trans, K)
    transition_mat <- transition_matrix(p_trans_valid, check_validity = FALSE)
    X_pred <- as.vector(transition_mat %*% X_t_prev)
    
    # Calculate total density (denominator)
    total_density <- sum(regime_densities * X_pred)
    
    # Calculate Fisher Information components
    # We need to sum over all regime pairs for the K-regime case
    info_star <- 0.0
    
    if (total_density > .Machine$double.eps) {  # Avoid division by zero
      
      # Calculate all pairwise likelihood differences
      # This generalizes the (d1-d2)² term from the original 2-regime case
      for (i in 1:K) {
        for (l in 1:K) {
          if (i != l) {  # Only consider different regimes
            
            # Likelihood difference between regimes i and l
            likelihood_diff <- regime_densities[i] - regime_densities[l]
            
            # Contribution to Fisher Information
            # Following the original formula: ((d1-d2)²/den)
            info_contribution <- (likelihood_diff^2) / total_density
            
            # Weight by the probability of being in regime i
            # This accounts for the regime-specific nature of transitions
            regime_weight <- X_pred[i]
            
            info_star <- info_star + regime_weight * info_contribution
          }
        }
      }
      
      # Apply the quadrature weight function correction
      # This accounts for the fact that we're integrating over a specific distribution
      # Formula from original: (2*π*σ²)^0.5 * exp((y-μ)²/(2*σ²))
      weight_correction <- sqrt(2 * pi * sigma_weights^2) * 
        exp((y_node - mu_weights)^2 / (2 * sigma_weights^2))
      
      info_star <- info_star * weight_correction
      
    } else {
      info_star <- 0.0
    }
    
    # Add weighted contribution to Fisher Information
    fisher_info <- fisher_info + info_star * current_weight
  }
  
  # Ensure Fisher Information is non-negative and finite
  if (!is.finite(fisher_info) || fisher_info < 0) {
    warning("Fisher Information calculation resulted in invalid value. Using fallback.")
    fisher_info <- 1.0  # Fallback to unit information
  }
  
  # Ensure we don't have zero information (would cause division by zero later)
  if (fisher_info < .Machine$double.eps) {
    fisher_info <- .Machine$double.eps
  }
  
  # Add attributes for debugging
  attr(fisher_info, "quadrature_info") <- list(
    n_nodes = n_nodes,
    node_range = range(nodes),
    weight_sum = sum(weights),
    mu_weights = mu_weights,
    sigma_weights = sigma_weights
  )
  
  attr(fisher_info, "regime_info") <- list(
    K = K,
    mu_range = range(mu),
    sigma2_range = range(sigma2),
    X_t_prev_entropy = -sum(X_t_prev * log(pmax(X_t_prev, .Machine$double.eps)))
  )
  
  return(fisher_info)
}

#' Calculate transition matrix predictive probabilities for Fisher Information
#'
#' @param p_trans Transition probabilities vector (length K*(K-1))
#' @param X_t_prev Previous filtered probabilities (length K)
#' @param K Number of regimes
#' @return Predictive probabilities (length K)
#' @details
#' Helper function to calculate predictive probabilities used in Fisher Information
#' calculation. This reconstructs the full transition matrix from the off-diagonal
#' probabilities and applies it to the previous filtered probabilities.
#'
#' @examples
#' # Internal helper function - typically not called directly
#' # X_pred <- calculate_predictive_probs(p_trans, X_t_prev, K)
calculate_predictive_probs <- function(p_trans, X_t_prev, K) {
  # Convert to valid transition probabilities and create matrix
  p_trans_valid <- convert_to_valid_probs(p_trans, K)
  transition_mat <- transition_matrix(p_trans_valid, check_validity = FALSE)
  
  # Calculate predictive probabilities
  X_pred <- as.vector(transition_mat %*% X_t_prev)
  
  # Ensure probabilities are valid
  X_pred <- pmax(X_pred, .Machine$double.eps)  # Avoid zeros
  X_pred <- X_pred / sum(X_pred)  # Normalize to sum to 1
  
  return(X_pred)
}

#' Validate Fisher Information calculation inputs
#'
#' @param mu Vector of regime means
#' @param sigma2 Vector of regime variances
#' @param X_t_prev Previous filtered probabilities
#' @param p_trans Transition probabilities
#' @param gh_setup Gauss-Hermite quadrature setup
#' @param K Number of regimes
#' @return TRUE if valid, otherwise throws an error
#' @details
#' Comprehensive validation for all inputs required for Fisher Information calculation.
#' Checks dimensions, ranges, and mathematical validity of all parameters.
validate_fisher_inputs <- function(mu, sigma2, X_t_prev, p_trans, gh_setup, K) {
  # Check basic score inputs
  validate_score_inputs(rep(1, K), 1, X_t_prev, p_trans, K)
  
  # Check quadrature setup
  validate_gauss_hermite_setup(gh_setup)
  
  # Check regime parameters
  if (has_invalid_values(mu)) {
    stop("mu contains invalid values (NA, NaN, or Inf)")
  }
  
  if (has_invalid_values(sigma2)) {
    stop("sigma2 contains invalid values (NA, NaN, or Inf)")
  }
  
  if (any(sigma2 <= 0)) {
    stop("All elements of sigma2 must be positive")
  }
  
  # Check that dimensions are consistent
  if (length(mu) != K || length(sigma2) != K) {
    stop("mu and sigma2 must have length K")
  }
  
  return(TRUE)
}

#' Calculate raw score vector for GAS model
#'
#' @param eta Vector of regime likelihoods at time t (length K)
#' @param tot_lik Total likelihood at time t (scalar)
#' @param X_t_prev Filtered probabilities from previous time step (length K)
#' @param p_trans Current transition probabilities (length K*(K-1))
#' @param K Number of regimes
#' @return Raw score vector (length K*(K-1))
#' @details 
#' Implements equation (13) from Bazzi et al. (2017) for K regimes:
#' r_t = (p(y_t|θ_i) - p(y_t|θ_j)) / p(y_t|θ) * g(f_t, θ, I_{t-1})
#' 
#' The raw score measures how much the observation at time t favors one regime
#' over another, weighted by the impact on transition probabilities.
#' 
#' For a K-regime model, we have K*(K-1) transition probabilities (off-diagonal
#' elements of the transition matrix), so the score vector has K*(K-1) elements.
#'
#' @examples
#' # Calculate raw score for a 3-regime model
#' K <- 3
#' eta <- c(0.1, 0.8, 0.1)  # Regime likelihoods
#' tot_lik <- sum(eta * c(0.3, 0.4, 0.3))  # Total likelihood
#' X_t_prev <- c(0.3, 0.4, 0.3)  # Previous filtered probabilities
#' p_trans <- rep(0.2, 6)  # Transition probabilities
#' 
#' raw_score <- calculate_raw_score_vector(eta, tot_lik, X_t_prev, p_trans, K)
calculate_raw_score_vector <- function(eta, tot_lik, X_t_prev, p_trans, K) {
  # Validate inputs
  if (!is.numeric(eta) || length(eta) != K) {
    stop("eta must be a numeric vector of length K")
  }
  
  if (!is.numeric(tot_lik) || length(tot_lik) != 1 || tot_lik <= 0) {
    stop("tot_lik must be a positive scalar")
  }
  
  if (!is.numeric(X_t_prev) || length(X_t_prev) != K) {
    stop("X_t_prev must be a numeric vector of length K")
  }
  
  if (!is.numeric(p_trans) || length(p_trans) != K*(K-1)) {
    stop("p_trans must be a numeric vector of length K*(K-1)")
  }
  
  if (!is.numeric(K) || K < 2 || K != as.integer(K)) {
    stop("K must be an integer >= 2")
  }
  
  # Check for valid probabilities
  if (any(X_t_prev < 0) || any(X_t_prev > 1) || abs(sum(X_t_prev) - 1) > 1e-10) {
    stop("X_t_prev must be valid probabilities that sum to 1")
  }
  
  if (any(p_trans < 0) || any(p_trans > 1)) {
    stop("p_trans must contain valid probabilities between 0 and 1")
  }
  
  # Initialize score vector
  n_transition <- K * (K - 1)
  raw_score <- numeric(n_transition)
  
  # Calculate likelihood differences and g-vector components
  # Following equation (13) and (14) from Bazzi et al. (2017)
  
  idx <- 1
  for (i in 1:K) {        # From regime i
    for (j in 1:K) {      # To regime j
      if (i != j) {       # Only off-diagonal transitions
        
        # Likelihood difference component (equation 13)
        # This measures how much the current observation favors regime i vs regime j
        likelihood_diff <- (eta[i] - eta[j]) / tot_lik
        
        # G-vector component (equation 14)
        # This measures the impact of filtered probabilities on transition probabilities
        # Following the original 2-regime implementation:
        # g[1] = X_t[from_regime] * p_trans[idx] * (1 - p_trans[idx])
        
        # The g-vector represents the derivative of the transition probability
        # with respect to the time-varying parameter f_{ij,t}
        g_component <- X_t_prev[i] * p_trans[idx] * (1 - p_trans[idx])
        
        # Raw score is the product of likelihood difference and g-component
        raw_score[idx] <- likelihood_diff * g_component
        
        idx <- idx + 1
      }
    }
  }
  
  # Add attributes for debugging and analysis
  attr(raw_score, "likelihood_diffs") <- {
    diffs <- matrix(0, K, K)
    idx <- 1
    for (i in 1:K) {
      for (j in 1:K) {
        if (i != j) {
          diffs[i, j] <- (eta[i] - eta[j]) / tot_lik
          idx <- idx + 1
        }
      }
    }
    diffs
  }
  
  attr(raw_score, "g_components") <- {
    g_comp <- numeric(n_transition)
    idx <- 1
    for (i in 1:K) {
      for (j in 1:K) {
        if (i != j) {
          g_comp[idx] <- X_t_prev[i] * p_trans[idx] * (1 - p_trans[idx])
          idx <- idx + 1
        }
      }
    }
    g_comp
  }
  
  return(raw_score)
}

#' Validate inputs for score vector calculation
#'
#' @param eta Vector of regime likelihoods
#' @param tot_lik Total likelihood
#' @param X_t_prev Previous filtered probabilities
#' @param p_trans Transition probabilities
#' @param K Number of regimes
#' @return TRUE if all inputs are valid, otherwise throws an error
#' @details
#' Helper function to validate all inputs for score vector calculations.
#' Ensures that dimensions are consistent and values are in valid ranges.
#'
#' @examples
#' # This function is typically called internally
#' # validate_score_inputs(eta, tot_lik, X_t_prev, p_trans, K)
validate_score_inputs <- function(eta, tot_lik, X_t_prev, p_trans, K) {
  # Check for NAs or infinite values
  if (has_invalid_values(eta)) {
    stop("eta contains invalid values (NA, NaN, or Inf)")
  }
  
  if (has_invalid_values(tot_lik)) {
    stop("tot_lik contains invalid values (NA, NaN, or Inf)")
  }
  
  if (has_invalid_values(X_t_prev)) {
    stop("X_t_prev contains invalid values (NA, NaN, or Inf)")
  }
  
  if (has_invalid_values(p_trans)) {
    stop("p_trans contains invalid values (NA, NaN, or Inf)")
  }
  
  # Check that eta values are non-negative (likelihoods)
  if (any(eta < 0)) {
    stop("All elements of eta must be non-negative (likelihoods)")
  }
  
  # Check consistency of dimensions
  n_transition_expected <- K * (K - 1)
  if (length(p_trans) != n_transition_expected) {
    stop(sprintf("Length of p_trans (%d) does not match expected K*(K-1) = %d", 
                 length(p_trans), n_transition_expected))
  }
  
  return(TRUE)
}



#' Calculate Fisher Information matrix using robust Gauss-Hermite quadrature
#'
#' @param mu Vector of regime means (length K)
#' @param sigma2 Vector of regime variances (length K)
#' @param X_t_prev Previous filtered probabilities (length K)
#' @param p_trans Current transition probabilities (length K*(K-1))
#' @param gh_setup Gauss-Hermite quadrature setup from setup_gauss_hermite_quadrature()
#' @param K Number of regimes
#' @param min_density_threshold Minimum density threshold for numerical stability (default: 1e-12)
#' @return Fisher Information value (scalar)
calculate_fisher_information_robust <- function(mu, sigma2, X_t_prev, p_trans, gh_setup, K, 
                                                min_density_threshold = 1e-12) {
  .Deprecated(
    msg = paste(
      "calculate_fisher_information_robust() is EXPERIMENTAL.",
      "It uses complex quadrature that can be computationally unstable.",
      "This function is kept for research purposes only."
    )
  )
  
  # Validate inputs (same as before)
  validate_fisher_inputs(mu, sigma2, X_t_prev, p_trans, gh_setup, K)
  
  # Extract quadrature nodes and weights
  nodes <- gh_setup$nodes
  weights <- gh_setup$weights
  n_nodes <- gh_setup$n_nodes
  
  # Initialize Fisher Information accumulator
  fisher_info <- 0.0
  
  # Parameters for the quadrature weight function
  mu_weights <- gh_setup$mu_quad
  sigma_weights <- gh_setup$sigma_quad
  
  # Pre-calculate transition matrix to avoid repeated computation
  p_trans_valid <- convert_to_valid_probs(p_trans, K)
  transition_mat <- transition_matrix(p_trans_valid, check_validity = FALSE)
  X_pred <- as.vector(transition_mat %*% X_t_prev)
  
  # Ensure predicted probabilities are valid
  X_pred <- pmax(X_pred, min_density_threshold)
  X_pred <- X_pred / sum(X_pred)  # Renormalize
  
  # Track numerical issues for diagnostics
  zero_density_count <- 0
  total_nodes <- 0
  
  # Numerical integration using Gauss-Hermite quadrature
  for (j in 1:n_nodes) {
    total_nodes <- total_nodes + 1
    
    # Current quadrature node (y-value)
    y_node <- nodes[j]
    current_weight <- weights[j]
    
    # Calculate regime densities at this node
    regime_densities <- numeric(K)
    for (k in 1:K) {
      regime_densities[k] <- dnorm(y_node, mu[k], sqrt(sigma2[k]))
    }
    
    # Check for invalid regime densities
    if (any(!is.finite(regime_densities)) || any(regime_densities < 0)) {
      zero_density_count <- zero_density_count + 1
      next  # Skip this node
    }
    
    # Calculate total density (denominator)
    total_density <- sum(regime_densities * X_pred)
    
    # Apply minimum threshold for numerical stability
    if (total_density < min_density_threshold) {
      zero_density_count <- zero_density_count + 1
      next  # Skip this node instead of using invalid density
    }
    
    # Calculate Fisher Information components
    info_star <- 0.0
    
    # Calculate all pairwise likelihood differences
    for (i in 1:K) {
      for (l in 1:K) {
        if (i != l) {  # Only consider different regimes
          
          # Likelihood difference between regimes i and l
          likelihood_diff <- regime_densities[i] - regime_densities[l]
          
          # Contribution to Fisher Information
          info_contribution <- (likelihood_diff^2) / total_density
          
          # Weight by the probability of being in regime i
          regime_weight <- X_pred[i]
          
          info_star <- info_star + regime_weight * info_contribution
        }
      }
    }
    
    # Apply the quadrature weight function correction
    weight_correction <- sqrt(2 * pi * sigma_weights^2) * 
      exp((y_node - mu_weights)^2 / (2 * sigma_weights^2))
    
    info_star <- info_star * weight_correction
    
    # Check for invalid contribution
    if (!is.finite(info_star) || info_star < 0) {
      zero_density_count <- zero_density_count + 1
      next  # Skip invalid contributions
    }
    
    # Add weighted contribution to Fisher Information
    fisher_info <- fisher_info + info_star * current_weight
  }
  
  # Diagnostic information
  valid_nodes <- total_nodes - zero_density_count
  if (zero_density_count > 0) {
    skip_percentage <- 100 * zero_density_count / total_nodes
    if (skip_percentage > 20) {  # Only warn if >20% of nodes skipped
      warning(paste("Skipped", zero_density_count, "out of", total_nodes, 
                    "quadrature nodes (", round(skip_percentage, 1), 
                    "%) due to numerical issues in Fisher Information calculation"))
    }
  }
  
  # Ensure Fisher Information is valid
  if (!is.finite(fisher_info) || fisher_info <= 0) {
    if (valid_nodes == 0) {
      warning("All quadrature nodes produced invalid densities. Using minimal Fisher Information.")
      fisher_info <- min_density_threshold
    } else {
      warning("Fisher Information calculation resulted in invalid value. Using fallback.")
      fisher_info <- 1.0  # Fallback to unit information
    }
  }
  
  # Ensure we don't have effectively zero information
  fisher_info <- max(fisher_info, min_density_threshold)
  
  # Add enhanced attributes for debugging
  attr(fisher_info, "quadrature_info") <- list(
    n_nodes = n_nodes,
    valid_nodes = valid_nodes,
    skipped_nodes = zero_density_count,
    skip_percentage = 100 * zero_density_count / total_nodes,
    node_range = range(nodes),
    weight_sum = sum(weights),
    mu_weights = mu_weights,
    sigma_weights = sigma_weights
  )
  
  attr(fisher_info, "regime_info") <- list(
    K = K,
    mu_range = range(mu),
    sigma2_range = range(sigma2),
    X_t_prev_entropy = -sum(X_t_prev * log(pmax(X_t_prev, min_density_threshold))),
    X_pred_entropy = -sum(X_pred * log(pmax(X_pred, min_density_threshold))),
    min_regime_density = min(regime_densities),
    max_regime_density = max(regime_densities)
  )
  
  return(fisher_info)
}



#' Calculate properly scaled GAS score vector with robust Fisher Information
#'
#' @param y_t Current observation at time t
#' @param mu Vector of regime means (length K)
#' @param sigma2 Vector of regime variances (length K)
#' @param X_t_lag Predicted probabilities before observation (length K)
#' @param X_t_prev Filtered probabilities from previous time step (length K)
#' @param p_trans Current transition probabilities (length K*(K-1))
#' @param gh_setup Gauss-Hermite quadrature setup from setup_gauss_hermite_quadrature()
#' @param K Number of regimes
#' @param scaling_method Method for scaling ("moore_penrose", "simple", "normalized", or "original")
#' @param min_density_threshold Minimum density threshold for numerical stability (default: 1e-12)
#' @return Properly scaled score vector (length K*(K-1))
calculate_gas_score_robust <- function(y_t, mu, sigma2, X_t_lag, X_t_prev, p_trans, 
                                       gh_setup, K, scaling_method = "simple",
                                       min_density_threshold = 1e-12) {
  .Deprecated(
    msg = paste(
      "calculate_gas_score_robust() is EXPERIMENTAL and NOT suitable for parameter estimation.",
      "It uses complex Moore-Penrose scaling that causes optimization instability.",
      "For estimation, use calculate_gas_score() with scaling_method='simple' instead.",
      "This function is kept for research/comparison purposes only."
    )
  )
  
  # Validate all inputs (same validation as original)
  if (!is.numeric(y_t) || length(y_t) != 1) {
    stop("y_t must be a numeric scalar")
  }
  
  if (has_invalid_values(y_t)) {
    stop("y_t contains invalid values (NA, NaN, or Inf)")
  }
  
  validate_fisher_inputs(mu, sigma2, X_t_prev, p_trans, gh_setup, K)
  
  if (!is.numeric(X_t_lag) || length(X_t_lag) != K) {
    stop("X_t_lag must be a numeric vector of length K")
  }
  
  if (any(X_t_lag < 0) || any(X_t_lag > 1) || abs(sum(X_t_lag) - 1) > 1e-10) {
    stop("X_t_lag must be valid probabilities that sum to 1")
  }
  
  scaling_method <- match.arg(scaling_method, c("moore_penrose", "simple", "normalized", "original"))
  
  # Step 1: Calculate regime likelihoods for the current observation
  eta <- numeric(K)
  for (k in 1:K) {
    eta[k] <- dnorm(y_t, mu[k], sqrt(sigma2[k]))
  }
  
  # Calculate total likelihood
  tot_lik <- sum(eta * X_t_lag)
  
  # Protect against numerical issues
  if (tot_lik <= min_density_threshold || !is.finite(tot_lik)) {
    warning("Total likelihood is too small or invalid. Using fallback calculation.")
    # In case of numerical issues, return zero score
    return(rep(0, K*(K-1)))
  }
  
  # Step 2: Calculate raw score vector
  raw_score <- calculate_raw_score_vector(eta, tot_lik, X_t_prev, p_trans, K)
  
  # If raw score is essentially zero, return it as-is
  if (all(abs(raw_score) < min_density_threshold)) {
    attr(raw_score, "calculation_info") <- list(
      method = "zero_score",
      y_t = y_t,
      tot_lik = tot_lik,
      fisher_info = NA,
      scaling_method = scaling_method
    )
    return(raw_score)
  }
  
  # Step 3: Calculate Fisher Information using ROBUST version
  fisher_info <- calculate_fisher_information_robust(
    mu, sigma2, X_t_prev, p_trans, gh_setup, K, min_density_threshold
  )
  
  # Step 4: Apply appropriate scaling method
  if (scaling_method == "original") {
    # Use the original 2-regime methodology
    g_vector <- attr(raw_score, "g_components")
    if (is.null(g_vector)) {
      # Recalculate g_vector if not available
      g_vector <- numeric(K*(K-1))
      idx <- 1
      for (i in 1:K) {
        for (j in 1:K) {
          if (i != j) {
            g_vector[idx] <- X_t_prev[i] * p_trans[idx] * (1 - p_trans[idx])
            idx <- idx + 1
          }
        }
      }
    }
    scaled_score <- calculate_original_scaling(raw_score, fisher_info, g_vector)
  } else {
    # Use Moore-Penrose or other modern scaling methods
    scaled_score <- apply_moore_penrose_scaling(raw_score, fisher_info, scaling_method)
  }
  
  # Step 5: Final validation and cleanup
  if (has_invalid_values(scaled_score)) {
    warning("Scaled score contains invalid values. Using conservative fallback.")
    # Conservative fallback: small scaled version of raw score
    score_norm <- sqrt(sum(raw_score^2))
    if (score_norm > 0) {
      scaled_score <- raw_score / score_norm * 0.01  # Very conservative scaling
    } else {
      scaled_score <- rep(0, K*(K-1))
    }
  }
  
  # Add comprehensive attributes for debugging and analysis
  attr(scaled_score, "calculation_info") <- list(
    y_t = y_t,
    tot_lik = tot_lik,
    fisher_info = fisher_info,
    scaling_method = scaling_method,
    raw_score_norm = sqrt(sum(raw_score^2)),
    scaled_score_norm = sqrt(sum(scaled_score^2)),
    regime_likelihoods = eta,
    most_likely_regime = which.max(eta),
    likelihood_spread = max(eta) - min(eta),
    min_density_threshold = min_density_threshold
  )
  
  # Add Fisher Information diagnostics
  fisher_quad_info <- attr(fisher_info, "quadrature_info")
  if (!is.null(fisher_quad_info)) {
    attr(scaled_score, "fisher_diagnostics") <- fisher_quad_info
  }
  
  attr(scaled_score, "regime_diagnostics") <- list(
    K = K,
    filtered_prob_entropy = -sum(X_t_prev * log(pmax(X_t_prev, min_density_threshold))),
    predicted_prob_entropy = -sum(X_t_lag * log(pmax(X_t_lag, min_density_threshold))),
    max_transition_prob = max(p_trans),
    min_transition_prob = min(p_trans)
  )
  
  return(scaled_score)
}
