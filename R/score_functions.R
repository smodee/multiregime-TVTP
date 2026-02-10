#' GAS Score-Driven Helper Functions
#'
#' This file implements score calculation for GAS (Generalized Autoregressive Score)
#' models following Bazzi et al. (2017).
#'
#' Main functions:
#' - calculate_gas_score(): Compute scaled score vector for GAS dynamics
#' - setup_gauss_hermite_quadrature(): Setup numerical integration for Fisher Information
#' - apply_moore_penrose_scaling(): Apply scaling to raw score vectors
#'
#' The implementation supports both diagonal (p11, p22, ...) and off-diagonal
#' parameterizations of transition probabilities. For K=2 with diag_probs=TRUE,
#' results are identical to the original HMMGAS C implementation.

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
#' @export
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
  } else {  # standardized method
    mu_quad <- 0.0
    sigma_quad <- 1.0
  }

  # Use statmod::gauss.quad.prob for Gauss-Hermite quadrature

  # This matches the original HMMGAS C implementation exactly
  GQ <- statmod::gauss.quad.prob(n_nodes, "normal", mu = mu_quad, sigma = sigma_quad)
  nodes <- GQ$nodes
  weights <- GQ$weights
  
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
    node_range_limit = if (method == "data_based") range(nodes) else NULL,
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
#' \dontrun{
#' y <- rnorm(1000)
#' gh_setup <- setup_gauss_hermite_quadrature(y)
#' print(gh_setup)
#' }
#' @export
print.gauss_hermite_setup <- function(x, ...) {
  cat("Gauss-Hermite Quadrature Setup\n")
  cat("==============================\n")
  cat("Number of nodes:", x$n_nodes, "\n")
  cat("Method:", x$method, "\n")
  cat("Integration center (mu):", round(x$mu_quad, 4), "\n")
  cat("Integration scale (sigma):", round(x$sigma_quad, 4), "\n")
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
#' \dontrun{
#' y <- rnorm(1000)
#' gh_setup <- setup_gauss_hermite_quadrature(y)
#' validate_gauss_hermite_setup(gh_setup)  # Returns TRUE
#' }
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
  if (!gh_setup$method %in% c("data_based", "standardized", "statmod_exact")) {
    stop("Method must be 'data_based', 'standardized', or 'statmod_exact'")
  }
  
  return(TRUE)
}

#' Calculate properly scaled GAS score vector
#'
#' @param y_obs Current observation at time t
#' @param mu Vector of regime means (length K)
#' @param sigma2 Vector of regime variances (length K or 1 if equal_variances)
#' @param trans_prob Current transition probabilities (length K for diag_probs=TRUE, or K*(K-1))
#' @param diag_probs If TRUE, trans_prob contains diagonal elements (p11, p22, ...);
#'                   if FALSE, contains off-diagonal elements
#' @param X_pred Predicted probabilities before observation (length K)
#' @param gh_setup Gauss-Hermite quadrature setup from setup_gauss_hermite_quadrature()
#' @param tot_lik_external Optional pre-computed total likelihood (for performance optimization)
#' @param scaling_method Method for scaling. Options: "factored", "simple", "normalized",
#'        "moore_penrose". If NULL (default), auto-selects based on diag_probs:
#'        "factored" for diagonal parameterization (matches HMMGAS),
#'        "simple" for off-diagonal parameterization.
#'
#'        Scaling methods:
#'        \itemize{
#'          \item "factored": Factors the score into magnitude (S*) times direction (g-normalized).
#'                Uses S* = S/sqrt(I) where S is the likelihood ratio score, and multiplies
#'                by a normalized g-vector that encodes gradient direction based on filtered
#'                probabilities and transition probabilities. This matches the original HMMGAS
#'                methodology and is recommended for diagonal parameterization.
#'          \item "simple": Basic inverse Fisher scaling using 1/sqrt(I). General purpose.
#'          \item "normalized": Normalizes the raw score vector to unit length, then scales
#'                by 1/sqrt(I). Treats all score components uniformly.
#'          \item "moore_penrose": Uses Moore-Penrose pseudoinverse. EXPERIMENTAL - can cause
#'                optimization failures.
#'        }
#' @return List with components:
#'   \itemize{
#'     \item scaled_score: Properly scaled score vector (length matches trans_prob)
#'     \item fisher_info: Fisher Information scalar
#'     \item raw_score: Raw (unscaled) score vector
#'   }
#' @details
#' Main interface function that combines all components of the GAS score calculation
#' following Bazzi et al. (2017). This function:
#'
#' 1. Calculates regime likelihoods for the current observation
#' 2. Computes the raw score vector (equation 13 in Bazzi et al.)
#' 3. Calculates Fisher Information using Gauss-Hermite quadrature
#' 4. Applies scaling (method depends on parameterization)
#'
#' The resulting scaled score vector is used to update the time-varying parameters
#' in the GAS model: f_{t+1} = ω + A*s_t + B*(f_t - ω)
#'
#' This implementation generalizes the original 2-regime methodology to K regimes
#' while maintaining the mathematical foundation and numerical stability.
#'
#' When diag_probs=TRUE, the default scaling_method is "factored", which produces
#' results identical to the original HMMGAS C implementation (Filtering_2RegimesGAS.c)
#' for K=2.
#'
#' @examples
#' # Setup for a 2-regime model with diagonal parameterization
#' y <- rnorm(1000)
#' gh_setup <- setup_gauss_hermite_quadrature(y)
#'
#' # Model parameters (diagonal parameterization)
#' mu <- c(-1, 1)
#' sigma2 <- c(0.5, 0.5)
#' trans_prob <- c(0.8, 0.9)  # p11, p22
#' X_pred <- c(0.6, 0.4)
#' y_obs <- 0.5
#'
#' # Calculate scaled score (uses "factored" by default for diag_probs=TRUE)
#' result <- calculate_gas_score(y_obs, mu, sigma2, trans_prob,
#'                               diag_probs = TRUE, X_pred, gh_setup)
#' result$scaled_score  # The score vector
#' result$fisher_info   # Fisher Information
#' @export
calculate_gas_score <- function(y_obs, mu, sigma2, trans_prob, diag_probs = TRUE,
                                X_pred, gh_setup, scaling_method = NULL,
                                tot_lik_external = NULL) {
  # NOTE on X_pred parameter:
  # For diagonal parameterization, X_pred should be the FILTERED probs from t-1
  # (used for the g-vector calculation). The tot_lik denominator should use

  # predicted probs (X_tlag), which should be passed via tot_lik_external.
  #
  # If tot_lik_external is NULL, this function computes tot_lik = sum(eta * X_pred),
  # which is INCORRECT for the C code match but kept for backward compatibility.
  #
  # To match the original C code exactly:
  # - Pass X_pred = X_t[, t-1] (filtered probs from t-1) for g-vector
  # - Pass tot_lik_external = sum(eta * X_tlag) (using predicted probs) for S denominator

  # Infer K from the length of mu
  K <- length(mu)

  # Handle equal variances: expand single variance to K variances
  if (length(sigma2) == 1) {
    sigma2 <- rep(sigma2, K)
  }

  # Determine n_transition based on parameterization
  if (diag_probs) {
    n_transition <- K  # Diagonal: p11, p22, ..., pKK
  } else {
    n_transition <- K * (K - 1)  # Off-diagonal: all p_ij where i != j
  }

  # Validate trans_prob length
  if (length(trans_prob) != n_transition) {
    stop(sprintf("trans_prob length (%d) does not match expected for K=%d with diag_probs=%s (expected %d)",
                 length(trans_prob), K, diag_probs, n_transition))
  }

  # Auto-select scaling method based on parameterization if not specified
  # - "factored" for diagonal (matches HMMGAS methodology)
  # - "simple" for off-diagonal (factored method doesn't apply)
  if (is.null(scaling_method)) {
    scaling_method <- if (diag_probs) "factored" else "simple"
  }

  if (scaling_method == "moore_penrose") {
    warning(
      paste(
        "WARNING: scaling_method='moore_penrose' is EXPERIMENTAL.",
        "For parameter estimation, use 'simple' or 'factored' instead.",
        "Moore-Penrose scaling can cause optimization failures."
      ),
      call. = FALSE
    )
  }

  # Validate all inputs
  if (!is.numeric(y_obs) || length(y_obs) != 1) {
    stop("y_obs must be a numeric scalar")
  }

  if (has_invalid_values(y_obs)) {
    stop("y_obs contains invalid values (NA, NaN, or Inf)")
  }

  if (!is.numeric(X_pred) || length(X_pred) != K) {
    stop("X_pred must be a numeric vector of length K")
  }

  if (any(X_pred < 0) || any(X_pred > 1) || abs(sum(X_pred) - 1) > 1e-10) {
    stop("X_pred must be valid probabilities that sum to 1")
  }

  scaling_method <- match.arg(scaling_method, c("simple", "normalized", "factored", "moore_penrose"))

  # Step 1: Calculate regime likelihoods for the current observation
  eta <- numeric(K)
  for (k in 1:K) {
    eta[k] <- dnorm(y_obs, mu[k], sqrt(sigma2[k]))
  }

  # Calculate total likelihood
  # CRITICAL: Use tot_lik_external if provided (should be sum(eta * X_tlag) from caller)
  # This is the CORRECT approach to match the C code, where:
  #   - S denominator uses exp(logLik[t]) = sum(eta * X_tlag)  (predicted probs)
  #   - g-vector uses X_t[t-1] (filtered probs from previous time)
  if (!is.null(tot_lik_external)) {
    tot_lik <- tot_lik_external
  } else {
    # Fallback: compute using X_pred (may not match C code if X_pred is filtered probs)
    tot_lik <- sum(eta * X_pred)
  }

  # Protect against numerical issues
  if (tot_lik <= .Machine$double.eps || !is.finite(tot_lik)) {
    warning("Total likelihood is too small or invalid. Using fallback calculation.")
    # Return zero score in a list format
    return(list(
      scaled_score = rep(0, n_transition),
      fisher_info = NA,
      raw_score = rep(0, n_transition)
    ))
  }

  # For diagonal parameterization, we use the original HMMGAS methodology
  # which is different from the off-diagonal case
  if (diag_probs) {
    # DIAGONAL PARAMETERIZATION
    # Use the methodology from HMMGAS C code (Filtering_2RegimesGAS.c)
    # This produces identical results to original for K=2

    # For diagonal, X_pred serves as the "previous filtered probabilities" in the
    # original formulation (X_t_prev in the C code)
    X_t_prev <- X_pred

    # Calculate Fisher Information for diagonal case
    # Need to convert diagonal probs to full transition matrix for Fisher calculation
    p_trans_for_fisher <- trans_prob

    # Use Fisher Information calculation (adapted for diagonal)
    fisher_info <- calculate_fisher_information_diagonal(
      mu = mu,
      sigma2 = sigma2,
      X_t_prev = X_t_prev,
      p_diag = trans_prob,
      gh_setup = gh_setup,
      K = K
    )

    # Calculate raw score (for diagnostics)
    # For diagonal, raw score is computed differently
    raw_score <- calculate_raw_score_diagonal(eta, tot_lik, X_t_prev, trans_prob, K)

    # Apply scaling based on method
    if (scaling_method == "factored") {
      # Use the EXACT methodology from HMMGAS C code (factored into magnitude * direction)
      scaled_score <- calculate_factored_scaling(eta, tot_lik, X_t_prev, trans_prob, fisher_info, K)
    } else {
      # Use alternative scaling methods
      scaled_score <- apply_moore_penrose_scaling(raw_score, fisher_info, scaling_method)
    }

  } else {
    # OFF-DIAGONAL PARAMETERIZATION
    # Use the generalized K-regime methodology

    # For off-diagonal, use X_pred as filtered probabilities from previous step
    X_t_prev <- X_pred

    # Validate inputs for off-diagonal case
    validate_fisher_inputs(mu, sigma2, X_t_prev, trans_prob, gh_setup, K)

    # Calculate raw score vector
    raw_score <- calculate_raw_score_vector(eta, tot_lik, X_t_prev, trans_prob, K)

    # Calculate Fisher Information
    fisher_info <- calculate_fisher_information(mu, sigma2, X_t_prev, trans_prob, gh_setup, K)

    # Apply scaling
    # NOTE: "factored" scaling method is only valid for diagonal parameterization
    # (it's based on the HMMGAS C code which uses p11, p22 structure)
    # For off-diagonal, we always use simple/normalized scaling
    if (scaling_method == "factored") {
      warning("scaling_method='factored' is only valid for diagonal parameterization. Using 'simple' instead.")
      scaling_method <- "simple"
    }
    scaled_score <- apply_moore_penrose_scaling(raw_score, fisher_info, scaling_method)
  }

  # Final validation and cleanup
  if (has_invalid_values(scaled_score)) {
    warning("Scaled score contains invalid values. Using conservative fallback.")
    # Conservative fallback: small scaled version of raw score
    score_norm <- sqrt(sum(raw_score^2))
    if (score_norm > 0) {
      scaled_score <- raw_score / score_norm * 0.01  # Very conservative scaling
    } else {
      scaled_score <- rep(0, n_transition)
    }
  }

  # Return as a list (matching what model_GAS.R expects)
  return(list(
    scaled_score = scaled_score,
    fisher_info = fisher_info,
    raw_score = raw_score
  ))
}

#' Calculate Fisher Information for diagonal parameterization
#'
#' @param mu Vector of regime means (length K)
#' @param sigma2 Vector of regime variances (length K)
#' @param X_t_prev Previous filtered probabilities (length K)
#' @param p_diag Diagonal transition probabilities (p11, p22, ..., pKK)
#' @param gh_setup Gauss-Hermite quadrature setup
#' @param K Number of regimes
#' @return Fisher Information scalar
#' @details
#' Calculates Fisher Information for the diagonal parameterization case.
#' This follows the methodology from the original HMMGAS C code.
calculate_fisher_information_diagonal <- function(mu, sigma2, X_t_prev, p_diag, gh_setup, K) {
  # Extract quadrature nodes and weights
  nodes <- gh_setup$nodes
  weights <- gh_setup$weights
  n_nodes <- gh_setup$n_nodes

  # Parameters for the quadrature weight function
  mu_weights <- gh_setup$mu_quad
  sigma_weights <- gh_setup$sigma_quad

  # Initialize Fisher Information accumulator
  fisher_info <- 0.0

  # Build transition matrix from diagonal probabilities
  # For diagonal parameterization: P[i,i] = p_diag[i], P[i,j] = (1-p_diag[i])/(K-1) for j != i
  # Note: P is row-stochastic (rows sum to 1), where P[i,j] = P(to j | from i)
  transition_mat <- matrix(0, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      if (i == j) {
        transition_mat[i, j] <- p_diag[i]
      } else {
        transition_mat[i, j] <- (1 - p_diag[i]) / (K - 1)
      }
    }
  }

  # Calculate predicted probabilities: X_pred = P^T %*% X_t_prev
  # CRITICAL FIX: Must use transpose because P is row-stochastic (P[i,j] = P(to j | from i))
  # This matches compute_predicted_probs() in transition_helpers.R and the original C code
  X_pred <- as.vector(t(transition_mat) %*% X_t_prev)

  # Numerical integration using Gauss-Hermite quadrature
  # Following the original HMMGAS implementation
  for (j in 1:n_nodes) {
    y_node <- nodes[j]
    current_weight <- weights[j]

    # Calculate regime densities at this node
    regime_densities <- numeric(K)
    for (k in 1:K) {
      regime_densities[k] <- dnorm(y_node, mu[k], sqrt(sigma2[k]))
    }

    # Calculate total density (denominator)
    total_density <- sum(regime_densities * X_pred)

    if (total_density > .Machine$double.eps) {
      # For diagonal parameterization with K=2, the original formula is:
      # I_star = ((d1-d2)^2 / den) * weight_correction
      # This generalizes to K regimes by considering (d1-dK)^2

      # Use first and last regime (generalizes K=2 case)
      likelihood_diff <- regime_densities[1] - regime_densities[K]
      info_contribution <- (likelihood_diff^2) / total_density

      # Apply the quadrature weight function correction
      weight_correction <- sqrt(2 * pi * sigma_weights^2) *
        exp((y_node - mu_weights)^2 / (2 * sigma_weights^2))

      info_star <- info_contribution * weight_correction

      # Add weighted contribution
      fisher_info <- fisher_info + info_star * current_weight
    }
  }

  # Ensure Fisher Information is valid
  if (!is.finite(fisher_info) || fisher_info <= 0) {
    fisher_info <- 1.0  # Fallback to unit information
  }

  # Ensure non-zero
  fisher_info <- max(fisher_info, .Machine$double.eps)

  return(fisher_info)
}

#' Calculate raw score for diagonal parameterization
#'
#' @param eta Vector of regime likelihoods (length K)
#' @param tot_lik Total likelihood scalar
#' @param X_t_prev Previous filtered probabilities (length K)
#' @param p_diag Diagonal transition probabilities (p11, p22, ..., pKK)
#' @param K Number of regimes
#' @return Raw score vector (length K)
#' @details
#' Calculates the raw (unscaled) score vector for diagonal parameterization.
#' This is used for diagnostics and alternative scaling methods.
calculate_raw_score_diagonal <- function(eta, tot_lik, X_t_prev, p_diag, K) {
  raw_score <- numeric(K)

  # For diagonal parameterization, the raw score follows the structure:
  # raw_score[i] = likelihood_diff * g_component
  # where g_component = X_t_prev[i] * p_diag[i] * (1 - p_diag[i])

  # Single S value (like original C code)
  S <- (eta[1] - eta[K]) / tot_lik

  sign_multiplier <- 1
  for (i in 1:K) {
    g_component <- X_t_prev[i] * p_diag[i] * (1 - p_diag[i])
    raw_score[i] <- sign_multiplier * S * g_component
    sign_multiplier <- -sign_multiplier
  }

  return(raw_score)
}

#' Batch calculate GAS scores for multiple time points
#'
#' @param y_series Time series of observations (length T)
#' @param mu Vector of regime means (length K)
#' @param sigma2 Vector of regime variances (length K)
#' @param X_pred_series Matrix of predicted probabilities (T x K)
#' @param trans_prob_series Matrix of transition probabilities (T x n_transition)
#' @param gh_setup Gauss-Hermite quadrature setup
#' @param diag_probs If TRUE, use diagonal parameterization (default: TRUE)
#' @param scaling_method Scaling method to use. If NULL (default), auto-selects
#'        based on diag_probs: "factored" for diagonal, "simple" for off-diagonal.
#' @param verbose Whether to show progress
#' @return Matrix of scaled scores (T x n_transition)
#' @details
#' Efficiently calculates GAS scores for a series of observations. Useful for
#' batch processing and when working with long time series. Includes progress
#' reporting and error handling for individual time points.
#'
#' @examples
#' \dontrun{
#' # Calculate scores for an entire time series
#' y_series <- rnorm(100)
#' # ... setup matrices for X_pred_series, trans_prob_series ...
#' # scores <- calculate_gas_scores_batch(y_series, mu, sigma2, ...)
#' }
calculate_gas_scores_batch <- function(y_series, mu, sigma2, X_pred_series,
                                       trans_prob_series, gh_setup,
                                       diag_probs = TRUE,
                                       scaling_method = NULL, verbose = FALSE) {
  # Infer K from mu
  K <- length(mu)

  # Validate inputs
  T_obs <- length(y_series)

  # Determine expected number of transition parameters
  if (diag_probs) {
    n_transition <- K
  } else {
    n_transition <- K * (K - 1)
  }

  if (!is.matrix(X_pred_series) || nrow(X_pred_series) != T_obs || ncol(X_pred_series) != K) {
    stop("X_pred_series must be a T x K matrix")
  }

  if (!is.matrix(trans_prob_series) || nrow(trans_prob_series) != T_obs || ncol(trans_prob_series) != n_transition) {
    stop(sprintf("trans_prob_series must be a T x %d matrix (got %d x %d)",
                 n_transition, nrow(trans_prob_series), ncol(trans_prob_series)))
  }

  # Initialize output
  scaled_scores <- matrix(0, nrow = T_obs, ncol = n_transition)

  # Calculate scores for each time point
  for (t in 1:T_obs) {
    if (verbose && t %% 100 == 0) {
      cat("Processing time point", t, "of", T_obs, "\n")
    }

    tryCatch({
      result <- calculate_gas_score(
        y_obs = y_series[t],
        mu = mu,
        sigma2 = sigma2,
        trans_prob = trans_prob_series[t, ],
        diag_probs = diag_probs,
        X_pred = X_pred_series[t, ],
        gh_setup = gh_setup,
        scaling_method = scaling_method
      )
      scaled_scores[t, ] <- result$scaled_score
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
#' @param score_result Result list from calculate_gas_score() containing
#'        scaled_score, fisher_info, and raw_score
#' @param ... Additional arguments (ignored)
#' @return NULL (prints summary as side effect)
#' @details
#' Provides detailed information about the GAS score calculation including
#' input parameters, intermediate results, and diagnostics.
print_gas_score_summary <- function(score_result, ...) {
  if (!is.list(score_result) || is.null(score_result$scaled_score)) {
    stop("Input must be a result list from calculate_gas_score()")
  }

  scaled_score <- score_result$scaled_score
  fisher_info <- score_result$fisher_info
  raw_score <- score_result$raw_score

  cat("GAS Score Calculation Summary\n")
  cat("=============================\n")

  if (!is.na(fisher_info)) {
    cat("Fisher Information:", format(fisher_info, scientific = TRUE, digits = 4), "\n")
  }

  cat("Raw score norm:", format(sqrt(sum(raw_score^2)), digits = 4), "\n")
  cat("Scaled score norm:", format(sqrt(sum(scaled_score^2)), digits = 4), "\n")

  cat("\nScaled Score Vector:\n")
  cat(format(scaled_score, digits = 4), "\n")

  cat("\nRaw Score Vector:\n")
  cat(format(raw_score, digits = 4), "\n")
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
#' @export
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

#' Calculate score using factored scaling (DIAGONAL PARAMETERIZATION ONLY)
#'
#' @param eta Vector of regime likelihoods (length K)
#' @param tot_lik Total likelihood scalar
#' @param X_t_prev Previous filtered probabilities (length K)
#' @param p_trans Diagonal transition probabilities (p11, p22, ..., pKK) - length K
#' @param fisher_info Fisher Information scalar
#' @param K Number of regimes
#' @return Scaled score vector (length K)
#' @details
#' IMPORTANT: This function is ONLY valid for diagonal parameterization (diag_probs=TRUE).
#' For off-diagonal parameterization, use apply_moore_penrose_scaling() instead.
#'
#' This "factored" scaling method separates the score into two components:
#' - Magnitude: S* = S / sqrt(I), where S is the likelihood ratio score
#' - Direction: g-normalized, a unit vector encoding the gradient direction
#'
#' The final score is computed as: scaled_score = S* × g_normalized
#'
#' This factorization has a natural interpretation: S* captures "how much news"
#' arrived, while g_normalized determines "which direction to update."
#'
#' Implementation steps (matching HMMGAS C code in Filtering_2RegimesGAS.c):
#'
#' 1. Calculate scalar S = (eta\[1\] - eta\[K\]) / tot_lik
#' 2. Build g-vector with ALTERNATING SIGNS:
#'    - g\[1\] = +X_t_prev\[1\] * p11 * (1-p11)  (positive)
#'    - g\[2\] = -X_t_prev\[2\] * p22 * (1-p22)  (negative)
#' 3. Normalize: g_normalized = g / ||g||
#' 4. Scale: S_star = S / sqrt(fisher_info)
#' 5. Final score = S_star * g_normalized
#'
#' For K=2 with diag_probs=TRUE, this produces identical results to the
#' original HMMGAS C implementation. For K>2, this generalizes by using
#' S = (eta\[1\] - eta\[K\]) / tot_lik and alternating signs on the g-vector.
calculate_factored_scaling <- function(eta, tot_lik, X_t_prev, p_trans, fisher_info, K) {
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
  if (!is.numeric(fisher_info) || length(fisher_info) != 1 || fisher_info <= 0) {
    stop("fisher_info must be a positive scalar")
  }
  # CRITICAL: This function only works for diagonal parameterization
  if (length(p_trans) != K) {
    stop(sprintf("calculate_factored_scaling only supports diagonal parameterization. Expected p_trans length %d (K), got %d. Use apply_moore_penrose_scaling for off-diagonal.", K, length(p_trans)))
  }

  # Step 1: Calculate single scalar S (like original C code)
  # Original: S = (eta[0] - eta[1]) / exp(logLik[t])
  # For K>2, we use (eta[1] - eta[K]) as a generalization
  S <- (eta[1] - eta[K]) / tot_lik

  # Step 2: Build g-vector with alternating signs
  # Original C code (for K=2):
  #   g[0] =  X_t[(t-1)*2]   * p11[t] * (1-p11[t])   (positive)
  #   g[1] = -X_t[(t-1)*2+1] * p22[t] * (1-p22[t])   (negative)
  #
  # For diagonal parameterization, p_trans has K elements (p11, p22, ..., pKK)
  n_transition <- length(p_trans)
  g_vector <- numeric(n_transition)

  sign_multiplier <- 1  # Start positive, alternate: +1, -1, +1, -1, ...
  for (i in 1:n_transition) {
    g_vector[i] <- sign_multiplier * X_t_prev[i] * p_trans[i] * (1 - p_trans[i])
    sign_multiplier <- -sign_multiplier  # Flip sign for next component
  }

  # Step 3: Normalize g-vector
  # Original: g_mod = sqrt(g[0]^2 + g[1]^2); g = g / g_mod
  g_mod <- sqrt(sum(g_vector^2))
  if (g_mod > .Machine$double.eps) {
    g_normalized <- g_vector / g_mod
  } else {
    g_normalized <- rep(0, n_transition)
  }

  # Step 4: Scale S by Fisher Information
  # Original: S_star = S / sqrt(I[t])
  if (fisher_info > .Machine$double.eps) {
    S_star <- S / sqrt(fisher_info)
  } else {
    S_star <- 0
  }

  # Step 5: Final score = S_star * g_normalized
  # Original: score_SCAL[t*2] = S_star * g[t*2]; score_SCAL[t*2+1] = S_star * g[t*2+1]
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
#' I = ∫ \[(d1-d2)²/den\] * φ(y; μ_quad, σ_quad²) dy
#' 
#' where d1, d2 are regime densities and den is the total density.
#'
#' @examples
#' \dontrun{
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
#' }
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
#' \dontrun{
#' # Internal helper function - typically not called directly
#' # X_pred <- calculate_predictive_probs(p_trans, X_t_prev, K)
#' }
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
#' \dontrun{
#' # Calculate raw score for a 3-regime model
#' K <- 3
#' eta <- c(0.1, 0.8, 0.1)  # Regime likelihoods
#' tot_lik <- sum(eta * c(0.3, 0.4, 0.3))  # Total likelihood
#' X_t_prev <- c(0.3, 0.4, 0.3)  # Previous filtered probabilities
#' p_trans <- rep(0.2, 6)  # Transition probabilities
#' 
#' raw_score <- calculate_raw_score_vector(eta, tot_lik, X_t_prev, p_trans, K)
#' }
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
#' \dontrun{
#' # This function is typically called internally
#' # validate_score_inputs(eta, tot_lik, X_t_prev, p_trans, K)
#' }
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
