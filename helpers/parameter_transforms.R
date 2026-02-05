#' Parameter Transformation Functions for Multi-Regime TVTP Models
#' 
#' This file provides robust parameter handling with attribute-based metadata
#' system that eliminates the need for fragile length-based inference.
#'
#' PARAMETER VECTOR STRUCTURE AND ATTRIBUTES
#' =========================================
#' 
#' All parameter vectors in this system carry metadata as R attributes:
#'
#' CORE ATTRIBUTES (required for all parameter vectors):
#' - K: Number of regimes (integer, >= 2)
#' - model_type: One of "constant", "tvp", "exogenous", "gas"
#' - diag_probs: TRUE for diagonal transition probabilities, FALSE for off-diagonal
#' - equal_variances: TRUE for single shared variance, FALSE for K separate variances
#'
#' STATE TRACKING ATTRIBUTES:
#' - parameterization: "natural" (original space) or "transformed" (unconstrained space)
#' - parameter_structure: Named list with start/end indices for each component
#' - original_length: Length before transformations (for validation)
#'
#' PARAMETER VECTOR LAYOUT:
#' ========================
#' 
#' All models follow this structure in natural space:
#' [mu_1, ..., mu_K, sigma2_components, transition_components, time_varying_components]
#'
#' Where:
#' - mu_components: Always K mean parameters
#' - sigma2_components: K parameters (equal_variances=FALSE) or 1 parameter (equal_variances=TRUE)
#' - transition_components: K parameters (diag_probs=TRUE) or K*(K-1) parameters (diag_probs=FALSE)
#' - time_varying_components: Depends on model_type
#'   * "constant": none
#'   * "tvp", "exogenous": A coefficients (matching transition_components length)
#'   * "gas": A coefficients + B coefficients (both matching transition_components length)
#'
#' DIAGONAL VS OFF-DIAGONAL TRANSITION PROBABILITIES:
#' ==================================================
#' 
#' diag_probs = TRUE (diagonal parameterization):
#'   - K parameters representing p11, p22, ..., pKK
#'   - Off-diagonal probabilities equal within each row: pij = (1-pii)/(K-1) for i≠j
#'   - Compatible with original simulation.R implementation for K=2
#'
#' diag_probs = FALSE (off-diagonal parameterization):  
#'   - K*(K-1) parameters representing all off-diagonal transition probabilities
#'   - Diagonal probabilities calculated as: pii = 1 - sum(pij for j≠i)
#'   - Default for new implementation, supports arbitrary transition structures
#'
#' EXAMPLES:
#' =========
#' 
#' # 2-regime constant model with diagonal probs (matches original simulation.R)
#' par <- c(-1, 1, 0.5, 0.5, 0.8, 0.9)  # mu1, mu2, sig1, sig2, p11, p22
#' par <- set_parameter_attributes(par, K=2, model_type="constant", 
#'                                 diag_probs=TRUE, equal_variances=FALSE)
#'
#' # 3-regime TVP model with off-diagonal probs and equal variances
#' par <- c(-2, 0, 2, 0.5, 0.1,0.1,0.1,0.1,0.1,0.1, 0.05,0.05,0.05,0.05,0.05,0.05)
#' par <- set_parameter_attributes(par, K=3, model_type="tvp",
#'                                 diag_probs=FALSE, equal_variances=TRUE)

# Load required helper functions
source("helpers/utility_functions.R")

#' Set parameter vector attributes with validation
#'
#' @param par Numeric parameter vector
#' @param K Number of regimes
#' @param model_type Model type ("constant", "tvp", "exogenous", "gas")
#' @param diag_probs TRUE for diagonal transition probs, FALSE for off-diagonal
#' @param equal_variances TRUE for shared variance, FALSE for separate variances
#' @param parameterization "natural" or "transformed" space
#' @return Parameter vector with proper attributes set
#' @details
#' This is the main function for creating properly attributed parameter vectors.
#' It validates the configuration and vector length, then sets all required attributes.
set_parameter_attributes <- function(par, K, model_type = c("constant", "tvp", "exogenous", "gas"),
                                     diag_probs = FALSE, equal_variances = FALSE, 
                                     parameterization = c("natural", "transformed")) {
  
  # Input validation
  if (!is.numeric(par) || length(par) == 0) {
    stop("Parameter vector must be a non-empty numeric vector")
  }
  if (!is.numeric(K) || K < 2 || K != round(K)) {
    stop("K must be an integer >= 2")
  }
  
  model_type <- match.arg(model_type)
  parameterization <- match.arg(parameterization)
  
  # Validate vector length against configuration
  expected_length <- calculate_expected_length(K, model_type, diag_probs, equal_variances)
  if (length(par) != expected_length) {
    stop(sprintf("Parameter vector length (%d) doesn't match expected length (%d) for configuration: K=%d, model_type='%s', diag_probs=%s, equal_variances=%s",
                 length(par), expected_length, K, model_type, diag_probs, equal_variances))
  }
  
  # Set core attributes
  attr(par, "K") <- as.integer(K)
  attr(par, "model_type") <- model_type
  attr(par, "diag_probs") <- as.logical(diag_probs)
  attr(par, "equal_variances") <- as.logical(equal_variances)
  attr(par, "parameterization") <- parameterization
  
  # Calculate and set parameter structure
  structure <- calculate_parameter_structure(K, model_type, diag_probs, equal_variances)
  attr(par, "parameter_structure") <- structure
  
  # Set metadata
  attr(par, "original_length") <- length(par)
  
  return(par)
}

#' Calculate expected parameter vector length for given configuration
#'
#' @param K Number of regimes  
#' @param model_type Model type
#' @param diag_probs Diagonal transition probabilities flag
#' @param equal_variances Equal variances flag
#' @return Expected parameter vector length
calculate_expected_length <- function(K, model_type, diag_probs, equal_variances) {
  
  # Mean parameters: always K
  n_mu <- K
  
  # Variance parameters: K or 1 depending on equal_variances
  n_sigma2 <- ifelse(equal_variances, 1, K)
  
  # Transition parameters: K (diagonal) or K*(K-1) (off-diagonal)
  n_trans <- ifelse(diag_probs, K, K * (K - 1))
  
  # Time-varying parameters depend on model type
  n_timevarying <- switch(model_type,
                          "constant" = 0,
                          "tvp" = n_trans,      # A coefficients matching transition structure
                          "exogenous" = n_trans, # A coefficients matching transition structure  
                          "gas" = 2 * n_trans   # A + B coefficients matching transition structure
  )
  
  return(n_mu + n_sigma2 + n_trans + n_timevarying)
}

#' Calculate parameter structure indices for given configuration
#'
#' @param K Number of regimes
#' @param model_type Model type  
#' @param diag_probs Diagonal transition probabilities flag
#' @param equal_variances Equal variances flag
#' @return Named list with start/end indices for each parameter component
calculate_parameter_structure <- function(K, model_type, diag_probs, equal_variances) {
  
  structure <- list()
  current_idx <- 1
  
  # Mean parameters
  structure$mu <- list(start = current_idx, end = current_idx + K - 1)
  current_idx <- current_idx + K
  
  # Variance parameters
  n_sigma2 <- ifelse(equal_variances, 1, K)
  structure$sigma2 <- list(start = current_idx, end = current_idx + n_sigma2 - 1)
  current_idx <- current_idx + n_sigma2
  
  # Transition parameters
  n_trans <- ifelse(diag_probs, K, K * (K - 1))
  structure$trans_prob <- list(start = current_idx, end = current_idx + n_trans - 1)
  current_idx <- current_idx + n_trans
  
  # Time-varying parameters
  if (model_type %in% c("tvp", "exogenous", "gas")) {
    structure$A <- list(start = current_idx, end = current_idx + n_trans - 1)
    current_idx <- current_idx + n_trans
  }
  
  if (model_type == "gas") {
    structure$B <- list(start = current_idx, end = current_idx + n_trans - 1)
    current_idx <- current_idx + n_trans
  }
  
  return(structure)
}

#' Extract parameter component using attributes
#'
#' @param par Parameter vector with attributes
#' @param component Component name ("mu", "sigma2", "trans_prob", "A", "B")
#' @return Extracted parameter component
#' @details
#' Uses the parameter_structure attribute to extract components safely.
#' Much more robust than the old index-based extraction.
extract_parameter_component <- function(par, component) {
  
  # Validate parameter vector has required attributes
  validate_parameter_attributes(par)
  
  structure <- attr(par, "parameter_structure")
  
  if (!component %in% names(structure)) {
    model_type <- attr(par, "model_type")
    stop(sprintf("Component '%s' not available for model type '%s'", component, model_type))
  }
  
  comp_info <- structure[[component]]
  return(par[comp_info$start:comp_info$end])
}

#' Validate parameter vector has required attributes
#'
#' @param par Parameter vector to validate
#' @return TRUE if valid, throws error if not
validate_parameter_attributes <- function(par) {
  
  required_attrs <- c("K", "model_type", "diag_probs", "equal_variances", 
                      "parameterization", "parameter_structure")
  
  missing_attrs <- required_attrs[!sapply(required_attrs, function(x) !is.null(attr(par, x)))]
  
  if (length(missing_attrs) > 0) {
    stop(sprintf("Parameter vector missing required attributes: %s", 
                 paste(missing_attrs, collapse = ", ")))
  }
  
  # Validate attribute types and ranges
  K <- attr(par, "K")
  if (!is.integer(K) || K < 2) {
    stop("Attribute 'K' must be an integer >= 2")
  }
  
  model_type <- attr(par, "model_type")
  if (!model_type %in% c("constant", "tvp", "exogenous", "gas")) {
    stop("Attribute 'model_type' must be one of: constant, tvp, exogenous, gas")
  }
  
  parameterization <- attr(par, "parameterization")
  if (!parameterization %in% c("natural", "transformed")) {
    stop("Attribute 'parameterization' must be 'natural' or 'transformed'")
  }
  
  return(TRUE)
}

#' Transform parameters from natural space to unconstrained space
#'
#' @param par Parameter vector with attributes in natural space
#' @return Parameter vector with attributes in transformed space
#' @details
#' Applies transformations suitable for unconstrained optimization:
#' - Log transformation for variances (ensures positivity)
#' - Logit transformation for probabilities and coefficients (ensures [0,1] range)
#' - Preserves all attributes, updates parameterization to "transformed"
transform_parameters <- function(par) {
  
  # Validate input
  validate_parameter_attributes(par)
  
  if (attr(par, "parameterization") != "natural") {
    stop("Can only transform parameters from natural space")
  }
  
  # Extract components using robust attribute-based extraction
  mu <- extract_parameter_component(par, "mu")
  sigma2 <- extract_parameter_component(par, "sigma2")
  trans_prob <- extract_parameter_component(par, "trans_prob")
  
  model_type <- attr(par, "model_type")
  
  # Validate bounds before transformation
  if (any(sigma2 <= 0)) {
    stop("Variance parameters must be positive for transformation")
  }
  if (any(trans_prob <= 0) || any(trans_prob >= 1)) {
    stop("Transition probabilities must be strictly between 0 and 1")
  }
  
  # Apply transformations
  mu_t <- mu  # Means unchanged
  sigma2_t <- log(sigma2)  # Log transformation
  trans_prob_t <- logit(trans_prob)  # Logit transformation
  
  # Build transformed vector
  par_t <- c(mu_t, sigma2_t, trans_prob_t)
  
  # Handle time-varying components if present
  if (model_type %in% c("tvp", "exogenous", "gas")) {
    A <- extract_parameter_component(par, "A")
    if (any(A <= 0) || any(A >= 1)) {
      stop("A coefficients must be strictly between 0 and 1")
    }
    A_t <- logit(A)
    par_t <- c(par_t, A_t)
  }
  
  if (model_type == "gas") {
    B <- extract_parameter_component(par, "B")
    if (any(B <= 0) || any(B >= 1)) {
      stop("B coefficients must be strictly between 0 and 1")
    }
    B_t <- logit(B)
    par_t <- c(par_t, B_t)
  }
  
  # Copy all attributes and update parameterization
  attributes(par_t) <- attributes(par)
  attr(par_t, "parameterization") <- "transformed"
  
  return(par_t)
}

#' Transform parameters from unconstrained space back to natural space
#'
#' @param par_t Parameter vector with attributes in transformed space
#' @return Parameter vector with attributes in natural space
untransform_parameters <- function(par_t) {
  
  # Validate input
  validate_parameter_attributes(par_t)
  
  if (attr(par_t, "parameterization") != "transformed") {
    stop("Can only untransform parameters from transformed space")
  }
  
  # Extract components  
  mu_t <- extract_parameter_component(par_t, "mu")
  sigma2_t <- extract_parameter_component(par_t, "sigma2")
  trans_prob_t <- extract_parameter_component(par_t, "trans_prob")
  
  model_type <- attr(par_t, "model_type")
  
  # Apply inverse transformations
  mu <- mu_t  # Means unchanged
  sigma2 <- exp(sigma2_t)  # Exp transformation
  trans_prob <- logistic(trans_prob_t)  # Logistic transformation
  
  # Build natural vector
  par <- c(mu, sigma2, trans_prob)
  
  # Handle time-varying components if present
  if (model_type %in% c("tvp", "exogenous", "gas")) {
    A_t <- extract_parameter_component(par_t, "A")
    A <- logistic(A_t)
    par <- c(par, A)
  }
  
  if (model_type == "gas") {
    B_t <- extract_parameter_component(par_t, "B")
    B <- logistic(B_t)
    par <- c(par, B)
  }
  
  # Copy all attributes and update parameterization
  attributes(par) <- attributes(par_t)
  attr(par, "parameterization") <- "natural"
  
  return(par)
}

#' Validate transformation round-trip consistency
#'
#' @param par Parameter vector with attributes in natural space
#' @param tolerance Numerical tolerance (default: 1e-10)
#' @return TRUE if consistent, throws error otherwise
validate_transformation_consistency <- function(par, tolerance = 1e-10) {
  
  # Round-trip transformation
  par_t <- transform_parameters(par)
  par_recovered <- untransform_parameters(par_t)
  
  # Compare values (attributes should be identical except parameterization)
  max_diff <- max(abs(par - par_recovered))
  
  if (max_diff > tolerance) {
    stop(sprintf("Transformation consistency failed. Max difference: %e > tolerance: %e",
                 max_diff, tolerance))
  }
  
  cat("✓ Transformation consistency check passed (max diff:", 
      formatC(max_diff, format = "e", digits = 2), ")\n")
  
  return(TRUE)
}

#' Print parameter vector structure information
#'
#' @param par Parameter vector with attributes
print_parameter_info <- function(par) {
  
  validate_parameter_attributes(par)
  
  K <- attr(par, "K")
  model_type <- attr(par, "model_type")
  diag_probs <- attr(par, "diag_probs")
  equal_variances <- attr(par, "equal_variances")
  parameterization <- attr(par, "parameterization")
  structure <- attr(par, "parameter_structure")
  
  cat("Parameter Vector Information\n")
  cat("============================\n")
  cat(sprintf("Length: %d\n", length(par)))
  cat(sprintf("K (regimes): %d\n", K))
  cat(sprintf("Model type: %s\n", model_type))
  cat(sprintf("Transition probabilities: %s\n", ifelse(diag_probs, "diagonal", "off-diagonal")))
  cat(sprintf("Variances: %s\n", ifelse(equal_variances, "equal", "separate")))
  cat(sprintf("Parameterization: %s\n", parameterization))
  cat("\nParameter structure:\n")
  
  for (comp_name in names(structure)) {
    comp_info <- structure[[comp_name]]
    comp_length <- comp_info$end - comp_info$start + 1
    cat(sprintf("  %s: indices %d-%d (length %d)\n", 
                comp_name, comp_info$start, comp_info$end, comp_length))
  }
  
  cat("\nParameter values:\n")
  for (comp_name in names(structure)) {
    comp_values <- extract_parameter_component(par, comp_name)
    cat(sprintf("  %s: [%s]\n", comp_name, paste(round(comp_values, 4), collapse = ", ")))
  }
}

# Convenience wrapper functions for backward compatibility with existing code

#' Extract mean parameters (backward compatibility)
mean_from_par <- function(par) {
  extract_parameter_component(par, "mu")
}

#' Extract variance parameters (backward compatibility) 
sigma2_from_par <- function(par) {
  extract_parameter_component(par, "sigma2")
}

#' Extract transition probabilities (backward compatibility)
transp_from_par <- function(par) {
  extract_parameter_component(par, "trans_prob")
}

#' Extract A coefficients (backward compatibility)
A_from_par <- function(par) {
  extract_parameter_component(par, "A")
}

#' Extract B coefficients (backward compatibility)
B_from_par <- function(par) {
  extract_parameter_component(par, "B")
}
