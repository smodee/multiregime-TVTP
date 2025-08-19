# Logit/logistic functions are needed from the utility_functions script
source("helpers/utility_functions.R")

#' Extract mean parameters from a parameter vector
#'
#' @param par Parameter vector
#' @param model_type Type of model ("constant", "tvp", "exogenous", "gas")
#' @return Vector of mean parameters (length K)
#' @details
#' For all model types, the mean parameters are the first K elements of the parameter vector.
#'
#' @examples
#' par <- c(1, 2, 3, 0.1, 0.2, 0.3)  # 3-regime constant model
#' mean_from_par(par, "constant")  # Returns c(1, 2, 3)
#' 
#' par_gas <- c(-1, 1, 0.5, 0.5, 0.3, 0.4, 0.1, 0.2, 0.8, 0.9)  # 2-regime GAS
#' mean_from_par(par_gas, "gas")  # Returns c(-1, 1)
mean_from_par <- function(par, model_type = c("constant", "tvp", "exogenous", "gas")) {
  # Validate inputs
  if (!is.numeric(par) || length(par) == 0) {
    stop("Parameter vector must be a non-empty numeric vector")
  }
  
  model_type <- match.arg(model_type)
  K <- count_regime(par, model_type)
  
  return(par[1:K])
}

#' Extract variance parameters from a parameter vector
#'
#' @param par Parameter vector
#' @param model_type Type of model ("constant", "tvp", "exogenous", "gas")
#' @return Vector of variance parameters (length K)
#' @details
#' For all model types, the variance parameters are elements (K+1) to (2K) of the parameter vector.
#'
#' @examples
#' par <- c(1, 2, 3, 0.1, 0.2, 0.3)  # 3-regime constant model
#' sigma2_from_par(par, "constant")  # Returns c(0.1, 0.2, 0.3)
#' 
#' par_gas <- c(-1, 1, 0.5, 0.5, 0.3, 0.4, 0.1, 0.2, 0.8, 0.9)  # 2-regime GAS
#' sigma2_from_par(par_gas, "gas")  # Returns c(0.5, 0.5)
sigma2_from_par <- function(par, model_type = c("constant", "tvp", "exogenous", "gas")) {
  # Validate inputs
  if (!is.numeric(par) || length(par) == 0) {
    stop("Parameter vector must be a non-empty numeric vector")
  }
  
  model_type <- match.arg(model_type)
  K <- count_regime(par, model_type)
  
  return(par[(K+1):(2*K)])
}

#' Extract transition probability parameters from a parameter vector
#'
#' @param par Parameter vector
#' @param model_type Type of model ("constant", "tvp", "exogenous", "gas")
#' @return Vector of transition probability parameters (length K*(K-1))
#' @details
#' For all models with transitions (including constant models), the transition 
#' probability parameters are elements (2K+1) to (2K + K*(K-1)) of the parameter vector.
#' For constant models, these are fixed transition probabilities.
#' For other models, these are initial/baseline transition probabilities.
#'
#' @examples
#' par_const <- c(1, 2, 0.1, 0.2, 0.3, 0.4)  # 2-regime constant model
#' transp_from_par(par_const, "constant")  # Returns c(0.3, 0.4)
#' 
#' par_gas <- c(-1, 1, 0.5, 0.5, 0.3, 0.4, 0.1, 0.2, 0.8, 0.9)  # 2-regime GAS
#' transp_from_par(par_gas, "gas")  # Returns c(0.3, 0.4)
transp_from_par <- function(par, model_type = c("constant", "tvp", "exogenous", "gas")) {
  # Validate inputs
  if (!is.numeric(par) || length(par) == 0) {
    stop("Parameter vector must be a non-empty numeric vector")
  }
  
  model_type <- match.arg(model_type)
  
  K <- count_regime(par, model_type)
  n_transition <- K * (K - 1)
  
  return(par[(2*K+1):(2*K + n_transition)])
}

#' Extract A coefficients from a parameter vector
#'
#' @param par Parameter vector
#' @param model_type Type of model ("tvp", "exogenous", "gas") - not applicable for "constant"
#' @return Vector of A coefficients (length K*(K-1))
#' @details
#' For models with time-varying transitions, the A coefficients control the
#' sensitivity of transitions to driving variables (observations, exogenous variables, or scores).
#' Not applicable for constant models.
#' 
#' - TVP: sensitivity to lagged observations
#' - Exogenous: sensitivity to exogenous variables  
#' - GAS: sensitivity to scaled scores
#'
#' @examples
#' par_tvp <- c(1, 2, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2)  # 2-regime TVP model
#' A_from_par(par_tvp, "tvp")  # Returns c(0.1, 0.2)
#' 
#' par_gas <- c(-1, 1, 0.5, 0.5, 0.3, 0.4, 0.1, 0.2, 0.8, 0.9)  # 2-regime GAS
#' A_from_par(par_gas, "gas")  # Returns c(0.1, 0.2)
A_from_par <- function(par, model_type = c("tvp", "exogenous", "gas")) {
  # Validate inputs
  if (!is.numeric(par) || length(par) == 0) {
    stop("Parameter vector must be a non-empty numeric vector")
  }
  
  model_type <- match.arg(model_type)
  
  if (model_type == "constant") {
    stop("A coefficients are not applicable for constant models")
  }
  
  K <- count_regime(par, model_type)
  n_transition <- K * (K - 1)
  
  # For TVP and exogenous: A coefficients are after transitions
  # For GAS: A coefficients are after transitions (B coefficients come after A)
  start_idx <- 2*K + n_transition + 1
  end_idx <- 2*K + 2*n_transition
  
  return(par[start_idx:end_idx])
}

#' Extract B coefficients from a parameter vector (GAS model only)
#'
#' @param par Parameter vector
#' @param model_type Must be "gas" - B coefficients only exist in GAS models
#' @return Vector of B coefficients (length K*(K-1))
#' @details
#' B coefficients control the persistence/memory in GAS models. They determine
#' how much the previous time-varying parameter value influences the current value.
#' Only applicable to GAS models.
#'
#' @examples
#' par_gas <- c(-1, 1, 0.5, 0.5, 0.3, 0.4, 0.1, 0.2, 0.8, 0.9)  # 2-regime GAS
#' B_from_par(par_gas, "gas")  # Returns c(0.8, 0.9)
B_from_par <- function(par, model_type = "gas") {
  # Validate inputs
  if (!is.numeric(par) || length(par) == 0) {
    stop("Parameter vector must be a non-empty numeric vector")
  }
  
  if (model_type != "gas") {
    stop("B coefficients are only applicable for GAS models")
  }
  
  K <- count_regime(par, model_type)
  n_transition <- K * (K - 1)
  
  # B coefficients are the last K*(K-1) elements in GAS models
  start_idx <- 2*K + 2*n_transition + 1
  end_idx <- 2*K + 3*n_transition
  
  return(par[start_idx:end_idx])
}

#' Count the number of regimes from parameter vector length
#'
#' @param par Parameter vector
#' @param model_type Type of model ("constant", "tvp", "exogenous", "gas")
#' @return Number of regimes (K)
#' @details
#' Determines the number of regimes based on the parameter vector length
#' and the specified model type. Each model type has a different parameter
#' structure:
#' 
#' - constant: K means + K variances + K*(K-1) transitions = K + K²
#' - tvp: K means + K variances + K*(K-1) transitions + K*(K-1) A = 2K²
#' - exogenous: K means + K variances + K*(K-1) transitions + K*(K-1) A = 2K²
#' - gas: K means + K variances + K*(K-1) transitions + K*(K-1) A + K*(K-1) B = 3K² - K
#'
#' @examples
#' # For a 3-regime constant model (12 parameters: 3 means + 3 variances + 6 transitions)
#' par_const <- c(1, 2, 3, 0.1, 0.2, 0.3, 0.3, 0.4, 0.2, 0.5, 0.1, 0.6)
#' count_regime(par_const, "constant")  # Returns 3
#' 
#' # For a 2-regime GAS model (14 parameters total)
#' par_gas <- c(-1, 1, 0.5, 0.5, 0.3, 0.4, 0.1, 0.2, 0.8, 0.9)
#' count_regime(par_gas, "gas")  # Returns 2
count_regime <- function(par, model_type = c("constant", "tvp", "exogenous", "gas")) {
  # Validate inputs
  if (!is.numeric(par) || length(par) == 0) {
    stop("Parameter vector must be a non-empty numeric vector")
  }
  
  if (has_invalid_values(par)) {
    stop("Parameter vector contains invalid values (NA, NaN, or Inf)")
  }
  
  # Match and validate model type
  model_type <- match.arg(model_type)
  
  # Calculate K based on model type
  K <- switch(model_type,
              "constant" = {
                # K means + K variances + K*(K-1) transitions = K + K²
                # So: K² + K - length(par) = 0
                # Using quadratic formula: K = (-1 + sqrt(1 + 4*length(par)))/2
                discriminant <- 1 + 4*length(par)
                if (discriminant < 0) {
                  stop("Invalid parameter vector length for constant model (discriminant < 0)")
                }
                (-1 + sqrt(discriminant))/2
              },
              "tvp" = {
                # K means + K variances + K*(K-1) transitions + K*(K-1) A = 2K²
                # So: 2K² = length(par) => K = sqrt(length(par)/2)
                sqrt(length(par)/2)
              },
              "exogenous" = {
                # K means + K variances + K*(K-1) transitions + K*(K-1) A = 2K²
                # So: 2K² = length(par) => K = sqrt(length(par)/2)
                sqrt(length(par)/2)
              },
              "gas" = {
                # K means + K variances + K*(K-1) transitions + K*(K-1) A + K*(K-1) B = 3K² - K
                # So: 3K² - K - length(par) = 0
                # Using quadratic formula: K = (1 + sqrt(1 + 12*length(par)))/6
                discriminant <- 1 + 12*length(par)
                if (discriminant < 0) {
                  stop("Invalid parameter vector length for GAS model (discriminant < 0)")
                }
                (1 + sqrt(discriminant))/6
              }
  )
  
  # Validate that K is a positive integer
  if (!is.finite(K) || K <= 0) {
    stop(paste("Invalid parameter vector length", length(par), 
               "for model type", model_type, "- calculated K =", K))
  }
  
  if (abs(K - round(K)) > 1e-10) {
    stop(paste("Invalid parameter vector length", length(par), 
               "for model type", model_type, 
               "- does not correspond to an integer number of regimes.",
               "Calculated K =", K))
  }
  
  K_int <- round(K)
  
  # Additional validation: verify the parameter count is correct
  expected_length <- switch(model_type,
                            "constant" = K_int + K_int^2,
                            "tvp" = 2 * K_int^2,
                            "exogenous" = 2 * K_int^2,
                            "gas" = 3 * K_int^2 - K_int
  )
  
  if (length(par) != expected_length) {
    stop(paste("Parameter vector length", length(par), 
               "does not match expected length", expected_length,
               "for", K_int, "regimes with model type", model_type))
  }
  
  return(K_int)
}

#' Validate parameter vector structure and values
#'
#' @param par Parameter vector to validate
#' @param model_type Type of model ("constant", "tvp", "exogenous", "gas")
#' @param K Optional: number of regimes (will be inferred if not provided)
#' @return TRUE if valid, otherwise throws an error with detailed explanation
#' @details
#' Performs comprehensive validation of parameter vector including:
#' - Correct length for the specified model type
#' - Positive variance parameters
#' - Valid probability ranges for transition parameters (if applicable)
#' - Valid coefficient ranges for A and B parameters (if applicable)
#'
#' @examples
#' # Valid parameter vector for 2-regime constant model
#' par_const <- c(1, 2, 0.1, 0.2, 0.3, 0.4)  # 2 means + 2 variances + 2 transitions
#' validate_parameter_vector(par_const, "constant")  # Returns TRUE
#' 
#' # Valid parameter vector for 2-regime GAS model
#' par_gas <- c(-1, 1, 0.5, 0.5, 0.3, 0.4, 0.1, 0.2, 0.8, 0.9)
#' validate_parameter_vector(par_gas, "gas")  # Returns TRUE
#'
#' @param par Parameter vector to validate
#' @param model_type Type of model ("constant", "tvp", "exogenous", "gas")
#' @param K Optional: number of regimes (will be inferred if not provided)
#' @return TRUE if valid, otherwise throws an error with detailed explanation
validate_parameter_vector <- function(par, model_type = c("constant", "tvp", "exogenous", "gas"), K = NULL) {
  # Validate basic inputs
  if (!is.numeric(par) || length(par) == 0) {
    stop("Parameter vector must be a non-empty numeric vector")
  }
  
  if (has_invalid_values(par)) {
    stop("Parameter vector contains invalid values (NA, NaN, or Inf)")
  }
  
  model_type <- match.arg(model_type)
  
  # Determine or validate K
  if (is.null(K)) {
    tryCatch({
      K <- count_regime(par, model_type)
    }, error = function(e) {
      stop("Invalid parameter vector length for model type '", model_type, "'. ", e$message)
    })
  } else {
    # Validate provided K
    if (!is.numeric(K) || length(K) != 1 || K != as.integer(K) || K < 1) {
      stop("K must be a positive integer")
    }
    K <- as.integer(K)
    
    # Check if parameter vector length matches expected length for this K
    expected_length <- switch(model_type,
                              "constant" = K + K^2,
                              "tvp" = 2 * K^2,
                              "exogenous" = 2 * K^2,
                              "gas" = 3 * K^2 - K
    )
    
    if (length(par) != expected_length) {
      stop("Parameter vector length (", length(par), ") does not match expected length (", 
           expected_length, ") for ", K, " regimes with model type '", model_type, "'")
    }
  }
  
  # Extract and validate variance parameters
  tryCatch({
    sigma2 <- sigma2_from_par(par, model_type)
    if (any(sigma2 <= 0)) {
      bad_variances <- sigma2[sigma2 <= 0]
      stop("All variance parameters must be positive. Found non-positive values: ", 
           paste(round(bad_variances, 6), collapse = ", "))
    }
  }, error = function(e) {
    if (grepl("must be positive", e$message)) {
      stop(e$message)  # Re-throw our own error
    } else {
      stop("Error in variance parameter extraction: ", e$message)
    }
  })
  
  # Validate transition probabilities (all models have them)
  tryCatch({
    trans_prob <- transp_from_par(par, model_type)
    
    # Check for values outside [0,1] range
    if (any(trans_prob < 0) || any(trans_prob > 1)) {
      invalid_probs <- trans_prob[trans_prob < 0 | trans_prob > 1]
      stop("Transition probabilities must be between 0 and 1. Found invalid values: ", 
           paste(round(invalid_probs, 6), collapse = ", "))
    }
  }, error = function(e) {
    if (grepl("must be between 0 and 1", e$message)) {
      stop(e$message)  # Re-throw our own error
    } else {
      stop("Error in transition probability extraction: ", e$message)
    }
  })
  
  # Additional validations for time-varying models
  if (model_type %in% c("tvp", "exogenous", "gas")) {
    
    # Validate A coefficients
    tryCatch({
      A_coeffs <- A_from_par(par, model_type)
      # Note: A coefficients can theoretically be any real number, but warn about extreme values
      if (any(abs(A_coeffs) > 10)) {
        extreme_A <- A_coeffs[abs(A_coeffs) > 10]
        warning("Some A coefficients are quite large (|A| > 10): ", 
                paste(round(extreme_A, 6), collapse = ", "), 
                ". This may cause numerical instability.")
      }
    }, error = function(e) {
      stop("Error in A coefficient extraction: ", e$message)
    })
    
    # Additional GAS-specific validation
    if (model_type == "gas") {
      tryCatch({
        B_coeffs <- B_from_par(par, model_type)
        
        # B coefficients should typically be between 0 and 1 for stability
        if (any(B_coeffs < 0) || any(B_coeffs > 1)) {
          invalid_B <- B_coeffs[B_coeffs < 0 | B_coeffs > 1]
          warning("B coefficients are typically between 0 and 1 for stability. Found: ", 
                  paste(round(invalid_B, 6), collapse = ", "))
        }
        
        # Check for explosive behavior: |A| + |B| should typically be < 1 for each parameter
        A_coeffs <- A_from_par(par, model_type)  # Get A coeffs again
        combined_coeffs <- abs(A_coeffs) + abs(B_coeffs)
        if (any(combined_coeffs > 1.5)) {
          explosive_idx <- which(combined_coeffs > 1.5)
          warning("Some parameter combinations |A| + |B| > 1.5 may cause explosive behavior. ",
                  "Indices: ", paste(explosive_idx, collapse = ", "))
        }
      }, error = function(e) {
        stop("Error in B coefficient extraction: ", e$message)
      })
    }
  }
  
  # If we reach here, all validations passed
  return(TRUE)
}

#' Transform parameters from natural space to unconstrained space
#'
#' @param par Parameters in their natural space
#' @param model_type Type of model ("constant", "tvp", "exogenous", "gas")
#' @return Parameters in unconstrained space suitable for optimization
#' @details
#' Applies transformations to make parameters suitable for unconstrained optimization:
#' - Log transformation for variances (ensures positivity)
#' - Logit transformation for transition probabilities (ensures [0,1] range)
#' - Logit transformation for A and B coefficients in time-varying models (ensures [0,1] range)
#' 
#' All models: means remain unchanged, variances → log(variances), transitions → logit(transitions)
#' Time-varying models additionally: A → logit(A), B → logit(B) (GAS only)
#'
#' @examples
#' # Transform 2-regime constant model parameters
#' par_const <- c(1, 2, 0.1, 0.2, 0.3, 0.4)
#' par_transformed <- transform_parameters(par_const, "constant")
#' 
#' # Transform 2-regime GAS model parameters  
#' par_gas <- c(-1, 1, 0.5, 0.5, 0.3, 0.4, 0.1, 0.2, 0.8, 0.9)
#' par_transformed <- transform_parameters(par_gas, "gas")
transform_parameters <- function(par, model_type = c("constant", "tvp", "exogenous", "gas")) {
  # Validate inputs
  if (!is.numeric(par) || length(par) == 0) {
    stop("Parameter vector must be a non-empty numeric vector")
  }
  
  model_type <- match.arg(model_type)
  
  # Get the number of regimes
  K <- count_regime(par, model_type)
  
  # Extract parameter components
  mu <- mean_from_par(par, model_type)
  sigma2 <- sigma2_from_par(par, model_type)
  
  # Check for invalid values before transformation
  if (any(sigma2 <= 0)) {
    stop("Variance parameters must be positive for transformation. Found: ", 
         paste(sigma2[sigma2 <= 0], collapse = ", "))
  }
  
  # Transform basic parameters (all models have these)
  mu_t <- mu  # Means remain unchanged
  sigma2_t <- log(sigma2)  # Log transformation for variances
  
  # Handle transition probabilities (all models have these)
  trans_prob <- transp_from_par(par, model_type)
  
  # Check transition probability bounds
  if (any(trans_prob <= 0) || any(trans_prob >= 1)) {
    stop("Transition probabilities must be strictly between 0 and 1 for transformation. Found: ",
         paste(trans_prob[trans_prob <= 0 | trans_prob >= 1], collapse = ", "))
  }
  
  trans_prob_t <- logit(trans_prob)  # Logit transformation for probabilities
  
  # Build transformed parameter vector based on model type
  if (model_type == "constant") {
    # Constant model: mu + sigma2 + transitions
    par_transformed <- c(mu_t, sigma2_t, trans_prob_t)
    
  } else if (model_type %in% c("tvp", "exogenous")) {
    # TVP/Exogenous models: mu + sigma2 + transitions + A
    A_coeffs <- A_from_par(par, model_type)
    
    # Check A coefficient bounds  
    if (any(A_coeffs < 0) || any(A_coeffs > 1)) {
      stop("A coefficients must be between 0 and 1 for transformation. Found: ",
           paste(A_coeffs[A_coeffs < 0 | A_coeffs > 1], collapse = ", "))
    }
    
    A_coeffs_t <- logit(A_coeffs)  # Logit transformation for A coefficients
    par_transformed <- c(mu_t, sigma2_t, trans_prob_t, A_coeffs_t)
    
  } else if (model_type == "gas") {
    # GAS model: mu + sigma2 + transitions + A + B
    A_coeffs <- A_from_par(par, model_type)
    B_coeffs <- B_from_par(par, model_type)
    
    # Check A coefficient bounds
    if (any(A_coeffs < 0) || any(A_coeffs > 1)) {
      stop("A coefficients must be between 0 and 1 for transformation. Found: ",
           paste(A_coeffs[A_coeffs < 0 | A_coeffs > 1], collapse = ", "))
    }
    
    # Check B coefficient bounds
    if (any(B_coeffs < 0) || any(B_coeffs > 1)) {
      stop("B coefficients must be between 0 and 1 for transformation. Found: ",
           paste(B_coeffs[B_coeffs < 0 | B_coeffs > 1], collapse = ", "))
    }
    
    A_coeffs_t <- logit(A_coeffs)  # Logit transformation for A coefficients
    B_coeffs_t <- logit(B_coeffs)  # Logit transformation for B coefficients
    par_transformed <- c(mu_t, sigma2_t, trans_prob_t, A_coeffs_t, B_coeffs_t)
  }
  
  # Add attributes to preserve model information
  attr(par_transformed, "model_type") <- model_type
  attr(par_transformed, "K") <- K
  attr(par_transformed, "original_length") <- length(par)
  
  return(par_transformed)
}

#' Transform parameters from unconstrained space back to natural space
#'
#' @param par_t Parameters in unconstrained space (from transform_parameters)
#' @param model_type Type of model ("constant", "tvp", "exogenous", "gas")
#' @return Parameters in their natural space
#' @details
#' Applies inverse transformations to convert parameters back to their natural bounds:
#' - Exp transformation for log-variances (ensures positivity)
#' - Logistic transformation for logit-probabilities (ensures [0,1] range)
#' - Logistic transformation for logit-coefficients (ensures [0,1] range)
#' 
#' This function is the inverse of transform_parameters().
#'
#' @examples
#' # Transform back from unconstrained space
#' par_const_orig <- c(1, 2, 0.1, 0.2, 0.3, 0.4)
#' par_transformed <- transform_parameters(par_const_orig, "constant")
#' par_recovered <- untransform_parameters(par_transformed, "constant")
#' # par_recovered should equal par_const_orig
untransform_parameters <- function(par_t, model_type = c("constant", "tvp", "exogenous", "gas")) {
  # Validate inputs
  if (!is.numeric(par_t) || length(par_t) == 0) {
    stop("Transformed parameter vector must be a non-empty numeric vector")
  }
  
  model_type <- match.arg(model_type)
  
  # Try to get K from attributes first, then infer from length
  K <- attr(par_t, "K")
  if (is.null(K)) {
    # Infer K from parameter vector length and model type
    expected_lengths <- switch(model_type,
                               "constant" = function(K) K + K^2,
                               "tvp" = function(K) 2 * K^2,
                               "exogenous" = function(K) 2 * K^2, 
                               "gas" = function(K) 3 * K^2 - K
    )
    
    # Find K that gives the correct length
    for (K_test in 1:10) {  # Reasonable upper bound
      if (expected_lengths(K_test) == length(par_t)) {
        K <- K_test
        break
      }
    }
    
    if (is.null(K)) {
      stop("Cannot determine number of regimes from parameter vector length ", 
           length(par_t), " for model type '", model_type, "'")
    }
  }
  
  # Calculate parameter structure
  n_transition <- K * (K - 1)
  
  # Extract and untransform parameter components based on model type
  if (model_type == "constant") {
    # Structure: [mu, log_sigma2, logit_trans]
    mu_t <- par_t[1:K]
    sigma2_t <- par_t[(K+1):(2*K)]
    trans_prob_t <- par_t[(2*K+1):(2*K + n_transition)]
    
    # Untransform
    mu <- mu_t  # Means unchanged
    sigma2 <- exp(sigma2_t)  # Exp transformation for variances
    trans_prob <- logistic(trans_prob_t)  # Logistic transformation for probabilities
    
    par_natural <- c(mu, sigma2, trans_prob)
    
  } else if (model_type %in% c("tvp", "exogenous")) {
    # Structure: [mu, log_sigma2, logit_trans, logit_A]
    mu_t <- par_t[1:K]
    sigma2_t <- par_t[(K+1):(2*K)]
    trans_prob_t <- par_t[(2*K+1):(2*K + n_transition)]
    A_coeffs_t <- par_t[(2*K + n_transition + 1):(2*K + 2*n_transition)]
    
    # Untransform
    mu <- mu_t  # Means unchanged
    sigma2 <- exp(sigma2_t)  # Exp transformation for variances
    trans_prob <- logistic(trans_prob_t)  # Logistic transformation for probabilities
    A_coeffs <- logistic(A_coeffs_t)  # Logistic transformation for A coefficients
    
    par_natural <- c(mu, sigma2, trans_prob, A_coeffs)
    
  } else if (model_type == "gas") {
    # Structure: [mu, log_sigma2, logit_trans, logit_A, logit_B]
    mu_t <- par_t[1:K]
    sigma2_t <- par_t[(K+1):(2*K)]
    trans_prob_t <- par_t[(2*K+1):(2*K + n_transition)]
    A_coeffs_t <- par_t[(2*K + n_transition + 1):(2*K + 2*n_transition)]
    B_coeffs_t <- par_t[(2*K + 2*n_transition + 1):(2*K + 3*n_transition)]
    
    # Untransform
    mu <- mu_t  # Means unchanged
    sigma2 <- exp(sigma2_t)  # Exp transformation for variances  
    trans_prob <- logistic(trans_prob_t)  # Logistic transformation for probabilities
    A_coeffs <- logistic(A_coeffs_t)  # Logistic transformation for A coefficients
    B_coeffs <- logistic(B_coeffs_t)  # Logistic transformation for B coefficients
    
    par_natural <- c(mu, sigma2, trans_prob, A_coeffs, B_coeffs)
  }
  
  # Add attributes to preserve information
  attr(par_natural, "model_type") <- model_type
  attr(par_natural, "K") <- K
  attr(par_natural, "transformed_length") <- length(par_t)
  
  return(par_natural)
}

#' Validate transformation round-trip consistency
#'
#' @param par Original parameters in natural space
#' @param model_type Type of model
#' @param tolerance Numerical tolerance for comparison (default: 1e-10)
#' @return TRUE if transformation is consistent, otherwise throws an error
#' @details
#' Tests that transform_parameters() followed by untransform_parameters()
#' recovers the original parameter values within numerical tolerance.
#' Useful for testing and debugging transformation functions.
#'
#' @examples
#' par_test <- c(1, 2, 0.1, 0.2, 0.3, 0.4)
#' validate_transformation_consistency(par_test, "constant")
validate_transformation_consistency <- function(par, model_type, tolerance = 1e-10) {
  # Transform and untransform
  par_transformed <- transform_parameters(par, model_type)
  par_recovered <- untransform_parameters(par_transformed, model_type)
  
  # Check consistency
  max_diff <- max(abs(par - par_recovered))
  
  if (max_diff > tolerance) {
    stop("Transformation consistency check failed. Maximum difference: ", 
         max_diff, " > tolerance: ", tolerance,
         "\nOriginal: ", paste(round(par, 8), collapse = ", "),
         "\nRecovered: ", paste(round(par_recovered, 8), collapse = ", "))
  }
  
  cat("✓ Transformation consistency check passed (max diff:", 
      formatC(max_diff, format = "e", digits = 2), ")\n")
  
  return(TRUE)
}

#' Compress variance parameters for equal variance constraint
#'
#' @param params Full parameter vector in natural space
#' @param K Number of regimes
#' @param model_type Type of model ("constant", "tvp", "exogenous", "gas")
#' @return Compressed parameter vector with single variance parameter
#' @details
#' Converts a parameter vector with K separate variance parameters to one with
#' a single shared variance parameter. Used when equal_variances = TRUE.
#' 
#' The compression uses the mean of the K variance values to create the single
#' shared variance parameter.
#'
#' @examples
#' # TVP model with 3 regimes: [mu1, mu2, mu3, σ²1, σ²2, σ²3, trans..., A...]
#' full_params <- c(-1, 0, 1, 0.5, 0.6, 0.4, rep(0.2, 6), rep(0.1, 6))
#' compressed <- compress_variances(full_params, K = 3, "tvp")
#' # Result: [mu1, mu2, mu3, σ²_avg, trans..., A...]
compress_variances <- function(params, K, model_type = c("constant", "tvp", "exogenous", "gas")) {
  model_type <- match.arg(model_type)
  
  if (length(params) == 0) {
    stop("Parameter vector cannot be empty")
  }
  
  # Extract components
  mu <- params[1:K]
  sigma2 <- params[(K+1):(2*K)]
  
  # Calculate shared variance (using mean of individual variances)
  sigma2_shared <- mean(sigma2)
  
  # Get remaining parameters (transition probs, A, B, etc.)
  if (length(params) > 2*K) {
    remaining <- params[(2*K+1):length(params)]
    compressed_params <- c(mu, sigma2_shared, remaining)
  } else {
    compressed_params <- c(mu, sigma2_shared)
  }
  
  # Add attributes to track the transformation
  attr(compressed_params, "compressed") <- TRUE
  attr(compressed_params, "original_K") <- K
  attr(compressed_params, "model_type") <- model_type
  attr(compressed_params, "original_length") <- length(params)
  
  return(compressed_params)
}

#' Expand single variance parameter to K regime-specific variances
#'
#' @param params Compressed parameter vector with single variance
#' @param K Number of regimes  
#' @param model_type Type of model ("constant", "tvp", "exogenous", "gas")
#' @return Full parameter vector with K identical variance parameters
#' @details
#' Converts a parameter vector with a single shared variance parameter back to
#' one with K identical variance parameters. Used to restore the expected
#' parameter structure after optimization with equal_variances = TRUE.
#'
#' @examples
#' # Compressed TVP: [mu1, mu2, mu3, σ²_shared, trans..., A...]
#' compressed <- c(-1, 0, 1, 0.5, rep(0.2, 6), rep(0.1, 6))
#' expanded <- expand_variances(compressed, K = 3, "tvp")  
#' # Result: [mu1, mu2, mu3, σ², σ², σ², trans..., A...]
expand_variances <- function(params, K, model_type = c("constant", "tvp", "exogenous", "gas")) {
  model_type <- match.arg(model_type)
  
  if (length(params) == 0) {
    stop("Parameter vector cannot be empty")
  }
  
  # Extract components
  mu <- params[1:K]
  sigma2_shared <- params[K+1]
  
  # Expand shared variance to K identical values
  sigma2_expanded <- rep(sigma2_shared, K)
  
  # Get remaining parameters
  if (length(params) > K+1) {
    remaining <- params[(K+2):length(params)]
    expanded_params <- c(mu, sigma2_expanded, remaining)
  } else {
    expanded_params <- c(mu, sigma2_expanded)
  }
  
  # Add attributes to track the transformation
  attr(expanded_params, "expanded") <- TRUE
  attr(expanded_params, "target_K") <- K
  attr(expanded_params, "model_type") <- model_type
  
  return(expanded_params)
}

#' Check if parameter vector represents compressed variances
#'
#' @param params Parameter vector to check
#' @return TRUE if parameters were compressed, FALSE otherwise
#' @details
#' Helper function to check if a parameter vector has been compressed
#' using compress_variances(). Useful for conditional logic.
is_compressed_variances <- function(params) {
  return(!is.null(attr(params, "compressed")) && attr(params, "compressed"))
}

#' Calculate expected parameter length for compressed vs full parameterization
#'
#' @param K Number of regimes
#' @param model_type Type of model ("constant", "tvp", "exogenous", "gas")
#' @param equal_variances Whether variances are constrained to be equal
#' @return Expected length of parameter vector
#' @details
#' Helper function to calculate the expected parameter vector length based on
#' the model type and variance constraint setting.
expected_param_length <- function(K, model_type = c("constant", "tvp", "exogenous", "gas"), 
                                  equal_variances = FALSE) {
  model_type <- match.arg(model_type)
  
  # Base lengths for unconstrained case
  base_lengths <- switch(model_type,
                         "constant" = K + K + K*(K-1),                    # mu + sigma2 + trans
                         "tvp" = K + K + K*(K-1) + K*(K-1),               # mu + sigma2 + trans + A  
                         "exogenous" = K + K + K*(K-1) + K*(K-1),         # mu + sigma2 + trans + A
                         "gas" = K + K + K*(K-1) + K*(K-1) + K*(K-1)      # mu + sigma2 + trans + A + B
  )
  
  # Adjust for variance constraint (saves K-1 parameters)
  if (equal_variances) {
    return(base_lengths - (K - 1))
  } else {
    return(base_lengths)
  }
}
