# Logit/logistic functions are needed from the utility_functions script
source("helpers/utility_functions.R")

#' Parameter transformation functions for regime switching models
#' 
#' These functions handle parameter transformations for optimization and estimation
#' in regime switching models with time-varying transition probabilities.
#' 
#' The naming convention is:
#' - transform_*: Transforms parameters from their natural space to the unconstrained space
#' - untransform_*: Transforms parameters from the unconstrained space back to their natural space

#' Extract mean parameters from a parameter vector
#'
#' @param par Parameter vector
#' @return Vector of mean parameters
#' @examples
#' par <- c(1, 2, 3, 0.1, 0.2, 0.3, 0.2, 0.3, 0.1, 0.2)
#' mean_from_par(par)  # Returns c(1, 2, 3)
mean_from_par <- function(par) {
  K <- count_regime(par)
  return(par[1:K])
}

#' Extract variance parameters from a parameter vector
#'
#' @param par Parameter vector
#' @return Vector of variance parameters
#' @examples
#' par <- c(1, 2, 3, 0.1, 0.2, 0.3, 0.2, 0.3, 0.1, 0.2)
#' sigma2_from_par(par)  # Returns c(0.1, 0.2, 0.3)
sigma2_from_par <- function(par) {
  K <- count_regime(par)
  return(par[(K+1):(2*K)])
}

#' Extract transition probability parameters from a parameter vector
#'
#' @param par Parameter vector
#' @return Vector of transition probability parameters
#' @examples
#' par <- c(1, 2, 3, 0.1, 0.2, 0.3, 0.2, 0.3, 0.1, 0.2)
#' transp_from_par(par)  # Returns c(0.2, 0.3, 0.1, 0.2) for a 3-state model
transp_from_par <- function(par) {
  K <- count_regime(par)
  return(par[(2*K+1):(K*(K+1))])
}

#' Extract coefficients for exogenous or autoregressive effects from a parameter vector
#'
#' @param par Parameter vector
#' @return Vector of A coefficients
#' @examples
#' par <- c(1, 2, 3, 0.1, 0.2, 0.3, 0.2, 0.3, 0.1, 0.2, 0.05, 0.1, 0.15, 0.2)
#' A_from_par(par)  # Returns c(0.05, 0.1, 0.15, 0.2) for a 3-state model
A_from_par <- function(par) {
  K <- count_regime(par)
  return(par[((K*(K+1))+1):(2*K^2)])
}

#' Determine the number of regimes from a parameter vector length
#'
#' @param par Parameter vector
#' @return Number of regimes (K)
#' @examples
#' par <- c(1, 2, 3, 0.1, 0.2, 0.3, 0.2, 0.3, 0.1, 0.2, 0.05, 0.1, 0.15, 0.2)
#' count_regime(par)  # Returns 3
count_regime <- function(par) {
  K <- sqrt(length(par)/2)
  
  if (K != as.integer(K)) {
    stop("The parameter vector length is not compatible with any valid number of regimes.")
  }
  
  return(as.integer(K))
}

#' Transform parameters from natural space to unconstrained space for TVP model
#'
#' @param par Parameters in their natural space
#' @return Parameters in unconstrained space
#' @details
#' Applies log transformation to variances and logit to transition probabilities
#' to make parameters suitable for unconstrained optimization.
transform_TVP <- function(par) {
  mu <- mean_from_par(par)
  sigma2 <- sigma2_from_par(par)
  init_trans <- transp_from_par(par)
  A <- A_from_par(par)
  
  # Check for invalid values
  if (any(sigma2 <= 0)) {
    stop("Variance parameters must be positive")
  }
  if (any(init_trans <= 0) || any(init_trans >= 1)) {
    stop("Transition probabilities must be between 0 and 1")
  }
  
  return(c(mu, log(sigma2), logit(init_trans), A))
}

#' Transform parameters from unconstrained space back to natural space for TVP model
#'
#' @param par_t Parameters in unconstrained space
#' @return Parameters in their natural space
#' @details
#' Applies exp transformation to log-variances and logistic to logit-probabilities
#' to convert back to naturally bounded parameters.
untransform_TVP <- function(par_t) {
  mu <- mean_from_par(par_t)
  log_sigma2 <- sigma2_from_par(par_t)
  logit_init_trans <- transp_from_par(par_t)
  A <- A_from_par(par_t)
  
  return(c(mu, exp(log_sigma2), logistic(logit_init_trans), A))
}

#' Transform parameters from natural space to unconstrained space for exogenous TVP model
#'
#' @param par Parameters in their natural space
#' @return Parameters in unconstrained space
#' @details
#' Applies log transformation to variances and logit to transition probabilities
#' to make parameters suitable for unconstrained optimization.
transform_TVPXExo <- function(par) {
  mu <- mean_from_par(par)
  sigma2 <- sigma2_from_par(par)
  init_trans <- transp_from_par(par)
  A <- A_from_par(par)
  
  # Check for invalid values
  if (any(sigma2 <= 0)) {
    stop("Variance parameters must be positive")
  }
  if (any(init_trans <= 0) || any(init_trans >= 1)) {
    stop("Transition probabilities must be between 0 and 1")
  }
  
  return(c(mu, log(sigma2), logit(init_trans), A))
}

#' Transform parameters from unconstrained space back to natural space for exogenous TVP model
#'
#' @param par_t Parameters in unconstrained space
#' @return Parameters in their natural space
#' @details
#' Applies exp transformation to log-variances and logistic to logit-probabilities
#' to convert back to naturally bounded parameters.
untransform_TVPXExo <- function(par_t) {
  mu <- mean_from_par(par_t)
  log_sigma2 <- sigma2_from_par(par_t)
  logit_init_trans <- transp_from_par(par_t)
  A <- A_from_par(par_t)
  
  return(c(mu, exp(log_sigma2), logistic(logit_init_trans), A))
}

# For backward compatibility - these will be deprecated in future versions
Rpar.trasf_TVPXExo <- function(par_t) {
  warning("This function is deprecated. Please use untransform_TVPXExo() instead.")
  return(untransform_TVPXExo(par_t))
}

Rpar.trasf.inv_TVPXExo <- function(par) {
  warning("This function is deprecated. Please use transform_TVPXExo() instead.")
  return(transform_TVPXExo(par))
}

Rpar.trasf_TVP <- function(par_t) {
  warning("This function is deprecated. Please use untransform_TVP() instead.")
  return(untransform_TVP(par_t))
}

Rpar.trasf.inv_TVP <- function(par) {
  warning("This function is deprecated. Please use transform_TVP() instead.")
  return(transform_TVP(par))
}

#' Check if parameter vector has valid structure
#'
#' @param par Parameter vector to check
#' @param K Number of regimes (optional, will be inferred if not provided)
#' @return TRUE if valid, otherwise throws an error with explanation
#' @examples
#' # Valid parameter vector for 3 regimes
#' par <- c(1, 2, 3, 0.1, 0.2, 0.3, 0.2, 0.3, 0.1, 0.2, 0.05, 0.1, 0.15, 0.2)
#' validate_parameter_vector(par)  # Returns TRUE
validate_parameter_vector <- function(par, K = NULL) {
  if (is.null(K)) {
    tryCatch({
      K <- count_regime(par)
    }, error = function(e) {
      stop("Invalid parameter vector length. Cannot determine number of regimes.")
    })
  }
  
  expected_length <- 2*K^2
  
  if (length(par) != expected_length) {
    stop(sprintf("Invalid parameter vector length. Expected %d parameters for %d regimes, got %d.", 
                 expected_length, K, length(par)))
  }
  
  # Check variance parameters
  sigma2 <- sigma2_from_par(par)
  if (any(sigma2 <= 0)) {
    stop("Variance parameters must be positive")
  }
  
  # Check transition probabilities
  trans_prob <- transp_from_par(par)
  if (any(trans_prob < 0) || any(trans_prob > 1)) {
    stop("Transition probabilities must be between 0 and 1")
  }
  
  # All checks passed
  return(TRUE)
}

#' Create initial parameter guesses for a regime switching model
#'
#' @param K Number of regimes
#' @param mu_range Range for mean parameters (default: c(-3, 3))
#' @param sigma2_range Range for variance parameters (default: c(0.1, 1))
#' @param trans_prob_range Range for transition probabilities (default: c(0.1, 0.4))
#' @param A_range Range for A coefficients (default: c(-0.2, 0.2))
#' @return A parameter vector with reasonable initial values
#' @examples
#' # Create initial parameters for a 3-regime model
#' create_initial_parameters(3)
create_initial_parameters <- function(K, 
                                      mu_range = c(-3, 3),
                                      sigma2_range = c(0.1, 1),
                                      trans_prob_range = c(0.1, 0.4),
                                      A_range = c(-0.2, 0.2)) {
  # Generate means that are well-separated
  mu_step <- diff(mu_range) / (K - 1)
  mu <- seq(mu_range[1], mu_range[2], length.out = K)
  
  # Generate increasing variances
  sigma2_step <- diff(sigma2_range) / (K - 1)
  sigma2 <- seq(sigma2_range[1], sigma2_range[2], length.out = K)
  
  # Generate transition probabilities
  trans_prob <- rep(mean(trans_prob_range), K*(K-1))
  
  # Generate A coefficients
  A <- rep(0, K*(K-1))  # Start with no effect
  
  # Combine all parameters
  par <- c(mu, sigma2, trans_prob, A)
  
  return(par)
}