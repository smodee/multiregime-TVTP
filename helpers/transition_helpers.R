#' Create a transition matrix from a vector of off-diagonal probabilities
#'
#' @param probs A vector of off-diagonal transition probabilities
#' @param check_validity Logical; whether to check that the resulting matrix is valid (default: TRUE)
#' @return A properly formatted Markov transition matrix
#'
#' @details 
#' This function accepts off-diagonal transition probabilities specified row by row.
#' For a K-state model, the length of probs should be K*(K-1).
#' 
#' For example, in a 2-state model, probs = c(p12, p21).
#' For a 3-state model, probs = c(p12, p13, p21, p23, p31, p32).
#'
#' The diagonal elements are calculated to ensure each row sums to 1.
#'
#' @examples
#' # 2x2 matrix
#' transition_matrix(c(0.3, 0.2))  # p12=0.3, p21=0.2
#'
#' # 3x3 matrix
#' transition_matrix(c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3))
transition_matrix <- function(probs, check_validity = TRUE) {
  # Check that input is numeric
  if (!is.numeric(probs)) {
    stop("Input probabilities must be numeric")
  }
  
  # Length of the input probabilities vector
  n <- length(probs)
  
  # Determine the size of the transition matrix
  m <- floor(0.5 + sqrt(0.25 + n))
  if (m * (m - 1) != n) {
    stop("Invalid length of probabilities vector. For K states, expected K*(K-1) probabilities.")
  }
  
  # Initialize an empty transition matrix
  t_mat <- matrix(0, nrow = m, ncol = m)
  
  # Fill the transition matrix row by row with off-diagonal elements
  idx <- 1
  for (i in 1:m) {
    for (j in 1:(m-1)) {
      if (j < i) {
        t_mat[i, j] <- probs[idx]
        idx <- idx + 1
      } else {
        t_mat[i, j+1] <- probs[idx]
        idx <- idx + 1
      }
    }
    
    # Calculate the diagonal element
    stay_prob <- 1 - sum(t_mat[i, ])
    t_mat[i, i] <- stay_prob
  }
  
  # Only check validity if requested and not in a performance-critical loop
  if (check_validity) {
    # Check for negative values
    if (any(t_mat < -1e-10)) {  # Allow for minor numerical imprecision
      neg_idx <- which(t_mat < 0, arr.ind = TRUE)
      neg_values <- t_mat[neg_idx]
      err_msg <- sprintf("Matrix contains negative values (min: %.6f). Off-diagonal elements may sum to >1.", 
                         min(t_mat))
      stop(err_msg)
    }
    
    # Check if each row sums to 1 (allowing for some numerical tolerance)
    row_sums <- rowSums(t_mat)
    if (!all(abs(row_sums - 1) < 1e-10)) {
      err_msg <- sprintf("Not all rows sum to 1. Row sums: %s", 
                         paste(format(row_sums, digits = 6), collapse = ", "))
      stop(err_msg)
    }
  }
  
  return(t_mat)
}

#' Calculate the stationary distribution of a Markov transition matrix
#'
#' @param probs Vector of off-diagonal probabilities or transition matrix
#' @param tol Tolerance for detecting eigenvalue 1 (default: 1e-8)
#' @param fallback_value What to return if no unique stationary distribution exists (default: NULL)
#' @return The stationary distribution as a vector
#'
#' @details
#' Calculates the stationary distribution π that satisfies π = π⋅P, where P is the transition matrix.
#' This is equivalent to finding the eigenvector corresponding to eigenvalue 1 of P^T.
#'
#' If fallback_value is provided and no unique stationary distribution exists, returns the fallback instead of stopping.
#'
#' @examples
#' # Using off-diagonal transition probabilities
#' stat_dist(c(0.3, 0.2))  # p12=0.3, p21=0.2
#'
#' # Using a transition matrix directly
#' P <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE)
#' stat_dist(P)
stat_dist <- function(probs, tol = 1e-8, fallback_value = NULL) {
  # Check if input is already a matrix
  if (is.matrix(probs)) {
    t_mat <- probs
    
    # Verify it's a valid stochastic matrix
    if (nrow(t_mat) != ncol(t_mat)) {
      stop("Input matrix must be square")
    }
    
    # Check row sums
    row_sums <- rowSums(t_mat)
    if (!all(abs(row_sums - 1) < 1e-10)) {
      warning("Input matrix rows don't sum to 1. This may not be a valid transition matrix.")
    }
  } else {
    # Convert vector to transition matrix
    t_mat <- transition_matrix(probs)
  }
  
  # Calculate the eigenvalues and eigenvectors of the transpose
  eigen_decomp <- eigen(t(t_mat))
  
  # Find the eigenvector corresponding to eigenvalue 1
  eigenvalue_1_idx <- which(abs(eigen_decomp$values - 1) < tol)
  
  if (length(eigenvalue_1_idx) != 1) {
    if (!is.null(fallback_value)) {
      warning("The transition matrix does not have a unique stationary distribution. Using fallback value.")
      return(fallback_value)
    } else {
      stop("The transition matrix does not have a unique stationary distribution.")
    }
  }
  
  # Extract the eigenvector and take the real part
  # (it should be real, but might have small imaginary parts due to numerical issues)
  stationary_vector <- Re(eigen_decomp$vectors[, eigenvalue_1_idx])
  
  # Normalize the vector to sum to 1
  stationary_distribution <- stationary_vector / sum(stationary_vector)
  
  return(stationary_distribution)
}

#' Calculate the stationary distribution using power iteration
#'
#' @param probs Vector of off-diagonal probabilities or transition matrix
#' @param max_iter Maximum number of iterations (default: 1000)
#' @param tol Convergence tolerance (default: 1e-10)
#' @param initial_dist Initial distribution (default: uniform)
#' @return The stationary distribution as a vector
#'
#' @details
#' An alternative method to calculate the stationary distribution through iterative application
#' of the transition matrix. This is sometimes more numerically stable than the eigenvalue method.
#'
#' @examples
#' # Using off-diagonal transition probabilities
#' stat_dist_power(c(0.3, 0.2))  # p12=0.3, p21=0.2
#'
#' # Using a transition matrix directly
#' P <- matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE)
#' stat_dist_power(P)
stat_dist_power <- function(probs, max_iter = 1000, tol = 1e-10, initial_dist = NULL) {
  # Check if input is already a matrix
  if (is.matrix(probs)) {
    t_mat <- probs
  } else {
    # Convert vector to transition matrix
    t_mat <- transition_matrix(probs)
  }
  
  n <- nrow(t_mat)
  
  # Initialize distribution
  if (is.null(initial_dist)) {
    dist <- rep(1/n, n)  # Uniform distribution
  } else {
    if (length(initial_dist) != n) {
      stop("Initial distribution length must match the number of states")
    }
    dist <- initial_dist / sum(initial_dist)  # Normalize
  }
  
  # Power iteration
  for (i in 1:max_iter) {
    new_dist <- dist %*% t_mat
    
    # Check convergence
    if (max(abs(new_dist - dist)) < tol) {
      return(as.vector(new_dist))
    }
    
    dist <- new_dist
  }
  
  warning("Power iteration did not converge within maximum iterations")
  return(as.vector(dist))
}

#' Convert a vector of transition probabilities to ensure they form a valid stochastic matrix
#'
#' @param probs Vector of off-diagonal probabilities to be validated
#' @param K Number of states
#' @param max_diag_sum Maximum sum of off-diagonal elements (default: 0.99)
#' @return A modified vector of probabilities that will form a valid stochastic matrix
#'
#' @details
#' This function ensures that the probabilities will result in a valid transition matrix
#' with non-negative elements and rows that sum to 1.
#'
#' @examples
#' # Fix probabilities that would lead to negative diagonal elements
#' convert_to_valid_probs(c(0.6, 0.6, 0.3, 0.3, 0.4, 0.4), K=3)
convert_to_valid_probs <- function(probs, K, max_diag_sum = 0.99) {
  # Check for NAs or NaNs in input
  if (any(is.na(probs))) {
    warning("NA values detected in probabilities, replacing with 0.5")
    probs[is.na(probs)] <- 0.5
  }
  
  # Check that we have the right number of probabilities
  expected_length <- K * (K - 1)
  if (length(probs) != expected_length) {
    stop(sprintf("Invalid length of probabilities vector. Expected %d, got %d.", 
                 expected_length, length(probs)))
  }
  
  # Initialize output
  valid_probs <- numeric(length(probs))
  idx <- 1
  
  for (i in 1:K) {
    # Calculate indices for this row's off-diagonal elements
    row_indices <- idx:(idx + K - 2)
    
    # Safety check for indexing
    if (length(row_indices) < 1 || max(row_indices) > length(probs)) {
      warning(paste("Invalid indices for row", i, "- Using defaults"))
      valid_probs[idx:(idx + K - 2)] <- rep(0.1, K - 1)
      idx <- idx + K - 1
      next
    }
    
    # Extract probabilities for this row
    row_probs <- probs[row_indices]
    
    # Check for invalid values
    if (any(is.na(row_probs)) || any(is.infinite(row_probs))) {
      warning(paste("Invalid values in row", i, "- Using defaults"))
      row_probs <- rep(0.1, length(row_indices))
    }
    
    # If sum > max_diag_sum, scale them down proportionally
    row_sum <- sum(row_probs)
    
    # Check if row_sum is valid
    if (is.na(row_sum) || is.infinite(row_sum)) {
      warning(paste("Invalid row sum for row", i, "- Using defaults"))
      row_probs <- rep(0.1, length(row_indices))
      row_sum <- sum(row_probs)
    }
    
    if (row_sum > max_diag_sum) {
      scale_factor <- max_diag_sum / row_sum
      row_probs <- row_probs * scale_factor
    }
    
    # Store the valid probabilities
    valid_probs[row_indices] <- row_probs
    idx <- idx + K - 1
  }
  
  return(valid_probs)
}