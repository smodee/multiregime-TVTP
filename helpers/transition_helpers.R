#' Transition Matrix Helper Functions for Multi-Regime TVTP Models
#'
#' This file provides functions for creating and manipulating transition matrices
#' with support for both diagonal and off-diagonal parameterizations.
#'
#' DIAGONAL VS OFF-DIAGONAL PARAMETERIZATIONS:
#' ===========================================
#'
#' DIAGONAL PARAMETERIZATION (diag_probs = TRUE) [DEFAULT]:
#' - Uses K parameters representing diagonal elements p_11, p_22, ..., p_KK
#' - Off-diagonal elements equal within each row: p_ij = (1-p_ii)/(K-1) for i≠j
#' - Compatible with original HMMGAS/simulation.R implementation
#' - For K=2: [p11, p22] -> matrix [[p11, 1-p11], [1-p22, p22]]
#' - For K=3: [p11, p22, p33] -> matrix with equal off-diagonal probabilities per row
#' - IMPORTANT FOR K>2: This parameterization assumes that when leaving a state,
#'   the probability of transitioning to any other state is uniform. This is a
#'   modeling simplification that reduces parameters but may not fit all data.
#'
#' OFF-DIAGONAL PARAMETERIZATION (diag_probs = FALSE):
#' - Uses K*(K-1) parameters representing all off-diagonal elements p_ij (i!=j)
#' - Diagonal elements calculated as: p_ii = 1 - sum(p_ij for j!=i)
#' - Allows arbitrary transition structures between states
#' - For K=2: [p12, p21] -> matrix [[p11, p12], [p21, p22]] where p11=1-p12, p22=1-p21
#' - More flexible but requires more parameters to estimate
#'
#' This dual system allows validation against the original implementation while
#' maintaining flexibility for K>2 regime models with complex transition dynamics.

# Load required helper functions
source("helpers/utility_functions.R")

#' Create a transition matrix from probability parameters
#'
#' @param probs Vector of transition probability parameters
#' @param diag_probs If TRUE, probs contains diagonal elements; if FALSE, off-diagonal elements
#' @param check_validity Logical; whether to validate the resulting matrix (default: TRUE)
#' @return A properly formatted Markov transition matrix
#' @details
#' This is the main function for creating transition matrices. It automatically
#' determines the parameterization and creates the appropriate matrix.
#' 
#' For off-diagonal parameterization (diag_probs=FALSE):
#' - Expects K*(K-1) parameters in row-major order excluding diagonal
#' - Example K=3: [p12, p13, p21, p23, p31, p32]
#' 
#' For diagonal parameterization (diag_probs=TRUE):
#' - Expects K parameters representing diagonal persistence probabilities
#' - Example K=3: [p11, p22, p33]
#'
#' @examples
#' # Off-diagonal parameterization (default)
#' transition_matrix(c(0.3, 0.2), diag_probs = FALSE)  # 2x2 matrix
#' 
#' # Diagonal parameterization (original simulation.R style)
#' transition_matrix(c(0.8, 0.9), diag_probs = TRUE)   # 2x2 matrix with p11=0.8, p22=0.9
transition_matrix <- function(probs, diag_probs = TRUE, check_validity = TRUE) {
  
  if (!is.numeric(probs) || length(probs) == 0) {
    stop("Input probabilities must be a non-empty numeric vector")
  }
  
  if (diag_probs) {
    return(transition_matrix_diagonal(probs, check_validity))
  } else {
    return(transition_matrix_offdiagonal(probs, check_validity))
  }
}

#' Create transition matrix from diagonal parameterization
#'
#' @param diag_probs Vector of diagonal transition probabilities [p11, p22, ..., pKK]
#' @param check_validity Logical; whether to validate the resulting matrix
#' @return K x K transition matrix
#' @details
#' Creates a transition matrix where:
#' - Diagonal elements are specified directly
#' - Off-diagonal elements are equal within each row: p_ij = (1-p_ii)/(K-1) for i≠j
#' 
#' This matches the original simulation.R parameterization for K=2 and generalizes
#' it to K>2 by assuming equal off-diagonal probabilities within each row.
transition_matrix_diagonal <- function(diag_probs, check_validity = TRUE) {
  
  K <- length(diag_probs)
  
  if (K < 2) {
    stop("Number of regimes must be >= 2")
  }
  
  # Validate diagonal probabilities are in valid range
  if (any(diag_probs <= 0) || any(diag_probs >= 1)) {
    stop("Diagonal probabilities must be strictly between 0 and 1")
  }
  
  # Initialize transition matrix
  t_mat <- matrix(0, nrow = K, ncol = K)
  
  # Set diagonal elements and equal off-diagonal elements for each row
  for (i in 1:K) {
    t_mat[i, i] <- diag_probs[i]
    
    # Calculate equal off-diagonal probability for this row
    off_diag_prob <- (1 - diag_probs[i]) / (K - 1)
    
    # Set all off-diagonal elements for this row
    for (j in 1:K) {
      if (i != j) {
        t_mat[i, j] <- off_diag_prob
      }
    }
  }
  
  if (check_validity) {
    validate_transition_matrix(t_mat)
  }
  
  return(t_mat)
}

#' Create transition matrix from off-diagonal parameterization  
#'
#' @param off_diag_probs Vector of off-diagonal transition probabilities
#' @param check_validity Logical; whether to validate the resulting matrix
#' @return K x K transition matrix
#' @details
#' Creates a transition matrix from K*(K-1) off-diagonal elements specified
#' in row-major order. Diagonal elements are calculated as p_ii = 1 - sum(p_ij for j≠i).
#' 
#' This is the default parameterization for the new implementation.
transition_matrix_offdiagonal <- function(off_diag_probs, check_validity = TRUE) {
  
  n <- length(off_diag_probs)
  
  # Determine matrix size: K*(K-1) = n, so K^2 - K - n = 0
  K <- (1 + sqrt(1 + 4*n)) / 2
  
  if (abs(K - round(K)) > 1e-10) {
    stop(sprintf("Invalid length %d for off-diagonal probabilities. Expected K*(K-1) for some integer K >= 2.", n))
  }
  
  K <- round(K)
  
  # Initialize transition matrix
  t_mat <- matrix(0, nrow = K, ncol = K)
  
  # Fill off-diagonal elements row by row
  idx <- 1
  for (i in 1:K) {
    for (j in 1:K) {
      if (i != j) {
        t_mat[i, j] <- off_diag_probs[idx]
        idx <- idx + 1
      }
    }
    
    # Calculate diagonal element to make row sum to 1
    t_mat[i, i] <- 1 - sum(t_mat[i, -i])
  }
  
  if (check_validity) {
    validate_transition_matrix(t_mat)
  }
  
  return(t_mat)
}

#' Convert between diagonal and off-diagonal parameterizations
#'
#' @param probs Vector of transition probability parameters
#' @param from_diag If TRUE, converting FROM diagonal TO off-diagonal; if FALSE, the reverse
#' @return Vector of probabilities in the target parameterization
#' @details
#' Converts between the two parameterization schemes. Useful for validation
#' studies where you need to compare results using different parameterizations.
#' 
#' Note: Converting from off-diagonal to diagonal may lose information if the
#' original off-diagonal probabilities were not equal within rows.
#'
#' @examples
#' # Convert diagonal to off-diagonal
#' diag_params <- c(0.8, 0.9)  # p11=0.8, p22=0.9
#' off_diag_params <- convert_parameterization(diag_params, from_diag = TRUE)
#' # Result: [0.2, 0.1] representing p12=0.2, p21=0.1
#' 
#' # Convert back (should recover original structure, not exact values)
#' recovered_diag <- convert_parameterization(off_diag_params, from_diag = FALSE)
convert_parameterization <- function(probs, from_diag = TRUE) {
  
  if (!is.numeric(probs) || length(probs) == 0) {
    stop("Input probabilities must be a non-empty numeric vector")
  }
  
  if (from_diag) {
    # Convert FROM diagonal TO off-diagonal
    K <- length(probs)
    
    if (K < 2) {
      stop("Need at least 2 diagonal probabilities")
    }
    
    # Create transition matrix from diagonal parameterization
    t_mat <- transition_matrix_diagonal(probs, check_validity = FALSE)
    
    # Extract off-diagonal elements in row-major order
    off_diag_probs <- numeric(K * (K - 1))
    idx <- 1
    
    for (i in 1:K) {
      for (j in 1:K) {
        if (i != j) {
          off_diag_probs[idx] <- t_mat[i, j]
          idx <- idx + 1
        }
      }
    }
    
    return(off_diag_probs)
    
  } else {
    # Convert FROM off-diagonal TO diagonal
    n <- length(probs)
    K <- (1 + sqrt(1 + 4*n)) / 2
    
    if (abs(K - round(K)) > 1e-10) {
      stop("Invalid length for off-diagonal probabilities")
    }
    
    K <- round(K)
    
    # Create transition matrix from off-diagonal parameterization
    t_mat <- transition_matrix_offdiagonal(probs, check_validity = FALSE)
    
    # Extract diagonal elements
    diag_probs <- diag(t_mat)
    
    return(diag_probs)
  }
}

#' Validate that a matrix is a proper stochastic transition matrix
#'
#' @param t_mat Transition matrix to validate
#' @param tol Numerical tolerance for validation checks
#' @return TRUE if valid, throws error otherwise
validate_transition_matrix <- function(t_mat, tol = 1e-10) {
  
  if (!is.matrix(t_mat)) {
    stop("Input must be a matrix")
  }
  
  if (nrow(t_mat) != ncol(t_mat)) {
    stop("Transition matrix must be square")
  }
  
  K <- nrow(t_mat)
  
  # Check for negative values
  if (any(t_mat < -tol)) {
    min_val <- min(t_mat)
    neg_positions <- which(t_mat < -tol, arr.ind = TRUE)
    stop(sprintf("Transition matrix contains negative values (min: %.6f) at positions: %s", 
                 min_val, paste(apply(neg_positions, 1, function(x) paste0("(", x[1], ",", x[2], ")")), collapse = ", ")))
  }
  
  # Check that all values are <= 1 
  if (any(t_mat > 1 + tol)) {
    max_val <- max(t_mat)
    large_positions <- which(t_mat > 1 + tol, arr.ind = TRUE)
    stop(sprintf("Transition matrix contains values > 1 (max: %.6f) at positions: %s",
                 max_val, paste(apply(large_positions, 1, function(x) paste0("(", x[1], ",", x[2], ")")), collapse = ", ")))
  }
  
  # Check that rows sum to 1
  row_sums <- rowSums(t_mat)
  bad_rows <- which(abs(row_sums - 1) > tol)
  
  if (length(bad_rows) > 0) {
    stop(sprintf("Rows %s do not sum to 1. Row sums: %s", 
                 paste(bad_rows, collapse = ", "),
                 paste(sprintf("%.6f", row_sums[bad_rows]), collapse = ", ")))
  }
  
  return(TRUE)
}

#' Convert probability parameters to ensure valid transition matrix
#'
#' @param probs Vector of transition probability parameters
#' @param diag_probs If TRUE, probs contains diagonal elements; if FALSE, off-diagonal elements  
#' @param max_constraint Maximum allowed value to ensure positive diagonal/off-diagonal elements
#' @return Vector of adjusted probabilities that will create a valid transition matrix
#' @details
#' Ensures that probability parameters will result in a valid stochastic matrix.
#' For diagonal parameterization, constrains values to (0,1).
#' For off-diagonal parameterization, ensures row sums don't exceed max_constraint.
convert_to_valid_probs <- function(probs, diag_probs = TRUE, max_constraint = 0.99) {
  
  if (!is.numeric(probs) || length(probs) == 0) {
    stop("Input probabilities must be a non-empty numeric vector")
  }
  
  # Handle NAs and infinite values
  if (any(is.na(probs) | is.infinite(probs))) {
    warning("Invalid values (NA/Inf) detected in probabilities, replacing with default values")
    probs[is.na(probs) | is.infinite(probs)] <- ifelse(diag_probs, 0.5, 0.1)
  }
  
  if (diag_probs) {
    # For diagonal parameterization: ensure all values in (0, 1)
    probs <- pmax(0.01, pmin(0.99, probs))
    return(probs)
    
  } else {
    # For off-diagonal parameterization: ensure row sums don't exceed max_constraint
    n <- length(probs)
    K <- (1 + sqrt(1 + 4*n)) / 2
    
    if (abs(K - round(K)) > 1e-10) {
      stop("Invalid length for off-diagonal probabilities")
    }
    
    K <- round(K)
    valid_probs <- numeric(length(probs))
    idx <- 1
    
    for (i in 1:K) {
      # Extract off-diagonal probabilities for this row
      row_indices <- idx:(idx + K - 2)
      row_probs <- probs[row_indices]
      
      # Ensure individual probabilities are in valid range
      row_probs <- pmax(0.001, pmin(0.99, row_probs))
      
      # If row sum exceeds max_constraint, scale proportionally
      row_sum <- sum(row_probs)
      if (row_sum > max_constraint) {
        scale_factor <- max_constraint / row_sum
        row_probs <- row_probs * scale_factor
      }
      
      valid_probs[row_indices] <- row_probs
      idx <- idx + K - 1
    }
    
    return(valid_probs)
  }
}

#' Calculate the stationary distribution of a transition matrix
#'
#' @param probs Vector of probability parameters OR transition matrix
#' @param diag_probs If probs is a vector, whether it uses diagonal parameterization
#' @param method Method for calculation ("eigen" or "power")
#' @param tol Tolerance for eigenvalue detection (default: 1e-8)
#' @param fallback_value Value to return if no unique stationary distribution exists
#' @return The stationary distribution as a vector
#' @details
#' Calculates the stationary distribution π such that π = π⋅P, where P is the transition matrix.
#' 
#' Two methods available:
#' - "eigen": Uses eigendecomposition (faster, default)
#' - "power": Uses power iteration (more stable for ill-conditioned matrices)
#'
#' @examples
#' # Using diagonal parameters
#' stat_dist(c(0.8, 0.9), diag_probs = TRUE)
#' 
#' # Using off-diagonal parameters  
#' stat_dist(c(0.2, 0.1), diag_probs = FALSE)
#' 
#' # Using transition matrix directly
#' P <- transition_matrix(c(0.8, 0.9), diag_probs = TRUE)
#' stat_dist(P)
stat_dist <- function(probs, diag_probs = TRUE, method = c("eigen", "power"), 
                      tol = 1e-8, fallback_value = NULL) {
  
  method <- match.arg(method)
  
  # Convert to transition matrix if needed
  if (is.matrix(probs)) {
    t_mat <- probs
    validate_transition_matrix(t_mat)
  } else {
    t_mat <- transition_matrix(probs, diag_probs = diag_probs, check_validity = TRUE)
  }
  
  K <- nrow(t_mat)
  
  if (method == "eigen") {
    # Eigenvalue method
    eigen_decomp <- eigen(t(t_mat))
    
    # Find eigenvalue closest to 1
    eigenvalue_1_idx <- which.min(abs(eigen_decomp$values - 1))
    
    if (abs(eigen_decomp$values[eigenvalue_1_idx] - 1) > tol) {
      if (!is.null(fallback_value)) {
        warning("No eigenvalue close to 1 found. Using fallback value.")
        return(fallback_value)
      } else {
        stop("The transition matrix does not have a stationary distribution (no eigenvalue ≈ 1).")
      }
    }
    
    # Extract and normalize eigenvector
    stationary_vector <- Re(eigen_decomp$vectors[, eigenvalue_1_idx])
    
    if (any(stationary_vector < 0)) {
      stationary_vector <- -stationary_vector  # Flip sign if needed
    }
    
    stationary_distribution <- stationary_vector / sum(stationary_vector)
    
  } else {
    # Power iteration method
    stationary_distribution <- stat_dist_power(t_mat, tol = tol)
  }
  
  return(stationary_distribution)
}

#' Calculate stationary distribution using power iteration
#'
#' @param t_mat Transition matrix
#' @param max_iter Maximum iterations (default: 1000) 
#' @param tol Convergence tolerance (default: 1e-10)
#' @param initial_dist Initial distribution (default: uniform)
#' @return Stationary distribution vector
stat_dist_power <- function(t_mat, max_iter = 1000, tol = 1e-10, initial_dist = NULL) {
  
  K <- nrow(t_mat)
  
  # Initialize distribution
  if (is.null(initial_dist)) {
    dist <- rep(1/K, K)
  } else {
    if (length(initial_dist) != K) {
      stop("Initial distribution length must match number of states")
    }
    dist <- initial_dist / sum(initial_dist)
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

#' Determine number of regimes from parameter vector length
#'
#' @param probs Vector of probability parameters
#' @param diag_probs Whether parameters use diagonal parameterization
#' @return Number of regimes K
infer_K_from_probs <- function(probs, diag_probs = TRUE) {
  
  n <- length(probs)
  
  if (diag_probs) {
    # For diagonal parameterization: K parameters
    return(n)
  } else {
    # For off-diagonal parameterization: K*(K-1) parameters
    K <- (1 + sqrt(1 + 4*n)) / 2
    
    if (abs(K - round(K)) > 1e-10) {
      stop(sprintf("Invalid parameter vector length %d for off-diagonal parameterization", n))
    }
    
    return(round(K))
  }
}

#' Create parameter vector with proper attributes for transition probabilities
#'
#' @param probs Numeric vector of probability values
#' @param diag_probs Whether these are diagonal or off-diagonal probabilities  
#' @return Parameter vector with attributes set
create_transition_params <- function(probs, diag_probs = TRUE) {
  
  if (!is.numeric(probs) || length(probs) == 0) {
    stop("Probabilities must be a non-empty numeric vector")
  }
  
  # Validate and clean probabilities
  probs <- convert_to_valid_probs(probs, diag_probs = diag_probs)
  
  # Set attributes for identification
  attr(probs, "diag_probs") <- diag_probs
  attr(probs, "K") <- infer_K_from_probs(probs, diag_probs)
  attr(probs, "param_type") <- "transition_probabilities"
  
  return(probs)
}

#' Print transition matrix in a readable format
#'
#' @param t_mat Transition matrix
#' @param digits Number of decimal places to display (default: 4)
print_transition_matrix <- function(t_mat, digits = 4) {
  
  K <- nrow(t_mat)
  
  cat("Transition Matrix (", K, "x", K, "):\n", sep = "")
  cat("========================\n")
  
  # Create column headers
  col_headers <- paste0("S", 1:K)
  row_headers <- paste0("S", 1:K)
  
  # Print with row and column labels
  cat(sprintf("%4s", ""), paste(sprintf("%8s", col_headers), collapse = ""), "\n")
  
  for (i in 1:K) {
    cat(sprintf("%4s", row_headers[i]))
    for (j in 1:K) {
      cat(sprintf("%8s", format(round(t_mat[i,j], digits), nsmall = digits)))
    }
    cat("\n")
  }
  
  # Print row sums as validation
  cat("\nRow sums: ")
  row_sums <- rowSums(t_mat)
  cat(paste(sprintf("%.6f", row_sums), collapse = ", "), "\n")
  
  # Calculate and print stationary distribution
  tryCatch({
    stat_dist_vec <- stat_dist(t_mat)
    cat("Stationary distribution: ")
    cat(paste(sprintf("%.4f", stat_dist_vec), collapse = ", "), "\n")
  }, error = function(e) {
    cat("Stationary distribution: Could not calculate\n")
  })
}
