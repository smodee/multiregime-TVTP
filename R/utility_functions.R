#' General utility functions for regime switching models
#' These functions provide general mathematical and data processing utilities
#' that are needed across different parts of the regime switching models.

#' Logit transformation (natural space -> unconstrained space)
#'
#' @param x Value in \[0,1\] to transform
#' @param eps Small value to prevent boundary issues (default: .Machine$double.eps)
#' @return Transformed value in (-∞,∞)
#' @examples
#' \dontrun{
#' logit(0.75)  # Returns 1.098612
#' }
logit <- function(x, eps = .Machine$double.eps) {
  # Add bounds checking to prevent errors
  x <- pmin(pmax(x, eps), 1 - eps)
  return(log(x/(1-x)))
}

#' Logistic transformation (unconstrained space -> natural space)
#'
#' @param x Value in (-∞,∞) to transform
#' @return Transformed value in (0,1)
#' @details
#' Standard logistic function used by TVP and Exogenous models.
#'
#' NOTE ON MODEL-SPECIFIC USAGE:
#' - TVP model: Uses standard logistic() - matches original HMMGAS C implementation
#' - Exogenous model: Uses standard logistic() - matches original HMMGAS C implementation
#' - GAS model: Uses logistic_clamped() - matches original HMMGAS C implementation
#'
#' The original HMMGAS package applies probability clamping \[1e-10, 1-1e-10\] only
#' in the GAS model, not in TVP or Exogenous models. This asymmetry exists in the
#' original C code and is preserved here for consistency.
#' @examples
#' \dontrun{
#' logistic(1.098612)  # Returns 0.75
#' }
logistic <- function(x) {
  return(1/(1+exp(-x)))
}

#' Clamped logistic transformation (unconstrained space -> bounded natural space)
#'
#' @param x Value in (-∞,∞) to transform
#' @param eps Small value for boundary clamping (default: 1e-10)
#' @return Transformed value in \[eps, 1-eps\]
#' @details
#' This function matches the original HMMGAS C implementation which uses:
#' p = eps + (1-2*eps)/(1+exp(-x))
#'
#' This maps the unconstrained value x to the interval \[eps, 1-eps\] rather than (0,1),
#' providing numerical stability when probabilities are used in subsequent calculations
#' involving log(), division by p*(1-p), etc.
#'
#' The default eps=1e-10 matches the original HMMGAS implementation.
#'
#' NOTE ON MODEL-SPECIFIC USAGE:
#' This clamped version is used ONLY by the GAS model. The original HMMGAS C code
#' applies this clamping specifically in GAS (Filtering_2RegimesGAS.c) but uses
#' standard logistic in TVP (Filtering_2RegimesTVP.c) and Exogenous models.
#' This is because the GAS score calculation involves terms like p*(1-p) in the
#' denominator, making numerical stability more critical.
#' @examples
#' \dontrun{
#' logistic_clamped(0)     # Returns 0.5
#' logistic_clamped(100)   # Returns 1-1e-10 (not exactly 1)
#' logistic_clamped(-100)  # Returns 1e-10 (not exactly 0)
#' }
logistic_clamped <- function(x, eps = 1e-10) {
  return(eps + (1 - 2*eps) / (1 + exp(-x)))
}

#' Softmax row transformation with reference category
#'
#' Converts K-1 unconstrained real parameters into a full row of K probabilities
#' using softmax normalization with the diagonal entry fixed at reference value 0.
#'
#' @param x Numeric vector of length K-1 (unconstrained softmax parameters for
#'   off-diagonal entries of one row)
#' @return A list with components:
#'   \describe{
#'     \item{offdiag}{Numeric vector of K-1 off-diagonal probabilities}
#'     \item{diag}{Scalar diagonal probability}
#'   }
#' @details
#' For K-1 unconstrained parameters x_1, ..., x_{K-1} and a fixed reference of 0
#' for the diagonal entry:
#' \itemize{
#'   \item p_j = exp(x_j) / (1 + sum(exp(x_k))) for off-diagonal entries
#'   \item p_diag = 1 / (1 + sum(exp(x_k))) for the diagonal entry
#' }
#'
#' All probabilities are guaranteed positive and sum to 1. Uses the max-subtraction
#' trick for numerical stability with large parameter values.
#'
#' @examples
#' \dontrun{
#' # Equal probabilities for K=3: all params = 0
#' softmax_row(c(0, 0))  # offdiag = c(1/3, 1/3), diag = 1/3
#'
#' # High persistence: negative off-diagonal params
#' softmax_row(c(-3, -3))  # Small off-diag, large diagonal
#' }
softmax_row <- function(x) {
  max_x <- max(c(x, 0))
  exp_x <- exp(x - max_x)
  exp_ref <- exp(-max_x)
  denom <- exp_ref + sum(exp_x)
  list(offdiag = exp_x / denom, diag = exp_ref / denom)
}

#' Safe logarithm to handle zeros and negative values
#'
#' @param x Value to take the logarithm of
#' @param min_value Minimum value to replace zeros/negative values (default: .Machine$double.eps)
#' @return Log of max(x, min_value)
#' @examples
#' \dontrun{
#' safe_log(0)      # Returns log of a small positive number instead of -Inf
#' safe_log(-5)     # Returns log of a small positive number instead of NaN
#' safe_log(c(2, 0, 3, -1))  # Handles vectors safely
#' }
safe_log <- function(x, min_value = .Machine$double.eps) {
  return(log(pmax(x, min_value)))
}

#' Safe exponentiation to prevent overflow
#'
#' @param x Value to exponentiate
#' @param max_value Maximum allowed value in exponent (default: 700)
#' @return exp(min(x, max_value))
#' @examples
#' \dontrun{
#' safe_exp(1000)  # Returns a large but finite number instead of Inf
#' }
safe_exp <- function(x, max_value = 700) {
  return(exp(pmin(x, max_value)))
}

#' Check if a vector contains any NA, NaN, Inf, or -Inf values
#'
#' @param x Vector to check
#' @return TRUE if any invalid values are found, FALSE otherwise
#' @examples
#' \dontrun{
#' has_invalid_values(c(1, 2, NA, 4))        # Returns TRUE
#' has_invalid_values(c(1, 2, NaN, 4))       # Returns TRUE
#' has_invalid_values(c(1, 2, Inf, 4))       # Returns TRUE
#' has_invalid_values(c(1, 2, -Inf, 4))      # Returns TRUE
#' has_invalid_values(c(1, 2, 3, 4))         # Returns FALSE
#' }
has_invalid_values <- function(x) {
  return(any(is.na(x) | is.nan(x) | is.infinite(x)))
}

#' Replace invalid values (NA, NaN, Inf, -Inf) with a default value
#'
#' @param x Vector to process
#' @param replacement Value to use for replacements (default: 0)
#' @return Vector with invalid values replaced
#' @examples
#' \dontrun{
#' replace_invalid_values(c(1, 2, NA, Inf, -Inf, 6))  # Returns c(1, 2, 0, 0, 0, 6)
#' }
replace_invalid_values <- function(x, replacement = 0) {
  x[is.na(x) | is.nan(x) | is.infinite(x)] <- replacement
  return(x)
}

#' Calculate moving average of a time series
#'
#' @param x Numeric vector
#' @param window_size Number of periods to include in the moving average
#' @param align Alignment ("left", "center", "right") (default: "center")
#' @return Vector of moving averages
#' @examples
#' \dontrun{
#' moving_average(1:10, 3)  # 3-period centered moving average
#' }
moving_average <- function(x, window_size, align = "center") {
  if (window_size < 1 || window_size > length(x)) {
    stop("Window size must be between 1 and the length of the vector")
  }
  
  n <- length(x)
  result <- numeric(n)
  
  if (align == "center") {
    # Center alignment (symmetric window)
    half_window <- floor(window_size/2)
    for (i in 1:n) {
      start <- max(1, i - half_window)
      end <- min(n, i + half_window)
      result[i] <- mean(x[start:end], na.rm = TRUE)
    }
  } else if (align == "left") {
    # Left alignment (use past values)
    for (i in 1:n) {
      start <- max(1, i - window_size + 1)
      end <- i
      result[i] <- mean(x[start:end], na.rm = TRUE)
    }
  } else if (align == "right") {
    # Right alignment (use future values)
    for (i in 1:n) {
      start <- i
      end <- min(n, i + window_size - 1)
      result[i] <- mean(x[start:end], na.rm = TRUE)
    }
  } else {
    stop("Alignment must be 'left', 'center', or 'right'")
  }
  
  return(result)
}

#' Standardize a vector (zero mean, unit variance)
#'
#' @param x Numeric vector to standardize
#' @param na.rm Whether to remove NA values (default: TRUE)
#' @return Standardized vector
#' @examples
#' \dontrun{
#' standardize(c(1, 2, 3, 4, 5))  # Returns values with mean 0 and SD 1
#' }
standardize <- function(x, na.rm = TRUE) {
  return((x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm))
}

#' Calculate median absolute deviation (robust measure of dispersion)
#'
#' @param x Numeric vector
#' @param na.rm Whether to remove NA values (default: TRUE)
#' @param constant Scaling factor for normal distribution (default: 1.4826)
#' @return Median absolute deviation
#' @examples
#' \dontrun{
#' mad_value(c(1, 2, 3, 100))  # Less affected by outliers than standard deviation
#' }
mad_value <- function(x, na.rm = TRUE, constant = 1.4826) {
  med <- median(x, na.rm = na.rm)
  return(constant * median(abs(x - med), na.rm = na.rm))
}

#' Create a progress bar for lengthy operations
#'
#' @param total Total number of iterations
#' @param title Title for the progress bar (default: "Progress")
#' @param width Width of the progress bar (default: 50)
#' @return Function that updates the progress bar
#' @examples
#' \dontrun{
#' n <- 100
#' pb <- progress_bar(n, "Processing")
#' for (i in 1:n) {
#'   # Do some work
#'   Sys.sleep(0.01)
#'   pb(i)  # Update progress bar
#' }
#' }
progress_bar <- function(total, title = "Progress", width = 50) {
  pb_format <- paste0("\r", title, ": [%-", width, "s] %3d%%")
  
  last_printed <- -1
  
  function(current) {
    # Calculate percentage
    pct <- round(current / total * 100)
    
    # Only update if percentage changed
    if (pct != last_printed) {
      # Calculate progress bar
      filled <- round(width * current / total)
      bar <- paste(rep("=", filled), collapse = "")
      
      # Print progress bar
      cat(sprintf(pb_format, bar, pct))
      
      # If 100%, add newline
      if (pct >= 100) {
        cat("\n")
      }
      
      # Update last printed
      last_printed <<- pct
    }
  }
}

#' Format execution time in a human-readable way
#'
#' @param seconds Time in seconds
#' @return String with formatted time
#' @examples
#' \dontrun{
#' format_time(65)  # Returns "1m 5s"
#' }
format_time <- function(seconds) {
  if (seconds < 60) {
    return(sprintf("%.1fs", seconds))
  } else if (seconds < 3600) {
    minutes <- floor(seconds / 60)
    secs <- seconds - minutes * 60
    return(sprintf("%dm %.0fs", minutes, secs))
  } else {
    hours <- floor(seconds / 3600)
    minutes <- floor((seconds - hours * 3600) / 60)
    return(sprintf("%dh %dm", hours, minutes))
  }
}

#' Generate diverse starting points for regime-switching model estimation
#'
#' @param y Observed time series data
#' @param K Number of regimes
#' @param model_type Type of model ("constant", "tvp", "exogenous", "gas")
#' @param n_starts Number of starting points to generate
#' @param diag_probs If TRUE, use diagonal transition probability parameterization
#' @param equal_variances If TRUE, use single shared variance parameter
#' @param seed Random seed for reproducibility (optional)
#' @return List of parameter vectors with proper attributes, each suitable for the specified model type
#' @details
#' Generates diverse starting points for maximum likelihood estimation to help
#' avoid local optima. Uses data-driven heuristics to create reasonable initial
#' guesses while adding random variation for diversity.
#' 
#' UPDATED to support both diagonal and off-diagonal transition probability
#' parameterizations, as well as equal variance constraints. All returned
#' parameter vectors have proper attributes set for robust handling.
#' 
#' Strategy:
#' - Means: Based on data quantiles with random noise
#' - Variances: Fractions of data variance (single if equal_variances=TRUE)
#' - Transition probabilities: Generated according to diag_probs setting
#' - A coefficients: Random values in \[0, 0.5\] for moderate sensitivity
#' - B coefficients: Random values in \[0.5, 0.99\] for persistence
#'
#' @export
#' @examples
#' # Generate starting points for diagonal parameterization (original style)
#' y <- rnorm(200)
#' starts_diag <- generate_starting_points(y, K=2, model_type="constant",
#'                                          n_starts=3, diag_probs=TRUE)
#'
#' # Generate starting points for off-diagonal parameterization (new style)
#' starts_offdiag <- generate_starting_points(y, K=3, model_type="tvp",
#'                                             n_starts=3, diag_probs=FALSE)
generate_starting_points <- function(y, K, model_type = c("constant", "tvp", "exogenous", "gas"), 
                                     n_starts, diag_probs = TRUE, equal_variances = FALSE, 
                                     seed = NULL) {
  # Validate inputs
  if (!is.numeric(y) || length(y) == 0) {
    stop("y must be a non-empty numeric vector")
  }
  
  if (!is.numeric(K) || K < 2 || K != as.integer(K)) {
    stop("K must be an integer >= 2")
  }
  
  if (!is.numeric(n_starts) || n_starts < 1 || n_starts != as.integer(n_starts)) {
    stop("n_starts must be a positive integer")
  }
  
  model_type <- match.arg(model_type)
  
  # Validate logical parameters
  if (!is.logical(diag_probs) || length(diag_probs) != 1) {
    stop("diag_probs must be a single logical value")
  }
  if (!is.logical(equal_variances) || length(equal_variances) != 1) {
    stop("equal_variances must be a single logical value")
  }
  
  # Set seed for reproducibility if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Clean data and calculate basic statistics
  y_clean <- y[!is.na(y) & !is.nan(y) & is.finite(y)]
  if (length(y_clean) == 0) {
    stop("No valid data points after removing NA/NaN/Inf values")
  }
  
  y_mean <- mean(y_clean)
  y_var <- var(y_clean)
  y_sd <- sqrt(y_var)
  
  # Ensure variance is positive
  if (y_var <= 0) {
    warning("Data has zero or negative variance. Using unit variance for starting points.")
    y_var <- 1.0
    y_sd <- 1.0
  }
  
  # Calculate number of transition probabilities based on parameterization
  if (diag_probs) {
    n_transition <- K  # Diagonal elements only
  } else {
    n_transition <- K * (K - 1)  # Off-diagonal elements
  }
  
  # Generate starting points
  #
  # Scaling strategy: start 1 uses sensible defaults (no noise). Subsequent
  # starts add randomisation that widens with n_starts so that investing in
  # more starts systematically explores a larger region of the parameter space.
  # The log-scaling formula  base * log(n_starts) / log(7)  equals the base
  # value when n_starts = 7 and grows slowly beyond that.

  starting_points <- vector("list", n_starts)

  for (i in 1:n_starts) {

    # ---- 1. Means (mu) ----
    # Quantile anchors spread regimes across the data range.
    # Start 1: pure anchors (no noise). Subsequent starts add Gaussian noise
    # whose SD scales with n_starts.
    quantile_positions <- seq(0.1, 0.9, length.out = K)
    base_means <- quantile(y_clean, probs = quantile_positions)

    if (i == 1 || n_starts == 1) {
      mu_start <- as.numeric(base_means)
    } else {
      noise_sd <- 0.5 * y_sd * log(n_starts) / log(7)
      mu_start <- as.numeric(base_means) + rnorm(K, mean = 0, sd = noise_sd)
    }

    # ---- 2. Variances (sigma2) ----
    # Sampled on the log scale so that the spread is multiplicative.
    # Start 1: data variance for every regime. Subsequent starts draw
    # log(sigma2) ~ Normal(log(y_var), log_sd) with log_sd scaling.
    n_sigma2 <- if (equal_variances) 1 else K

    if (i == 1 || n_starts == 1) {
      sigma2_start <- rep(y_var, n_sigma2)
    } else {
      log_sd <- 0.5 * log(n_starts) / log(7)
      sigma2_start <- exp(rnorm(n_sigma2, mean = log(y_var), sd = log_sd))
    }

    # ---- 3. Transition probabilities ----
    if (diag_probs) {
      # Diagonal parameterisation: probabilities in (0, 1).
      # Start 1: 0.8 (sensible persistent default). Subsequent starts use
      # Beta distribution centered at 0.8 whose concentration decreases
      # (i.e. spread increases) with n_starts.
      if (i == 1 || n_starts == 1) {
        trans_prob_start <- rep(0.8, K)
      } else {
        concentration <- 7 * log(7) / log(n_starts)
        a_beta <- concentration * 0.8
        b_beta <- concentration * 0.2
        trans_prob_start <- rbeta(K, shape1 = a_beta, shape2 = b_beta)
        # Clamp to avoid extreme values that cause numerical issues
        trans_prob_start <- pmin(pmax(trans_prob_start, 0.01), 0.99)
      }

    } else {
      # Off-diagonal logit parameterisation: probabilities in (0, 1).
      # Small values ensure row sums don't exceed 1 (diagonal stays positive).
      # Start 1: sensible default (0.1). Subsequent starts use Beta distribution
      # centered at 0.15 whose concentration decreases with n_starts.
      max_per_entry <- 0.9 / (K - 1)  # Leave at least 10% for diagonal
      if (i == 1 || n_starts == 1) {
        trans_prob_start <- rep(0.1, n_transition)
      } else {
        center <- min(0.15, max_per_entry * 0.5)
        concentration <- 7 * log(7) / log(n_starts)
        a_beta <- concentration * center
        b_beta <- concentration * (1 - center)
        trans_prob_start <- rbeta(n_transition, shape1 = a_beta, shape2 = b_beta)
        # Clamp to valid range
        trans_prob_start <- pmin(pmax(trans_prob_start, 0.01), min(0.3, max_per_entry))
      }
    }

    # ---- 4. Build parameter vector by model type ----
    if (model_type == "constant") {
      # Constant model: [mu, sigma2, trans_prob]
      param_vector <- c(mu_start, sigma2_start, trans_prob_start)

    } else if (model_type %in% c("tvp", "exogenous")) {
      # TVP/Exogenous models: [mu, sigma2, trans_prob, A]
      # A coefficients are unbounded. First starting point is 0 (matching original),
      # subsequent points use normal distribution with log-scaled SD for exploration
      if (i == 1 || n_starts == 1) {
        A_start <- rep(0, n_transition)
      } else {
        # SD scales logarithmically: SD=1 when n_starts=7
        sd_A <- log(n_starts) / log(7)
        A_start <- rnorm(n_transition, mean = 0, sd = sd_A)
      }
      param_vector <- c(mu_start, sigma2_start, trans_prob_start, A_start)

    } else if (model_type == "gas") {
      # GAS model: [mu, sigma2, trans_prob, A, B]
      # A coefficients: same log-scaled strategy as TVP
      if (i == 1 || n_starts == 1) {
        A_start <- rep(0, n_transition)
      } else {
        sd_A <- log(n_starts) / log(7)
        A_start <- rnorm(n_transition, mean = 0, sd = sd_A)
      }
      # B coefficients (persistence): Beta distribution mirroring diagonal
      # transition probs, centered at 0.8.
      if (i == 1 || n_starts == 1) {
        B_start <- rep(0.8, n_transition)
      } else {
        concentration <- 7 * log(7) / log(n_starts)
        a_beta <- concentration * 0.8
        b_beta <- concentration * 0.2
        B_start <- rbeta(n_transition, shape1 = a_beta, shape2 = b_beta)
        B_start <- pmin(pmax(B_start, 0.01), 0.99)
      }
      param_vector <- c(mu_start, sigma2_start, trans_prob_start, A_start, B_start)
    }
    
    # 5. Set proper attributes using the new parameter system
    tryCatch({
      param_vector <- set_parameter_attributes(
        par = param_vector,
        K = K,
        model_type = model_type,
        diag_probs = diag_probs,
        equal_variances = equal_variances,
        parameterization = "natural"
      )
      
      starting_points[[i]] <- param_vector
      
    }, error = function(e) {
      # If validation fails, create a safe fallback
      warning(paste("Starting point", i, "validation failed, using fallback:", e$message))
      
      # Create conservative fallback with proper structure
      mu_fallback <- seq(y_mean - y_sd, y_mean + y_sd, length.out = K)
      
      if (equal_variances) {
        sigma2_fallback <- y_var
      } else {
        sigma2_fallback <- rep(y_var, K)
      }
      
      if (diag_probs) {
        trans_prob_fallback <- rep(0.8, K)  # High persistence
      } else {
        trans_prob_fallback <- rep(0.1, n_transition)  # Small off-diagonal probs
      }
      
      if (model_type == "constant") {
        fallback_vector <- c(mu_fallback, sigma2_fallback, trans_prob_fallback)
      } else if (model_type %in% c("tvp", "exogenous")) {
        A_fallback <- rep(0, n_transition)  # Match original: A=0 as neutral starting point
        fallback_vector <- c(mu_fallback, sigma2_fallback, trans_prob_fallback, A_fallback)
      } else if (model_type == "gas") {
        A_fallback <- rep(0, n_transition)  # Match original: A=0 as neutral starting point
        B_fallback <- rep(0.9, n_transition)
        fallback_vector <- c(mu_fallback, sigma2_fallback, trans_prob_fallback, A_fallback, B_fallback)
      }
      
      # Set attributes for fallback vector
      fallback_vector <- set_parameter_attributes(
        par = fallback_vector,
        K = K,
        model_type = model_type,
        diag_probs = diag_probs,
        equal_variances = equal_variances,
        parameterization = "natural"
      )
      
      starting_points[[i]] <<- fallback_vector
    })
  }
  
  # Add metadata about the generation process
  attr(starting_points, "generation_info") <- list(
    K = K,
    model_type = model_type,
    diag_probs = diag_probs,
    equal_variances = equal_variances,
    n_starts = n_starts,
    data_characteristics = list(
      n_obs = length(y_clean),
      mean = y_mean,
      variance = y_var,
      range = range(y_clean)
    ),
    generation_seed = seed
  )
  
  return(starting_points)
}


# --- Early Stopping Utilities for Multi-Start Optimization ---

#' Create early stopping configuration
#'
#' Builds a configuration list for the early stopping mechanism used during
#' multi-start optimization. Internal function, not exported.
#'
#' @param enabled Logical; whether early stopping is active.
#' @param patience Integer; number of evaluations without improvement before
#'   stopping. Default 500.
#' @param rel_tol Numeric; minimum relative improvement to count as progress.
#'   Default 1e-8.
#' @param max_objective Numeric; objective values above this trigger immediate
#'   stopping. Default 1e10.
#' @param max_evals Integer; maximum total evaluations per start. Default 50000.
#' @return A list with the configuration values.
#' @keywords internal
create_early_stop_config <- function(enabled = FALSE,
                                     patience = 500L,
                                     rel_tol = 1e-8,
                                     max_objective = 1e10,
                                     max_evals = 50000L) {
  list(
    enabled = enabled,
    patience = as.integer(patience),
    rel_tol = rel_tol,
    max_objective = max_objective,
    max_evals = as.integer(max_evals)
  )
}

#' Wrap an objective function with early stopping tracking
#'
#' Creates a closure around the base objective function that monitors evaluation
#' history and returns a large penalty value when divergence is detected.
#' The tracker environment is mutable and shared between the wrapper and caller.
#'
#' @param base_objective_fn The original objective function (takes a parameter
#'   vector, returns a scalar).
#' @param early_stop_config Configuration list from
#'   \code{create_early_stop_config}.
#' @return A list with two elements:
#'   \describe{
#'     \item{objective}{The wrapped objective function.}
#'     \item{tracker}{An environment with fields: \code{eval_count},
#'       \code{best_value}, \code{evals_since_improvement},
#'       \code{early_stopped}, \code{stop_reason}.}
#'   }
#' @keywords internal
create_early_stop_objective <- function(base_objective_fn, early_stop_config) {
  tracker <- new.env(parent = emptyenv())
  tracker$eval_count <- 0L
  tracker$best_value <- Inf
  tracker$evals_since_improvement <- 0L
  tracker$early_stopped <- FALSE
  tracker$stop_reason <- ""

  objective_fn <- function(par_t) {
    # If already stopped, return penalty immediately
    if (tracker$early_stopped) {
      return(1e20)
    }

    tracker$eval_count <- tracker$eval_count + 1L

    # Call the real objective
    value <- base_objective_fn(par_t)

    # Non-finite values get penalty (existing behavior, now tracked)
    if (!is.finite(value)) {
      return(1e20)
    }

    # Criterion: explosion
    if (value > early_stop_config$max_objective) {
      tracker$early_stopped <- TRUE
      tracker$stop_reason <- "objective_explosion"
      return(1e20)
    }

    # Criterion: stagnation
    if (is.infinite(tracker$best_value)) {
      # First finite evaluation always counts as improvement
      improvement <- 1
    } else {
      improvement <- (tracker$best_value - value) /
        (abs(tracker$best_value) + 1e-10)
    }
    if (improvement > early_stop_config$rel_tol) {
      tracker$best_value <- value
      tracker$evals_since_improvement <- 0L
    } else {
      tracker$evals_since_improvement <- tracker$evals_since_improvement + 1L
    }

    if (tracker$evals_since_improvement > early_stop_config$patience) {
      tracker$early_stopped <- TRUE
      tracker$stop_reason <- "stagnation"
      return(1e20)
    }

    # Criterion: budget exceeded
    if (tracker$eval_count > early_stop_config$max_evals) {
      tracker$early_stopped <- TRUE
      tracker$stop_reason <- "max_evals_exceeded"
      return(1e20)
    }

    return(value)
  }

  list(objective = objective_fn, tracker = tracker)
}
