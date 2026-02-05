#' General utility functions for regime switching models
#' These functions provide general mathematical and data processing utilities
#' that are needed across different parts of the regime switching models.

#' Logit transformation (natural space -> unconstrained space)
#'
#' @param x Value in [0,1] to transform
#' @param eps Small value to prevent boundary issues (default: .Machine$double.eps)
#' @return Transformed value in (-∞,∞)
#' @examples
#' logit(0.75)  # Returns 1.098612
logit <- function(x, eps = .Machine$double.eps) {
  # Add bounds checking to prevent errors
  x <- pmin(pmax(x, eps), 1 - eps)
  return(log(x/(1-x)))
}

#' Logistic transformation (unconstrained space -> natural space)
#'
#' @param x Value in (-∞,∞) to transform
#' @return Transformed value in (0,1)
#' @examples
#' logistic(1.098612)  # Returns 0.75
logistic <- function(x) {
  return(1/(1+exp(-x)))
}

#' Safe logarithm to handle zeros and negative values
#'
#' @param x Value to take the logarithm of
#' @param min_value Minimum value to replace zeros/negative values (default: .Machine$double.eps)
#' @return Log of max(x, min_value)
#' @examples
#' safe_log(0)      # Returns log of a small positive number instead of -Inf
#' safe_log(-5)     # Returns log of a small positive number instead of NaN
#' safe_log(c(2, 0, 3, -1))  # Handles vectors safely
safe_log <- function(x, min_value = .Machine$double.eps) {
  return(log(pmax(x, min_value)))
}

#' Safe exponentiation to prevent overflow
#'
#' @param x Value to exponentiate
#' @param max_value Maximum allowed value in exponent (default: 700)
#' @return exp(min(x, max_value))
#' @examples
#' safe_exp(1000)  # Returns a large but finite number instead of Inf
safe_exp <- function(x, max_value = 700) {
  return(exp(pmin(x, max_value)))
}

#' Check if a vector contains any NA, NaN, Inf, or -Inf values
#'
#' @param x Vector to check
#' @return TRUE if any invalid values are found, FALSE otherwise
#' @examples
#' has_invalid_values(c(1, 2, NA, 4))        # Returns TRUE
#' has_invalid_values(c(1, 2, NaN, 4))       # Returns TRUE
#' has_invalid_values(c(1, 2, Inf, 4))       # Returns TRUE
#' has_invalid_values(c(1, 2, -Inf, 4))      # Returns TRUE
#' has_invalid_values(c(1, 2, 3, 4))         # Returns FALSE
has_invalid_values <- function(x) {
  return(any(is.na(x) | is.nan(x) | is.infinite(x)))
}

#' Replace invalid values (NA, NaN, Inf, -Inf) with a default value
#'
#' @param x Vector to process
#' @param replacement Value to use for replacements (default: 0)
#' @return Vector with invalid values replaced
#' @examples
#' replace_invalid_values(c(1, 2, NA, Inf, -Inf, 6))  # Returns c(1, 2, 0, 0, 0, 6)
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
#' moving_average(1:10, 3)  # 3-period centered moving average
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
#' standardize(c(1, 2, 3, 4, 5))  # Returns values with mean 0 and SD 1
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
#' mad_value(c(1, 2, 3, 100))  # Less affected by outliers than standard deviation
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
#' n <- 100
#' pb <- progress_bar(n, "Processing")
#' for (i in 1:n) {
#'   # Do some work
#'   Sys.sleep(0.01)
#'   pb(i)  # Update progress bar
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
#' format_time(65)  # Returns "1m 5s"
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
#' - A coefficients: Random values in [0, 0.5] for moderate sensitivity
#' - B coefficients: Random values in [0.5, 0.99] for persistence
#'
#' @examples
#' # Generate starting points for diagonal parameterization (original style)
#' y <- rnorm(1000)
#' starts_diag <- generate_starting_points(y, K=2, model_type="constant", 
#'                                          n_starts=5, diag_probs=TRUE)
#' 
#' # Generate starting points for off-diagonal parameterization (new style)
#' starts_offdiag <- generate_starting_points(y, K=3, model_type="tvp", 
#'                                             n_starts=5, diag_probs=FALSE)
generate_starting_points <- function(y, K, model_type = c("constant", "tvp", "exogenous", "gas"), 
                                     n_starts, diag_probs = FALSE, equal_variances = FALSE, 
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
  starting_points <- vector("list", n_starts)
  
  for (i in 1:n_starts) {
    
    # 1. Generate means using quantiles with noise
    # Use quantiles to spread regimes across data range
    quantile_positions <- seq(0.1, 0.9, length.out = K)
    base_means <- quantile(y_clean, probs = quantile_positions)
    
    # Add random noise (±0.5 standard deviations)
    noise_scale <- 0.5 * y_sd
    mu_start <- base_means + runif(K, -noise_scale, noise_scale)
    
    # 2. Generate variances based on equal_variances setting
    if (equal_variances) {
      # Single shared variance parameter
      var_fraction <- runif(1, 0.5, 1.5)  # 50% to 150% of data variance
      sigma2_start <- var_fraction * y_var
    } else {
      # Separate variance for each regime
      var_fractions <- runif(K, 0.3, 1.5)  # 30% to 150% of data variance
      sigma2_start <- var_fractions * y_var
    }
    
    # 3. Generate transition probabilities based on parameterization
    if (diag_probs) {
      # Generate diagonal persistence probabilities
      # Higher values (0.5-0.95) for persistence, with some variation
      trans_prob_start <- runif(K, 0.5, 0.95)
      
    } else {
      # Generate off-diagonal transition probabilities
      # Start with base probability around 1/(K+1), add noise
      base_prob <- 1.0 / (K + 1)
      noise_range <- 0.1  # ±10% variation
      
      # Generate raw probabilities with noise
      raw_trans_probs <- pmax(0.01, pmin(0.8, 
                                         base_prob + runif(n_transition, -noise_range, noise_range)))
      
      # Ensure they form a valid stochastic matrix using new helper
      trans_prob_start <- convert_to_valid_probs(raw_trans_probs, diag_probs = FALSE)
    }
    
    # 4. Build parameter vector based on model type
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
      # A coefficients are unbounded. First starting point is 0 (matching original),
      # subsequent points use normal distribution with log-scaled SD for exploration
      if (i == 1 || n_starts == 1) {
        A_start <- rep(0, n_transition)
      } else {
        # SD scales logarithmically: SD=1 when n_starts=7
        sd_A <- log(n_starts) / log(7)
        A_start <- rnorm(n_transition, mean = 0, sd = sd_A)
      }
      B_start <- runif(n_transition, 0.5, 0.99) # High persistence
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
        trans_prob_fallback <- rep(0.2, n_transition)
        trans_prob_fallback <- convert_to_valid_probs(trans_prob_fallback, diag_probs = FALSE)
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
