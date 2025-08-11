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
