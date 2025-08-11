#' Post-Estimation Analysis Functions for Regime-Switching Models
#' 
#' This file contains functions for analyzing and summarizing results from
#' regime-switching models with time-varying transition probabilities.
#' These functions work with output from all implemented model types:
#' - Constant transition probability models
#' - Time-varying probability (TVP) models  
#' - Exogenous variable-driven models
#' - Generalized Autoregressive Score (GAS) models
#'
#' All functions expect filtered probabilities in the standard format:
#' a T×K matrix where T is the number of time periods and K is the number
#' of regimes, with each row representing probabilities at time t.
#'
#' Main functions:
#' - calculate_persistence(): Regime persistence and transition metrics
#' - [Future]: Model comparison utilities, forecasting metrics, etc.
#'
#' Dependencies: None (base R only)

#' Helper function to calculate regime persistence metrics
#'
#' @param filtered_probs Matrix of filtered probabilities (T x K)
#' @return List with persistence metrics
#' @details
#' Calculates various metrics about regime persistence and transitions
#' from the filtered probability sequences.
calculate_persistence <- function(filtered_probs) {
  if (!is.matrix(filtered_probs)) {
    stop("filtered_probs must be a matrix")
  }
  
  T_obs <- nrow(filtered_probs)
  K <- ncol(filtered_probs)
  
  # Get most likely regime at each time point
  regime_sequence <- apply(filtered_probs, 1, which.max)
  
  # Count transitions
  transitions <- sum(diff(regime_sequence) != 0)
  transition_rate <- transitions / (T_obs - 1)
  
  # Calculate average durations in each regime
  rle_result <- rle(regime_sequence)
  regime_durations <- split(rle_result$lengths, rle_result$values)
  
  avg_durations <- numeric(K)
  for (k in 1:K) {
    if (k %in% names(regime_durations)) {
      avg_durations[k] <- mean(regime_durations[[as.character(k)]])
    } else {
      avg_durations[k] <- 0
    }
  }
  
  # Calculate regime occupancy percentages
  regime_counts <- table(factor(regime_sequence, levels = 1:K))
  regime_percentages <- as.numeric(regime_counts) / T_obs * 100
  
  return(list(
    num_transitions = transitions,
    transition_rate = transition_rate,
    avg_durations = avg_durations,
    regime_percentages = regime_percentages,
    regime_sequence = regime_sequence
  ))
}
