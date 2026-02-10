#' Parallel Computing Utilities
#' 
#' This file provides utilities for conservative parallel processing across
#' all simulation studies. Functions handle core detection, allocation strategy,
#' and graceful fallback when parallel packages are unavailable.
#'
#' Main functions:
#' - get_conservative_cores(): Detect cores with conservative reservation
#' - get_parallel_allocation(): Determine optimal worker and core allocation
#'
#' Dependencies: parallel package (graceful fallback if unavailable)

#' Get conservative core count for parallel processing
#'
#' @param max_cores User-specified maximum cores (NULL for auto-detection)
#' @param reserve_cores Number of cores to reserve for system (default: 2)
#' @param min_cores Minimum cores to use even if detection fails (default: 1)
#' @return Conservative number of cores to use safely
#' @details
#' Implements a conservative approach to core allocation that:
#' - Reserves cores for system use to avoid oversubscription
#' - Handles edge cases where core detection fails
#' - Never exceeds available cores
#' - Provides sensible fallbacks
#' 
#' The reserve_cores parameter helps avoid the "200% CPU load" warnings
#' that occur when parallel workers compete for limited resources.
#'
#' @export
#' @examples
#' # Auto-detect with 2 cores reserved for system
#' cores <- get_conservative_cores()
#'
#' # Manual limit with 1 core reserved
#' cores <- get_conservative_cores(max_cores = 4, reserve_cores = 1)
#'
#' # Force single core (useful for debugging)
#' cores <- get_conservative_cores(max_cores = 1, reserve_cores = 0)
get_conservative_cores <- function(max_cores = NULL, reserve_cores = 2, min_cores = 1) {
  
  # Check if parallel package is available
  if (!requireNamespace("parallel", quietly = TRUE)) {
    if (!is.null(max_cores) && max_cores > 1) {
      warning("parallel package not available, using single core")
    }
    return(1)
  }
  
  # Detect available cores
  detected_cores <- tryCatch({
    parallel::detectCores()
  }, error = function(e) {
    NA
  })
  
  # Handle detection failure
  if (is.na(detected_cores) || detected_cores < 1) {
    warning("Could not detect CPU cores, defaulting to single core")
    return(1)
  }
  
  # Calculate conservative core count
  if (is.null(max_cores)) {
    # Auto mode: use detected cores minus reserve
    conservative_cores <- max(min_cores, detected_cores - reserve_cores)
  } else {
    # Manual mode: respect user limit but still be conservative
    conservative_cores <- max(min_cores, min(max_cores, detected_cores - reserve_cores))
  }
  
  # Final safety check: never exceed actual detected cores
  conservative_cores <- max(min_cores, min(conservative_cores, detected_cores))
  
  return(conservative_cores)
}

#' Calculate optimal parallel allocation for simulation studies
#'
#' @param total_cores Total cores available (from get_conservative_cores)
#' @param num_repetitions Number of simulation repetitions to run
#' @param n_starts Number of starting points per estimation
#' @return List with sim_workers and cores_per_sim allocation
#' @details
#' Determines the optimal allocation between simulation-level and estimation-level
#' parallelization based on available resources and workload characteristics.
#' 
#' Strategy selection logic:
#' - Calculates potential throughput for both simulation-level and estimation-level parallelization
#' - Chooses simulation-level when num_repetitions provides better parallelization opportunity
#' - Chooses estimation-level when n_starts provides better parallelization opportunity
#' - Always avoids nested parallelization (coordinated mode)
#' 
#' The function ensures efficient core utilization while maintaining simple,
#' single-level parallelization to avoid resource conflicts.
#'
#' @export
#' @examples
#' # Typical usage in simulation study
#' total_cores <- get_conservative_cores(reserve_cores = 2)
#' allocation <- get_parallel_allocation(total_cores, num_repetitions = 50, n_starts = 3)
#'
#' # Result interpretation:
#' # allocation$sim_workers = 1, cores_per_sim = 8  -> Estimation-level parallelization
#' # allocation$sim_workers = 8, cores_per_sim = 1  -> Simulation-level parallelization
#' # allocation$sim_workers = 1, cores_per_sim = 1  -> Sequential (no parallelization)
get_parallel_allocation <- function(total_cores, num_repetitions, n_starts) {
  
  # Input validation
  if (!is.numeric(total_cores) || total_cores < 1) {
    stop("total_cores must be a positive integer")
  }
  if (!is.numeric(num_repetitions) || num_repetitions < 1) {
    stop("num_repetitions must be a positive integer")
  }
  if (!is.numeric(n_starts) || n_starts < 1) {
    stop("n_starts must be a positive integer")
  }
  
  total_cores <- as.integer(total_cores)
  num_repetitions <- as.integer(num_repetitions)
  n_starts <- as.integer(n_starts)
  
  # Handle single core case
  if (total_cores == 1) {
    return(list(
      sim_workers = 1L,
      cores_per_sim = 1L
    ))
  }
  
  # Handle single repetition case - must use estimation-level
  if (num_repetitions == 1) {
    return(list(
      sim_workers = 1L,
      cores_per_sim = total_cores
    ))
  }
  
  # Calculate potential throughput for each strategy
  sim_throughput <- min(total_cores, num_repetitions)  # Max simulations in parallel
  est_throughput <- min(total_cores, n_starts)        # Max estimation parallelization
  
  # Choose strategy based on which provides better parallelization opportunity
  if (sim_throughput >= est_throughput && num_repetitions > 1) {
    # Prefer simulation-level parallelization
    # Use as many simulation workers as we have cores (up to num_repetitions)
    sim_workers <- min(total_cores, num_repetitions)
    return(list(
      sim_workers = sim_workers,
      cores_per_sim = 1L  # Each simulation runs sequentially
    ))
  } else {
    # Use estimation-level parallelization
    # Single simulation worker using all available cores for estimation
    return(list(
      sim_workers = 1L,
      cores_per_sim = total_cores
    ))
  }
}

#' Check if coordinated parallelization will be used
#'
#' @param allocation Result from get_parallel_allocation()
#' @return TRUE if coordinated, FALSE if estimation-only or sequential
#' @details
#' Helper function to check parallelization strategy from allocation result.
#' Useful for conditional logic in simulation functions.
#'
#' @examples
#' \dontrun{
#' allocation <- get_parallel_allocation(6, 50, 3)
#' if (is_coordinated_parallel(allocation)) {
#'   # Set up simulation-level parallelization
#' } else {
#'   # Use estimation-only parallelization
#' }
#' }
is_coordinated_parallel <- function(allocation) {
  if (!is.list(allocation) || !all(c("sim_workers", "cores_per_sim") %in% names(allocation))) {
    stop("allocation must be a list with sim_workers and cores_per_sim elements")
  }
  
  return(allocation$sim_workers > 1)
}

#' Check if any parallelization will be used
#'
#' @param allocation Result from get_parallel_allocation()
#' @return TRUE if parallel, FALSE if sequential
#' @details
#' Helper function to check if any parallel processing will occur.
#'
#' @examples
#' \dontrun{
#' allocation <- get_parallel_allocation(1, 10, 3)
#' if (is_parallel(allocation)) {
#'   # Some form of parallelization enabled
#' } else {
#'   # Fully sequential execution
#' }
#' }
is_parallel <- function(allocation) {
  if (!is.list(allocation) || !all(c("sim_workers", "cores_per_sim") %in% names(allocation))) {
    stop("allocation must be a list with sim_workers and cores_per_sim elements")
  }
  
  return(allocation$sim_workers > 1 || allocation$cores_per_sim > 1)
}

#' Get human-readable description of parallelization strategy
#'
#' @param allocation Result from get_parallel_allocation()
#' @return Character string describing the strategy
#' @details
#' Provides a readable summary of the parallelization approach for logging
#' and user feedback.
#'
#' @examples
#' \dontrun{
#' allocation <- get_parallel_allocation(4, 20, 3)
#' cat("Strategy:", describe_parallel_strategy(allocation))
#' }
describe_parallel_strategy <- function(allocation) {
  if (!is.list(allocation) || !all(c("sim_workers", "cores_per_sim") %in% names(allocation))) {
    stop("allocation must be a list with sim_workers and cores_per_sim elements")
  }
  
  sim_workers <- allocation$sim_workers
  cores_per_sim <- allocation$cores_per_sim
  
  if (sim_workers == 1 && cores_per_sim == 1) {
    return("Sequential (no parallelization)")
  } else if (sim_workers == 1 && cores_per_sim > 1) {
    return(paste0("Estimation-only (", cores_per_sim, " cores per estimation)"))
  } else if (sim_workers > 1) {
    return(paste0("Coordinated (", sim_workers, " sim workers, ", 
                  cores_per_sim, " cores per simulation)"))
  } else {
    return("Unknown strategy")
  }
}

#' Diagnostic function to check parallel setup
#'
#' @param max_cores Optional maximum cores to test (NULL for auto-detection)
#' @param reserve_cores Number of cores to reserve (default: 2)
#' @param example_repetitions Example number of repetitions for allocation test
#' @param example_starts Example number of starting points for allocation test
#' @return Invisible list with diagnostic information
#' @details
#' Useful diagnostic function to understand how parallel allocation will work
#' on the current system. Helps debug parallelization issues.
#'
#' @export
#' @examples
#' # Check your system's parallel capabilities
#' check_parallel_setup()
#'
#' # Test with specific parameters
#' check_parallel_setup(max_cores = 4, reserve_cores = 1)
check_parallel_setup <- function(max_cores = NULL, reserve_cores = 2, 
                                 example_repetitions = 50, example_starts = 3) {
  
  cat("=== PARALLEL SETUP DIAGNOSTICS ===\n")
  
  # Core detection
  if (requireNamespace("parallel", quietly = TRUE)) {
    detected <- tryCatch(parallel::detectCores(), error = function(e) NA)
    cat("Detected cores:", ifelse(is.na(detected), "detection failed", detected), "\n")
  } else {
    cat("parallel package: NOT AVAILABLE\n")
    detected <- NA
  }
  
  # Conservative allocation
  conservative <- get_conservative_cores(max_cores, reserve_cores)
  cat("Conservative cores (reserve", reserve_cores, "):", conservative, "\n")
  
  # Example allocation
  allocation <- get_parallel_allocation(conservative, example_repetitions, example_starts)
  cat("\nExample allocation (", example_repetitions, " reps, ", example_starts, " starts):\n", sep = "")
  cat("Strategy:", describe_parallel_strategy(allocation), "\n")
  cat("Simulation workers:", allocation$sim_workers, "\n")
  cat("Cores per simulation:", allocation$cores_per_sim, "\n")
  
  # Efficiency analysis
  total_used <- allocation$sim_workers * allocation$cores_per_sim
  efficiency <- total_used / conservative * 100
  cat("Core utilization:", total_used, "/", conservative, "(", round(efficiency, 1), "%)\n")
  
  # Warnings
  if (!is.na(detected) && conservative >= detected) {
    cat("\nWARNING: Conservative allocation uses all detected cores\n")
    cat("Consider increasing reserve_cores parameter\n")
  }
  
  if (allocation$cores_per_sim > example_starts) {
    cat("\nNOTE: cores_per_sim (", allocation$cores_per_sim, 
        ") > n_starts (", example_starts, ")\n", sep = "")
    cat("Estimation parallelization will be limited by n_starts\n")
  }
  
  cat("===================================\n")
  
  # Return diagnostic info invisibly
  invisible(list(
    detected_cores = detected,
    conservative_cores = conservative,
    allocation = allocation,
    efficiency = efficiency
  ))
}
