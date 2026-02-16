# Helper functions for parameter recovery tests
# Auto-sourced by testthat before any test file runs.

# --- Constants ----------------------------------------------------------------

RECOVERY_THRESHOLDS <- list(
  mu         = list(method = "relative", tol = 0.15),
  sigma2     = list(method = "relative", tol = 0.20),
  A          = list(method = "absolute", tol = 0.20),
  B          = list(method = "absolute", tol = 0.15),
  trans_prob = list(method = "absolute", tol = 0.15)
)

RECOVERY_NEAR_ZERO_THRESHOLD <- 0.05
RECOVERY_NEAR_ZERO_ABS_TOL   <- 0.15

RECOVERY_SIM_BURNIN <- 200L
RECOVERY_FILTER_B   <- 100L
RECOVERY_FILTER_C   <- 10L
RECOVERY_N_STARTS   <- 10L

# --- Permutation helpers ------------------------------------------------------

recovery_all_permutations <- function(K) {
  if (K == 2L) {
    list(c(1L, 2L), c(2L, 1L))
  } else if (K == 3L) {
    list(
      c(1L, 2L, 3L), c(1L, 3L, 2L), c(2L, 1L, 3L),
      c(2L, 3L, 1L), c(3L, 1L, 2L), c(3L, 2L, 1L)
    )
  } else {
    stop("Permutations only implemented for K = 2 and K = 3")
  }
}

recovery_extract_offdiag_row_major <- function(mat) {
  K <- nrow(mat)
  result <- numeric(K * (K - 1L))
  idx <- 1L
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      if (i != j) {
        result[idx] <- mat[i, j]
        idx <- idx + 1L
      }
    }
  }
  result
}

recovery_build_matrix_from_offdiag <- function(offdiag, K) {
  mat <- matrix(0, K, K)
  idx <- 1L
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      if (i != j) {
        mat[i, j] <- offdiag[idx]
        idx <- idx + 1L
      }
    }
  }
  mat
}

# --- Label switching resolution -----------------------------------------------

recovery_resolve_label_switching <- function(mu_true, est_result, K,
                                             model_type, diag_probs,
                                             equal_variances) {
  perms <- recovery_all_permutations(K)
  mu_est <- est_result$mu_est

  best_perm <- NULL
  best_cost <- Inf
  for (perm in perms) {
    cost <- sum(abs(mu_est[perm] - mu_true))
    if (cost < best_cost) {
      best_cost <- cost
      best_perm <- perm
    }
  }

  aligned <- list()
  aligned$mu <- mu_est[best_perm]

  if (equal_variances) {
    aligned$sigma2 <- est_result$sigma2_est
  } else {
    aligned$sigma2 <- est_result$sigma2_est[best_perm]
  }

  trans_est <- est_result$init_trans_est
  if (diag_probs) {
    aligned$trans_prob <- trans_est[best_perm]
  } else {
    t_mat <- transition_matrix(trans_est, diag_probs = FALSE,
                               check_validity = FALSE)
    t_mat_perm <- t_mat[best_perm, best_perm]
    aligned$trans_prob <- recovery_extract_offdiag_row_major(t_mat_perm)
  }

  if (model_type %in% c("tvp", "exogenous", "gas")) {
    A_est <- est_result$A_est
    if (diag_probs) {
      aligned$A <- A_est[best_perm]
    } else {
      A_mat <- recovery_build_matrix_from_offdiag(A_est, K)
      A_mat_perm <- A_mat[best_perm, best_perm]
      aligned$A <- recovery_extract_offdiag_row_major(A_mat_perm)
    }
  }

  if (model_type == "gas") {
    B_est <- est_result$B_est
    if (diag_probs) {
      aligned$B <- B_est[best_perm]
    } else {
      B_mat <- recovery_build_matrix_from_offdiag(B_est, K)
      B_mat_perm <- B_mat[best_perm, best_perm]
      aligned$B <- recovery_extract_offdiag_row_major(B_mat_perm)
    }
  }

  aligned$permutation <- best_perm
  aligned
}

# --- Accuracy checking --------------------------------------------------------

recovery_check_accuracy <- function(true_vals, est_vals,
                                    method = "relative", tol = 0.10) {
  if (length(true_vals) != length(est_vals)) {
    return(list(pass = FALSE, errors = NA,
                detail = "Length mismatch"))
  }

  n <- length(true_vals)
  errors   <- numeric(n)
  pass_vec <- logical(n)

  for (i in seq_len(n)) {
    if (method == "relative" && abs(true_vals[i]) < RECOVERY_NEAR_ZERO_THRESHOLD) {
      errors[i]   <- abs(est_vals[i] - true_vals[i])
      pass_vec[i] <- errors[i] < RECOVERY_NEAR_ZERO_ABS_TOL
    } else if (method == "relative") {
      errors[i]   <- abs(est_vals[i] - true_vals[i]) / abs(true_vals[i])
      pass_vec[i] <- errors[i] < tol
    } else {
      errors[i]   <- abs(est_vals[i] - true_vals[i])
      pass_vec[i] <- errors[i] < tol
    }
  }

  detail <- paste0(
    "true=[", paste(round(true_vals, 4), collapse = ", "), "] ",
    "est=[", paste(round(est_vals, 4), collapse = ", "), "] ",
    "err=[", paste(round(errors, 4), collapse = ", "), "] ",
    "tol=", tol, " method=", method, " pass=", all(pass_vec)
  )

  list(pass = all(pass_vec), errors = errors, detail = detail)
}

# --- Scenario builder ---------------------------------------------------------

recovery_build_scenario <- function(level, model_type) {

  attr_model_type <- if (model_type == "exo") "exogenous" else model_type

  if (level == 1L) {
    K <- 2L; diag_probs <- TRUE; equal_variances <- TRUE; N <- 1000L
    pars <- list(
      constant = c(-1, 1, 0.5, 0.8, 0.9),
      tvp      = c(-1, 1, 0.5, 0.8, 0.9, 0.1, -0.1),
      exo      = c(-1, 1, 0.5, 0.8, 0.9, 0.2, -0.2),
      gas      = c(-1, 1, 0.5, 0.8, 0.9, 0.1, -0.1, 0.9, 0.85)
    )
  } else if (level == 2L) {
    K <- 2L; diag_probs <- FALSE; equal_variances <- FALSE; N <- 1000L
    pars <- list(
      constant = c(-1, 1, 0.3, 0.7, 0.2, 0.15),
      tvp      = c(-1, 1, 0.3, 0.7, 0.2, 0.15, 0.1, -0.1),
      exo      = c(-1, 1, 0.3, 0.7, 0.2, 0.15, 0.15, -0.15),
      gas      = c(-1, 1, 0.3, 0.7, 0.2, 0.15, 0.05, -0.05, 0.9, 0.85)
    )
  } else if (level == 3L) {
    K <- 3L; diag_probs <- FALSE; equal_variances <- TRUE; N <- 1500L
    base <- c(-2, 0, 2, 0.5, 0.08, 0.08, 0.10, 0.10, 0.06, 0.06)
    pars <- list(
      constant = base,
      tvp      = c(base, 0.05, -0.03, 0.04, -0.04, 0.03, -0.05),
      exo      = c(base, 0.08, -0.04, 0.05, -0.06, 0.04, -0.07),
      gas      = c(base, rep(0.03, 6), rep(0.85, 6))
    )
  } else if (level == 4L) {
    K <- 3L; diag_probs <- FALSE; equal_variances <- FALSE; N <- 2000L
    base <- c(-2, 0, 2, 0.3, 0.5, 0.8, 0.08, 0.08, 0.10, 0.10, 0.06, 0.06)
    pars <- list(
      constant = base,
      tvp      = c(base, 0.05, -0.03, 0.04, -0.04, 0.03, -0.05),
      exo      = c(base, 0.08, -0.04, 0.05, -0.06, 0.04, -0.07),
      gas      = c(base, rep(0.03, 6), rep(0.85, 6))
    )
  }

  par <- set_parameter_attributes(pars[[model_type]], K = K,
                                  model_type = attr_model_type,
                                  diag_probs = diag_probs,
                                  equal_variances = equal_variances)

  list(
    true_par        = par,
    K               = K,
    model_type      = model_type,
    attr_model_type = attr_model_type,
    diag_probs      = diag_probs,
    equal_variances = equal_variances,
    N               = N,
    burn_in         = RECOVERY_SIM_BURNIN,
    level           = level
  )
}

# --- Single recovery runner ---------------------------------------------------

recovery_run_single <- function(scenario, seed = 42L) {

  par <- scenario$true_par
  K   <- scenario$K
  N   <- scenario$N

  set.seed(seed)

  # Simulate data (single run)
  if (scenario$model_type == "constant") {
    sim_data <- dataConstCD(1L, N, par, burn_in = scenario$burn_in)
  } else if (scenario$model_type == "tvp") {
    sim_data <- dataTVPCD(1L, N, par, burn_in = scenario$burn_in)
  } else if (scenario$model_type == "exo") {
    X_Exo_full <- rnorm(N + scenario$burn_in)
    sim_data <- dataTVPXExoCD(1L, N, par, X_Exo_full,
                              burn_in = scenario$burn_in)
  } else if (scenario$model_type == "gas") {
    sim_data <- dataGASCD(1L, N, par, burn_in = scenario$burn_in)
  }

  y <- as.numeric(sim_data[1L, ])

  # Estimate
  est_args <- list(
    y               = y,
    K               = K,
    diag_probs      = scenario$diag_probs,
    equal_variances = scenario$equal_variances,
    n_starts        = RECOVERY_N_STARTS,
    C               = RECOVERY_FILTER_C,
    parallel        = TRUE,
    seed            = seed + 1000L,
    verbose         = 0L
  )

  if (scenario$model_type == "gas") {
    est_args$B_burnin     <- RECOVERY_FILTER_B
    est_args$use_fallback <- FALSE
    est_result <- do.call(estimate_gas_model, est_args)
  } else if (scenario$model_type == "exo") {
    est_args$B     <- RECOVERY_FILTER_B
    est_args$X_Exo <- X_Exo_full[1:N]
    est_result <- do.call(estimate_exo_model, est_args)
  } else if (scenario$model_type == "tvp") {
    est_args$B <- RECOVERY_FILTER_B
    est_result <- do.call(estimate_tvp_model, est_args)
  } else {
    est_args$B <- RECOVERY_FILTER_B
    est_result <- do.call(estimate_constant_model, est_args)
  }

  # True parameters
  mu_true     <- extract_parameter_component(par, "mu")
  sigma2_true <- extract_parameter_component(par, "sigma2")
  trans_true  <- extract_parameter_component(par, "trans_prob")

  # Resolve label switching
  aligned <- recovery_resolve_label_switching(
    mu_true, est_result, K,
    model_type      = scenario$attr_model_type,
    diag_probs      = scenario$diag_probs,
    equal_variances = scenario$equal_variances
  )

  # Accuracy checks
  acc <- list()
  acc$mu     <- recovery_check_accuracy(mu_true, aligned$mu,
                  RECOVERY_THRESHOLDS$mu$method, RECOVERY_THRESHOLDS$mu$tol)
  acc$sigma2 <- recovery_check_accuracy(sigma2_true, aligned$sigma2,
                  RECOVERY_THRESHOLDS$sigma2$method,
                  RECOVERY_THRESHOLDS$sigma2$tol)
  acc$trans_prob <- recovery_check_accuracy(trans_true, aligned$trans_prob,
                     RECOVERY_THRESHOLDS$trans_prob$method,
                     RECOVERY_THRESHOLDS$trans_prob$tol)

  if (scenario$model_type %in% c("tvp", "exo", "gas")) {
    A_true <- extract_parameter_component(par, "A")
    acc$A <- recovery_check_accuracy(A_true, aligned$A,
               RECOVERY_THRESHOLDS$A$method, RECOVERY_THRESHOLDS$A$tol)
  }

  if (scenario$model_type == "gas") {
    B_true <- extract_parameter_component(par, "B")
    acc$B <- recovery_check_accuracy(B_true, aligned$B,
               RECOVERY_THRESHOLDS$B$method, RECOVERY_THRESHOLDS$B$tol)
  }

  list(
    est_result = est_result,
    aligned    = aligned,
    accuracy   = acc,
    converged  = (est_result$diagnostics$convergence_code == 0)
  )
}
