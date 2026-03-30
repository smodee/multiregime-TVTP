#' Check if C filtering backend is available
#'
#' @return Logical scalar: TRUE if C backend is available, FALSE otherwise.
#' @export
cpp_available <- function() {
  is.loaded("C_filtering_Const", PACKAGE = "multiregimeTVTP")
}

# R wrappers for C filtering functions (called via .Call)

Cfiltering_Const <- function(mu, sigma2, trans_prob, y, n_burnin, n_cutoff, diag_probs) {
  .Call("C_filtering_Const", as.double(mu), as.double(sigma2), as.double(trans_prob),
        as.double(y), as.integer(n_burnin), as.integer(n_cutoff), as.logical(diag_probs),
        PACKAGE = "multiregimeTVTP")
}

Cfiltering_TVP <- function(mu, sigma2, init_trans, A, y, n_burnin, n_cutoff, diag_probs) {
  .Call("C_filtering_TVP", as.double(mu), as.double(sigma2), as.double(init_trans),
        as.double(A), as.double(y), as.integer(n_burnin), as.integer(n_cutoff),
        as.logical(diag_probs), PACKAGE = "multiregimeTVTP")
}

Cfiltering_Exo <- function(mu, sigma2, init_trans, A, y, X_Exo, n_burnin, n_cutoff, diag_probs) {
  .Call("C_filtering_Exo", as.double(mu), as.double(sigma2), as.double(init_trans),
        as.double(A), as.double(y), as.double(X_Exo), as.integer(n_burnin),
        as.integer(n_cutoff), as.logical(diag_probs), PACKAGE = "multiregimeTVTP")
}

Cfiltering_GAS <- function(mu, sigma2, init_trans, A, B, y, n_burnin, n_cutoff,
                            diag_probs, gh_nodes, gh_weights, mu_quad, sigma_quad,
                            use_fallback, A_threshold) {
  .Call("C_filtering_GAS", as.double(mu), as.double(sigma2), as.double(init_trans),
        as.double(A), as.double(B), as.double(y), as.integer(n_burnin),
        as.integer(n_cutoff), as.logical(diag_probs), as.double(gh_nodes),
        as.double(gh_weights), as.double(mu_quad), as.double(sigma_quad),
        as.logical(use_fallback), as.double(A_threshold), PACKAGE = "multiregimeTVTP")
}
