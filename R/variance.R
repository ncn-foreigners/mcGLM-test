# ---------------------------------------------------------------------------
# Sandwich variance estimators for each method
# ---------------------------------------------------------------------------

#' Sandwich variance for the naive estimator
#'
#' V_naive = A^{-1} C A^{-1} / n, where
#' A = (1/n) sum dot_mu(eta_i) xi_hat xi_hat',
#' C = (1/n) sum eps_i^2 xi_hat xi_hat'.
#' @keywords internal
vcov_naive <- function(psi, y, z_hat, x, family) {
  fam    <- get_link_funs(family)
  xi_hat <- cbind(z_hat, x)
  n      <- length(y)
  eta    <- as.numeric(xi_hat %*% psi)
  mu_val <- fam$mu(eta)
  w      <- fam$mu_dot(eta)
  eps    <- y - mu_val

  A <- crossprod(xi_hat * w, xi_hat) / n
  C <- crossprod(xi_hat * eps, xi_hat * eps) / n
  A_inv <- solve(A)
  A_inv %*% C %*% A_inv / n
}

#' Sandwich variance for BCA/BCM estimators (binary)
#'
#' Under drifting regime, same asymptotic variance as naive: A^{-1} C A^{-1} / n
#' evaluated at the corrected estimate.
#' @keywords internal
vcov_bc_bin <- function(psi_bc, y, z_hat, x, family) {
  vcov_naive(psi_bc, y, z_hat, x, family)
}

#' Sandwich variance for corrected-score estimator (binary)
#'
#' V_cs = J^{-1} S (J^{-1})' / n, where
#' J = (1/n) sum d phi_i / d psi',
#' S = (1/n) sum phi_i phi_i'.
#' @keywords internal
vcov_cs_bin <- function(psi, y, z_hat, x, family, p01, p10, pi_z) {
  fam    <- get_link_funs(family)
  xi_hat <- cbind(z_hat, x)
  n      <- length(y)
  p      <- length(psi)

  eta_tilde <- as.numeric(xi_hat %*% psi)
  resid     <- y - fam$mu(eta_tilde)
  m_mat     <- compute_m_bin(psi, x, fam$mu, p01, p10, pi_z)

  # phi_i = xi_hat_i * resid_i - m_i(psi)
  phi_mat <- xi_hat * resid - m_mat   # n x p

  S <- crossprod(phi_mat) / n

  # J = -(I + M) where I = Ihat, M = Mhat
  I_hat <- compute_Ihat(psi, z_hat, x, fam$mu_dot)
  M_hat <- compute_Mhat_bin(psi, x, fam$mu, p01, p10, pi_z)
  J     <- -(I_hat + M_hat)
  J_inv <- solve(J)

  J_inv %*% S %*% t(J_inv) / n
}


# ---- Multicategory variance estimators ----

#' Sandwich variance for naive estimator (multicategory)
#' @keywords internal
vcov_naive_multi <- function(psi, y, z_hat, x, K, family) {
  fam <- get_link_funs(family)
  n   <- length(y)
  s   <- K - 1
  r   <- ncol(x)

  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) d_hat[, k] <- as.numeric(z_hat == k)
  xi_hat <- cbind(d_hat, x)

  gamma <- c(0, psi[seq_len(s)])
  alpha <- psi[(s + 1):(s + r)]
  eta_base  <- as.numeric(x %*% alpha)
  eta_tilde <- eta_base + gamma[z_hat + 1]

  w   <- fam$mu_dot(eta_tilde)
  eps <- y - fam$mu(eta_tilde)

  A <- crossprod(xi_hat * w, xi_hat) / n
  C <- crossprod(xi_hat * eps, xi_hat * eps) / n
  A_inv <- solve(A)
  A_inv %*% C %*% A_inv / n
}

#' @keywords internal
vcov_bc_multi <- function(psi_bc, y, z_hat, x, K, family) {
  vcov_naive_multi(psi_bc, y, z_hat, x, K, family)
}

#' Sandwich variance for corrected-score estimator (multicategory)
#' @keywords internal
vcov_cs_multi <- function(psi, y, z_hat, x, K, family, Pi, pi_z) {
  fam <- get_link_funs(family)
  n   <- length(y)
  s   <- K - 1
  r   <- ncol(x)
  p   <- s + r

  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) d_hat[, k] <- as.numeric(z_hat == k)
  xi_hat <- cbind(d_hat, x)

  gamma <- c(0, psi[seq_len(s)])
  alpha <- psi[(s + 1):(s + r)]
  eta_base  <- as.numeric(x %*% alpha)
  eta_tilde <- eta_base + gamma[z_hat + 1]

  resid <- y - fam$mu(eta_tilde)
  m_mat <- compute_m_multi(psi, x, K, fam$mu, Pi, pi_z)

  phi_mat <- xi_hat * resid - m_mat
  S <- crossprod(phi_mat) / n

  I_hat <- compute_Ihat_multi(psi, z_hat, x, K, fam$mu_dot)
  M_hat <- compute_Mhat_multi(psi, x, K, fam$mu, Pi, pi_z)
  J     <- -(I_hat + M_hat)
  J_inv <- solve(J)

  J_inv %*% S %*% t(J_inv) / n
}
