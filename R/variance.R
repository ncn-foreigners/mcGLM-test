# ---------------------------------------------------------------------------
# Sandwich variance estimators for each method
# ---------------------------------------------------------------------------

#' Sandwich variance for the naive estimator
#'
#' V_naive = A^(-1) C A^(-1) / N, where
#' A = (1/N) sum wt_i * dot_mu(eta_i) * xi_hat * xi_hat',
#' C = (1/N) sum wt_i * eps_i^2 * xi_hat * xi_hat'.
#' @keywords internal
vcov_naive <- function(psi, y, z_hat, x, family, wt = NULL) {
  fam    <- get_link_funs(family)
  xi_hat <- cbind(z_hat, x)
  n      <- length(y)
  N      <- if (is.null(wt)) n else sum(wt)
  eta    <- as.numeric(xi_hat %*% psi)
  mu_val <- fam$mu(eta)
  w      <- fam$mu_dot(eta)
  eps    <- y - mu_val

  if (is.null(wt)) {
    A <- crossprod(xi_hat * w, xi_hat) / N
    C <- crossprod(xi_hat * eps, xi_hat * eps) / N
  } else {
    A <- crossprod(xi_hat * (wt * w), xi_hat) / N
    C <- crossprod(xi_hat * (wt * eps), xi_hat * eps) / N
  }
  A_inv <- solve(A)
  A_inv %*% C %*% A_inv / N
}

#' Sandwich variance for BCA/BCM estimators (binary)
#'
#' Under drifting regime (\code{corrected = FALSE}), uses the naive sandwich
#' A^(-1) C A^(-1) / N evaluated at the corrected estimate.
#'
#' Under fixed misclassification (\code{corrected = TRUE}), uses the
#' influence-function approach via delta method.
#'
#' @param psi_naive The naive estimator (linearization point).
#' @param psi_bc The corrected estimator (BCA or BCM).
#' @keywords internal
vcov_bc_bin <- function(psi_bc, y, z_hat, x, family,
                        p01 = NULL, p10 = NULL, pi_z = NULL,
                        psi_naive = NULL,
                        type = c("bca", "bcm"), corrected = FALSE,
                        wt = NULL) {
  if (!corrected) return(vcov_naive(psi_bc, y, z_hat, x, family, wt = wt))

  type   <- match.arg(type)
  fam    <- get_link_funs(family)
  xi_hat <- cbind(z_hat, x)
  n      <- length(y)
  N      <- if (is.null(wt)) n else sum(wt)
  p      <- length(psi_bc)

  # Evaluate residuals and Fisher info at psi_naive (the linearization point)
  if (is.null(psi_naive)) psi_naive <- psi_bc
  eta    <- as.numeric(xi_hat %*% psi_naive)
  w      <- fam$mu_dot(eta)
  resid  <- y - fam$mu(eta)

  if (is.null(wt)) {
    A_hat <- crossprod(xi_hat * w, xi_hat) / N
  } else {
    A_hat <- crossprod(xi_hat * (wt * w), xi_hat) / N
  }
  A_inv <- solve(A_hat)
  M_hat <- compute_Mhat_bin(psi_naive, x, fam$mu, p01, p10, pi_z, wt = wt)
  m_mat <- compute_m_bin(psi_naive, x, fam$mu, p01, p10, pi_z)  # n x p
  if (is.null(wt)) {
    m_bar <- colMeans(m_mat)
  } else {
    m_bar <- colSums(wt * m_mat) / N
  }

  # Jacobian of the correction map and the "H" matrix
  if (type == "bca") {
    G     <- diag(p) - A_inv %*% M_hat
    H_inv <- A_inv
  } else {
    H     <- A_hat + M_hat
    H_inv <- solve(H)
    G     <- diag(p) - H_inv %*% M_hat
  }

  # Per-observation influence function (as rows of n x p matrix):
  score_mat    <- xi_hat * resid                       # n x p (rows = u_i^T)
  centered_m   <- sweep(m_mat, 2, m_bar)               # n x p (rows = w_i^T)

  L1 <- A_inv %*% t(G)    # p x p: maps u_i -> first IF component
  L2 <- t(H_inv)          # p x p: maps w_i -> second IF component

  IF_mat <- score_mat %*% L1 - centered_m %*% L2       # n x p

  if (is.null(wt)) {
    crossprod(IF_mat) / N^2
  } else {
    crossprod(IF_mat * wt, IF_mat) / N^2
  }
}

#' Sandwich variance for corrected-score estimator (binary)
#'
#' V_cs = J^{-1} S (J^{-1})' / N, where
#' J = weighted (1/N) sum d phi_i / d psi',
#' S = weighted (1/N) sum phi_i phi_i'.
#' @keywords internal
vcov_cs_bin <- function(psi, y, z_hat, x, family, p01, p10, pi_z,
                        wt = NULL) {
  fam    <- get_link_funs(family)
  xi_hat <- cbind(z_hat, x)
  n      <- length(y)
  N      <- if (is.null(wt)) n else sum(wt)

  eta_tilde <- as.numeric(xi_hat %*% psi)
  resid     <- y - fam$mu(eta_tilde)
  m_mat     <- compute_m_bin(psi, x, fam$mu, p01, p10, pi_z)

  # phi_i = xi_hat_i * resid_i - m_i(psi)
  phi_mat <- xi_hat * resid - m_mat   # n x p

  if (is.null(wt)) {
    S <- crossprod(phi_mat) / N
  } else {
    S <- crossprod(phi_mat * wt, phi_mat) / N
  }

  # J = -(I + M) where I = Ihat, M = Mhat
  I_hat <- compute_Ihat(psi, z_hat, x, fam$mu_dot, wt = wt)
  M_hat <- compute_Mhat_bin(psi, x, fam$mu, p01, p10, pi_z, wt = wt)
  J     <- -(I_hat + M_hat)
  J_inv <- solve(J)

  J_inv %*% S %*% t(J_inv) / N
}


# ---- Multicategory variance estimators ----

#' Sandwich variance for naive estimator (multicategory)
#' @keywords internal
vcov_naive_multi <- function(psi, y, z_hat, x, K, family, wt = NULL) {
  fam <- get_link_funs(family)
  n   <- length(y)
  N   <- if (is.null(wt)) n else sum(wt)
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

  if (is.null(wt)) {
    A <- crossprod(xi_hat * w, xi_hat) / N
    C <- crossprod(xi_hat * eps, xi_hat * eps) / N
  } else {
    A <- crossprod(xi_hat * (wt * w), xi_hat) / N
    C <- crossprod(xi_hat * (wt * eps), xi_hat * eps) / N
  }
  A_inv <- solve(A)
  A_inv %*% C %*% A_inv / N
}

#' Sandwich variance for BCA/BCM estimators (multicategory)
#' @keywords internal
vcov_bc_multi <- function(psi_bc, y, z_hat, x, K, family,
                          Pi = NULL, pi_z = NULL,
                          psi_naive = NULL,
                          type = c("bca", "bcm"), corrected = FALSE,
                          wt = NULL) {
  if (!corrected) return(vcov_naive_multi(psi_bc, y, z_hat, x, K, family,
                                           wt = wt))

  type <- match.arg(type)
  fam  <- get_link_funs(family)
  n    <- length(y)
  N    <- if (is.null(wt)) n else sum(wt)
  s    <- K - 1
  r    <- ncol(x)
  p    <- s + r

  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) d_hat[, k] <- as.numeric(z_hat == k)
  xi_hat <- cbind(d_hat, x)

  if (is.null(psi_naive)) psi_naive <- psi_bc

  gamma <- c(0, psi_naive[seq_len(s)])
  alpha <- psi_naive[(s + 1):(s + r)]
  eta_base  <- as.numeric(x %*% alpha)
  eta_tilde <- eta_base + gamma[z_hat + 1]

  w     <- fam$mu_dot(eta_tilde)
  resid <- y - fam$mu(eta_tilde)

  if (is.null(wt)) {
    A_hat <- crossprod(xi_hat * w, xi_hat) / N
  } else {
    A_hat <- crossprod(xi_hat * (wt * w), xi_hat) / N
  }
  A_inv <- solve(A_hat)
  M_hat <- compute_Mhat_multi(psi_naive, x, K, fam$mu, Pi, pi_z, wt = wt)
  m_mat <- compute_m_multi(psi_naive, x, K, fam$mu, Pi, pi_z)
  if (is.null(wt)) {
    m_bar <- colMeans(m_mat)
  } else {
    m_bar <- colSums(wt * m_mat) / N
  }

  if (type == "bca") {
    G     <- diag(p) - A_inv %*% M_hat
    H_inv <- A_inv
  } else {
    H     <- A_hat + M_hat
    H_inv <- solve(H)
    G     <- diag(p) - H_inv %*% M_hat
  }

  score_mat  <- xi_hat * resid
  centered_m <- sweep(m_mat, 2, m_bar)

  L1 <- A_inv %*% t(G)
  L2 <- t(H_inv)

  IF_mat <- score_mat %*% L1 - centered_m %*% L2

  if (is.null(wt)) {
    crossprod(IF_mat) / N^2
  } else {
    crossprod(IF_mat * wt, IF_mat) / N^2
  }
}

#' Sandwich variance for corrected-score estimator (multicategory)
#' @keywords internal
vcov_cs_multi <- function(psi, y, z_hat, x, K, family, Pi, pi_z,
                          wt = NULL) {
  fam <- get_link_funs(family)
  n   <- length(y)
  N   <- if (is.null(wt)) n else sum(wt)
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
  if (is.null(wt)) {
    S <- crossprod(phi_mat) / N
  } else {
    S <- crossprod(phi_mat * wt, phi_mat) / N
  }

  I_hat <- compute_Ihat_multi(psi, z_hat, x, K, fam$mu_dot, wt = wt)
  M_hat <- compute_Mhat_multi(psi, x, K, fam$mu, Pi, pi_z, wt = wt)
  J     <- -(I_hat + M_hat)
  J_inv <- solve(J)

  J_inv %*% S %*% t(J_inv) / N
}
