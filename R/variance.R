# ---------------------------------------------------------------------------
# Sandwich variance estimators for each method
# ---------------------------------------------------------------------------

#' Sandwich variance for the naive estimator
#' @keywords internal
vcov_naive <- function(psi, y, xi_hat, family, wt = NULL) {
  fam    <- get_link_funs(family)
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
#' Avoids materializing the full n x p influence function matrix by computing
#' the sandwich components from cross-products of score and drift matrices.
#' @keywords internal
vcov_bc_bin <- function(psi_bc, y, xi_hat, x, family,
                        p01 = NULL, p10 = NULL, pi_z = NULL,
                        psi_naive = NULL,
                        type = c("bca", "bcm"), corrected = FALSE,
                        wt = NULL) {
  if (!corrected) return(vcov_naive(psi_bc, y, xi_hat, family, wt = wt))

  type   <- match.arg(type)
  fam    <- get_link_funs(family)
  n      <- length(y)
  N      <- if (is.null(wt)) n else sum(wt)
  p      <- length(psi_bc)

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
  M_hat <- compute_Mhat_bin(psi_naive, x, fam$mu, fam$mu_dot, p01, p10, pi_z,
                            wt = wt)
  m_mat <- compute_m_bin(psi_naive, x, fam$mu, p01, p10, pi_z)
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

  # Compute V = crossprod(IF_mat * sqrt(wt)) / N^2 without materializing IF_mat
  # IF_mat = score_mat %*% L1 - centered_m %*% L2
  # V = L1' S_ss L1 + L2' S_mm L2 - L1' S_sm L2 - L2' S_ms L1
  score_mat  <- xi_hat * resid                 # n x p
  centered_m <- sweep(m_mat, 2, m_bar)         # n x p
  L1 <- A_inv %*% t(G)
  L2 <- t(H_inv)

  if (is.null(wt)) {
    S_ss <- crossprod(score_mat) / N^2
    S_mm <- crossprod(centered_m) / N^2
    S_sm <- crossprod(score_mat, centered_m) / N^2
  } else {
    S_ss <- crossprod(score_mat * wt, score_mat) / N^2
    S_mm <- crossprod(centered_m * wt, centered_m) / N^2
    S_sm <- crossprod(score_mat * wt, centered_m) / N^2
  }

  t(L1) %*% S_ss %*% L1 + t(L2) %*% S_mm %*% L2 -
    t(L1) %*% S_sm %*% L2 - t(L2) %*% t(S_sm) %*% L1
}

#' Sandwich variance for corrected-score estimator (binary)
#' @keywords internal
vcov_cs_bin <- function(psi, y, xi_hat, x, family, p01, p10, pi_z,
                        wt = NULL) {
  fam    <- get_link_funs(family)
  n      <- length(y)
  N      <- if (is.null(wt)) n else sum(wt)

  eta_tilde <- as.numeric(xi_hat %*% psi)
  resid     <- y - fam$mu(eta_tilde)
  m_mat     <- compute_m_bin(psi, x, fam$mu, p01, p10, pi_z)

  phi_mat <- xi_hat * resid - m_mat

  if (is.null(wt)) {
    S <- crossprod(phi_mat) / N
  } else {
    S <- crossprod(phi_mat * wt, phi_mat) / N
  }

  I_hat <- compute_Ihat(psi, xi_hat, fam$mu_dot, wt = wt)
  M_hat <- compute_Mhat_bin(psi, x, fam$mu, fam$mu_dot, p01, p10, pi_z,
                            wt = wt)
  J     <- -(I_hat + M_hat)
  J_inv <- solve(J)

  J_inv %*% S %*% t(J_inv) / N
}


# ---- Multicategory variance estimators ----

#' Sandwich variance for naive estimator (multicategory)
#' @keywords internal
vcov_naive_multi <- function(psi, y, xi_hat, z_hat, x, K, family, wt = NULL) {
  fam <- get_link_funs(family)
  n   <- length(y)
  N   <- if (is.null(wt)) n else sum(wt)
  s   <- K - 1
  r   <- ncol(x)

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
vcov_bc_multi <- function(psi_bc, y, xi_hat, z_hat, x, K, family,
                          Pi = NULL, pi_z = NULL,
                          psi_naive = NULL,
                          type = c("bca", "bcm"), corrected = FALSE,
                          wt = NULL) {
  if (!corrected) return(vcov_naive_multi(psi_bc, y, xi_hat, z_hat, x, K,
                                           family, wt = wt))

  type <- match.arg(type)
  fam  <- get_link_funs(family)
  n    <- length(y)
  N    <- if (is.null(wt)) n else sum(wt)
  s    <- K - 1
  r    <- ncol(x)
  p    <- s + r

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

  if (is.null(wt)) {
    S_ss <- crossprod(score_mat) / N^2
    S_mm <- crossprod(centered_m) / N^2
    S_sm <- crossprod(score_mat, centered_m) / N^2
  } else {
    S_ss <- crossprod(score_mat * wt, score_mat) / N^2
    S_mm <- crossprod(centered_m * wt, centered_m) / N^2
    S_sm <- crossprod(score_mat * wt, centered_m) / N^2
  }

  t(L1) %*% S_ss %*% L1 + t(L2) %*% S_mm %*% L2 -
    t(L1) %*% S_sm %*% L2 - t(L2) %*% t(S_sm) %*% L1
}

#' Sandwich variance for corrected-score estimator (multicategory)
#' @keywords internal
vcov_cs_multi <- function(psi, y, xi_hat, z_hat, x, K, family, Pi, pi_z,
                          wt = NULL) {
  fam <- get_link_funs(family)
  n   <- length(y)
  N   <- if (is.null(wt)) n else sum(wt)
  s   <- K - 1
  r   <- ncol(x)

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

  I_hat <- compute_Ihat_multi(psi, xi_hat, z_hat, K, fam$mu_dot, wt = wt)
  M_hat <- compute_Mhat_multi(psi, x, K, fam$mu, Pi, pi_z, wt = wt)
  J     <- -(I_hat + M_hat)
  J_inv <- solve(J)

  J_inv %*% S %*% t(J_inv) / N
}
