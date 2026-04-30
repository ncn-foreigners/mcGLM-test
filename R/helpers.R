# ---------------------------------------------------------------------------
# Internal helper functions for computing drift, Fisher info, Jacobian
# ---------------------------------------------------------------------------

# ---- Binary misclassification helpers ----

#' Compute the toggle gap delta_i(psi) for binary misclassification
#' delta_i = mu(gamma + alpha'x_i) - mu(alpha'x_i)
#' @keywords internal
compute_delta <- function(psi, x, mu_fun) {
  gamma <- psi[1]
  alpha <- psi[-1]
  eta0 <- as.numeric(x %*% alpha)          # alpha'x_i
  mu_fun(gamma + eta0) - mu_fun(eta0)      # n-vector
}

#' Compute per-observation drift m_hat(psi) for binary misclassification
#'
#' Returns n x p matrix where row i is m_i(psi).
#' m_i = ( -c1 * delta_i,  -c2 * delta_i * x_i )
#' @keywords internal
compute_m_bin <- function(psi, x, mu_fun, p01, p10, pi_z) {
  n <- nrow(x)
  p <- length(psi)
  delta <- compute_delta(psi, x, mu_fun)   # n-vector
  c1 <- p01 * (1 - pi_z)
  c2 <- p01 * (1 - pi_z) - p10 * pi_z
  m <- matrix(0, n, p)
  m[, 1]  <- -c1 * delta
  m[, -1] <- -c2 * delta * x
  m
}

#' Compute m_hat(psi) = colMeans(m_i)
#' @keywords internal
compute_mhat_bin <- function(psi, x, mu_fun, p01, p10, pi_z) {
  colMeans(compute_m_bin(psi, x, mu_fun, p01, p10, pi_z))
}

#' Compute I_hat(psi) = (1/n) sum dot_mu(tilde_eta_i) * xi_hat * xi_hat'
#' @keywords internal
compute_Ihat <- function(psi, z_hat, x, mu_dot_fun) {
  xi_hat <- cbind(z_hat, x)
  eta_tilde <- as.numeric(xi_hat %*% psi)
  w <- mu_dot_fun(eta_tilde)               # n-vector of weights
  crossprod(xi_hat * w, xi_hat) / nrow(x)
}

#' Compute M_hat(psi) = d m_hat / d psi' for binary misclassification
#'
#' M_hat is the Jacobian of m_hat w.r.t. psi (p x p matrix).
#' Uses numerical differentiation.
#' @keywords internal
compute_Mhat_bin <- function(psi, x, mu_fun, p01, p10, pi_z) {
  p <- length(psi)
  M <- matrix(0, p, p)
  h <- 1e-7
  m0 <- compute_mhat_bin(psi, x, mu_fun, p01, p10, pi_z)
  for (j in seq_len(p)) {
    psi_h <- psi
    psi_h[j] <- psi_h[j] + h
    m1 <- compute_mhat_bin(psi_h, x, mu_fun, p01, p10, pi_z)
    M[, j] <- (m1 - m0) / h
  }
  M
}


# ---- Multicategory misclassification helpers ----

#' Compute per-observation drift m_i(psi) for multicategory misclassification
#'
#' @param psi Parameter vector (gamma_1,...,gamma_(K-1), alpha).
#' @param x  n x r covariate matrix.
#' @param K  Number of categories (including baseline 0).
#' @param mu_fun Inverse link function.
#' @param Pi K x K misclassification matrix (columns = true, rows = proxy).
#' @param pi_z K-vector of class prevalences.
#' @return n x p matrix of per-observation drift contributions.
#' @keywords internal
compute_m_multi <- function(psi, x, K, mu_fun, Pi, pi_z) {
  n  <- nrow(x)
  s  <- K - 1          # number of gamma parameters
  r  <- ncol(x)
  p  <- s + r

  gamma <- c(0, psi[seq_len(s)])           # gamma_0 = 0, gamma_1,...,gamma_{K-1}
  alpha <- psi[(s + 1):p]

  eta_base <- as.numeric(x %*% alpha)      # n-vector: alpha'x_i

  # mu_{i,ell} for each category ell = 0,...,K-1: n x K matrix
  mu_mat <- sapply(seq_len(K), function(ell) {
    mu_fun(eta_base + gamma[ell])
  })  # n x K

  m <- matrix(0, n, p)

  # gamma_k block (k = 1,...,K-1)
  for (k in seq_len(s)) {
    for (ell in seq_len(K)) {
      # pi_z[ell] * Pi[k+1, ell] * (mu_{i,ell} - mu_{i,k})
      # Note: k+1 because R is 1-indexed; proxy row k corresponds to category k
      # but category indices in Pi are 1:K mapping to 0:(K-1)
      coeff <- pi_z[ell] * Pi[k + 1, ell]
      m[, k] <- m[, k] + coeff * (mu_mat[, ell] - mu_mat[, k + 1])
    }
  }

  # alpha block
  for (ell in seq_len(K)) {
    for (j in seq_len(K)) {
      coeff <- pi_z[ell] * Pi[j, ell]
      diff_mu <- mu_mat[, ell] - mu_mat[, j]  # n-vector
      m[, (s + 1):p] <- m[, (s + 1):p] + coeff * diff_mu * x
    }
  }

  m
}

#' @keywords internal
compute_mhat_multi <- function(psi, x, K, mu_fun, Pi, pi_z) {
  colMeans(compute_m_multi(psi, x, K, mu_fun, Pi, pi_z))
}

#' Compute I_hat(psi) for multicategory case
#'
#' I_hat = (1/n) sum dot_mu(eta_i) * xi_hat * xi_hat'
#' @keywords internal
compute_Ihat_multi <- function(psi, z_hat, x, K, mu_dot_fun) {
  n <- nrow(x)
  s <- K - 1
  r <- ncol(x)
  gamma <- c(0, psi[seq_len(s)])
  alpha <- psi[(s + 1):(s + r)]
  eta_base <- as.numeric(x %*% alpha)

  # dummy encoding of z_hat
  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) {
    d_hat[, k] <- as.numeric(z_hat == k)
  }
  xi_hat <- cbind(d_hat, x)

  # linear predictor for each observation using its proxy category
  eta_tilde <- eta_base + gamma[z_hat + 1]  # +1 for R indexing
  w <- mu_dot_fun(eta_tilde)
  crossprod(xi_hat * w, xi_hat) / n
}

#' Numerical Jacobian M_hat for multicategory drift
#' @keywords internal
compute_Mhat_multi <- function(psi, x, K, mu_fun, Pi, pi_z) {
  p <- length(psi)
  M <- matrix(0, p, p)
  h <- 1e-7
  m0 <- compute_mhat_multi(psi, x, K, mu_fun, Pi, pi_z)
  for (j in seq_len(p)) {
    psi_h <- psi
    psi_h[j] <- psi_h[j] + h
    m1 <- compute_mhat_multi(psi_h, x, K, mu_fun, Pi, pi_z)
    M[, j] <- (m1 - m0) / h
  }
  M
}
