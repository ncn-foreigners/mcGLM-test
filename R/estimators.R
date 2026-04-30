# ---------------------------------------------------------------------------
# Estimators: naive, BCA, BCM, CS — binary and multicategory
# ---------------------------------------------------------------------------

# ========================== BINARY CASE ==================================

#' Fit the naive GLM estimator (binary misclassification)
#' @keywords internal
fit_naive_bin <- function(y, z_hat, x, family) {
  fam <- get_link_funs(family)
  xi_hat <- cbind(z_hat, x)
  dat <- data.frame(y = y, xi_hat)
  fit <- stats::glm(y ~ . - 1, data = dat, family = fam$family)
  list(
    coefficients = unname(stats::coef(fit)),
    fitted       = stats::fitted(fit),
    glm_fit      = fit
  )
}

#' Additive bias correction (BCA) for binary misclassification
#'
#' psi_bca = psi_hat - I_hat(psi_hat)^(-1) * m_hat(psi_hat)
#'
#' With \code{iterate = TRUE}, uses a fixed-point iteration
#' psi^(k+1) = psi_hat - I_hat(psi_hat)^(-1) * m_hat(psi^(k)),
#' keeping I_hat fixed at the naive estimate and only updating the drift
#' m_hat. This removes higher-order bias from evaluating m at the wrong point.
#' @keywords internal
fit_bca_bin <- function(psi_naive, y, z_hat, x, family, p01, p10, pi_z,
                        iterate = FALSE, max_iter = 50, tol = 1e-8) {
  fam <- get_link_funs(family)
  # I_hat fixed at the naive estimate (well-conditioned)
  I_hat_inv <- solve(compute_Ihat(psi_naive, z_hat, x, fam$mu_dot))

  m_hat <- compute_mhat_bin(psi_naive, x, fam$mu, p01, p10, pi_z)
  psi   <- psi_naive - I_hat_inv %*% m_hat
  if (!iterate) return(as.numeric(psi))

  for (iter in seq_len(max_iter)) {
    m_hat   <- compute_mhat_bin(psi, x, fam$mu, p01, p10, pi_z)
    psi_new <- psi_naive - I_hat_inv %*% m_hat
    if (max(abs(psi_new - psi)) < tol) break
    psi <- as.numeric(psi_new)
  }
  as.numeric(psi)
}

#' Multiplicative bias correction (BCM) for binary misclassification
#'
#' psi_bcm = psi_hat - (I_hat + M_hat)^(-1) * m_hat
#'
#' With \code{iterate = TRUE}, iterates Newton steps on the corrected-score
#' equation Phi_n(psi) = U_hat(psi) - m_hat(psi) = 0. The one-step BCM
#' is exactly one Newton step starting from psi_hat; iterating to convergence
#' yields the corrected-score (CS) estimator.
#' @keywords internal
fit_bcm_bin <- function(psi_naive, y, z_hat, x, family, p01, p10, pi_z,
                        iterate = FALSE, max_iter = 50, tol = 1e-8) {
  fam <- get_link_funs(family)
  xi_hat <- cbind(z_hat, x)
  n <- length(y)

  psi <- psi_naive
  for (iter in seq_len(if (iterate) max_iter else 1L)) {
    eta_tilde <- as.numeric(xi_hat %*% psi)
    resid <- y - fam$mu(eta_tilde)
    U_hat <- colMeans(xi_hat * resid)
    m_hat <- compute_mhat_bin(psi, x, fam$mu, p01, p10, pi_z)
    Phi   <- U_hat - m_hat

    I_hat <- compute_Ihat(psi, z_hat, x, fam$mu_dot)
    M_hat <- compute_Mhat_bin(psi, x, fam$mu, p01, p10, pi_z)
    step  <- solve(I_hat + M_hat, Phi)

    psi_new <- psi + step
    if (iterate && max(abs(psi_new - psi)) < tol) { psi <- psi_new; break }
    psi <- psi_new
  }
  psi
}

#' Corrected-score estimator for binary misclassification
#'
#' Solves (1/n) sum phi_i(psi) = 0, where
#' phi_i = xi_hat_i * (Y_i - mu(tilde_eta_i)) - m_i(psi)
#' @keywords internal
fit_cs_bin <- function(psi_init, y, z_hat, x, family, p01, p10, pi_z) {
  fam <- get_link_funs(family)
  n   <- length(y)
  xi_hat <- cbind(z_hat, x)

  # corrected score function: returns p-vector that should be zero
  phi_mean <- function(psi) {
    eta_tilde <- as.numeric(xi_hat %*% psi)
    resid     <- y - fam$mu(eta_tilde)
    score     <- colMeans(xi_hat * resid)
    m_hat     <- compute_mhat_bin(psi, x, fam$mu, p01, p10, pi_z)
    score - m_hat
  }

  # Jacobian of the corrected score
  phi_jac <- function(psi) {
    I_hat <- compute_Ihat(psi, z_hat, x, fam$mu_dot)
    M_hat <- compute_Mhat_bin(psi, x, fam$mu, p01, p10, pi_z)
    -(I_hat + M_hat)
  }

  sol <- nleqslv::nleqslv(psi_init, phi_mean, jac = phi_jac,
                           control = list(maxit = 500, ftol = 1e-12))
  if (sol$termcd > 2) {
    warning("CS solver did not converge (termcd = ", sol$termcd, ")")
  }
  sol$x
}


# ========================= MULTICATEGORY CASE ============================

#' Fit the naive GLM estimator (multicategory misclassification)
#' @keywords internal
fit_naive_multi <- function(y, z_hat, x, K, family) {
  fam <- get_link_funs(family)
  s <- K - 1
  # dummy-encode z_hat with baseline = 0
  d_hat <- matrix(0, length(y), s)
  for (k in seq_len(s)) {
    d_hat[, k] <- as.numeric(z_hat == k)
  }
  xi_hat <- cbind(d_hat, x)
  dat <- data.frame(y = y, xi_hat)
  fit <- stats::glm(y ~ . - 1, data = dat, family = fam$family)
  list(
    coefficients = unname(stats::coef(fit)),
    fitted       = stats::fitted(fit),
    glm_fit      = fit
  )
}

#' BCA for multicategory misclassification
#' @keywords internal
fit_bca_multi <- function(psi_naive, y, z_hat, x, K, family, Pi, pi_z,
                          iterate = FALSE, max_iter = 50, tol = 1e-8) {
  fam <- get_link_funs(family)
  I_hat_inv <- solve(compute_Ihat_multi(psi_naive, z_hat, x, K, fam$mu_dot))

  m_hat <- compute_mhat_multi(psi_naive, x, K, fam$mu, Pi, pi_z)
  psi   <- psi_naive - I_hat_inv %*% m_hat
  if (!iterate) return(as.numeric(psi))

  for (iter in seq_len(max_iter)) {
    m_hat   <- compute_mhat_multi(psi, x, K, fam$mu, Pi, pi_z)
    psi_new <- psi_naive - I_hat_inv %*% m_hat
    if (max(abs(psi_new - psi)) < tol) break
    psi <- as.numeric(psi_new)
  }
  as.numeric(psi)
}

#' BCM for multicategory misclassification
#' @keywords internal
fit_bcm_multi <- function(psi_naive, y, z_hat, x, K, family, Pi, pi_z,
                          iterate = FALSE, max_iter = 50, tol = 1e-8) {
  fam <- get_link_funs(family)
  n   <- length(y)
  s   <- K - 1
  r   <- ncol(x)

  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) d_hat[, k] <- as.numeric(z_hat == k)
  xi_hat <- cbind(d_hat, x)

  psi <- psi_naive
  for (iter in seq_len(if (iterate) max_iter else 1L)) {
    gamma <- c(0, psi[seq_len(s)])
    alpha <- psi[(s + 1):(s + r)]
    eta_base  <- as.numeric(x %*% alpha)
    eta_tilde <- eta_base + gamma[z_hat + 1]
    resid <- y - fam$mu(eta_tilde)
    U_hat <- colMeans(xi_hat * resid)
    m_hat <- compute_mhat_multi(psi, x, K, fam$mu, Pi, pi_z)
    Phi   <- U_hat - m_hat

    I_hat <- compute_Ihat_multi(psi, z_hat, x, K, fam$mu_dot)
    M_hat <- compute_Mhat_multi(psi, x, K, fam$mu, Pi, pi_z)
    step  <- solve(I_hat + M_hat, Phi)

    psi_new <- psi + step
    if (iterate && max(abs(psi_new - psi)) < tol) { psi <- psi_new; break }
    psi <- psi_new
  }
  psi
}

#' Corrected-score estimator for multicategory misclassification
#' @keywords internal
fit_cs_multi <- function(psi_init, y, z_hat, x, K, family, Pi, pi_z) {
  fam <- get_link_funs(family)
  n   <- length(y)
  s   <- K - 1
  r   <- ncol(x)

  # dummy-encode z_hat
  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) {
    d_hat[, k] <- as.numeric(z_hat == k)
  }
  xi_hat <- cbind(d_hat, x)

  phi_mean <- function(psi) {
    gamma <- c(0, psi[seq_len(s)])
    alpha <- psi[(s + 1):(s + r)]
    eta_base <- as.numeric(x %*% alpha)
    eta_tilde <- eta_base + gamma[z_hat + 1]
    resid <- y - fam$mu(eta_tilde)
    score <- colMeans(xi_hat * resid)
    m_hat <- compute_mhat_multi(psi, x, K, fam$mu, Pi, pi_z)
    score - m_hat
  }

  phi_jac <- function(psi) {
    I_hat <- compute_Ihat_multi(psi, z_hat, x, K, fam$mu_dot)
    M_hat <- compute_Mhat_multi(psi, x, K, fam$mu, Pi, pi_z)
    -(I_hat + M_hat)
  }

  sol <- nleqslv::nleqslv(psi_init, phi_mean, jac = phi_jac,
                           control = list(maxit = 500, ftol = 1e-12))
  if (sol$termcd > 2) {
    warning("CS solver did not converge (termcd = ", sol$termcd, ")")
  }
  sol$x
}
