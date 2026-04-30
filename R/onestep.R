# ---------------------------------------------------------------------------
# One-step joint estimation via TMB mixture likelihood
# ---------------------------------------------------------------------------

#' Compute mixture weights from known misclassification rates (binary)
#'
#' Returns a length-4 vector: omega[j*2 + l] = P(z_hat=j | Z=l) * P(Z=l)
#' for j in {0,1}, l in {0,1}.
#' @keywords internal
compute_omega_bin <- function(p01, p10, pi_z) {
  # Pi[j,l] = P(z_hat=j | Z=l)
  # Pi[0,0] = 1-p01,  Pi[1,0] = p01
  # Pi[0,1] = p10,    Pi[1,1] = 1-p10
  c(
    (1 - p01) * (1 - pi_z),  # omega[0,0] = P(z_hat=0|Z=0)*P(Z=0)
    p10 * pi_z,               # omega[0,1] = P(z_hat=0|Z=1)*P(Z=1)
    p01 * (1 - pi_z),         # omega[1,0] = P(z_hat=1|Z=0)*P(Z=0)
    (1 - p10) * pi_z          # omega[1,1] = P(z_hat=1|Z=1)*P(Z=1)
  )
}

#' Compute mixture weights from known misclassification matrix (multicategory)
#'
#' Returns a length K*K vector: omega[j*K + l] = Pi[j,l] * pi_z[l]
#' @keywords internal
compute_omega_multi <- function(Pi, pi_z, K) {
  omega <- numeric(K * K)
  for (j in seq_len(K)) {
    for (l in seq_len(K)) {
      omega[(j - 1) * K + l] <- Pi[j, l] * pi_z[l]
    }
  }
  omega
}

#' One-step joint estimator for binary misclassification
#'
#' Maximizes the integrated likelihood that marginalizes over the latent true
#' label via a mixture, using TMB's automatic differentiation.
#'
#' @param y Numeric response vector.
#' @param z_hat Integer proxy covariate (0/1).
#' @param x Numeric matrix of correctly observed covariates (n x r), with intercept.
#' @param family Family string or object.
#' @param p01 False positive rate (required if weights = "fixed").
#' @param p10 False negative rate (required if weights = "fixed").
#' @param pi_z Prevalence P(Z=1) (required if weights = "fixed").
#' @param weights Character: "fixed" uses known rates, "estimated" estimates via softmax.
#' @param homoskedastic Logical, for Gaussian only.
#' @param optim_control List of control parameters for nlminb.
#' @param wt Optional frequency weights vector (length n).
#' @return List with coefficients, vcov, loglik, convergence.
#' @keywords internal
fit_onestep_bin <- function(y, z_hat, x, family,
                            p01 = NULL, p10 = NULL, pi_z = NULL,
                            weights = "fixed", homoskedastic = TRUE,
                            optim_control = list(), wt = NULL) {
  fam <- get_link_funs(family)
  n   <- length(y)
  r   <- ncol(x)
  K   <- 2L
  s   <- 1L
  d   <- s + r  # gamma, alpha_0, ..., alpha_{r-1}

  # Determine distribution code
  dist_code <- switch(fam$family$family,
    gaussian = 1L,
    poisson  = 2L,
    binomial = 3L,
    stop("Unsupported family for one-step: ", fam$family$family)
  )

  # Compute omega
  weights_fixed <- as.integer(weights == "fixed")
  if (weights_fixed) {
    stopifnot(!is.null(p01), !is.null(p10), !is.null(pi_z))
    omega_data <- compute_omega_bin(p01, p10, pi_z)
  } else {
    omega_data <- rep(0.25, 4)  # placeholder
  }

  # Frequency weights
  if (is.null(wt)) wt_data <- rep(1.0, n) else wt_data <- as.numeric(wt)

  # Starting values from naive GLM
  xi_hat <- cbind(z_hat, x)
  naive_fit <- stats::glm(y ~ . - 1,
    data = data.frame(y = y, xi_hat),
    family = fam$family,
    weights = wt_data
  )
  b_init <- unname(stats::coef(naive_fit))

  # Build theta vector
  theta_init <- b_init
  if (!weights_fixed) {
    # K*K - 1 = 3 softmax parameters; initialize near uniform
    theta_init <- c(theta_init, rep(0, K * K - 1))
  }
  if (dist_code == 1L) {
    # Gaussian: log(sigma)
    resid_sd <- stats::sd(stats::residuals(naive_fit))
    if (homoskedastic) {
      theta_init <- c(theta_init, log(resid_sd))
    } else {
      theta_init <- c(theta_init, log(resid_sd), log(resid_sd))
    }
  }

  # TMB data list
  data_list <- list(
    Y             = y,
    X             = x,
    z_hat         = as.integer(z_hat),
    dist_code     = dist_code,
    weights_fixed = weights_fixed,
    omega_data    = omega_data,
    homoskedastic = as.integer(homoskedastic),
    K             = K,
    wt            = wt_data
  )

  obj <- TMB::MakeADFun(
    data       = data_list,
    parameters = list(theta = theta_init),
    DLL        = "mcGLM",
    silent     = TRUE
  )

  # Optimize
  ctrl <- list(iter.max = 500, abs.tol = 1e-12, rel.tol = 1e-10)
  ctrl[names(optim_control)] <- optim_control

  opt <- tryCatch({
    stats::nlminb(
      start     = obj$par,
      objective = obj$fn,
      gradient  = obj$gr,
      control   = ctrl
    )
  }, error = function(e) {
    warning("nlminb failed, trying BFGS fallback")
    stats::optim(
      par    = obj$par,
      fn     = obj$fn,
      gr     = obj$gr,
      method = "BFGS",
      control = list(maxit = 500, reltol = 1e-10)
    )
  })

  if (!is.null(opt$convergence) && opt$convergence != 0) {
    warning("One-step optimization may not have converged (code: ", opt$convergence, ")")
  }

  # Extract coefficients and variance
  b_hat <- opt$par[1:d]
  V     <- vcov_onestep(obj, opt, d)

  nms <- c("gamma", paste0("alpha", seq_len(r) - 1))
  names(b_hat)     <- nms
  colnames(V)      <- nms
  rownames(V)      <- nms

  list(
    coefficients = b_hat,
    vcov         = V,
    loglik       = -opt$objective,
    convergence  = opt$convergence
  )
}


#' One-step joint estimator for multicategory misclassification
#' @keywords internal
fit_onestep_multi <- function(y, z_hat, x, K, family,
                              Pi = NULL, pi_z = NULL,
                              weights = "fixed", homoskedastic = TRUE,
                              optim_control = list(), wt = NULL) {
  fam <- get_link_funs(family)
  n   <- length(y)
  r   <- ncol(x)
  s   <- K - 1
  d   <- s + r

  dist_code <- switch(fam$family$family,
    gaussian = 1L,
    poisson  = 2L,
    binomial = 3L,
    stop("Unsupported family for one-step: ", fam$family$family)
  )

  weights_fixed <- as.integer(weights == "fixed")
  if (weights_fixed) {
    stopifnot(!is.null(Pi), !is.null(pi_z))
    omega_data <- compute_omega_multi(Pi, pi_z, K)
  } else {
    omega_data <- rep(1 / (K * K), K * K)
  }

  # Frequency weights
  if (is.null(wt)) wt_data <- rep(1.0, n) else wt_data <- as.numeric(wt)

  # Starting values from naive GLM
  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) d_hat[, k] <- as.numeric(z_hat == k)
  xi_hat <- cbind(d_hat, x)
  naive_fit <- stats::glm(y ~ . - 1,
    data = data.frame(y = y, xi_hat),
    family = fam$family,
    weights = wt_data
  )
  b_init <- unname(stats::coef(naive_fit))

  theta_init <- b_init
  if (!weights_fixed) {
    theta_init <- c(theta_init, rep(0, K * K - 1))
  }
  if (dist_code == 1L) {
    resid_sd <- stats::sd(stats::residuals(naive_fit))
    if (homoskedastic) {
      theta_init <- c(theta_init, log(resid_sd))
    } else {
      theta_init <- c(theta_init, log(resid_sd), log(resid_sd))
    }
  }

  data_list <- list(
    Y             = y,
    X             = x,
    z_hat         = as.integer(z_hat),
    dist_code     = dist_code,
    weights_fixed = weights_fixed,
    omega_data    = omega_data,
    homoskedastic = as.integer(homoskedastic),
    K             = K,
    wt            = wt_data
  )

  obj <- TMB::MakeADFun(
    data       = data_list,
    parameters = list(theta = theta_init),
    DLL        = "mcGLM",
    silent     = TRUE
  )

  ctrl <- list(iter.max = 500, abs.tol = 1e-12, rel.tol = 1e-10)
  ctrl[names(optim_control)] <- optim_control

  opt <- tryCatch({
    stats::nlminb(
      start     = obj$par,
      objective = obj$fn,
      gradient  = obj$gr,
      control   = ctrl
    )
  }, error = function(e) {
    warning("nlminb failed, trying BFGS fallback")
    stats::optim(
      par    = obj$par,
      fn     = obj$fn,
      gr     = obj$gr,
      method = "BFGS",
      control = list(maxit = 500, reltol = 1e-10)
    )
  })

  if (!is.null(opt$convergence) && opt$convergence != 0) {
    warning("One-step optimization may not have converged (code: ", opt$convergence, ")")
  }

  b_hat <- opt$par[1:d]
  V     <- vcov_onestep(obj, opt, d)

  nms <- c(paste0("gamma", seq_len(s)), paste0("alpha", seq_len(r) - 1))
  names(b_hat)     <- nms
  colnames(V)      <- nms
  rownames(V)      <- nms

  list(
    coefficients = b_hat,
    vcov         = V,
    loglik       = -opt$objective,
    convergence  = opt$convergence
  )
}


#' Compute variance-covariance matrix from TMB Hessian
#'
#' Uses the inverse Hessian at the optimum, with fallbacks following MLBC's
#' pattern: cholesky -> solve -> generalized inverse -> diagonal.
#' @param obj TMB objective (from MakeADFun).
#' @param opt Optimization result (from nlminb/optim).
#' @param d Number of regression coefficients to extract.
#' @return d x d variance-covariance matrix.
#' @keywords internal
vcov_onestep <- function(obj, opt, d) {
  H <- tryCatch({
    obj$he(opt$par)
  }, error = function(e) {
    warning("Hessian computation failed, using numerical approximation")
    if (requireNamespace("numDeriv", quietly = TRUE)) {
      numDeriv::hessian(obj$fn, opt$par)
    } else {
      warning("numDeriv not available; returning diagonal variance")
      return(diag(1e-4, length(opt$par)))
    }
  })

  if (any(!is.finite(H))) {
    warning("Hessian contains non-finite values")
    H[!is.finite(H)] <- 0
  }

  # Regularize if not positive definite
  eig_vals <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  min_eig  <- min(eig_vals)
  if (min_eig <= 1e-12) {
    warning("Hessian not positive definite (min eigenvalue: ", min_eig,
            "), adding regularization")
    H <- H + diag(abs(min_eig) + 1e-8, nrow(H))
  }

  # Invert with fallbacks
  V_full <- tryCatch({
    chol_H <- chol(H)
    chol2inv(chol_H)
  }, error = function(e) {
    tryCatch({
      solve(H)
    }, error = function(e2) {
      if (requireNamespace("MASS", quietly = TRUE)) {
        MASS::ginv(H)
      } else {
        warning("All covariance computations failed, returning diagonal")
        diag(1e-4, nrow(H))
      }
    })
  })

  V <- V_full[1:d, 1:d, drop = FALSE]

  if (any(!is.finite(V))) {
    warning("Covariance matrix contains non-finite values")
    V[!is.finite(V)] <- 0
    diag(V) <- pmax(diag(V), 1e-8)
  }

  V
}
