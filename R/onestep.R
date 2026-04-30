# ---------------------------------------------------------------------------
# One-step joint estimation via TMB mixture likelihood
# ---------------------------------------------------------------------------

#' Compute mixture weights from known misclassification rates (binary)
#'
#' Returns a length-4 vector: omega[j*2 + l] = P(z_hat=j | Z=l) * P(Z=l)
#' for j in {0,1}, l in {0,1}.
#' @keywords internal
compute_omega_bin <- function(p01, p10, pi_z) {
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

#' One-step joint estimator for binary misclassification (GLM)
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
  d   <- s + r

  dist_code <- switch(fam$family$family,
    gaussian = 1L, poisson = 2L, binomial = 3L,
    stop("Unsupported family for one-step: ", fam$family$family)
  )

  weights_fixed <- as.integer(weights == "fixed")
  if (weights_fixed) {
    stopifnot(!is.null(p01), !is.null(p10), !is.null(pi_z))
    omega_data <- compute_omega_bin(p01, p10, pi_z)
  } else {
    omega_data <- rep(0.25, 4)
  }

  if (is.null(wt)) wt_data <- rep(1.0, n) else wt_data <- as.numeric(wt)

  xi_hat <- cbind(z_hat, x)
  naive_fit <- stats::glm(y ~ . - 1,
    data = data.frame(y = y, xi_hat),
    family = fam$family, weights = wt_data
  )
  b_init <- unname(stats::coef(naive_fit))

  theta_init <- b_init
  if (!weights_fixed) theta_init <- c(theta_init, rep(0, K * K - 1))
  if (dist_code == 1L) {
    resid_sd <- stats::sd(stats::residuals(naive_fit))
    if (homoskedastic) {
      theta_init <- c(theta_init, log(resid_sd))
    } else {
      theta_init <- c(theta_init, log(resid_sd), log(resid_sd))
    }
  }

  data_list <- list(
    model_type    = 0L,
    X             = x,
    z_hat         = as.integer(z_hat),
    weights_fixed = weights_fixed,
    omega_data    = omega_data,
    K             = K,
    wt            = wt_data,
    Y             = y,
    dist_code     = dist_code,
    homoskedastic = as.integer(homoskedastic)
  )

  obj <- TMB::MakeADFun(
    data = data_list, parameters = list(theta = theta_init),
    DLL = "mcGLM", silent = TRUE
  )

  opt <- run_nlminb(obj, optim_control)

  b_hat <- opt$par[1:d]
  V     <- vcov_onestep(obj, opt, d)
  nms   <- c("gamma", paste0("alpha", seq_len(r) - 1))
  names(b_hat) <- nms
  colnames(V) <- rownames(V) <- nms

  list(coefficients = b_hat, vcov = V,
       loglik = -opt$objective, convergence = opt$convergence)
}


#' One-step joint estimator for multicategory misclassification (GLM)
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
    gaussian = 1L, poisson = 2L, binomial = 3L,
    stop("Unsupported family for one-step: ", fam$family$family)
  )

  weights_fixed <- as.integer(weights == "fixed")
  if (weights_fixed) {
    stopifnot(!is.null(Pi), !is.null(pi_z))
    omega_data <- compute_omega_multi(Pi, pi_z, K)
  } else {
    omega_data <- rep(1 / (K * K), K * K)
  }

  if (is.null(wt)) wt_data <- rep(1.0, n) else wt_data <- as.numeric(wt)

  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) d_hat[, k] <- as.numeric(z_hat == k)
  xi_hat <- cbind(d_hat, x)
  naive_fit <- stats::glm(y ~ . - 1,
    data = data.frame(y = y, xi_hat),
    family = fam$family, weights = wt_data
  )
  b_init <- unname(stats::coef(naive_fit))

  theta_init <- b_init
  if (!weights_fixed) theta_init <- c(theta_init, rep(0, K * K - 1))
  if (dist_code == 1L) {
    resid_sd <- stats::sd(stats::residuals(naive_fit))
    if (homoskedastic) {
      theta_init <- c(theta_init, log(resid_sd))
    } else {
      theta_init <- c(theta_init, log(resid_sd), log(resid_sd))
    }
  }

  data_list <- list(
    model_type    = 0L,
    X             = x,
    z_hat         = as.integer(z_hat),
    weights_fixed = weights_fixed,
    omega_data    = omega_data,
    K             = K,
    wt            = wt_data,
    Y             = y,
    dist_code     = dist_code,
    homoskedastic = as.integer(homoskedastic)
  )

  obj <- TMB::MakeADFun(
    data = data_list, parameters = list(theta = theta_init),
    DLL = "mcGLM", silent = TRUE
  )

  opt <- run_nlminb(obj, optim_control)

  b_hat <- opt$par[1:d]
  V     <- vcov_onestep(obj, opt, d)
  nms   <- c(paste0("gamma", seq_len(s)), paste0("alpha", seq_len(r) - 1))
  names(b_hat) <- nms
  colnames(V) <- rownames(V) <- nms

  list(coefficients = b_hat, vcov = V,
       loglik = -opt$objective, convergence = opt$convergence)
}


# ======================== MULTINOMIAL ONE-STEP =============================

#' One-step joint estimator for multinomial logistic with misclassified covariate
#'
#' Maximizes the integrated multinomial likelihood that marginalizes over the
#' latent true label Z via a mixture.
#'
#' The model is:
#' \deqn{P(Y=j | Z=l, x) = \exp(\eta_{j,l}) / \sum_{j'} \exp(\eta_{j',l})}
#' where \eqn{\eta_{j,l} = \gamma_{j,l} + \alpha_j' x} with \eqn{\gamma_{j,0}=0}
#' (baseline Z) and \eqn{\eta_{0,l}=0} (baseline Y).
#'
#' @param y Integer vector of categorical responses in \{0, ..., J-1\}.
#' @param z_hat Integer vector of observed proxy covariate in \{0, ..., K-1\}.
#' @param x Numeric matrix of correctly observed covariates (n x r), including
#'   intercept if desired.
#' @param J Number of response categories.
#' @param K Number of categories for the misclassified covariate.
#' @param p01,p10,pi_z Binary misclassification parameters (for K=2).
#' @param Pi K x K misclassification matrix (for K >= 2).
#' @param pi_z Prevalence vector (length K).
#' @param weights Character: \code{"fixed"} or \code{"estimated"}.
#' @param optim_control List of control parameters for \code{nlminb}.
#' @param wt Optional frequency weights (length n).
#' @return List with coefficients, vcov, loglik, convergence.
#' @keywords internal
fit_onestep_multinomial <- function(y, z_hat, x, J, K,
                                    p01 = NULL, p10 = NULL, pi_z = NULL,
                                    Pi = NULL,
                                    weights = "fixed",
                                    optim_control = list(), wt = NULL) {
  n <- length(y)
  r <- ncol(x)
  s <- K - 1
  block_size <- s + r
  d <- (J - 1) * block_size   # total regression parameters

  is_binary_z <- (K == 2L)

  weights_fixed <- as.integer(weights == "fixed")
  if (weights_fixed) {
    if (is_binary_z) {
      stopifnot(!is.null(p01), !is.null(p10), !is.null(pi_z))
      omega_data <- compute_omega_bin(p01, p10, pi_z)
    } else {
      stopifnot(!is.null(Pi), !is.null(pi_z))
      omega_data <- compute_omega_multi(Pi, pi_z, K)
    }
  } else {
    omega_data <- rep(1 / (K * K), K * K)
  }

  if (is.null(wt)) wt_data <- rep(1.0, n) else wt_data <- as.numeric(wt)

  # --- Starting values ---
  # Fit naive multinomial via category-specific binomial regressions
  theta_init <- numeric(d)

  # Dummy-encode z_hat (baseline = 0)
  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) d_hat[, k] <- as.numeric(z_hat == k)
  xi_naive <- cbind(d_hat, x)

  for (jj in seq_len(J - 1)) {
    y_bin <- as.numeric(y == jj)
    tryCatch({
      fit_j <- stats::glm(y_bin ~ . - 1,
        data = data.frame(y_bin = y_bin, xi_naive),
        family = stats::binomial(), weights = wt_data
      )
      theta_init[((jj - 1) * block_size + 1):(jj * block_size)] <-
        unname(stats::coef(fit_j))
    }, error = function(e) {
      # leave as zeros
    })
  }

  if (!weights_fixed) {
    theta_init <- c(theta_init, rep(0, K * K - 1))
  }

  data_list <- list(
    model_type    = 1L,
    X             = x,
    z_hat         = as.integer(z_hat),
    weights_fixed = weights_fixed,
    omega_data    = omega_data,
    K             = K,
    wt            = wt_data,
    Y_int         = as.integer(y),
    J             = as.integer(J)
  )

  obj <- TMB::MakeADFun(
    data = data_list, parameters = list(theta = theta_init),
    DLL = "mcGLM", silent = TRUE
  )

  opt <- run_nlminb(obj, optim_control)

  b_hat <- opt$par[1:d]
  V     <- vcov_onestep(obj, opt, d)

  # Parameter names: for each response category j=1,...,J-1
  nms <- character(d)
  for (jj in seq_len(J - 1)) {
    offset <- (jj - 1) * block_size
    if (s > 0) {
      if (s == 1) {
        nms[offset + 1] <- paste0("y", jj, ":gamma")
      } else {
        nms[offset + seq_len(s)] <- paste0("y", jj, ":gamma", seq_len(s))
      }
    }
    nms[offset + s + seq_len(r)] <- paste0("y", jj, ":alpha", seq_len(r) - 1)
  }
  names(b_hat) <- nms
  colnames(V) <- rownames(V) <- nms

  list(coefficients = b_hat, vcov = V,
       loglik = -opt$objective, convergence = opt$convergence)
}


#' Fit naive multinomial logistic regression with proxy covariate
#'
#' Fits a standard multinomial logistic model treating the proxy z_hat as
#' the true covariate, using nnet::multinom if available, otherwise
#' category-specific binomial GLMs.
#'
#' @param y Integer response vector in \{0,...,J-1\}.
#' @param z_hat Integer proxy covariate vector.
#' @param x Covariate matrix (n x r).
#' @param J Number of response categories.
#' @param K Number of Z categories.
#' @param wt Optional frequency weights.
#' @return List with coefficients vector (ordered by response category blocks).
#' @keywords internal
fit_naive_multinomial <- function(y, z_hat, x, J, K, wt = NULL) {
  n <- length(y)
  r <- ncol(x)
  s <- K - 1
  block_size <- s + r
  d <- (J - 1) * block_size

  d_hat <- matrix(0, n, s)
  for (k in seq_len(s)) d_hat[, k] <- as.numeric(z_hat == k)
  xi <- cbind(d_hat, x)

  if (is.null(wt)) wt_data <- rep(1, n) else wt_data <- wt

  # Try nnet::multinom first
  if (requireNamespace("nnet", quietly = TRUE)) {
    dat <- data.frame(y = factor(y, levels = 0:(J - 1)), xi)
    fit <- nnet::multinom(y ~ . - 1, data = dat, weights = wt_data, trace = FALSE)
    # nnet returns a (J-1) x p matrix of coefficients
    coef_mat <- stats::coef(fit)
    if (!is.matrix(coef_mat)) coef_mat <- matrix(coef_mat, nrow = 1)
    psi <- as.numeric(t(coef_mat))  # stack by row = by response category
    return(list(coefficients = psi, multinom_fit = fit))
  }

  # Fallback: category-specific binomial GLMs
  psi <- numeric(d)
  for (jj in seq_len(J - 1)) {
    y_bin <- as.numeric(y == jj)
    tryCatch({
      fit_j <- stats::glm(y_bin ~ . - 1,
        data = data.frame(y_bin = y_bin, xi),
        family = stats::binomial(), weights = wt_data
      )
      psi[((jj - 1) * block_size + 1):(jj * block_size)] <-
        unname(stats::coef(fit_j))
    }, error = function(e) NULL)
  }
  list(coefficients = psi, multinom_fit = NULL)
}


# ======================== SHARED UTILITIES ================================

#' Run nlminb with BFGS fallback
#' @keywords internal
run_nlminb <- function(obj, optim_control = list()) {
  ctrl <- list(iter.max = 500, abs.tol = 1e-12, rel.tol = 1e-10)
  ctrl[names(optim_control)] <- optim_control

  opt <- tryCatch({
    stats::nlminb(
      start = obj$par, objective = obj$fn, gradient = obj$gr,
      control = ctrl
    )
  }, error = function(e) {
    warning("nlminb failed, trying BFGS fallback")
    stats::optim(
      par = obj$par, fn = obj$fn, gr = obj$gr,
      method = "BFGS",
      control = list(maxit = 500, reltol = 1e-10)
    )
  })

  if (!is.null(opt$convergence) && opt$convergence != 0) {
    warning("Optimization may not have converged (code: ", opt$convergence, ")")
  }
  opt
}


#' Compute variance-covariance matrix from TMB Hessian
#'
#' Uses the inverse Hessian at the optimum, with fallbacks.
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

  eig_vals <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  min_eig  <- min(eig_vals)
  if (min_eig <= 1e-12) {
    warning("Hessian not positive definite (min eigenvalue: ", min_eig,
            "), adding regularization")
    H <- H + diag(abs(min_eig) + 1e-8, nrow(H))
  }

  V_full <- tryCatch({
    chol2inv(chol(H))
  }, error = function(e) {
    tryCatch(solve(H), error = function(e2) {
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
