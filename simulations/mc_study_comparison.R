# ==========================================================================
# Comparative Monte Carlo simulation study
# mcGLM methods vs MC-SIMEX vs One-step vs Regression Calibration
#
# Methods:
#   Oracle  — GLM on true Z (upper bound)
#   Naive   — GLM on proxy Z_hat
#   BCA     — additive bias correction (1-step)
#   BCM     — multiplicative bias correction (1-step, corrected SE)
#   BCM-iter— iterated BCM (converges to CS)
#   CS      — corrected-score estimator
#   MC-SIMEX— Küchenhoff et al. (2006) via mcGLM::mcsimex() [C++ core]
#   Onestep — joint mixture-likelihood via TMB
#   RC      — regression calibration (Bayes' rule substitution)
#
# 20 scenarios: 14 binary + 6 multicategory
# B = 500 replications per scenario
#
# Usage:
#   Rscript simulations/mc_study_comparison.R              # sequential
#   Rscript simulations/mc_study_comparison.R --ncores=8   # parallel
#
# Install the package first:
#   remotes::install_github("ncn-foreigners/mcGLM-test")
# ==========================================================================

library(mcGLM)
library(parallel)

# ---- Parse command-line args ----
args <- commandArgs(trailingOnly = TRUE)
ncores_arg <- grep("^--ncores=", args, value = TRUE)
NCORES <- if (length(ncores_arg) > 0) {
  as.integer(sub("^--ncores=", "", ncores_arg))
} else {
  1L
}
cat(sprintf("Using %d core(s)\n", NCORES))

# ---- Import internal mcGLM functions via ::: ----
get_link_funs       <- mcGLM:::get_link_funs
build_xi_hat        <- mcGLM:::build_xi_hat
fit_naive_bin       <- mcGLM:::fit_naive_bin
fit_bca_bin         <- mcGLM:::fit_bca_bin
fit_bcm_bin         <- mcGLM:::fit_bcm_bin
fit_cs_bin          <- mcGLM:::fit_cs_bin
fit_naive_multi     <- mcGLM:::fit_naive_multi
fit_bca_multi       <- mcGLM:::fit_bca_multi
fit_bcm_multi       <- mcGLM:::fit_bcm_multi
fit_cs_multi        <- mcGLM:::fit_cs_multi
vcov_naive          <- mcGLM:::vcov_naive
vcov_bc_bin         <- mcGLM:::vcov_bc_bin
vcov_cs_bin         <- mcGLM:::vcov_cs_bin
vcov_naive_multi    <- mcGLM:::vcov_naive_multi
vcov_bc_multi       <- mcGLM:::vcov_bc_multi
vcov_cs_multi       <- mcGLM:::vcov_cs_multi
compute_mhat_bin    <- mcGLM:::compute_mhat_bin
compute_Mhat_bin    <- mcGLM:::compute_Mhat_bin
compute_mhat_multi  <- mcGLM:::compute_mhat_multi
compute_Mhat_multi  <- mcGLM:::compute_Mhat_multi

# ---- Regression calibration helpers ----

#' Regression calibration for binary misclassification
rc_fit_bin <- function(y, z_hat, x, family, p01, p10, pi_z) {
  prob_z1_given_1 <- (1 - p10) * pi_z / ((1 - p10) * pi_z + p01 * (1 - pi_z))
  prob_z1_given_0 <- p10 * pi_z / (p10 * pi_z + (1 - p01) * (1 - pi_z))
  z_cal <- ifelse(z_hat == 1, prob_z1_given_1, prob_z1_given_0)
  xi_cal <- cbind(z_cal, x)
  dat <- data.frame(y = y, xi_cal)
  fam <- get_link_funs(family)
  fit <- stats::glm(y ~ . - 1, data = dat, family = fam$family)
  list(coefficients = unname(stats::coef(fit)), glm_fit = fit)
}

#' Regression calibration for multicategory misclassification
rc_fit_multi <- function(y, z_hat, x, K, family, Pi, pi_z) {
  n <- length(y); s <- K - 1
  post_mat <- matrix(0, K, K)
  for (j in seq_len(K)) {
    numer <- Pi[j, ] * pi_z
    post_mat[j, ] <- numer / sum(numer)
  }
  d_cal <- matrix(0, n, s)
  for (i in seq_len(n)) {
    j <- z_hat[i] + 1
    for (k in seq_len(s)) d_cal[i, k] <- post_mat[j, k + 1]
  }
  xi_cal <- cbind(d_cal, x)
  dat <- data.frame(y = y, xi_cal)
  fam <- get_link_funs(family)
  fit <- stats::glm(y ~ . - 1, data = dat, family = fam$family)
  list(coefficients = unname(stats::coef(fit)), glm_fit = fit)
}


# ---- Oracle fit ----

oracle_fit_bin <- function(y, z, x, family) {
  fam <- get_link_funs(family)
  dat <- data.frame(y = y, cbind(z, x))
  fit <- stats::glm(y ~ . - 1, data = dat, family = fam$family)
  list(coefficients = unname(stats::coef(fit)),
       se = unname(summary(fit)$coefficients[, "Std. Error"]))
}

oracle_fit_multi <- function(y, z, x, K, family) {
  fam <- get_link_funs(family)
  xi <- build_xi_hat(z, x, K)
  dat <- data.frame(y = y, xi)
  fit <- stats::glm(y ~ . - 1, data = dat, family = fam$family)
  list(coefficients = unname(stats::coef(fit)),
       se = unname(summary(fit)$coefficients[, "Std. Error"]))
}


# ---- Single-replication engines ----

run_one_comparison_bin <- function(n, psi0, pi_z, p01, p10, family,
                                   regime = "known", n_v = 500,
                                   B_simex = 100) {
  p <- length(psi0)
  safe   <- function(expr) tryCatch(expr, error = function(e) rep(NA, p))
  safe_v <- function(expr) tryCatch(sqrt(pmax(diag(expr), 0)),
                                     error = function(e) rep(NA, p))
  safe_l <- function(expr) tryCatch(expr,
                                     error = function(e) list(
                                       coefficients = rep(NA, p),
                                       se = rep(NA, p)))

  # Generate data
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, pi_z)
  fam <- get_link_funs(family)
  eta <- psi0[1] * z + as.numeric(x %*% psi0[-1])
  mu_val <- fam$mu(eta)
  y <- if (family == "poisson") rpois(n, mu_val) else rbinom(n, 1, mu_val)

  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  # Misclassification params (known or estimated)
  if (regime == "validation") {
    val_idx <- sample(seq_len(n), min(n_v, n))
    tab <- table(factor(z_hat[val_idx], 0:1), factor(z[val_idx], 0:1))
    Pi_hat <- prop.table(tab, margin = 2)
    p01_use <- Pi_hat["1", "0"]
    p10_use <- Pi_hat["0", "1"]
    pi_z_use <- mean(z[val_idx])
  } else {
    p01_use <- p01; p10_use <- p10; pi_z_use <- pi_z
  }

  # Pre-compute xi_hat (used by all methods)
  xi_hat <- build_xi_hat(z_hat, x, 2L)

  # Pi matrix for mcsimex and onestep
  Pi_mat <- matrix(c(1 - p01_use, p01_use, p10_use, 1 - p10_use), 2, 2)

  # Oracle
  t0 <- proc.time()[3]
  oracle <- oracle_fit_bin(y, z, x, family)
  t_oracle <- proc.time()[3] - t0

  # Naive
  t0 <- proc.time()[3]
  naive_obj <- fit_naive_bin(y, xi_hat, family)
  psi_n <- naive_obj$coefficients
  se_n <- safe_v(vcov_naive(psi_n, y, xi_hat, family))
  t_naive <- proc.time()[3] - t0

  # BCA
  t0 <- proc.time()[3]
  psi_bca <- safe(fit_bca_bin(psi_n, y, xi_hat, x, family,
                               p01_use, p10_use, pi_z_use))
  se_bca <- safe_v(vcov_bc_bin(psi_bca, y, xi_hat, x, family,
                                p01_use, p10_use, pi_z_use, psi_naive = psi_n,
                                type = "bca", corrected = TRUE))
  t_bca <- proc.time()[3] - t0

  # BCM
  t0 <- proc.time()[3]
  psi_bcm <- safe(fit_bcm_bin(psi_n, y, xi_hat, x, family,
                               p01_use, p10_use, pi_z_use))
  se_bcm <- safe_v(vcov_bc_bin(psi_bcm, y, xi_hat, x, family,
                                p01_use, p10_use, pi_z_use, psi_naive = psi_n,
                                type = "bcm", corrected = TRUE))
  t_bcm <- proc.time()[3] - t0

  # BCM-iter
  t0 <- proc.time()[3]
  psi_bcmI <- safe(fit_bcm_bin(psi_n, y, xi_hat, x, family,
                                p01_use, p10_use, pi_z_use, iterate = TRUE))
  se_bcmI <- safe_v(vcov_cs_bin(psi_bcmI, y, xi_hat, x, family,
                                 p01_use, p10_use, pi_z_use))
  t_bcmI <- proc.time()[3] - t0

  # CS
  t0 <- proc.time()[3]
  psi_cs <- safe(fit_cs_bin(psi_n, y, xi_hat, x, family,
                             p01_use, p10_use, pi_z_use))
  se_cs <- safe_v(vcov_cs_bin(psi_cs, y, xi_hat, x, family,
                               p01_use, p10_use, pi_z_use))
  t_cs <- proc.time()[3] - t0

  # MC-SIMEX (C++ core via mcGLM::mcsimex)
  t0 <- proc.time()[3]
  simex_res <- safe_l({
    fit_sx <- mcsimex(y, z_hat, x, family = family, Pi = Pi_mat,
                      B = B_simex, jackknife = TRUE,
                      seed = sample.int(.Machine$integer.max, 1))
    se_sx <- if (!is.null(fit_sx$vcov)) sqrt(pmax(diag(fit_sx$vcov), 0))
             else rep(NA, p)
    list(coefficients = unname(fit_sx$coefficients), se = se_sx)
  })
  psi_simex <- simex_res$coefficients
  se_simex <- simex_res$se
  t_simex <- proc.time()[3] - t0

  # Onestep (TMB joint mixture likelihood)
  t0 <- proc.time()[3]
  onestep_res <- safe_l({
    fit_os <- mcGLM:::fit_onestep_bin(y, z_hat, x, family,
                                       p01_use, p10_use, pi_z_use,
                                       engine = "tmb")
    list(coefficients = unname(fit_os$coefficients),
         se = sqrt(pmax(diag(fit_os$vcov), 0)))
  })
  psi_onestep <- onestep_res$coefficients
  se_onestep <- onestep_res$se
  t_onestep <- proc.time()[3] - t0

  # Regression Calibration
  t0 <- proc.time()[3]
  rc_res <- safe_l(rc_fit_bin(y, z_hat, x, family, p01_use, p10_use, pi_z_use))
  psi_rc <- rc_res$coefficients
  se_rc <- rep(NA, p)
  t_rc <- proc.time()[3] - t0

  list(
    oracle = oracle$coefficients, naive = psi_n, bca = psi_bca,
    bcm = psi_bcm, bcm_iter = psi_bcmI, cs = psi_cs,
    mcsimex = psi_simex, onestep = psi_onestep, rc = psi_rc,
    se_oracle = oracle$se, se_naive = se_n, se_bca = se_bca,
    se_bcm = se_bcm, se_bcm_iter = se_bcmI, se_cs = se_cs,
    se_mcsimex = se_simex, se_onestep = se_onestep, se_rc = se_rc,
    t_oracle = t_oracle, t_naive = t_naive, t_bca = t_bca,
    t_bcm = t_bcm, t_bcm_iter = t_bcmI, t_cs = t_cs,
    t_mcsimex = t_simex, t_onestep = t_onestep, t_rc = t_rc
  )
}


run_one_comparison_multi <- function(n, psi0, K, pi_z, Pi, family,
                                      regime = "known", n_v = 2000,
                                      B_simex = 100) {
  p <- length(psi0)
  s <- K - 1; r <- 2
  safe   <- function(expr) tryCatch(expr, error = function(e) rep(NA, p))
  safe_v <- function(expr) tryCatch(sqrt(pmax(diag(expr), 0)),
                                     error = function(e) rep(NA, p))
  safe_l <- function(expr) tryCatch(expr,
                                     error = function(e) list(
                                       coefficients = rep(NA, p),
                                       se = rep(NA, p)))

  x <- cbind(1, rnorm(n))
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  fam <- get_link_funs(family)
  gamma_full <- c(0, psi0[seq_len(s)])
  alpha <- psi0[(s + 1):(s + r)]
  eta <- gamma_full[z + 1] + as.numeric(x %*% alpha)
  mu_val <- fam$mu(eta)
  y <- if (family == "poisson") rpois(n, mu_val) else rbinom(n, 1, mu_val)

  z_hat <- integer(n)
  for (i in seq_len(n)) z_hat[i] <- sample(0:(K - 1), 1, prob = Pi[, z[i] + 1])

  if (regime == "validation") {
    val_idx <- sample(seq_len(n), min(n_v, n))
    tab <- table(factor(z_hat[val_idx], 0:(K-1)), factor(z[val_idx], 0:(K-1)))
    Pi_use <- prop.table(as.matrix(tab), margin = 2)
    pi_z_use <- as.numeric(prop.table(table(factor(z[val_idx], 0:(K-1)))))
  } else {
    Pi_use <- Pi; pi_z_use <- pi_z
  }

  xi_hat <- build_xi_hat(z_hat, x, K)

  # Oracle
  t0 <- proc.time()[3]
  oracle <- oracle_fit_multi(y, z, x, K, family)
  t_oracle <- proc.time()[3] - t0

  # Naive
  t0 <- proc.time()[3]
  naive_obj <- fit_naive_multi(y, xi_hat, family)
  psi_n <- naive_obj$coefficients
  se_n <- safe_v(vcov_naive_multi(psi_n, y, xi_hat, z_hat, x, K, family))
  t_naive <- proc.time()[3] - t0

  # BCA
  t0 <- proc.time()[3]
  psi_bca <- safe(fit_bca_multi(psi_n, y, xi_hat, z_hat, x, K, family,
                                 Pi_use, pi_z_use))
  se_bca <- safe_v(vcov_bc_multi(psi_bca, y, xi_hat, z_hat, x, K, family,
                                  Pi_use, pi_z_use, psi_naive = psi_n,
                                  type = "bca", corrected = TRUE))
  t_bca <- proc.time()[3] - t0

  # BCM
  t0 <- proc.time()[3]
  psi_bcm <- safe(fit_bcm_multi(psi_n, y, xi_hat, z_hat, x, K, family,
                                 Pi_use, pi_z_use))
  se_bcm <- safe_v(vcov_bc_multi(psi_bcm, y, xi_hat, z_hat, x, K, family,
                                  Pi_use, pi_z_use, psi_naive = psi_n,
                                  type = "bcm", corrected = TRUE))
  t_bcm <- proc.time()[3] - t0

  # BCM-iter
  t0 <- proc.time()[3]
  psi_bcmI <- safe(fit_bcm_multi(psi_n, y, xi_hat, z_hat, x, K, family,
                                  Pi_use, pi_z_use, iterate = TRUE))
  se_bcmI <- safe_v(vcov_cs_multi(psi_bcmI, y, xi_hat, z_hat, x, K, family,
                                   Pi_use, pi_z_use))
  t_bcmI <- proc.time()[3] - t0

  # CS
  t0 <- proc.time()[3]
  psi_cs <- safe(fit_cs_multi(psi_n, y, xi_hat, z_hat, x, K, family,
                               Pi_use, pi_z_use))
  se_cs <- safe_v(vcov_cs_multi(psi_cs, y, xi_hat, z_hat, x, K, family,
                                 Pi_use, pi_z_use))
  t_cs <- proc.time()[3] - t0

  # MC-SIMEX (C++ core via mcGLM::mcsimex)
  t0 <- proc.time()[3]
  simex_res <- safe_l({
    fit_sx <- mcsimex(y, z_hat, x, family = family, Pi = Pi_use, K = K,
                      B = B_simex, jackknife = TRUE,
                      seed = sample.int(.Machine$integer.max, 1))
    se_sx <- if (!is.null(fit_sx$vcov)) sqrt(pmax(diag(fit_sx$vcov), 0))
             else rep(NA, p)
    list(coefficients = unname(fit_sx$coefficients), se = se_sx)
  })
  psi_simex <- simex_res$coefficients
  se_simex <- simex_res$se
  t_simex <- proc.time()[3] - t0

  # Onestep (TMB joint mixture likelihood)
  t0 <- proc.time()[3]
  onestep_res <- safe_l({
    fit_os <- mcGLM:::fit_onestep_multi(y, z_hat, x, K, family,
                                         Pi_use, pi_z_use, engine = "tmb")
    list(coefficients = unname(fit_os$coefficients),
         se = sqrt(pmax(diag(fit_os$vcov), 0)))
  })
  psi_onestep <- onestep_res$coefficients
  se_onestep <- onestep_res$se
  t_onestep <- proc.time()[3] - t0

  # RC
  t0 <- proc.time()[3]
  rc_res <- safe_l(rc_fit_multi(y, z_hat, x, K, family, Pi_use, pi_z_use))
  psi_rc <- rc_res$coefficients
  se_rc <- rep(NA, p)
  t_rc <- proc.time()[3] - t0

  list(
    oracle = oracle$coefficients, naive = psi_n, bca = psi_bca,
    bcm = psi_bcm, bcm_iter = psi_bcmI, cs = psi_cs,
    mcsimex = psi_simex, onestep = psi_onestep, rc = psi_rc,
    se_oracle = oracle$se, se_naive = se_n, se_bca = se_bca,
    se_bcm = se_bcm, se_bcm_iter = se_bcmI, se_cs = se_cs,
    se_mcsimex = se_simex, se_onestep = se_onestep, se_rc = se_rc,
    t_oracle = t_oracle, t_naive = t_naive, t_bca = t_bca,
    t_bcm = t_bcm, t_bcm_iter = t_bcmI, t_cs = t_cs,
    t_mcsimex = t_simex, t_onestep = t_onestep, t_rc = t_rc
  )
}


# ---- Summary function ----

summarize_comparison <- function(results, psi0, param_names) {
  methods <- c("oracle", "naive", "bca", "bcm", "bcm_iter", "cs",
               "mcsimex", "onestep", "rc")
  method_labels <- c("Oracle", "Naive", "BCA", "BCM", "BCM-iter", "CS",
                     "MC-SIMEX", "Onestep", "RC")
  p <- length(psi0)
  z_crit <- qnorm(0.975)
  out <- data.frame()

  for (m_idx in seq_along(methods)) {
    meth <- methods[m_idx]
    lab  <- method_labels[m_idx]
    est_mat <- do.call(rbind, lapply(results, `[[`, meth))
    se_mat  <- do.call(rbind, lapply(results, `[[`, paste0("se_", meth)))
    t_vec   <- sapply(results, `[[`, paste0("t_", meth))

    ok_est <- complete.cases(est_mat)
    ok_se  <- if (!all(is.na(se_mat))) complete.cases(se_mat) else rep(FALSE, nrow(est_mat))
    ok <- ok_est & ok_se
    B_ok <- sum(ok)
    B_est <- sum(ok_est)

    for (j in seq_len(p)) {
      est_j <- est_mat[ok_est, j]
      bias <- mean(est_j) - psi0[j]
      emp_sd <- sd(est_j)
      rmse <- sqrt(mean((est_j - psi0[j])^2))

      if (B_ok > 0) {
        se_j <- se_mat[ok, j]
        mean_se <- mean(se_j)
        lo <- est_mat[ok, j] - z_crit * se_j
        hi <- est_mat[ok, j] + z_crit * se_j
        coverage <- mean(lo <= psi0[j] & psi0[j] <= hi)
      } else {
        mean_se <- NA
        coverage <- NA
      }

      out <- rbind(out, data.frame(
        Method = lab, Parameter = param_names[j], True = psi0[j],
        Bias = round(bias, 5), Emp.SD = round(emp_sd, 5),
        Mean.SE = round(mean_se, 5), RMSE = round(rmse, 5),
        Coverage = round(coverage, 4),
        Med.Time = round(median(t_vec, na.rm = TRUE), 4),
        B = B_est, B.SE = B_ok,
        stringsAsFactors = FALSE
      ))
    }
  }
  out
}


# ==========================================================================
# Scenario definitions
# ==========================================================================

scenarios_bin <- list(
  B1  = list(family="poisson",  n=10000, p01=0.10, p10=0.15, pi_z=0.40, regime="known",
             label="B1: Poisson, n=10k, moderate miscl."),
  B2  = list(family="binomial", n=10000, p01=0.10, p10=0.15, pi_z=0.40, regime="known",
             label="B2: Logistic, n=10k, moderate miscl."),
  B3  = list(family="poisson",  n=500,   p01=0.10, p10=0.15, pi_z=0.40, regime="known",
             label="B3: Poisson, n=500, small sample"),
  B4  = list(family="poisson",  n=2000,  p01=0.10, p10=0.15, pi_z=0.40, regime="known",
             label="B4: Poisson, n=2k, moderate sample"),
  B5  = list(family="binomial", n=500,   p01=0.10, p10=0.15, pi_z=0.40, regime="known",
             label="B5: Logistic, n=500, small sample"),
  B6  = list(family="binomial", n=2000,  p01=0.10, p10=0.15, pi_z=0.40, regime="known",
             label="B6: Logistic, n=2k, moderate sample"),
  B7  = list(family="poisson",  n=10000, p01=0.30, p10=0.30, pi_z=0.40, regime="known",
             label="B7: Poisson, n=10k, heavy symmetric miscl."),
  B8  = list(family="binomial", n=10000, p01=0.30, p10=0.30, pi_z=0.40, regime="known",
             label="B8: Logistic, n=10k, heavy symmetric miscl."),
  B9  = list(family="poisson",  n=10000, p01=0.05, p10=0.05, pi_z=0.40, regime="known",
             label="B9: Poisson, n=10k, mild miscl."),
  B10 = list(family="poisson",  n=10000, p01=0.10, p10=0.15, pi_z=0.05, regime="known",
             label="B10: Poisson, n=10k, rare category"),
  B11 = list(family="poisson",  n=10000, p01=0.40, p10=0.10, pi_z=0.40, regime="known",
             label="B11: Poisson, n=10k, asymmetric miscl."),
  B12 = list(family="poisson",  n=10000, p01=0.10, p10=0.15, pi_z=0.40, regime="validation",
             label="B12: Poisson, n=10k, validation regime"),
  B13 = list(family="binomial", n=10000, p01=0.10, p10=0.15, pi_z=0.40, regime="validation",
             label="B13: Logistic, n=10k, validation regime"),
  B14 = list(family="poisson",  n=10000, p01=0.30, p10=0.30, pi_z=0.40, regime="validation",
             label="B14: Poisson, n=10k, validation + heavy miscl.")
)

psi0_bin  <- c(0.8, -0.5, 0.7)
names_bin <- c("gamma(0.8)", "alpha0(-0.5)", "alpha1(0.7)")

scenarios_multi <- list(
  M1 = list(family="poisson",  n=20000, p_miscl=0.20, regime="known",
            label="M1: Poisson, K=4, n=20k, 20% miscl."),
  M2 = list(family="binomial", n=20000, p_miscl=0.20, regime="known",
            label="M2: Logistic, K=4, n=20k, 20% miscl."),
  M3 = list(family="poisson",  n=5000,  p_miscl=0.20, regime="known",
            label="M3: Poisson, K=4, n=5k, small sample"),
  M4 = list(family="poisson",  n=20000, p_miscl=0.40, regime="known",
            label="M4: Poisson, K=4, n=20k, heavy 40% miscl."),
  M5 = list(family="poisson",  n=20000, p_miscl=0.10, regime="known",
            label="M5: Poisson, K=4, n=20k, mild 10% miscl."),
  M6 = list(family="poisson",  n=20000, p_miscl=0.20, regime="validation",
            label="M6: Poisson, K=4, n=20k, validation regime")
)

K_multi    <- 4
psi0_multi <- c(1.0, -0.9, 0.2, 0.8, -0.7)
pi_z_multi <- c(0.50, 0.25, 0.15, 0.10)
names_multi <- c("gamma1(1.0)", "gamma2(-0.9)", "gamma3(0.2)",
                  "alpha0(0.8)", "alpha1(-0.7)")

make_Pi <- function(K, p_miscl) {
  Pi <- matrix(p_miscl / (K - 1), K, K)
  diag(Pi) <- 1 - p_miscl
  Pi
}


# ==========================================================================
# Main simulation loop
# ==========================================================================

B <- 500
#B <- 10  # uncomment for quick testing
set.seed(2026)
all_tables <- list()

run_scenario <- function(run_fn, B, ncores, ...) {
  dots <- list(...)
  if (ncores > 1L) {
    seeds <- sample.int(.Machine$integer.max, B)
    mclapply(seq_len(B), function(b) {
      set.seed(seeds[b])
      do.call(run_fn, dots)
    }, mc.cores = ncores)
  } else {
    results <- vector("list", B)
    for (b in seq_len(B)) {
      if (b %% 100 == 0) cat(sprintf("  rep %d/%d\n", b, B))
      results[[b]] <- do.call(run_fn, dots)
    }
    results
  }
}

# ---- Binary scenarios ----
for (sc_name in names(scenarios_bin)) {
  sc <- scenarios_bin[[sc_name]]
  cat(sprintf("\n=== %s ===\n", sc$label))

  t0 <- proc.time()
  results <- run_scenario(
    run_one_comparison_bin, B = B, ncores = NCORES,
    n = sc$n, psi0 = psi0_bin, pi_z = sc$pi_z,
    p01 = sc$p01, p10 = sc$p10, family = sc$family,
    regime = sc$regime
  )
  t1 <- proc.time()
  cat(sprintf("  Elapsed: %.1f sec\n", (t1 - t0)[3]))

  tab <- summarize_comparison(results, psi0_bin, names_bin)
  tab <- cbind(Scenario = sc_name, tab)
  all_tables[[sc_name]] <- tab
  print(tab[tab$Parameter == "gamma(0.8)", ], row.names = FALSE)
}

# ---- Multicategory scenarios ----
for (sc_name in names(scenarios_multi)) {
  sc <- scenarios_multi[[sc_name]]
  cat(sprintf("\n=== %s ===\n", sc$label))

  Pi <- make_Pi(K_multi, sc$p_miscl)

  t0 <- proc.time()
  results <- run_scenario(
    run_one_comparison_multi, B = B, ncores = NCORES,
    n = sc$n, psi0 = psi0_multi, K = K_multi,
    pi_z = pi_z_multi, Pi = Pi, family = sc$family,
    regime = sc$regime
  )
  t1 <- proc.time()
  cat(sprintf("  Elapsed: %.1f sec\n", (t1 - t0)[3]))

  tab <- summarize_comparison(results, psi0_multi, names_multi)
  tab <- cbind(Scenario = sc_name, tab)
  all_tables[[sc_name]] <- tab
  print(tab[tab$Parameter == "gamma1(1.0)", ], row.names = FALSE)
}

# ---- Save ----
all_results <- do.call(rbind, all_tables)
write.csv(all_results, file = "mc_comparison_results.csv", row.names = FALSE)
cat("\nAll results saved to mc_comparison_results.csv\n")
