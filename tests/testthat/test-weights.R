# --------------------------------------------------------------------------
# Tests for frequency weights (freq_weights parameter)
# --------------------------------------------------------------------------

# Helper: generate binary data and return both original and
# a "collapsed" version with frequency weights
gen_weighted_binary <- function(seed = 77) {
  set.seed(seed)
  n <- 1000
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  psi0 <- c(0.8, -0.5, 0.7)
  eta <- psi0[1] * z + x %*% psi0[-1]
  y <- rpois(n, exp(eta))

  p01 <- 0.10; p10 <- 0.15; pi_z <- 0.4
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  list(y = y, z_hat = z_hat, x = x, p01 = p01, p10 = p10, pi_z = pi_z,
       psi0 = psi0, n = n)
}

# ---- Input validation ----

test_that("mcglm errors on wrong-length freq_weights", {
  d <- gen_weighted_binary()
  expect_error(
    mcglm(d$y, d$z_hat, d$x, family = "poisson", method = "naive",
          freq_weights = rep(1, 5)),
    "freq_weights must have length"
  )
})

test_that("mcglm errors on non-positive freq_weights", {
  d <- gen_weighted_binary()
  wt <- rep(1, d$n)
  wt[1] <- 0
  expect_error(
    mcglm(d$y, d$z_hat, d$x, family = "poisson", method = "naive",
          freq_weights = wt),
    "positive"
  )
  wt[1] <- -1
  expect_error(
    mcglm(d$y, d$z_hat, d$x, family = "poisson", method = "naive",
          freq_weights = wt),
    "positive"
  )
})

# ---- Unit weights = no weights ----

test_that("unit freq_weights produce identical results to no weights (binary)", {
  d <- gen_weighted_binary()
  fit_no_wt <- mcglm(d$y, d$z_hat, d$x, family = "poisson",
                      method = c("naive", "bca", "bcm", "cs"),
                      p01 = d$p01, p10 = d$p10, pi_z = d$pi_z)
  fit_unit  <- mcglm(d$y, d$z_hat, d$x, family = "poisson",
                      method = c("naive", "bca", "bcm", "cs"),
                      p01 = d$p01, p10 = d$p10, pi_z = d$pi_z,
                      freq_weights = rep(1, d$n))

  for (m in c("naive", "bca", "bcm", "cs")) {
    expect_equal(coef(fit_unit, m), coef(fit_no_wt, m), tolerance = 1e-10,
                 info = paste("method:", m))
  }
})

# ---- Duplicated data equivalence ----

test_that("freq_weights=2 on all obs equals fitting on doubled data (binary naive)", {
  set.seed(88)
  n <- 200
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.5 * z_hat + x %*% c(-0.3, 0.4)))

  # Fit with weights = 2
  fit_wt <- mcglm(y, z_hat, x, family = "poisson", method = "naive",
                   freq_weights = rep(2, n))

  # Fit on duplicated data
  y2 <- rep(y, 2)
  z_hat2 <- rep(z_hat, 2)
  x2 <- rbind(x, x)
  fit_dup <- mcglm(y2, z_hat2, x2, family = "poisson", method = "naive")

  expect_equal(coef(fit_wt, "naive"), coef(fit_dup, "naive"), tolerance = 1e-8)
})

test_that("freq_weights with duplication equivalence for BCA/BCM/CS (binary)", {
  set.seed(89)
  n <- 300
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  p01 <- 0.1; p10 <- 0.15; pi_z <- 0.4
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
  eta <- 0.8 * z + x %*% c(-0.5, 0.7)
  y <- rpois(n, exp(eta))

  # Double first 50 observations via weights
  wt <- rep(1, n)
  wt[1:50] <- 2

  fit_wt <- mcglm(y, z_hat, x, family = "poisson",
                   method = c("naive", "bca", "bcm", "cs"),
                   p01 = p01, p10 = p10, pi_z = pi_z,
                   freq_weights = wt)

  # Duplicate those rows in the data
  idx <- c(1:n, 1:50)
  fit_dup <- mcglm(y[idx], z_hat[idx], x[idx, ], family = "poisson",
                    method = c("naive", "bca", "bcm", "cs"),
                    p01 = p01, p10 = p10, pi_z = pi_z)

  for (m in c("naive", "bca", "bcm", "cs")) {
    expect_equal(unname(coef(fit_wt, m)), unname(coef(fit_dup, m)),
                 tolerance = 1e-6, info = paste("method:", m))
  }
})

# ---- Multicategory with weights ----

test_that("unit freq_weights produce identical results to no weights (multi)", {
  set.seed(90)
  n <- 500
  K <- 3
  pi_z <- c(0.5, 0.3, 0.2)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  Pi <- matrix(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8),
               nrow = K, byrow = TRUE)
  z_hat <- integer(n)
  for (i in seq_len(n)) z_hat[i] <- sample(0:(K - 1), 1, prob = Pi[, z[i] + 1])
  x <- cbind(1, rnorm(n))
  gamma <- c(0, 1.0, -0.9)
  alpha <- c(0.8, -0.7)
  eta <- as.numeric(x %*% alpha) + gamma[z + 1]
  y <- rpois(n, exp(eta))

  fit_no_wt <- mcglm(y, z_hat, x, family = "poisson",
                      method = c("naive", "bca", "bcm", "cs"),
                      Pi = Pi, pi_z = pi_z)
  fit_unit  <- mcglm(y, z_hat, x, family = "poisson",
                      method = c("naive", "bca", "bcm", "cs"),
                      Pi = Pi, pi_z = pi_z,
                      freq_weights = rep(1, n))

  for (m in c("naive", "bca", "bcm", "cs")) {
    expect_equal(coef(fit_unit, m), coef(fit_no_wt, m), tolerance = 1e-10,
                 info = paste("multi method:", m))
  }
})

test_that("freq_weights duplication equivalence for multicategory", {
  set.seed(91)
  n <- 300
  K <- 3
  pi_z <- c(0.5, 0.3, 0.2)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  Pi <- matrix(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8),
               nrow = K, byrow = TRUE)
  z_hat <- integer(n)
  for (i in seq_len(n)) z_hat[i] <- sample(0:(K - 1), 1, prob = Pi[, z[i] + 1])
  x <- cbind(1, rnorm(n))
  gamma <- c(0, 1.0, -0.9)
  alpha <- c(0.8, -0.7)
  eta <- as.numeric(x %*% alpha) + gamma[z + 1]
  y <- rpois(n, exp(eta))

  # Weight first 40 observations by 3
  wt <- rep(1, n)
  wt[1:40] <- 3

  fit_wt <- mcglm(y, z_hat, x, family = "poisson",
                   method = c("naive", "bca", "bcm", "cs"),
                   Pi = Pi, pi_z = pi_z, freq_weights = wt)

  idx <- c(1:n, rep(1:40, 2))
  fit_dup <- mcglm(y[idx], z_hat[idx], x[idx, ], family = "poisson",
                    method = c("naive", "bca", "bcm", "cs"),
                    Pi = Pi, pi_z = pi_z)

  for (m in c("naive", "bca", "bcm", "cs")) {
    expect_equal(unname(coef(fit_wt, m)), unname(coef(fit_dup, m)),
                 tolerance = 1e-5, info = paste("multi method:", m))
  }
})

# ---- Weighted helpers directly ----

test_that("compute_mhat_bin with unit weights matches unweighted", {
  set.seed(100)
  n <- 200
  x <- cbind(1, rnorm(n))
  psi <- c(0.8, -0.5, 0.7)
  fam <- mcGLM:::get_link_funs("poisson")

  m1 <- mcGLM:::compute_mhat_bin(psi, x, fam$mu, 0.1, 0.15, 0.4)
  m2 <- mcGLM:::compute_mhat_bin(psi, x, fam$mu, 0.1, 0.15, 0.4,
                                  wt = rep(1, n))
  expect_equal(m1, m2, tolerance = 1e-12)
})

test_that("compute_Ihat with unit weights matches unweighted", {
  set.seed(101)
  n <- 200
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  psi <- c(0.8, -0.5, 0.7)
  fam <- mcGLM:::get_link_funs("poisson")

  I1 <- mcGLM:::compute_Ihat(psi, z_hat, x, fam$mu_dot)
  I2 <- mcGLM:::compute_Ihat(psi, z_hat, x, fam$mu_dot, wt = rep(1, n))
  expect_equal(I1, I2, tolerance = 1e-12)
})

test_that("compute_Mhat_bin with unit weights matches unweighted", {
  set.seed(102)
  n <- 200
  x <- cbind(1, rnorm(n))
  psi <- c(0.8, -0.5, 0.7)
  fam <- mcGLM:::get_link_funs("poisson")

  M1 <- mcGLM:::compute_Mhat_bin(psi, x, fam$mu, 0.1, 0.15, 0.4)
  M2 <- mcGLM:::compute_Mhat_bin(psi, x, fam$mu, 0.1, 0.15, 0.4,
                                  wt = rep(1, n))
  expect_equal(M1, M2, tolerance = 1e-10)
})

# ---- Variance estimators with weights ----

test_that("vcov_naive with unit weights matches unweighted", {
  d <- gen_weighted_binary()
  naive <- mcGLM:::fit_naive_bin(d$y, d$z_hat, d$x, "poisson")
  psi <- naive$coefficients

  V1 <- mcGLM:::vcov_naive(psi, d$y, d$z_hat, d$x, "poisson")
  V2 <- mcGLM:::vcov_naive(psi, d$y, d$z_hat, d$x, "poisson",
                            wt = rep(1, d$n))
  expect_equal(V1, V2, tolerance = 1e-12)
})

test_that("vcov_cs_bin with unit weights matches unweighted", {
  d <- gen_weighted_binary()
  naive <- mcGLM:::fit_naive_bin(d$y, d$z_hat, d$x, "poisson")
  psi_cs <- mcGLM:::fit_cs_bin(naive$coefficients, d$y, d$z_hat, d$x,
                                "poisson", d$p01, d$p10, d$pi_z)

  V1 <- mcGLM:::vcov_cs_bin(psi_cs, d$y, d$z_hat, d$x, "poisson",
                             d$p01, d$p10, d$pi_z)
  V2 <- mcGLM:::vcov_cs_bin(psi_cs, d$y, d$z_hat, d$x, "poisson",
                             d$p01, d$p10, d$pi_z, wt = rep(1, d$n))
  expect_equal(V1, V2, tolerance = 1e-12)
})

test_that("vcov_bc_bin corrected with unit weights matches unweighted", {
  d <- gen_weighted_binary()
  naive <- mcGLM:::fit_naive_bin(d$y, d$z_hat, d$x, "poisson")
  psi_naive <- naive$coefficients
  psi_bca <- mcGLM:::fit_bca_bin(psi_naive, d$y, d$z_hat, d$x, "poisson",
                                  d$p01, d$p10, d$pi_z)

  V1 <- mcGLM:::vcov_bc_bin(psi_bca, d$y, d$z_hat, d$x, "poisson",
                             d$p01, d$p10, d$pi_z, psi_naive = psi_naive,
                             type = "bca", corrected = TRUE)
  V2 <- mcGLM:::vcov_bc_bin(psi_bca, d$y, d$z_hat, d$x, "poisson",
                             d$p01, d$p10, d$pi_z, psi_naive = psi_naive,
                             type = "bca", corrected = TRUE,
                             wt = rep(1, d$n))
  expect_equal(V1, V2, tolerance = 1e-10)
})

# ---- Weighted variance duplication check ----

test_that("vcov_naive with weights=2 matches duplicated data", {
  set.seed(103)
  n <- 200
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.5 * z_hat + x %*% c(-0.3, 0.4)))

  naive_wt <- mcGLM:::fit_naive_bin(y, z_hat, x, "poisson", wt = rep(2, n))
  V_wt <- mcGLM:::vcov_naive(naive_wt$coefficients, y, z_hat, x, "poisson",
                              wt = rep(2, n))

  y2 <- rep(y, 2)
  z2 <- rep(z_hat, 2)
  x2 <- rbind(x, x)
  naive_dup <- mcGLM:::fit_naive_bin(y2, z2, x2, "poisson")
  V_dup <- mcGLM:::vcov_naive(naive_dup$coefficients, y2, z2, x2, "poisson")

  expect_equal(V_wt, V_dup, tolerance = 1e-6)
})

# ---- Output stores freq_weights ----

test_that("mcglm output includes freq_weights", {
  d <- gen_weighted_binary()
  fit <- mcglm(d$y, d$z_hat, d$x, family = "poisson", method = "naive",
               freq_weights = rep(2, d$n))
  expect_equal(fit$freq_weights, rep(2, d$n))

  fit_no <- mcglm(d$y, d$z_hat, d$x, family = "poisson", method = "naive")
  expect_null(fit_no$freq_weights)
})

# ---- Non-integer weights ----

test_that("fractional freq_weights are accepted and produce finite results", {
  d <- gen_weighted_binary()
  wt <- runif(d$n, 0.5, 3.0)
  fit <- mcglm(d$y, d$z_hat, d$x, family = "poisson",
               method = c("naive", "bca", "bcm", "cs"),
               p01 = d$p01, p10 = d$p10, pi_z = d$pi_z,
               freq_weights = wt)

  for (m in c("naive", "bca", "bcm", "cs")) {
    expect_true(all(is.finite(coef(fit, m))), info = paste("method:", m))
  }
})

# ---- Binomial family with weights ----

test_that("freq_weights work with binomial family", {
  set.seed(104)
  n <- 500
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  psi0 <- c(0.8, -0.5, 0.7)
  eta <- psi0[1] * z + x %*% psi0[-1]
  p <- 1 / (1 + exp(-eta))
  y <- rbinom(n, 1, p)

  p01 <- 0.10; p10 <- 0.15; pi_z <- 0.4
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  # Unit weights = no weights
  fit1 <- mcglm(y, z_hat, x, family = "binomial",
                method = c("naive", "bca", "bcm", "cs"),
                p01 = p01, p10 = p10, pi_z = pi_z)
  fit2 <- mcglm(y, z_hat, x, family = "binomial",
                method = c("naive", "bca", "bcm", "cs"),
                p01 = p01, p10 = p10, pi_z = pi_z,
                freq_weights = rep(1, n))

  for (m in c("naive", "bca", "bcm", "cs")) {
    expect_equal(coef(fit1, m), coef(fit2, m), tolerance = 1e-10,
                 info = paste("binomial method:", m))
  }
})
