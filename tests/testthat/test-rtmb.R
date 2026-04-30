# --------------------------------------------------------------------------
# Tests for RTMB engine: verify it matches TMB results
# --------------------------------------------------------------------------

skip_if_no_rtmb <- function() {
  skip_if_not_installed("RTMB")
}

skip_if_no_tmb_dll <- function() {
  skip_if_not(
    tryCatch({ TMB::config(DLL = "mcGLM"); TRUE }, error = function(e) FALSE),
    "TMB DLL 'mcGLM' not compiled"
  )
}

# --- GLM: TMB vs RTMB produce same coefficients ---

test_that("RTMB matches TMB for binary Poisson onestep", {
  skip_if_no_rtmb()
  skip_if_no_tmb_dll()

  set.seed(200)
  n <- 1000
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  eta <- 0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]
  y <- rpois(n, exp(eta))
  p01 <- 0.10; p10 <- 0.15; pi_z <- 0.4
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  fit_tmb  <- mcglm(y, z_hat, x, family = "poisson",
                     method = c("naive", "onestep"),
                     p01 = p01, p10 = p10, pi_z = pi_z,
                     engine = "tmb")
  fit_rtmb <- mcglm(y, z_hat, x, family = "poisson",
                     method = c("naive", "onestep"),
                     p01 = p01, p10 = p10, pi_z = pi_z,
                     engine = "rtmb")

  expect_equal(coef(fit_tmb, "onestep"), coef(fit_rtmb, "onestep"),
               tolerance = 1e-4)
  expect_equal(fit_tmb$loglik_onestep, fit_rtmb$loglik_onestep,
               tolerance = 1e-2)
})

test_that("RTMB matches TMB for binary Binomial onestep", {
  skip_if_no_rtmb()
  skip_if_no_tmb_dll()

  set.seed(201)
  n <- 1000
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  eta <- 0.5 * z - 0.3 * x[, 1] + 0.4 * x[, 2]
  y <- rbinom(n, 1, plogis(eta))
  p01 <- 0.08; p10 <- 0.12; pi_z <- 0.4
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  fit_tmb  <- mcglm(y, z_hat, x, family = "binomial",
                     method = c("naive", "onestep"),
                     p01 = p01, p10 = p10, pi_z = pi_z,
                     engine = "tmb")
  fit_rtmb <- mcglm(y, z_hat, x, family = "binomial",
                     method = c("naive", "onestep"),
                     p01 = p01, p10 = p10, pi_z = pi_z,
                     engine = "rtmb")

  expect_equal(coef(fit_tmb, "onestep"), coef(fit_rtmb, "onestep"),
               tolerance = 1e-4)
})

test_that("RTMB matches TMB for multicategory Z Poisson onestep", {
  skip_if_no_rtmb()
  skip_if_no_tmb_dll()

  set.seed(202)
  n <- 1000
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

  fit_tmb  <- mcglm(y, z_hat, x, family = "poisson",
                     method = c("naive", "onestep"),
                     Pi = Pi, pi_z = pi_z, engine = "tmb")
  fit_rtmb <- mcglm(y, z_hat, x, family = "poisson",
                     method = c("naive", "onestep"),
                     Pi = Pi, pi_z = pi_z, engine = "rtmb")

  expect_equal(coef(fit_tmb, "onestep"), coef(fit_rtmb, "onestep"),
               tolerance = 1e-4)
})

# --- Multinomial: TMB vs RTMB ---

test_that("RTMB matches TMB for multinomial onestep (binary Z)", {
  skip_if_no_rtmb()
  skip_if_no_tmb_dll()

  set.seed(203)
  n <- 1000; J <- 3; K <- 2
  z <- rbinom(n, 1, 0.4)
  x <- cbind(1, rnorm(n))
  eta <- matrix(0, n, J)
  eta[, 2] <- 0.8 * z + x %*% c(-0.5, 0.7)
  eta[, 3] <- -0.5 * z + x %*% c(0.3, -0.4)
  prob <- exp(eta) / rowSums(exp(eta))
  y <- integer(n)
  for (i in 1:n) y[i] <- sample(0:(J - 1), 1, prob = prob[i, ])
  p01 <- 0.10; p10 <- 0.15; pi_z <- 0.4
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  fit_tmb  <- mcglm(y, z_hat, x, family = "multinomial",
                     method = c("naive", "onestep"),
                     p01 = p01, p10 = p10, pi_z = pi_z,
                     engine = "tmb")
  fit_rtmb <- mcglm(y, z_hat, x, family = "multinomial",
                     method = c("naive", "onestep"),
                     p01 = p01, p10 = p10, pi_z = pi_z,
                     engine = "rtmb")

  expect_equal(coef(fit_tmb, "onestep"), coef(fit_rtmb, "onestep"),
               tolerance = 1e-4)
  expect_equal(fit_tmb$loglik_onestep, fit_rtmb$loglik_onestep,
               tolerance = 1e-2)
})

test_that("RTMB matches TMB for multinomial onestep (multicategory Z)", {
  skip_if_no_rtmb()
  skip_if_no_tmb_dll()

  set.seed(204)
  n <- 1500; J <- 3; K <- 3
  pi_z <- c(0.5, 0.3, 0.2)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  Pi <- matrix(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8),
               nrow = K, byrow = TRUE)
  z_hat <- integer(n)
  for (i in 1:n) z_hat[i] <- sample(0:(K - 1), 1, prob = Pi[, z[i] + 1])
  x <- cbind(1, rnorm(n))
  eta <- matrix(0, n, J)
  gamma_true <- matrix(c(0.8, -0.3, -0.5, 0.6), nrow = 2, byrow = TRUE)
  alpha_true <- matrix(c(-0.5, 0.7, 0.3, -0.4), nrow = 2, byrow = TRUE)
  for (jj in 1:2) {
    gvec <- c(0, gamma_true[jj, ])
    eta[, jj + 1] <- gvec[z + 1] + x %*% alpha_true[jj, ]
  }
  prob <- exp(eta) / rowSums(exp(eta))
  y <- integer(n)
  for (i in 1:n) y[i] <- sample(0:(J - 1), 1, prob = prob[i, ])

  fit_tmb  <- mcglm(y, z_hat, x, family = "multinomial",
                     method = c("naive", "onestep"),
                     Pi = Pi, pi_z = pi_z, engine = "tmb")
  fit_rtmb <- mcglm(y, z_hat, x, family = "multinomial",
                     method = c("naive", "onestep"),
                     Pi = Pi, pi_z = pi_z, engine = "rtmb")

  expect_equal(coef(fit_tmb, "onestep"), coef(fit_rtmb, "onestep"),
               tolerance = 1e-3)
})

# --- RTMB-only tests (don't need TMB DLL) ---

test_that("RTMB onestep works standalone for Poisson", {
  skip_if_no_rtmb()

  set.seed(205)
  n <- 500
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.5 * z_hat + x %*% c(-0.3, 0.4)))

  fit <- mcglm(y, z_hat, x, family = "poisson",
               method = c("naive", "onestep"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.5,
               engine = "rtmb")

  expect_true("onestep" %in% names(fit$coefficients))
  expect_true(all(is.finite(coef(fit, "onestep"))))
  expect_true(is.finite(fit$loglik_onestep))
})

test_that("RTMB onestep with freq_weights matches unweighted", {
  skip_if_no_rtmb()

  set.seed(206)
  n <- 300
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.5 * z_hat + x %*% c(-0.3, 0.4)))

  fit1 <- mcglm(y, z_hat, x, family = "poisson",
                method = c("naive", "onestep"),
                p01 = 0.1, p10 = 0.1, pi_z = 0.5,
                engine = "rtmb")
  fit2 <- mcglm(y, z_hat, x, family = "poisson",
                method = c("naive", "onestep"),
                p01 = 0.1, p10 = 0.1, pi_z = 0.5,
                engine = "rtmb", freq_weights = rep(1, n))

  expect_equal(coef(fit1, "onestep"), coef(fit2, "onestep"), tolerance = 1e-6)
})

