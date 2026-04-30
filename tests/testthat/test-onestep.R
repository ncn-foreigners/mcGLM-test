# --------------------------------------------------------------------------
# Tests for one-step (joint estimation) and omega helpers
# --------------------------------------------------------------------------

# --- compute_omega_bin ---

test_that("compute_omega_bin sums to 1", {
  omega <- mcGLM:::compute_omega_bin(0.1, 0.15, 0.4)
  expect_equal(sum(omega), 1, tolerance = 1e-12)
})

test_that("compute_omega_bin has correct individual entries", {
  p01 <- 0.10; p10 <- 0.15; pi_z <- 0.4
  omega <- mcGLM:::compute_omega_bin(p01, p10, pi_z)

  expect_equal(omega[1], (1 - p01) * (1 - pi_z))  # P(z_hat=0|Z=0)*P(Z=0)
  expect_equal(omega[2], p10 * pi_z)                # P(z_hat=0|Z=1)*P(Z=1)
  expect_equal(omega[3], p01 * (1 - pi_z))          # P(z_hat=1|Z=0)*P(Z=0)
  expect_equal(omega[4], (1 - p10) * pi_z)          # P(z_hat=1|Z=1)*P(Z=1)
})

test_that("compute_omega_bin with zero error gives degenerate weights", {
  omega <- mcGLM:::compute_omega_bin(0, 0, 0.3)
  # No misclassification: omega[0,0] = P(Z=0), omega[1,1] = P(Z=1), rest 0
  expect_equal(omega[1], 0.7)   # P(Z=0)
  expect_equal(omega[2], 0)     # P(z_hat=0|Z=1)*P(Z=1)
  expect_equal(omega[3], 0)     # P(z_hat=1|Z=0)*P(Z=0)
  expect_equal(omega[4], 0.3)   # P(Z=1)
})

# --- compute_omega_multi ---

test_that("compute_omega_multi sums to 1", {
  K <- 3
  Pi <- matrix(c(0.8, 0.1, 0.1,
                 0.1, 0.8, 0.1,
                 0.1, 0.1, 0.8), nrow = K, byrow = TRUE)
  pi_z <- c(0.5, 0.3, 0.2)
  omega <- mcGLM:::compute_omega_multi(Pi, pi_z, K)

  expect_length(omega, K * K)
  expect_equal(sum(omega), 1, tolerance = 1e-12)
})

test_that("compute_omega_multi entries match Pi * pi_z", {
  K <- 3
  Pi <- matrix(c(0.7, 0.2, 0.1,
                 0.1, 0.7, 0.2,
                 0.2, 0.1, 0.7), nrow = K, byrow = TRUE)
  pi_z <- c(0.5, 0.3, 0.2)
  omega <- mcGLM:::compute_omega_multi(Pi, pi_z, K)

  for (j in 1:K) {
    for (l in 1:K) {
      idx <- (j - 1) * K + l
      expect_equal(omega[idx], Pi[j, l] * pi_z[l])
    }
  }
})

# --- fit_onestep_bin ---

test_that("fit_onestep_bin returns correct structure (poisson)", {
  skip_if_not(
    tryCatch({ TMB::config(DLL = "mcGLM"); TRUE }, error = function(e) FALSE),
    "TMB DLL 'mcGLM' not compiled"
  )

  set.seed(50)
  n <- 1000
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  eta <- 0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]
  y <- rpois(n, exp(eta))
  p01 <- 0.10; p10 <- 0.15; pi_z <- 0.4
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  res <- mcGLM:::fit_onestep_bin(y, z_hat, x, "poisson", p01, p10, pi_z)

  expect_named(res, c("coefficients", "vcov", "loglik", "convergence"))
  expect_length(res$coefficients, 3)
  expect_equal(dim(res$vcov), c(3, 3))
  expect_true(is.finite(res$loglik))
  expect_named(res$coefficients, c("gamma", "alpha0", "alpha1"))
})

test_that("fit_onestep_bin via mcglm interface", {
  skip_if_not(
    tryCatch({ TMB::config(DLL = "mcGLM"); TRUE }, error = function(e) FALSE),
    "TMB DLL 'mcGLM' not compiled"
  )

  set.seed(51)
  n <- 1000
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  eta <- 0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]
  y <- rpois(n, exp(eta))
  p01 <- 0.10; p10 <- 0.15; pi_z <- 0.4
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  fit <- mcglm(y, z_hat, x, family = "poisson",
               method = c("naive", "onestep"),
               p01 = p01, p10 = p10, pi_z = pi_z)

  expect_true("onestep" %in% names(fit$coefficients))
  expect_false(is.null(fit$vcov_onestep))
  expect_false(is.null(fit$loglik_onestep))
})
