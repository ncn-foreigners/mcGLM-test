# --------------------------------------------------------------------------
# Tests for the main mcglm() interface and S3 methods
# --------------------------------------------------------------------------

# --- Input validation ---

test_that("mcglm errors on mismatched dimensions", {
  expect_error(
    mcglm(y = 1:10, z_hat = 1:5, x = matrix(1, 10, 1), family = "poisson",
          method = "naive")
  )
})

test_that("mcglm errors when correction params missing", {
  expect_error(
    mcglm(y = rpois(10, 1), z_hat = rbinom(10, 1, 0.5),
          x = cbind(1, rnorm(10)), family = "poisson",
          method = c("naive", "bca")),
    "supply"
  )
})

# --- Binary interface ---

test_that("mcglm binary returns correct structure", {
  set.seed(1)
  n <- 500
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  eta <- 0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]
  y <- rpois(n, exp(eta))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)

  fit <- mcglm(y, z_hat, x, family = "poisson",
               method = c("naive", "bca", "bcm", "cs"),
               p01 = 0.10, p10 = 0.15, pi_z = 0.4)

  expect_s3_class(fit, "mcglm")
  expect_named(fit$coefficients, c("naive", "bca", "bcm", "cs"))
  expect_equal(fit$K, 2L)
  expect_equal(fit$n, n)
  expect_equal(fit$p, 3)

  # Each coefficient vector should have names
  for (nm in names(fit$coefficients)) {
    expect_named(fit$coefficients[[nm]])
    expect_length(fit$coefficients[[nm]], 3)
  }
})

test_that("mcglm auto-detects binary case", {
  set.seed(2)
  n <- 200
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.5 * z_hat + x %*% c(-0.3, 0.4)))

  fit <- mcglm(y, z_hat, x, family = "poisson", method = "naive")
  expect_equal(fit$K, 2L)
})

# --- Multicategory interface ---

test_that("mcglm multicategory returns correct structure", {
  set.seed(3)
  n <- 500
  K <- 3
  pi_z <- c(0.5, 0.3, 0.2)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  Pi <- matrix(c(0.8, 0.1, 0.1,
                 0.1, 0.8, 0.1,
                 0.1, 0.1, 0.8), nrow = K, byrow = TRUE)
  z_hat <- integer(n)
  for (i in seq_len(n)) z_hat[i] <- sample(0:(K - 1), 1, prob = Pi[, z[i] + 1])

  x <- cbind(1, rnorm(n))
  gamma <- c(0, 1.0, -0.9)
  alpha <- c(0.8, -0.7)
  eta <- as.numeric(x %*% alpha) + gamma[z + 1]
  y <- rpois(n, exp(eta))

  fit <- mcglm(y, z_hat, x, family = "poisson",
               method = c("naive", "bca", "bcm", "cs"),
               Pi = Pi, pi_z = pi_z)

  expect_s3_class(fit, "mcglm")
  expect_equal(fit$K, K)
  expect_equal(fit$p, K - 1 + ncol(x))
  expect_named(fit$coefficients, c("naive", "bca", "bcm", "cs"))
})

# --- Subset methods ---

test_that("mcglm accepts subset of methods", {
  set.seed(4)
  n <- 200
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.5 * z_hat + x %*% c(-0.3, 0.4)))

  fit <- mcglm(y, z_hat, x, family = "poisson",
               method = c("naive", "cs"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.5)

  expect_named(fit$coefficients, c("naive", "cs"))
})

# --- S3 methods ---

test_that("print.mcglm runs without error", {
  set.seed(5)
  n <- 200
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.3 * z_hat + x %*% c(-0.2, 0.4)))

  fit <- mcglm(y, z_hat, x, family = "poisson",
               method = c("naive", "bca"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.5)

  expect_output(print(fit), "Bias-corrected GLM")
})

test_that("summary.mcglm runs without error", {
  set.seed(5)
  n <- 200
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.3 * z_hat + x %*% c(-0.2, 0.4)))

  fit <- mcglm(y, z_hat, x, family = "poisson",
               method = c("naive", "bca"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.5)

  expect_output(summary(fit), "Bias correction")
})

test_that("coef.mcglm returns all or selected method", {
  set.seed(6)
  n <- 200
  x <- cbind(1, rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.3 * z_hat + x %*% c(-0.2, 0.4)))

  fit <- mcglm(y, z_hat, x, family = "poisson",
               method = c("naive", "bca", "bcm"),
               p01 = 0.1, p10 = 0.1, pi_z = 0.5)

  all_coefs <- coef(fit)
  expect_type(all_coefs, "list")
  expect_named(all_coefs, c("naive", "bca", "bcm"))

  bca_coef <- coef(fit, method = "bca")
  expect_true(is.numeric(bca_coef))
  expect_length(bca_coef, 3)
})

# --- Iterate flag ---

test_that("iterate flag changes BCA/BCM results", {
  set.seed(7)
  n <- 1000
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  eta <- 0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]
  y <- rpois(n, exp(eta))
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, 0.10)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, 0.15)

  fit1 <- mcglm(y, z_hat, x, family = "poisson",
                method = c("naive", "bca", "bcm"),
                p01 = 0.10, p10 = 0.15, pi_z = 0.4,
                iterate = FALSE)
  fit2 <- mcglm(y, z_hat, x, family = "poisson",
                method = c("naive", "bca", "bcm"),
                p01 = 0.10, p10 = 0.15, pi_z = 0.4,
                iterate = TRUE)

  # Naive should be identical
  expect_equal(fit1$coefficients$naive, fit2$coefficients$naive)
  # BCA/BCM should differ
  expect_false(isTRUE(all.equal(fit1$coefficients$bca, fit2$coefficients$bca)))
})

# --- Parameter naming ---

test_that("binary parameter names are gamma, alpha0, alpha1, ...", {
  set.seed(8)
  n <- 200
  x <- cbind(1, rnorm(n), rnorm(n))
  z_hat <- rbinom(n, 1, 0.5)
  y <- rpois(n, exp(0.3 * z_hat + x %*% c(-0.2, 0.4, 0.1)))

  fit <- mcglm(y, z_hat, x, family = "poisson", method = "naive")
  expect_equal(names(fit$coefficients$naive),
               c("gamma", "alpha0", "alpha1", "alpha2"))
})

test_that("multicategory parameter names are gamma1, gamma2, ..., alpha0, ...", {
  set.seed(9)
  n <- 300
  K <- 4
  z_hat <- sample(0:(K - 1), n, replace = TRUE)
  x <- cbind(1, rnorm(n))
  y <- rpois(n, 2)

  Pi <- matrix(0.05, K, K)
  diag(Pi) <- 0.85
  pi_z <- rep(1 / K, K)

  fit <- mcglm(y, z_hat, x, family = "poisson",
               method = "naive", Pi = Pi, pi_z = pi_z)
  expect_equal(names(fit$coefficients$naive),
               c("gamma1", "gamma2", "gamma3", "alpha0", "alpha1"))
})
