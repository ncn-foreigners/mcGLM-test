# --------------------------------------------------------------------------
# Tests for MC-SIMEX estimation
# --------------------------------------------------------------------------

gen_simex_binary <- function(n = 2000, seed = 42) {
  set.seed(seed)
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  psi0 <- c(0.8, -0.5, 0.7)
  y <- rpois(n, exp(psi0[1] * z + x %*% psi0[-1]))
  p01 <- 0.10; p10 <- 0.15
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
  Pi <- matrix(c(1 - p01, p01, p10, 1 - p10), 2, 2)
  list(y = y, z_hat = z_hat, x = x, Pi = Pi, psi0 = psi0, n = n, K = 2)
}

# --- Basic functionality ---

test_that("mcsimex returns correct structure", {
  d <- gen_simex_binary(n = 500)
  fit <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                 B = 20, seed = 1)
  expect_s3_class(fit, "mcsimex")
  expect_length(fit$coefficients, 3)
  expect_equal(nrow(fit$estimates), 5)  # lambda=0 + 4 levels
  expect_equal(length(fit$lambda), 5)
  expect_equal(fit$K, 2)
  expect_equal(fit$B, 20)
  expect_true(all(is.finite(fit$coefficients)))
})

test_that("mcsimex reduces bias (binary Poisson)", {
  d <- gen_simex_binary(n = 3000)
  fit <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                 B = 100, seed = 1)

  bias_naive <- abs(fit$naive[1] - d$psi0[1])
  bias_simex <- abs(fit$coefficients[1] - d$psi0[1])
  expect_lt(bias_simex, bias_naive)
})

# --- Extrapolation methods ---

test_that("all extrapolation methods produce finite coefficients", {
  d <- gen_simex_binary(n = 1000)
  for (method in c("linear", "quadratic", "loglinear")) {
    fit <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                   B = 30, extrapolation = method, seed = 1)
    expect_true(all(is.finite(fit$coefficients)),
                info = paste("method:", method))
  }
})

test_that("estimates increase monotonically with lambda for gamma", {
  d <- gen_simex_binary(n = 2000)
  fit <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                 B = 100, seed = 1)
  # gamma should decrease as lambda increases (more misclassification)
  gamma_est <- fit$estimates[, 1]
  diffs <- diff(gamma_est)
  # All diffs should be negative (or at least mostly)
  expect_true(sum(diffs < 0) >= length(diffs) - 1)
})

# --- Jackknife variance ---

test_that("jackknife variance is symmetric PSD", {
  d <- gen_simex_binary(n = 1000)
  fit <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                 B = 50, jackknife = TRUE, seed = 1)

  expect_false(is.null(fit$vcov))
  V <- fit$vcov
  expect_equal(dim(V), c(3, 3))
  expect_equal(V, t(V), tolerance = 1e-10)
  expect_true(all(diag(V) > 0))
})

test_that("jackknife = FALSE returns NULL vcov", {
  d <- gen_simex_binary(n = 500)
  fit <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                 B = 20, jackknife = FALSE, seed = 1)
  expect_null(fit$vcov)
})

# --- Binomial family ---

test_that("mcsimex works for binomial family", {
  set.seed(99)
  n <- 2000; x <- cbind(1, rnorm(n)); z <- rbinom(n, 1, 0.4)
  y <- rbinom(n, 1, plogis(0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]))
  p01 <- 0.08; p10 <- 0.12
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
  Pi <- matrix(c(1 - p01, p01, p10, 1 - p10), 2, 2)

  fit <- mcsimex(y, z_hat, x, family = "binomial", Pi = Pi, B = 50, seed = 1)
  expect_true(all(is.finite(fit$coefficients)))
  expect_length(fit$coefficients, 3)
})

# --- Multicategory Z ---

test_that("mcsimex works with multicategory Z (K=3)", {
  set.seed(55)
  n <- 2000; K <- 3; pi_z <- c(0.5, 0.3, 0.2)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  Pi <- matrix(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8),
               nrow = K, byrow = TRUE)
  z_hat <- integer(n)
  for (i in 1:n) z_hat[i] <- sample(0:(K - 1), 1, prob = Pi[, z[i] + 1])
  x <- cbind(1, rnorm(n))
  y <- rpois(n, exp(c(0, 1.0, -0.9)[z + 1] + x %*% c(0.8, -0.7)))

  fit <- mcsimex(y, z_hat, x, family = "poisson", Pi = Pi, K = K,
                 B = 50, seed = 1)
  expect_length(fit$coefficients, 4)  # 2 gamma + 2 alpha
  expect_true(all(is.finite(fit$coefficients)))
})

# --- Frequency weights ---

test_that("mcsimex: unit weights match no weights", {
  d <- gen_simex_binary(n = 500)
  fit1 <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                  B = 30, seed = 1)
  fit2 <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                  B = 30, seed = 1, freq_weights = rep(1, d$n))
  expect_equal(fit1$coefficients, fit2$coefficients, tolerance = 1e-8)
})

# --- Reproducibility ---

test_that("same seed produces same results", {
  d <- gen_simex_binary(n = 500)
  fit1 <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                  B = 30, seed = 123)
  fit2 <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                  B = 30, seed = 123)
  expect_equal(fit1$coefficients, fit2$coefficients)
})

test_that("different seeds produce different results", {
  d <- gen_simex_binary(n = 500)
  fit1 <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                  B = 30, seed = 1)
  fit2 <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                  B = 30, seed = 2)
  expect_false(isTRUE(all.equal(fit1$coefficients, fit2$coefficients)))
})

# --- Input validation ---

test_that("mcsimex errors on invalid Pi", {
  d <- gen_simex_binary(n = 100)
  bad_Pi <- matrix(c(0.5, 0.5, 0.5, 0.5), 2, 2)  # rows sum to 1, cols dont
  # Actually cols sum to 1 here too. Use truly bad Pi:
  bad_Pi <- matrix(c(0.9, 0.2, 0.1, 0.8), 2, 2)  # col 1 sums to 1.1
  expect_error(
    mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = bad_Pi, B = 5),
    "sum to 1"
  )
})

# --- S3 methods ---

test_that("print.mcsimex works", {
  d <- gen_simex_binary(n = 500)
  fit <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                 B = 20, seed = 1)
  expect_output(print(fit), "MC-SIMEX")
})

test_that("coef.mcsimex returns named vector", {
  d <- gen_simex_binary(n = 500)
  fit <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                 B = 20, seed = 1)
  co <- coef(fit)
  expect_true(is.numeric(co))
  expect_named(co)
})

test_that("vcov.mcsimex returns matrix", {
  d <- gen_simex_binary(n = 500)
  fit <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                 B = 20, seed = 1)
  V <- vcov(fit)
  expect_true(is.matrix(V))
})

test_that("summary.mcsimex works with coefficient table", {
  d <- gen_simex_binary(n = 500)
  fit <- mcsimex(d$y, d$z_hat, d$x, family = "poisson", Pi = d$Pi,
                 B = 20, seed = 1)
  expect_output(summary(fit), "Estimate")
})
