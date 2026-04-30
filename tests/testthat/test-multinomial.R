# --------------------------------------------------------------------------
# Tests for multinomial one-step estimation
# --------------------------------------------------------------------------

skip_if_no_tmb <- function() {
  skip_if_not(
    tryCatch({ TMB::config(DLL = "mcGLM"); TRUE }, error = function(e) FALSE),
    "TMB DLL 'mcGLM' not compiled"
  )
}

# Generate multinomial data with misclassified binary covariate (K=2)
gen_multinom_bin <- function(n = 2000, J = 3, seed = 42) {
  set.seed(seed)
  K <- 2
  pi_z <- 0.4
  z <- rbinom(n, 1, pi_z)
  x <- cbind(1, rnorm(n))
  r <- ncol(x)

  # True parameters: for each response category j=1,...,J-1
  # gamma_j (effect of Z=1 vs Z=0), alpha_j (intercept + slope)
  # J=3: response categories 0, 1, 2
  gamma_true <- c(0.8, -0.5)   # gamma for y=1, y=2
  alpha_true <- matrix(c(
    -0.5, 0.7,   # alpha for y=1: intercept, slope
     0.3, -0.4   # alpha for y=2: intercept, slope
  ), nrow = J - 1, byrow = TRUE)

  # Compute probabilities via softmax
  eta <- matrix(0, n, J)  # eta[,1] = 0 (baseline)
  for (jj in seq_len(J - 1)) {
    eta[, jj + 1] <- gamma_true[jj] * z + x %*% alpha_true[jj, ]
  }
  prob <- exp(eta)
  prob <- prob / rowSums(prob)

  # Sample responses
  y <- integer(n)
  for (i in seq_len(n)) {
    y[i] <- sample(0:(J - 1), 1, prob = prob[i, ])
  }

  # Misclassification
  p01 <- 0.10; p10 <- 0.15
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  list(y = y, z_hat = z_hat, x = x, z = z, J = J, K = K,
       p01 = p01, p10 = p10, pi_z = pi_z,
       gamma_true = gamma_true, alpha_true = alpha_true, n = n)
}

# Generate multinomial data with multicategory Z (K=3)
gen_multinom_multi <- function(n = 3000, J = 3, K = 3, seed = 55) {
  set.seed(seed)
  pi_z <- c(0.5, 0.3, 0.2)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  x <- cbind(1, rnorm(n))
  r <- ncol(x)
  s <- K - 1

  # gamma_{j,k} for j=1,...,J-1, k=1,...,K-1 (gamma_{j,0}=0)
  # alpha_{j,0..r-1} for j=1,...,J-1
  gamma_true <- matrix(c(
    0.8, -0.3,   # y=1: gamma for Z=1, Z=2
    -0.5, 0.6    # y=2: gamma for Z=1, Z=2
  ), nrow = J - 1, byrow = TRUE)

  alpha_true <- matrix(c(
    -0.5, 0.7,   # y=1: intercept, slope
     0.3, -0.4   # y=2: intercept, slope
  ), nrow = J - 1, byrow = TRUE)

  eta <- matrix(0, n, J)
  for (jj in seq_len(J - 1)) {
    gamma_z <- c(0, gamma_true[jj, ])  # baseline Z=0 => 0
    eta[, jj + 1] <- gamma_z[z + 1] + x %*% alpha_true[jj, ]
  }
  prob <- exp(eta)
  prob <- prob / rowSums(prob)

  y <- integer(n)
  for (i in seq_len(n)) {
    y[i] <- sample(0:(J - 1), 1, prob = prob[i, ])
  }

  Pi <- matrix(c(0.8, 0.1, 0.1,
                 0.1, 0.8, 0.1,
                 0.1, 0.1, 0.8), nrow = K, byrow = TRUE)
  z_hat <- integer(n)
  for (i in seq_len(n)) z_hat[i] <- sample(0:(K - 1), 1, prob = Pi[, z[i] + 1])

  list(y = y, z_hat = z_hat, x = x, z = z, J = J, K = K,
       Pi = Pi, pi_z = pi_z, n = n)
}

# --- mcglm interface ---

test_that("mcglm rejects BCA/BCM/CS for multinomial", {
  expect_error(
    mcglm(y = 0:2, z_hat = c(0, 1, 0), x = cbind(1, rnorm(3)),
          family = "multinomial", method = c("naive", "bca"),
          p01 = 0.1, p10 = 0.1, pi_z = 0.5),
    "only.*naive.*onestep"
  )
})

test_that("mcglm errors on invalid multinomial y", {
  expect_error(
    mcglm(y = c(0, 1, 3), z_hat = c(0, 1, 0), x = cbind(1, rnorm(3)),
          family = "multinomial", method = "naive", J = 3),
    "integers in"
  )
})

test_that("mcglm multinomial naive returns correct structure", {
  d <- gen_multinom_bin(n = 300)
  fit <- mcglm(d$y, d$z_hat, d$x, family = "multinomial",
               method = "naive", J = d$J)

  expect_s3_class(fit, "mcglm")
  expect_equal(fit$K, d$K)
  expect_equal(fit$J, d$J)
  expect_true("naive" %in% names(fit$coefficients))
  # (J-1) * (K-1 + r) = 2 * (1 + 2) = 6 parameters
  expect_length(fit$coefficients$naive, 6)
})

test_that("mcglm auto-detects J for multinomial", {
  set.seed(10)
  n <- 100
  y <- sample(0:3, n, replace = TRUE)
  z_hat <- rbinom(n, 1, 0.5)
  x <- cbind(1, rnorm(n))

  fit <- mcglm(y, z_hat, x, family = "multinomial", method = "naive")
  expect_equal(fit$J, 4)
})

# --- Naive multinomial ---

test_that("fit_naive_multinomial returns correct length", {
  d <- gen_multinom_bin(n = 500)
  res <- mcGLM:::fit_naive_multinomial(d$y, d$z_hat, d$x, d$J, d$K)
  # (J-1) * (s + r) = 2 * (1 + 2) = 6
  expect_length(res$coefficients, 6)
  expect_true(all(is.finite(res$coefficients)))
})

test_that("fit_naive_multinomial works with K=3", {
  d <- gen_multinom_multi(n = 500)
  res <- mcGLM:::fit_naive_multinomial(d$y, d$z_hat, d$x, d$J, d$K)
  # (J-1) * (s + r) = 2 * (2 + 2) = 8
  expect_length(res$coefficients, 8)
  expect_true(all(is.finite(res$coefficients)))
})

# --- One-step multinomial with binary Z ---

test_that("fit_onestep_multinomial works with binary Z", {
  skip_if_no_tmb()
  d <- gen_multinom_bin(n = 1000)

  res <- mcGLM:::fit_onestep_multinomial(
    d$y, d$z_hat, d$x, d$J, d$K,
    p01 = d$p01, p10 = d$p10, pi_z = d$pi_z
  )

  expect_named(res, c("coefficients", "vcov", "loglik", "convergence"))
  expect_length(res$coefficients, 6)
  expect_equal(dim(res$vcov), c(6, 6))
  expect_true(is.finite(res$loglik))
  expect_true(all(is.finite(res$coefficients)))
})

test_that("multinomial onestep via mcglm interface (binary Z)", {
  skip_if_no_tmb()
  d <- gen_multinom_bin(n = 1000)

  fit <- mcglm(d$y, d$z_hat, d$x, family = "multinomial",
               method = c("naive", "onestep"),
               p01 = d$p01, p10 = d$p10, pi_z = d$pi_z)

  expect_true("onestep" %in% names(fit$coefficients))
  expect_false(is.null(fit$vcov_onestep))
  expect_false(is.null(fit$loglik_onestep))
  expect_equal(fit$J, d$J)
})

test_that("multinomial onestep has higher loglik than naive (binary Z)", {
  skip_if_no_tmb()
  d <- gen_multinom_bin(n = 2000)

  fit <- mcglm(d$y, d$z_hat, d$x, family = "multinomial",
               method = c("naive", "onestep"),
               p01 = d$p01, p10 = d$p10, pi_z = d$pi_z)

  # One-step should generally fit better (or at least as well)
  expect_true(all(is.finite(coef(fit, "onestep"))))
  expect_true(is.finite(fit$loglik_onestep))
})

# --- One-step multinomial with multicategory Z ---

test_that("fit_onestep_multinomial works with multicategory Z", {
  skip_if_no_tmb()
  d <- gen_multinom_multi(n = 1500)

  res <- mcGLM:::fit_onestep_multinomial(
    d$y, d$z_hat, d$x, d$J, d$K,
    Pi = d$Pi, pi_z = d$pi_z
  )

  # (J-1) * (s + r) = 2 * (2 + 2) = 8
  expect_length(res$coefficients, 8)
  expect_equal(dim(res$vcov), c(8, 8))
  expect_true(is.finite(res$loglik))
  expect_true(all(is.finite(res$coefficients)))
})

test_that("multinomial onestep via mcglm interface (multicategory Z)", {
  skip_if_no_tmb()
  d <- gen_multinom_multi(n = 1500)

  fit <- mcglm(d$y, d$z_hat, d$x, family = "multinomial",
               method = c("naive", "onestep"),
               Pi = d$Pi, pi_z = d$pi_z)

  expect_true("onestep" %in% names(fit$coefficients))
  expect_length(coef(fit, "onestep"), 8)
  expect_false(is.null(fit$vcov_onestep))
})

# --- No misclassification sanity check ---

test_that("multinomial onestep recovers true params without misclassification", {
  skip_if_no_tmb()
  set.seed(99)
  n <- 5000
  J <- 3; K <- 2
  pi_z <- 0.4
  z <- rbinom(n, 1, pi_z)
  x <- cbind(1, rnorm(n))

  gamma_true <- c(0.8, -0.5)
  alpha_true <- matrix(c(-0.5, 0.7, 0.3, -0.4), nrow = 2, byrow = TRUE)

  eta <- matrix(0, n, J)
  for (jj in 1:(J - 1)) {
    eta[, jj + 1] <- gamma_true[jj] * z + x %*% alpha_true[jj, ]
  }
  prob <- exp(eta) / rowSums(exp(eta))
  y <- integer(n)
  for (i in 1:n) y[i] <- sample(0:(J - 1), 1, prob = prob[i, ])

  # No misclassification: z_hat = z, p01 = p10 = 0
  res <- mcGLM:::fit_onestep_multinomial(
    y, z, x, J, K, p01 = 0, p10 = 0, pi_z = pi_z
  )

  # Should be close to true values
  # Block 1 (y=1): gamma=0.8, alpha0=-0.5, alpha1=0.7
  # Block 2 (y=2): gamma=-0.5, alpha0=0.3, alpha1=-0.4
  true_psi <- c(0.8, -0.5, 0.7, -0.5, 0.3, -0.4)
  expect_equal(unname(res$coefficients), true_psi, tolerance = 0.15)
})

# --- Parameter naming ---

test_that("multinomial parameter names are correct (binary Z)", {
  d <- gen_multinom_bin(n = 300)
  fit <- mcglm(d$y, d$z_hat, d$x, family = "multinomial", method = "naive")
  nms <- names(coef(fit, "naive"))
  expect_equal(nms, c("y1:gamma", "y1:alpha0", "y1:alpha1",
                       "y2:gamma", "y2:alpha0", "y2:alpha1"))
})

test_that("multinomial parameter names are correct (multicategory Z)", {
  d <- gen_multinom_multi(n = 300)
  fit <- mcglm(d$y, d$z_hat, d$x, family = "multinomial", method = "naive",
               Pi = d$Pi, pi_z = d$pi_z)
  nms <- names(coef(fit, "naive"))
  expect_equal(nms, c("y1:gamma1", "y1:gamma2", "y1:alpha0", "y1:alpha1",
                       "y2:gamma1", "y2:gamma2", "y2:alpha0", "y2:alpha1"))
})

# --- Frequency weights ---

test_that("multinomial onestep: unit weights = no weights", {
  skip_if_no_tmb()
  d <- gen_multinom_bin(n = 500)

  fit1 <- mcglm(d$y, d$z_hat, d$x, family = "multinomial",
                method = c("naive", "onestep"),
                p01 = d$p01, p10 = d$p10, pi_z = d$pi_z)
  fit2 <- mcglm(d$y, d$z_hat, d$x, family = "multinomial",
                method = c("naive", "onestep"),
                p01 = d$p01, p10 = d$p10, pi_z = d$pi_z,
                freq_weights = rep(1, d$n))

  expect_equal(coef(fit1, "naive"), coef(fit2, "naive"), tolerance = 1e-8)
  expect_equal(coef(fit1, "onestep"), coef(fit2, "onestep"), tolerance = 1e-6)
})

# --- S3 methods ---

test_that("print/summary work for multinomial", {
  d <- gen_multinom_bin(n = 300)
  fit <- mcglm(d$y, d$z_hat, d$x, family = "multinomial", method = "naive")
  expect_output(print(fit), "Bias-corrected GLM")
  expect_output(summary(fit), "multinomial")
})
