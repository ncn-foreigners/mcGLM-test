# --------------------------------------------------------------------------
# Tests for binary BCA / BCM / CS estimators
# --------------------------------------------------------------------------

gen_binary_data <- function(n = 5000, seed = 42) {
  set.seed(seed)
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  gamma0 <- 0.8; alpha <- c(-0.5, 0.7)
  psi0 <- c(gamma0, alpha)
  eta <- gamma0 * z + x %*% alpha
  y <- rpois(n, exp(eta))

  p01 <- 0.10; p10 <- 0.15; pi_z <- 0.4
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)

  xi_hat <- cbind(z_hat, x)
  list(y = y, z_hat = z_hat, x = x, xi_hat = xi_hat,
       p01 = p01, p10 = p10, pi_z = pi_z, psi0 = psi0, n = n, family = "poisson")
}

test_that("fit_naive_bin returns correct structure", {
  d <- gen_binary_data()
  res <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, d$family)
  expect_named(res, c("coefficients", "fitted", "glm_fit"))
  expect_length(res$coefficients, 3)
  expect_length(res$fitted, d$n)
  expect_s3_class(res$glm_fit, "glm")
})

test_that("BCA reduces bias relative to naive (binary poisson)", {
  d <- gen_binary_data()
  naive <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, d$family)
  psi_naive <- naive$coefficients
  psi_bca <- mcGLM:::fit_bca_bin(psi_naive, d$y, d$xi_hat, d$x, d$family,
                                  d$p01, d$p10, d$pi_z)
  expect_lt(abs(psi_bca[1] - d$psi0[1]), abs(psi_naive[1] - d$psi0[1]))
})

test_that("BCA returns numeric vector of correct length", {
  d <- gen_binary_data()
  naive <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, d$family)
  psi_bca <- mcGLM:::fit_bca_bin(naive$coefficients, d$y, d$xi_hat, d$x,
                                  d$family, d$p01, d$p10, d$pi_z)
  expect_true(is.numeric(psi_bca))
  expect_length(psi_bca, 3)
})

test_that("iterated BCA converges", {
  d <- gen_binary_data()
  naive <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, d$family)
  psi_bca_1 <- mcGLM:::fit_bca_bin(naive$coefficients, d$y, d$xi_hat, d$x,
                                    d$family, d$p01, d$p10, d$pi_z,
                                    iterate = FALSE)
  psi_bca_it <- mcGLM:::fit_bca_bin(naive$coefficients, d$y, d$xi_hat, d$x,
                                     d$family, d$p01, d$p10, d$pi_z,
                                     iterate = TRUE)
  expect_false(isTRUE(all.equal(psi_bca_1, psi_bca_it, tolerance = 1e-6)))
  expect_true(is.numeric(psi_bca_it))
})

test_that("BCM reduces bias relative to naive (binary poisson)", {
  d <- gen_binary_data()
  naive <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, d$family)
  psi_naive <- naive$coefficients
  psi_bcm <- mcGLM:::fit_bcm_bin(psi_naive, d$y, d$xi_hat, d$x, d$family,
                                  d$p01, d$p10, d$pi_z)
  expect_lt(abs(psi_bcm[1] - d$psi0[1]), abs(psi_naive[1] - d$psi0[1]))
})

test_that("iterated BCM converges to CS", {
  d <- gen_binary_data()
  naive <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, d$family)
  psi_naive <- naive$coefficients
  psi_bcm_it <- mcGLM:::fit_bcm_bin(psi_naive, d$y, d$xi_hat, d$x, d$family,
                                     d$p01, d$p10, d$pi_z, iterate = TRUE)
  psi_cs <- mcGLM:::fit_cs_bin(psi_naive, d$y, d$xi_hat, d$x, d$family,
                                d$p01, d$p10, d$pi_z)
  expect_equal(unname(psi_bcm_it), unname(psi_cs), tolerance = 1e-6)
})

test_that("CS reduces bias relative to naive (binary poisson)", {
  d <- gen_binary_data()
  naive <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, d$family)
  psi_naive <- naive$coefficients
  psi_cs <- mcGLM:::fit_cs_bin(psi_naive, d$y, d$xi_hat, d$x, d$family,
                                d$p01, d$p10, d$pi_z)
  expect_lt(abs(psi_cs[1] - d$psi0[1]), abs(psi_naive[1] - d$psi0[1]))
})

test_that("CS corrected-score equation is near zero at solution", {
  d <- gen_binary_data()
  naive <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, d$family)
  psi_cs <- mcGLM:::fit_cs_bin(naive$coefficients, d$y, d$xi_hat, d$x,
                                d$family, d$p01, d$p10, d$pi_z)
  fam <- mcGLM:::get_link_funs(d$family)
  eta_tilde <- as.numeric(d$xi_hat %*% psi_cs)
  resid <- d$y - fam$mu(eta_tilde)
  score <- colMeans(d$xi_hat * resid)
  m_hat <- mcGLM:::compute_mhat_bin(psi_cs, d$x, fam$mu, d$p01, d$p10, d$pi_z)
  expect_equal(unname(score - m_hat), rep(0, length(psi_cs)), tolerance = 1e-8)
})

test_that("BCA/BCM/CS work for binomial family (binary)", {
  set.seed(99)
  n <- 3000; x <- cbind(1, rnorm(n)); z <- rbinom(n, 1, 0.4)
  psi0 <- c(0.8, -0.5, 0.7)
  eta <- psi0[1] * z + x %*% psi0[-1]
  y <- rbinom(n, 1, plogis(eta))
  p01 <- 0.10; p10 <- 0.15; pi_z <- 0.4
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
  xi_hat <- cbind(z_hat, x)

  naive <- mcGLM:::fit_naive_bin(y, xi_hat, "binomial")
  psi_naive <- naive$coefficients
  psi_bca <- mcGLM:::fit_bca_bin(psi_naive, y, xi_hat, x, "binomial",
                                  p01, p10, pi_z)
  psi_bcm <- mcGLM:::fit_bcm_bin(psi_naive, y, xi_hat, x, "binomial",
                                  p01, p10, pi_z)
  psi_cs  <- mcGLM:::fit_cs_bin(psi_naive, y, xi_hat, x, "binomial",
                                 p01, p10, pi_z)
  expect_true(all(is.finite(c(psi_bca, psi_bcm, psi_cs))))
})

test_that("with zero misclassification BCA/BCM/CS match naive", {
  set.seed(42); n <- 2000; x <- cbind(1, rnorm(n)); z <- rbinom(n, 1, 0.4)
  eta <- 0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]
  y <- rpois(n, exp(eta)); z_hat <- z
  xi_hat <- cbind(z_hat, x)

  naive <- mcGLM:::fit_naive_bin(y, xi_hat, "poisson")
  psi_naive <- naive$coefficients
  psi_bca <- mcGLM:::fit_bca_bin(psi_naive, y, xi_hat, x, "poisson", 0, 0, 0.4)
  psi_bcm <- mcGLM:::fit_bcm_bin(psi_naive, y, xi_hat, x, "poisson", 0, 0, 0.4)
  expect_equal(unname(psi_bca), unname(psi_naive), tolerance = 1e-8)
  expect_equal(unname(psi_bcm), unname(psi_naive), tolerance = 1e-8)
})
