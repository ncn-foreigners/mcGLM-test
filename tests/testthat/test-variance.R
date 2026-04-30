# --------------------------------------------------------------------------
# Tests for sandwich variance estimators
# --------------------------------------------------------------------------

gen_bin_for_vcov <- function(n = 2000, seed = 42) {
  set.seed(seed)
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  psi0 <- c(0.8, -0.5, 0.7)
  eta <- psi0[1] * z + x %*% psi0[-1]
  y <- rpois(n, exp(eta))
  p01 <- 0.10; p10 <- 0.15; pi_z <- 0.4
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
  xi_hat <- cbind(z_hat, x)
  list(y = y, z_hat = z_hat, x = x, xi_hat = xi_hat,
       p01 = p01, p10 = p10, pi_z = pi_z, psi0 = psi0, n = n)
}

test_that("vcov_naive is symmetric PSD", {
  d <- gen_bin_for_vcov()
  naive <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, "poisson")
  V <- mcGLM:::vcov_naive(naive$coefficients, d$y, d$xi_hat, "poisson")
  expect_equal(dim(V), c(3, 3))
  expect_equal(V, t(V), tolerance = 1e-12)
  expect_true(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values >= -1e-12))
})

test_that("vcov_naive diagonal entries are positive", {
  d <- gen_bin_for_vcov()
  naive <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, "poisson")
  V <- mcGLM:::vcov_naive(naive$coefficients, d$y, d$xi_hat, "poisson")
  expect_true(all(diag(V) > 0))
})

test_that("vcov_cs_bin is symmetric PSD", {
  d <- gen_bin_for_vcov()
  naive <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, "poisson")
  psi_cs <- mcGLM:::fit_cs_bin(naive$coefficients, d$y, d$xi_hat, d$x,
                                "poisson", d$p01, d$p10, d$pi_z)
  V <- mcGLM:::vcov_cs_bin(psi_cs, d$y, d$xi_hat, d$x, "poisson",
                            d$p01, d$p10, d$pi_z)
  expect_equal(dim(V), c(3, 3))
  expect_equal(V, t(V), tolerance = 1e-10)
  expect_true(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values >= -1e-12))
})

test_that("vcov_bc_bin uncorrected equals vcov_naive at same point", {
  d <- gen_bin_for_vcov()
  naive <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, "poisson")
  psi <- naive$coefficients
  V_naive <- mcGLM:::vcov_naive(psi, d$y, d$xi_hat, "poisson")
  V_bc    <- mcGLM:::vcov_bc_bin(psi, d$y, d$xi_hat, d$x, "poisson",
                                  type = "bca", corrected = FALSE)
  expect_equal(V_bc, V_naive, tolerance = 1e-12)
})

test_that("vcov_bc_bin corrected is symmetric PSD for BCA", {
  d <- gen_bin_for_vcov()
  naive <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, "poisson")
  psi_naive <- naive$coefficients
  psi_bca <- mcGLM:::fit_bca_bin(psi_naive, d$y, d$xi_hat, d$x, "poisson",
                                  d$p01, d$p10, d$pi_z)
  V <- mcGLM:::vcov_bc_bin(psi_bca, d$y, d$xi_hat, d$x, "poisson",
                            d$p01, d$p10, d$pi_z, psi_naive = psi_naive,
                            type = "bca", corrected = TRUE)
  expect_equal(dim(V), c(3, 3))
  expect_equal(V, t(V), tolerance = 1e-10)
  expect_true(all(diag(V) > 0))
})

test_that("vcov_bc_bin corrected is symmetric PSD for BCM", {
  d <- gen_bin_for_vcov()
  naive <- mcGLM:::fit_naive_bin(d$y, d$xi_hat, "poisson")
  psi_naive <- naive$coefficients
  psi_bcm <- mcGLM:::fit_bcm_bin(psi_naive, d$y, d$xi_hat, d$x, "poisson",
                                  d$p01, d$p10, d$pi_z)
  V <- mcGLM:::vcov_bc_bin(psi_bcm, d$y, d$xi_hat, d$x, "poisson",
                            d$p01, d$p10, d$pi_z, psi_naive = psi_naive,
                            type = "bcm", corrected = TRUE)
  expect_equal(dim(V), c(3, 3))
  expect_equal(V, t(V), tolerance = 1e-10)
  expect_true(all(diag(V) > 0))
})

# --- Multicategory ---

gen_multi_for_vcov <- function(n = 2000, K = 3, seed = 123) {
  set.seed(seed)
  pi_z <- c(0.5, 0.3, 0.2)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  Pi <- matrix(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8),
               nrow = K, byrow = TRUE)
  z_hat <- integer(n)
  for (i in seq_len(n)) z_hat[i] <- sample(0:(K - 1), 1, prob = Pi[, z[i] + 1])
  x <- cbind(1, rnorm(n)); psi0 <- c(1.0, -0.9, 0.8, -0.7)
  gamma <- c(0, psi0[1:2]); alpha <- psi0[3:4]
  y <- rpois(n, exp(as.numeric(x %*% alpha) + gamma[z + 1]))
  xi_hat <- mcGLM:::build_xi_hat(z_hat, x, K)
  list(y = y, z_hat = z_hat, x = x, xi_hat = xi_hat,
       K = K, Pi = Pi, pi_z = pi_z, psi0 = psi0, n = n)
}

test_that("vcov_naive_multi is symmetric PSD", {
  d <- gen_multi_for_vcov()
  naive <- mcGLM:::fit_naive_multi(d$y, d$xi_hat, "poisson")
  V <- mcGLM:::vcov_naive_multi(naive$coefficients, d$y, d$xi_hat,
                                 d$z_hat, d$x, d$K, "poisson")
  p <- length(d$psi0)
  expect_equal(dim(V), c(p, p))
  expect_equal(V, t(V), tolerance = 1e-12)
  expect_true(all(diag(V) > 0))
})

test_that("vcov_cs_multi is symmetric PSD", {
  d <- gen_multi_for_vcov()
  naive <- mcGLM:::fit_naive_multi(d$y, d$xi_hat, "poisson")
  psi_cs <- mcGLM:::fit_cs_multi(naive$coefficients, d$y, d$xi_hat, d$z_hat,
                                  d$x, d$K, "poisson", d$Pi, d$pi_z)
  V <- mcGLM:::vcov_cs_multi(psi_cs, d$y, d$xi_hat, d$z_hat, d$x,
                              d$K, "poisson", d$Pi, d$pi_z)
  p <- length(d$psi0)
  expect_equal(dim(V), c(p, p))
  expect_equal(V, t(V), tolerance = 1e-10)
  expect_true(all(diag(V) > 0))
})
