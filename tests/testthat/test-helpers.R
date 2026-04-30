# --------------------------------------------------------------------------
# Tests for internal helper functions (compute_delta, compute_m_bin, etc.)
# --------------------------------------------------------------------------

# Shared fixture: binary setup
setup_binary <- function(n = 200, seed = 42) {
  set.seed(seed)
  x <- cbind(1, rnorm(n))
  z <- rbinom(n, 1, 0.4)
  p01 <- 0.10; p10 <- 0.15; pi_z <- 0.4
  z_hat <- z
  z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
  z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
  psi <- c(0.8, -0.5, 0.7)  # gamma, alpha0, alpha1
  list(x = x, z = z, z_hat = z_hat, psi = psi, p01 = p01, p10 = p10,
       pi_z = pi_z, n = n)
}

# --- compute_delta ---

test_that("compute_delta returns correct toggle gap", {
  d <- setup_binary()
  mu_fun <- exp  # poisson inverse link
  delta <- mcGLM:::compute_delta(d$psi, d$x, mu_fun)

  gamma <- d$psi[1]
  alpha <- d$psi[-1]
  eta0 <- as.numeric(d$x %*% alpha)
  expected <- mu_fun(gamma + eta0) - mu_fun(eta0)

  expect_equal(delta, expected)
  expect_length(delta, d$n)
})

test_that("compute_delta is zero when gamma is zero", {
  d <- setup_binary()
  psi0 <- c(0, d$psi[-1])
  delta <- mcGLM:::compute_delta(psi0, d$x, exp)
  expect_true(all(delta == 0))
})

# --- compute_m_bin ---

test_that("compute_m_bin has correct dimensions", {
  d <- setup_binary()
  fam <- mcGLM:::get_link_funs("poisson")
  m <- mcGLM:::compute_m_bin(d$psi, d$x, fam$mu, d$p01, d$p10, d$pi_z)
  expect_equal(dim(m), c(d$n, length(d$psi)))
})

test_that("compute_m_bin uses correct constants c1, c2", {
  d <- setup_binary()
  fam <- mcGLM:::get_link_funs("poisson")
  m <- mcGLM:::compute_m_bin(d$psi, d$x, fam$mu, d$p01, d$p10, d$pi_z)
  delta <- mcGLM:::compute_delta(d$psi, d$x, fam$mu)

  c1 <- d$p01 * (1 - d$pi_z)
  c2 <- d$p01 * (1 - d$pi_z) - d$p10 * d$pi_z

  expect_equal(m[, 1], -c1 * delta)
  expect_equal(m[, 2], -c2 * delta * d$x[, 1])
  expect_equal(m[, 3], -c2 * delta * d$x[, 2])
})

test_that("compute_mhat_bin returns column means of compute_m_bin", {
  d <- setup_binary()
  fam <- mcGLM:::get_link_funs("poisson")
  m <- mcGLM:::compute_m_bin(d$psi, d$x, fam$mu, d$p01, d$p10, d$pi_z)
  mhat <- mcGLM:::compute_mhat_bin(d$psi, d$x, fam$mu, d$p01, d$p10, d$pi_z)
  expect_equal(mhat, colMeans(m))
})

# --- compute_Ihat ---

test_that("compute_Ihat is symmetric positive semi-definite", {
  d <- setup_binary()
  fam <- mcGLM:::get_link_funs("poisson")
  I <- mcGLM:::compute_Ihat(d$psi, d$z_hat, d$x, fam$mu_dot)

  expect_equal(I, t(I), tolerance = 1e-12)
  eigs <- eigen(I, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-10))
})

test_that("compute_Ihat has correct dimensions", {
  d <- setup_binary()
  fam <- mcGLM:::get_link_funs("poisson")
  I <- mcGLM:::compute_Ihat(d$psi, d$z_hat, d$x, fam$mu_dot)
  p <- length(d$psi)
  expect_equal(dim(I), c(p, p))
})

# --- compute_Mhat_bin ---

test_that("compute_Mhat_bin approximates analytical Jacobian", {
  d <- setup_binary()
  fam <- mcGLM:::get_link_funs("poisson")
  M <- mcGLM:::compute_Mhat_bin(d$psi, d$x, fam$mu, d$p01, d$p10, d$pi_z)
  expect_equal(dim(M), c(length(d$psi), length(d$psi)))

  # Cross-check with numDeriv if available
  skip_if_not_installed("numDeriv")
  M_num <- numDeriv::jacobian(
    function(p) mcGLM:::compute_mhat_bin(p, d$x, fam$mu, d$p01, d$p10, d$pi_z),
    d$psi
  )
  expect_equal(M, M_num, tolerance = 1e-4)
})

# --- Multicategory helpers ---

setup_multi <- function(n = 300, K = 3, seed = 123) {
  set.seed(seed)
  pi_z <- c(0.5, 0.3, 0.2)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  Pi <- matrix(c(
    0.8, 0.1, 0.1,
    0.1, 0.8, 0.1,
    0.1, 0.1, 0.8
  ), nrow = K, byrow = TRUE)
  z_hat <- integer(n)
  for (i in seq_len(n)) z_hat[i] <- sample(0:(K - 1), 1, prob = Pi[, z[i] + 1])
  x <- cbind(1, rnorm(n))
  psi <- c(0.5, -0.3, 0.8, -0.7)  # gamma1, gamma2, alpha0, alpha1
  list(x = x, z = z, z_hat = z_hat, psi = psi, Pi = Pi, pi_z = pi_z,
       K = K, n = n)
}

test_that("compute_m_multi has correct dimensions", {
  d <- setup_multi()
  fam <- mcGLM:::get_link_funs("poisson")
  m <- mcGLM:::compute_m_multi(d$psi, d$x, d$K, fam$mu, d$Pi, d$pi_z)
  expect_equal(dim(m), c(d$n, length(d$psi)))
})

test_that("compute_mhat_multi returns column means", {
  d <- setup_multi()
  fam <- mcGLM:::get_link_funs("poisson")
  m <- mcGLM:::compute_m_multi(d$psi, d$x, d$K, fam$mu, d$Pi, d$pi_z)
  mhat <- mcGLM:::compute_mhat_multi(d$psi, d$x, d$K, fam$mu, d$Pi, d$pi_z)
  expect_equal(mhat, colMeans(m))
})

test_that("compute_m_multi drift is zero when Pi is identity", {
  d <- setup_multi()
  fam <- mcGLM:::get_link_funs("poisson")
  Pi_id <- diag(d$K)
  m <- mcGLM:::compute_m_multi(d$psi, d$x, d$K, fam$mu, Pi_id, d$pi_z)
  # With identity Pi, Pi_{k,l}*(mu_l - mu_k) sums to zero for each k
  mhat <- colMeans(m)
  expect_equal(mhat, rep(0, length(d$psi)), tolerance = 1e-10)
})

test_that("compute_Ihat_multi is symmetric PSD", {
  d <- setup_multi()
  fam <- mcGLM:::get_link_funs("poisson")
  I <- mcGLM:::compute_Ihat_multi(d$psi, d$z_hat, d$x, d$K, fam$mu_dot)
  expect_equal(I, t(I), tolerance = 1e-12)
  eigs <- eigen(I, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-10))
})

test_that("compute_Mhat_multi matches numDeriv jacobian", {
  d <- setup_multi()
  fam <- mcGLM:::get_link_funs("poisson")
  M <- mcGLM:::compute_Mhat_multi(d$psi, d$x, d$K, fam$mu, d$Pi, d$pi_z)

  skip_if_not_installed("numDeriv")
  M_num <- numDeriv::jacobian(
    function(p) mcGLM:::compute_mhat_multi(p, d$x, d$K, fam$mu, d$Pi, d$pi_z),
    d$psi
  )
  expect_equal(M, M_num, tolerance = 1e-4)
})
