# --------------------------------------------------------------------------
# Tests for multicategory BCA / BCM / CS estimators
# --------------------------------------------------------------------------

gen_multi_data <- function(n = 5000, K = 3, seed = 123) {
  set.seed(seed)
  pi_z <- c(0.5, 0.3, 0.2)
  z <- sample(0:(K - 1), n, replace = TRUE, prob = pi_z)
  Pi <- matrix(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8),
               nrow = K, byrow = TRUE)
  z_hat <- integer(n)
  for (i in seq_len(n)) z_hat[i] <- sample(0:(K - 1), 1, prob = Pi[, z[i] + 1])
  x <- cbind(1, rnorm(n))
  psi0 <- c(1.0, -0.9, 0.8, -0.7)
  gamma <- c(0, psi0[1:2]); alpha <- psi0[3:4]
  eta <- as.numeric(x %*% alpha) + gamma[z + 1]
  y <- rpois(n, exp(eta))
  xi_hat <- mcGLM:::build_xi_hat(z_hat, x, K)
  list(y = y, z_hat = z_hat, x = x, xi_hat = xi_hat,
       K = K, Pi = Pi, pi_z = pi_z, psi0 = psi0, n = n, family = "poisson")
}

test_that("fit_naive_multi returns correct structure", {
  d <- gen_multi_data()
  res <- mcGLM:::fit_naive_multi(d$y, d$xi_hat, d$family)
  expect_named(res, c("coefficients", "fitted", "glm_fit"))
  expect_length(res$coefficients, length(d$psi0))
})

test_that("BCA reduces bias for multicategory poisson", {
  d <- gen_multi_data()
  naive <- mcGLM:::fit_naive_multi(d$y, d$xi_hat, d$family)
  psi_naive <- naive$coefficients
  psi_bca <- mcGLM:::fit_bca_multi(psi_naive, d$y, d$xi_hat, d$z_hat, d$x,
                                    d$K, d$family, d$Pi, d$pi_z)
  expect_lt(abs(psi_bca[1] - d$psi0[1]), abs(psi_naive[1] - d$psi0[1]))
})

test_that("BCM reduces bias for multicategory poisson", {
  d <- gen_multi_data()
  naive <- mcGLM:::fit_naive_multi(d$y, d$xi_hat, d$family)
  psi_naive <- naive$coefficients
  psi_bcm <- mcGLM:::fit_bcm_multi(psi_naive, d$y, d$xi_hat, d$z_hat, d$x,
                                    d$K, d$family, d$Pi, d$pi_z)
  expect_lt(abs(psi_bcm[1] - d$psi0[1]), abs(psi_naive[1] - d$psi0[1]))
})

test_that("iterated BCM converges to CS for multicategory", {
  d <- gen_multi_data(n = 3000)
  naive <- mcGLM:::fit_naive_multi(d$y, d$xi_hat, d$family)
  psi_naive <- naive$coefficients
  psi_bcm_it <- mcGLM:::fit_bcm_multi(psi_naive, d$y, d$xi_hat, d$z_hat, d$x,
                                       d$K, d$family, d$Pi, d$pi_z, iterate = TRUE)
  psi_cs <- mcGLM:::fit_cs_multi(psi_naive, d$y, d$xi_hat, d$z_hat, d$x,
                                  d$K, d$family, d$Pi, d$pi_z)
  expect_equal(psi_bcm_it, psi_cs, tolerance = 1e-5)
})

test_that("CS reduces bias for multicategory poisson", {
  d <- gen_multi_data()
  naive <- mcGLM:::fit_naive_multi(d$y, d$xi_hat, d$family)
  psi_cs <- mcGLM:::fit_cs_multi(naive$coefficients, d$y, d$xi_hat, d$z_hat,
                                  d$x, d$K, d$family, d$Pi, d$pi_z)
  expect_lt(abs(psi_cs[1] - d$psi0[1]), abs(naive$coefficients[1] - d$psi0[1]))
})

test_that("CS corrected-score equation is near zero (multicategory)", {
  d <- gen_multi_data()
  naive <- mcGLM:::fit_naive_multi(d$y, d$xi_hat, d$family)
  psi_cs <- mcGLM:::fit_cs_multi(naive$coefficients, d$y, d$xi_hat, d$z_hat,
                                  d$x, d$K, d$family, d$Pi, d$pi_z)
  fam <- mcGLM:::get_link_funs(d$family)
  s <- d$K - 1; r <- ncol(d$x)
  gamma <- c(0, psi_cs[1:s]); alpha <- psi_cs[(s + 1):(s + r)]
  eta_tilde <- as.numeric(d$x %*% alpha) + gamma[d$z_hat + 1]
  resid <- d$y - fam$mu(eta_tilde)
  score <- colMeans(d$xi_hat * resid)
  m_hat <- mcGLM:::compute_mhat_multi(psi_cs, d$x, d$K, fam$mu, d$Pi, d$pi_z)
  expect_equal(score - m_hat, rep(0, length(psi_cs)), tolerance = 1e-7)
})

test_that("with identity Pi all estimators match naive (multicategory)", {
  set.seed(123); n <- 2000; pi_z <- c(0.5, 0.3, 0.2)
  z <- sample(0:2, n, replace = TRUE, prob = pi_z)
  x <- cbind(1, rnorm(n)); psi0 <- c(1.0, -0.9, 0.8, -0.7)
  gamma <- c(0, psi0[1:2]); alpha <- psi0[3:4]
  y <- rpois(n, exp(as.numeric(x %*% alpha) + gamma[z + 1]))
  Pi_id <- diag(3)
  xi_hat <- mcGLM:::build_xi_hat(z, x, 3)

  naive <- mcGLM:::fit_naive_multi(y, xi_hat, "poisson")
  psi_naive <- naive$coefficients
  psi_bca <- mcGLM:::fit_bca_multi(psi_naive, y, xi_hat, z, x, 3, "poisson",
                                    Pi_id, pi_z)
  expect_equal(psi_bca, psi_naive, tolerance = 1e-7)
})
