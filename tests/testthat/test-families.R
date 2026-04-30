test_that("get_link_funs returns correct inverse link for poisson", {
  funs <- mcGLM:::get_link_funs("poisson")
  expect_equal(funs$mu(0), 1)
  expect_equal(funs$mu(1), exp(1))
  expect_equal(funs$mu(-2), exp(-2))
})

test_that("get_link_funs returns correct inverse link for binomial", {
  funs <- mcGLM:::get_link_funs("binomial")
  expect_equal(funs$mu(0), 0.5)
  expect_true(funs$mu(10) > 0.999)
  expect_true(funs$mu(-10) < 0.001)
})

test_that("get_link_funs returns correct inverse link for gaussian", {
  funs <- mcGLM:::get_link_funs("gaussian")
  expect_equal(funs$mu(3.7), 3.7)
})

test_that("get_link_funs mu_dot matches numerical derivative", {
  funs <- mcGLM:::get_link_funs("poisson")
  eta <- 1.5
  h <- 1e-6
  numerical <- (funs$mu(eta + h) - funs$mu(eta - h)) / (2 * h)
  expect_equal(funs$mu_dot(eta), numerical, tolerance = 1e-5)
})

test_that("get_link_funs mu_ddot matches numerical second derivative", {
  funs <- mcGLM:::get_link_funs("poisson")
  eta <- 1.0
  h <- 1e-5
  numerical <- (funs$mu_dot(eta + h) - funs$mu_dot(eta - h)) / (2 * h)
  expect_equal(funs$mu_ddot(eta), numerical, tolerance = 1e-4)
})

test_that("get_link_funs accepts family objects", {
  funs <- mcGLM:::get_link_funs(stats::poisson())
  expect_equal(funs$mu(0), 1)
})

test_that("get_link_funs errors on unsupported family", {
  expect_error(mcGLM:::get_link_funs("gamma"), "Unsupported")
})

test_that("get_link_funs returns dlog for supported families", {
  for (fam in c("poisson", "binomial", "gaussian")) {
    funs <- mcGLM:::get_link_funs(fam)
    expect_true(is.function(funs$dlog))
  }
})

test_that("dlog returns correct log-densities", {
  funs_p <- mcGLM:::get_link_funs("poisson")
  expect_equal(funs_p$dlog(3, 2), dpois(3, 2, log = TRUE))

  funs_b <- mcGLM:::get_link_funs("binomial")
  expect_equal(funs_b$dlog(1, 0.7), dbinom(1, 1, 0.7, log = TRUE))

  funs_g <- mcGLM:::get_link_funs("gaussian")
  expect_equal(funs_g$dlog(1.5, 2.0, 0.5), dnorm(1.5, 2.0, 0.5, log = TRUE))
})
