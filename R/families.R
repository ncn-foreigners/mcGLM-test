#' Get inverse-link, its derivative, and second derivative for a GLM family
#'
#' @param family A \code{\link[stats]{family}} object or character string
#'   ("poisson", "binomial", "gaussian").
#' @return A list with components \code{mu} (inverse link), \code{mu_dot}
#'   (first derivative), and \code{mu_ddot} (second derivative).
#' @keywords internal
get_link_funs <- function(family) {
  if (is.character(family)) {
    family <- switch(family,
      poisson  = stats::poisson(),
      binomial = stats::binomial(),
      gaussian = stats::gaussian(),
      stop("Unsupported family string: ", family)
    )
  }

  mu     <- family$linkinv
  mu_eta <- family$mu.eta          # d mu / d eta

  # second derivative via numerical differentiation
  mu_ddot <- function(eta) {
    h <- 1e-7 * pmax(abs(eta), 1)
    (mu_eta(eta + h) - mu_eta(eta - h)) / (2 * h)
  }

  # log-density function (for one-step method diagnostics)
  dlog <- switch(family$family,
    gaussian = function(y, mu, sigma) stats::dnorm(y, mu, sigma, log = TRUE),
    poisson  = function(y, mu, ...) stats::dpois(y, mu, log = TRUE),
    binomial = function(y, mu, ...) stats::dbinom(y, 1, mu, log = TRUE),
    NULL
  )

  list(mu = mu, mu_dot = mu_eta, mu_ddot = mu_ddot, dlog = dlog, family = family)
}
