#' Bias-Corrected GLM Estimation with Misclassified Covariates
#'
#' Fits a GLM where one covariate is observed with misclassification error and
#' applies bias corrections. Supports binary and multicategory misclassified
#' covariates and four estimation methods: naive, BCA (additive bias correction),
#' BCM (multiplicative bias correction), and CS (corrected score).
#'
#' @param y Numeric vector of responses (length n).
#' @param z_hat Integer vector of observed (proxy) covariate values. For the
#'   binary case, values in \{0, 1\}. For multicategory, values in
#'   \{0, 1, ..., K-1\} with baseline category 0.
#' @param x Numeric matrix of correctly observed covariates (n x r).
#'   An intercept column should be included if desired.
#' @param family A \code{\link[stats]{family}} object or one of
#'   \code{"poisson"}, \code{"binomial"}, \code{"gaussian"}.
#' @param method Character vector of methods to compute. Any subset of
#'   \code{c("naive", "bca", "bcm", "cs")}. Default is all four.
#' @param p01 False-positive rate P(Z_hat = 1 | Z = 0). Required for binary
#'   misclassification when \code{method} includes \code{"bca"}, \code{"bcm"},
#'   or \code{"cs"}.
#' @param p10 False-negative rate P(Z_hat = 0 | Z = 1). Required for binary
#'   misclassification when using corrections.
#' @param pi_z For binary case: scalar prevalence P(Z = 1). For multicategory:
#'   K-vector of class prevalences.
#' @param Pi K x K misclassification matrix for the multicategory case. Columns
#'   correspond to true categories, rows to proxy categories:
#'   \code{Pi[j, l] = P(Z_hat = j-1 | Z = l-1)}.
#' @param K Number of categories for multicategory case. If \code{NULL}
#'   (default), automatically determined: if all \code{z_hat} are in
#'   0, 1 and no \code{Pi} is supplied, binary; otherwise multicategory.
#' @param iterate Logical. If \code{TRUE}, BCA/BCM corrections are iterated
#'   until convergence. Iterated BCM converges to the corrected-score (CS)
#'   estimator. Default is \code{FALSE} (one-step correction).
#'
#' @return An object of class \code{"mcglm"}, a list with components:
#'   \describe{
#'     \item{coefficients}{Named list of coefficient vectors, one per method.}
#'     \item{naive_fit}{The \code{glm} object from the naive fit.}
#'     \item{method}{Methods used.}
#'     \item{family}{The GLM family.}
#'     \item{K}{Number of categories.}
#'     \item{n, p}{Sample size and number of parameters.}
#'   }
#'
#' @examples
#' set.seed(42)
#' n <- 5000
#' x <- cbind(1, rnorm(n))
#' z <- rbinom(n, 1, 0.4)
#' eta <- 0.8 * z - 0.5 * x[, 1] + 0.7 * x[, 2]
#' y <- rpois(n, exp(eta))
#'
#' # introduce misclassification
#' p01 <- 0.10; p10 <- 0.15
#' z_hat <- z
#' z_hat[z == 0] <- rbinom(sum(z == 0), 1, p01)
#' z_hat[z == 1] <- 1 - rbinom(sum(z == 1), 1, p10)
#'
#' fit <- mcglm(y, z_hat, x, family = "poisson",
#'              p01 = p01, p10 = p10, pi_z = 0.4)
#' fit
#'
#' @export
mcglm <- function(y, z_hat, x, family = "poisson",
                  method = c("naive", "bca", "bcm", "cs"),
                  p01 = NULL, p10 = NULL, pi_z = NULL,
                  Pi = NULL, K = NULL,
                  iterate = FALSE) {

  # --- input validation ---
  y     <- as.numeric(y)
  z_hat <- as.integer(z_hat)
  x     <- as.matrix(x)
  n     <- length(y)
  method <- match.arg(method, c("naive", "bca", "bcm", "cs"),
                      several.ok = TRUE)

  stopifnot(nrow(x) == n, length(z_hat) == n)

  # Determine binary vs multicategory
  if (is.null(K)) {
    unique_z <- sort(unique(z_hat))
    if (is.null(Pi) && all(unique_z %in% c(0L, 1L))) {
      K <- 2L
    } else {
      K <- length(unique_z)
    }
  }
  is_binary <- (K == 2L)

  # --- validate misclassification parameters ---
  needs_correction <- any(method %in% c("bca", "bcm", "cs"))
  if (needs_correction) {
    if (is_binary) {
      if (is.null(p01) || is.null(p10) || is.null(pi_z))
        stop("For binary corrections, supply p01, p10, and pi_z.")
    } else {
      if (is.null(Pi) || is.null(pi_z))
        stop("For multicategory corrections, supply Pi and pi_z.")
      Pi   <- as.matrix(Pi)
      pi_z <- as.numeric(pi_z)
      stopifnot(nrow(Pi) == K, ncol(Pi) == K, length(pi_z) == K)
    }
  }

  # --- fit ---
  results <- list()

  if (is_binary) {
    # ---- binary path ----
    naive <- fit_naive_bin(y, z_hat, x, family)
    psi   <- naive$coefficients
    results$naive <- psi
    p_total <- length(psi)

    if ("bca" %in% method)
      results$bca <- fit_bca_bin(psi, y, z_hat, x, family, p01, p10, pi_z,
                                 iterate = iterate)
    if ("bcm" %in% method)
      results$bcm <- fit_bcm_bin(psi, y, z_hat, x, family, p01, p10, pi_z,
                                 iterate = iterate)
    if ("cs"  %in% method)
      results$cs  <- fit_cs_bin(psi, y, z_hat, x, family, p01, p10, pi_z)

  } else {
    # ---- multicategory path ----
    naive <- fit_naive_multi(y, z_hat, x, K, family)
    psi   <- naive$coefficients
    results$naive <- psi
    p_total <- length(psi)

    if ("bca" %in% method)
      results$bca <- fit_bca_multi(psi, y, z_hat, x, K, family, Pi, pi_z,
                                   iterate = iterate)
    if ("bcm" %in% method)
      results$bcm <- fit_bcm_multi(psi, y, z_hat, x, K, family, Pi, pi_z,
                                   iterate = iterate)
    if ("cs"  %in% method)
      results$cs  <- fit_cs_multi(psi, y, z_hat, x, K, family, Pi, pi_z)
  }

  # --- parameter names ---
  if (is_binary) {
    nms <- c("gamma", paste0("alpha", seq_len(ncol(x)) - 1))
  } else {
    nms <- c(paste0("gamma", seq_len(K - 1)),
             paste0("alpha", seq_len(ncol(x)) - 1))
  }
  results <- lapply(results, function(v) { names(v) <- nms; v })

  structure(
    list(
      coefficients = results,
      naive_fit    = naive$glm_fit,
      method       = method,
      family       = family,
      K            = K,
      n            = n,
      p            = p_total
    ),
    class = "mcglm"
  )
}


#' @export
print.mcglm <- function(x, digits = 4, ...) {
  cat("Bias-corrected GLM with misclassified covariate\n")
  cat(sprintf("  n = %d, p = %d, K = %d\n\n", x$n, x$p, x$K))

  coefs <- do.call(cbind, x$coefficients)
  colnames(coefs) <- toupper(names(x$coefficients))
  print(round(coefs, digits))
  invisible(x)
}

#' @export
coef.mcglm <- function(object, method = NULL, ...) {
  if (is.null(method)) return(object$coefficients)
  object$coefficients[[method]]
}

#' @export
summary.mcglm <- function(object, ...) {
  cat("Bias-corrected GLM with misclassified covariate\n")
  cat(sprintf("  Family: %s\n", if (is.character(object$family))
    object$family else object$family$family))
  cat(sprintf("  n = %d, p = %d, K = %d\n", object$n, object$p, object$K))
  cat(sprintf("  Methods: %s\n\n", paste(object$method, collapse = ", ")))

  coefs <- do.call(cbind, object$coefficients)
  colnames(coefs) <- toupper(names(object$coefficients))
  cat("Coefficients:\n")
  print(round(coefs, 6))

  # bias relative to naive
  if (length(object$coefficients) > 1) {
    cat("\nBias correction (difference from naive):\n")
    naive <- object$coefficients$naive
    diffs <- sapply(object$coefficients[-1], function(v) v - naive)
    if (!is.matrix(diffs)) diffs <- matrix(diffs, nrow = 1)
    rownames(diffs) <- names(naive)
    print(round(diffs, 6))
  }

  invisible(object)
}
