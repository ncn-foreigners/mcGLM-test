# mcGLM 0.4.0

## New features

* Added `mcsimex()`: **MC-SIMEX** (Monte Carlo Simulation Extrapolation)
  estimator for GLMs with misclassified covariates, following Küchenhoff,
  Mwalili & Lesaffre (2006). The simulation step (matrix exponentiation via
  eigendecomposition, multinomial resampling, and IRLS GLM fitting) is
  implemented entirely in C++ using RcppEigen for high performance.
  Supports Poisson, Binomial, and Gaussian families with binary or
  multicategory misclassification. Extrapolation methods: linear,
  quadratic, loglinear. Jackknife variance estimation included.

# mcGLM 0.3.0

## New features

* Added one-step joint estimation for **multinomial logistic** response models
  (`family = "multinomial"`). The model supports categorical responses
  `Y ∈ {0,...,J-1}` with misclassified covariates (binary or multicategory Z).
  The TMB mixture likelihood marginalizes over the latent true Z, estimating
  response-category-specific coefficients `(gamma_{j,k}, alpha_j)` for each
  non-baseline response category.

* Naive multinomial estimation via `nnet::multinom` (if available) or
  category-specific binomial GLMs as fallback.

* BCA/BCM/CS methods are not yet available for multinomial — requesting them
  produces an informative error.

# mcGLM 0.2.0

## New features

* Added frequency weights (`freq_weights` parameter) to all estimation methods
  (naive, BCA, BCM, CS, onestep) for both binary and multicategory
  misclassification. Frequency weights allow working with aggregated or
  tabulated data where each row represents multiple identical observations.
  Variance estimators are also weight-aware.

* Added comprehensive test suite (`testthat`) covering families, helpers,
  binary and multicategory estimators, variance estimators, one-step estimation,
  and frequency weights (182 tests total).

# mcGLM 0.1.0

* Initial release with naive, BCA, BCM, and CS estimators for GLMs with
  misclassified covariates.
* Support for binary and multicategory misclassification.
* One-step joint estimation via TMB mixture likelihood.
* Sandwich variance estimators for all methods.
* Iterated BCA/BCM corrections.
