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
