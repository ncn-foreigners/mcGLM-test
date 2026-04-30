// src/mcsimex_core.cpp
// High-performance MC-SIMEX simulation step in C++
//
// Handles: matrix exponentiation via eigendecomposition, multinomial
// sampling of misclassified data, and IRLS GLM fitting — all in C++
// to avoid the R overhead of B × |lambda| model refits.

#include <RcppEigen.h>
#include <random>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

// --------------------------------------------------------------------------
// Matrix exponentiation via eigendecomposition: Pi^lambda
// --------------------------------------------------------------------------
// Pi is a K×K misclassification matrix. We compute Pi^lambda using
// the eigendecomposition Pi = V * diag(d) * V^{-1}, so
// Pi^lambda = V * diag(d^lambda) * V^{-1}.
// For real eigenvalues only (which holds for stochastic matrices
// satisfying the diagonal dominance condition).

MatrixXd mat_power(const MatrixXd& Pi, double lambda) {
  int K = Pi.rows();
  if (lambda == 0.0) return MatrixXd::Identity(K, K);
  if (lambda == 1.0) return Pi;

  Eigen::EigenSolver<MatrixXd> es(Pi);
  // Use real parts (valid for well-conditioned stochastic matrices)
  MatrixXd V = es.eigenvectors().real();
  VectorXd d = es.eigenvalues().real();

  // d^lambda (element-wise)
  VectorXd d_pow(K);
  for (int k = 0; k < K; k++) {
    d_pow(k) = std::pow(std::abs(d(k)), lambda);
    if (d(k) < 0) d_pow(k) = -d_pow(k);  // preserve sign
  }

  MatrixXd V_inv = V.inverse();
  return V * d_pow.asDiagonal() * V_inv;
}


// --------------------------------------------------------------------------
// Multinomial sampling: resample z_hat according to Pi^lambda
// --------------------------------------------------------------------------
// For each observation i with z_hat[i] = k, sample a new category
// from the k-th column of Pi^lambda (the conditional distribution of
// the true category given the observed proxy).

VectorXi resample_z(const VectorXi& z_hat, const MatrixXd& Pi_lam,
                     std::mt19937& rng) {
  int n = z_hat.size();
  int K = Pi_lam.rows();
  VectorXi z_new(n);
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  for (int i = 0; i < n; i++) {
    int k = z_hat(i);
    double u = unif(rng);
    double cumprob = 0.0;
    int cat = K - 1;  // default to last category
    for (int j = 0; j < K; j++) {
      cumprob += Pi_lam(j, k);
      if (u <= cumprob) { cat = j; break; }
    }
    z_new(i) = cat;
  }
  return z_new;
}


// --------------------------------------------------------------------------
// Build design matrix xi_hat = [dummies(z), x] for a given z vector
// --------------------------------------------------------------------------

MatrixXd build_design(const VectorXi& z, const MatrixXd& x, int K) {
  int n = z.size();
  int s = K - 1;
  int r = x.cols();
  int p = s + r;
  MatrixXd xi(n, p);

  // Dummy columns (baseline = 0)
  for (int k = 0; k < s; k++) {
    for (int i = 0; i < n; i++) {
      xi(i, k) = (z(i) == k + 1) ? 1.0 : 0.0;
    }
  }
  // Covariate columns
  xi.rightCols(r) = x;
  return xi;
}


// --------------------------------------------------------------------------
// IRLS GLM fitting (Poisson log link, Binomial logit, Gaussian identity)
// --------------------------------------------------------------------------
// Returns the p-vector of coefficients.
// dist_code: 1=Gaussian, 2=Poisson, 3=Binomial
// wt: frequency weights (length n, all positive)

// Inverse link functions
inline double mu_gaussian(double eta) { return eta; }
inline double mu_poisson(double eta) { return std::exp(eta); }
inline double mu_binomial(double eta) {
  return 1.0 / (1.0 + std::exp(-eta));
}

// d mu / d eta
inline double mu_dot_gaussian(double) { return 1.0; }
inline double mu_dot_poisson(double eta) { return std::exp(eta); }
inline double mu_dot_binomial(double eta) {
  double p = 1.0 / (1.0 + std::exp(-eta));
  return p * (1.0 - p);
}

// Variance function V(mu)
inline double var_gaussian(double) { return 1.0; }
inline double var_poisson(double mu) { return mu; }
inline double var_binomial(double mu) { return mu * (1.0 - mu); }

VectorXd irls_fit(const VectorXd& y, const MatrixXd& X,
                   int dist_code, const VectorXd& wt,
                   int max_iter = 25, double tol = 1e-8) {
  int n = X.rows();
  int p = X.cols();

  // Function pointers for the family
  double (*mu_fn)(double);
  double (*mu_dot_fn)(double);
  double (*var_fn)(double);

  switch (dist_code) {
    case 1: mu_fn = mu_gaussian; mu_dot_fn = mu_dot_gaussian;
            var_fn = var_gaussian; break;
    case 2: mu_fn = mu_poisson;  mu_dot_fn = mu_dot_poisson;
            var_fn = var_poisson;  break;
    case 3: mu_fn = mu_binomial; mu_dot_fn = mu_dot_binomial;
            var_fn = var_binomial; break;
    default: Rcpp::stop("Unknown dist_code");
  }

  // Initialize beta
  VectorXd beta = VectorXd::Zero(p);
  // Warm start: for Poisson, initialize intercept to log(mean(y))
  if (dist_code == 2) {
    double ymean = 0;
    double wtsum = 0;
    for (int i = 0; i < n; i++) { ymean += wt(i) * y(i); wtsum += wt(i); }
    ymean /= wtsum;
    if (ymean > 0) beta(0) = std::log(ymean);
  }
  // For Binomial, initialize to logit(mean(y))
  if (dist_code == 3) {
    double ymean = 0;
    double wtsum = 0;
    for (int i = 0; i < n; i++) { ymean += wt(i) * y(i); wtsum += wt(i); }
    ymean /= wtsum;
    ymean = std::max(0.01, std::min(0.99, ymean));
    beta(0) = std::log(ymean / (1.0 - ymean));
  }

  VectorXd eta(n), mu(n), w_irls(n), z_irls(n);

  for (int iter = 0; iter < max_iter; iter++) {
    // Compute eta, mu, working weights and working response
    eta = X * beta;
    for (int i = 0; i < n; i++) {
      mu(i) = mu_fn(eta(i));
      double md = mu_dot_fn(eta(i));
      double v = var_fn(mu(i));
      // IRLS weight: wt_i * (d mu/d eta)^2 / V(mu)
      w_irls(i) = wt(i) * md * md / std::max(v, 1e-10);
      // Working response: eta + (y - mu) / (d mu / d eta)
      z_irls(i) = eta(i) + (y(i) - mu(i)) / std::max(std::abs(md), 1e-10);
    }

    // Weighted least squares: (X' W X) beta = X' W z
    // Form X'WX and X'Wz efficiently
    MatrixXd XtWX(p, p);
    VectorXd XtWz(p);

    // Use Eigen's efficient operations
    MatrixXd Xw = X.array().colwise() * w_irls.array();  // X * diag(w)
    XtWX.noalias() = Xw.transpose() * X;
    XtWz.noalias() = Xw.transpose() * z_irls;

    VectorXd beta_new = XtWX.ldlt().solve(XtWz);

    double change = (beta_new - beta).cwiseAbs().maxCoeff();
    beta = beta_new;
    if (change < tol) break;
  }

  return beta;
}


// --------------------------------------------------------------------------
// Main MC-SIMEX simulation step (exported to R)
// --------------------------------------------------------------------------
// Returns a matrix of dimension (B * n_lambda) x p, where each row
// is a coefficient vector from one simulation replicate at one lambda level.
// The R side then averages over B and does the extrapolation.

// [[Rcpp::export]]
Rcpp::NumericMatrix mcsimex_sim_cpp(
    Rcpp::NumericVector y_r,
    Rcpp::IntegerVector z_hat_r,
    Rcpp::NumericMatrix x_r,
    Rcpp::NumericMatrix Pi_r,
    int K,
    int dist_code,
    Rcpp::NumericVector lambda_r,
    int B,
    Rcpp::NumericVector wt_r,
    unsigned int seed) {

  int n = y_r.size();
  int r = x_r.ncol();
  int s = K - 1;
  int p = s + r;
  int n_lambda = lambda_r.size();

  // Map R objects to Eigen
  Eigen::Map<VectorXd> y(y_r.begin(), n);
  Eigen::Map<VectorXi> z_hat(z_hat_r.begin(), n);
  Eigen::Map<MatrixXd> x(x_r.begin(), n, r);
  Eigen::Map<MatrixXd> Pi(Pi_r.begin(), K, K);
  Eigen::Map<VectorXd> wt(wt_r.begin(), n);

  // Pre-compute Pi^lambda for each lambda level
  std::vector<MatrixXd> Pi_powers(n_lambda);
  for (int l = 0; l < n_lambda; l++) {
    Pi_powers[l] = mat_power(Pi, lambda_r[l]);
  }

  // Output: (n_lambda * B) rows x p cols
  Rcpp::NumericMatrix result(n_lambda * B, p);

  // RNG
  std::mt19937 rng(seed);

  for (int l = 0; l < n_lambda; l++) {
    const MatrixXd& Pi_lam = Pi_powers[l];

    for (int b = 0; b < B; b++) {
      // Resample z according to Pi^lambda
      VectorXi z_sim = resample_z(z_hat, Pi_lam, rng);

      // Build design matrix
      MatrixXd xi = build_design(z_sim, x, K);

      // Fit GLM via IRLS
      VectorXd beta = irls_fit(y, xi, dist_code, wt);

      // Store result
      int row = l * B + b;
      for (int j = 0; j < p; j++) {
        result(row, j) = beta(j);
      }
    }

    Rcpp::checkUserInterrupt();  // Allow user to cancel between lambda levels
  }

  return result;
}
