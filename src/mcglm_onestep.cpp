// src/mcglm_onestep.cpp
// TMB template for one-step joint estimation of GLMs with misclassified covariates
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {

  // ---- Data ----
  DATA_VECTOR(Y);             // response (n)
  DATA_MATRIX(X);             // correctly observed covariates incl intercept (n x r)
  DATA_IVECTOR(z_hat);        // observed proxy covariate (n), values in {0,...,K-1}
  DATA_INTEGER(dist_code);    // 1=Gaussian, 2=Poisson, 3=Binomial
  DATA_INTEGER(weights_fixed);// 1=use omega_data, 0=estimate via softmax
  DATA_VECTOR(omega_data);    // pre-computed mixture weights (K*K for multi, 4 for binary)
  DATA_INTEGER(homoskedastic);// for Gaussian: 1=common sigma, 0=per-component sigma
  DATA_INTEGER(K);            // number of categories (2=binary)
  DATA_VECTOR(wt);            // frequency weights (n)

  // ---- Parameters ----
  PARAMETER_VECTOR(theta);

  int n = Y.size();
  int r = X.cols();
  int s = K - 1;      // number of dummy variables (baseline = 0)
  int d = s + r;       // total regression coefficients: gamma_1,...,gamma_{K-1}, alpha_0,...,alpha_{r-1}

  // 1) Unpack regression coefficients
  // Layout: gamma_1, ..., gamma_{K-1}, alpha_0, ..., alpha_{r-1}
  vector<Type> b = theta.segment(0, d);
  vector<Type> gamma_vec(K);
  gamma_vec(0) = Type(0);  // baseline
  for (int k = 1; k < K; k++) {
    gamma_vec(k) = b(k - 1);
  }
  vector<Type> alpha = b.segment(s, r);

  // 2) Mixture weights: omega[j * K + l] = P(z_hat = j | Z = l) * P(Z = l)
  int n_omega = K * K;
  vector<Type> omega(n_omega);

  int idx_after_b = d;
  if (weights_fixed == 1) {
    for (int q = 0; q < n_omega; q++) {
      omega(q) = omega_data(q);
    }
  } else {
    // Estimate via softmax: K*K - 1 free parameters
    int n_free = n_omega - 1;
    vector<Type> vraw = theta.segment(idx_after_b, n_free);
    idx_after_b += n_free;
    // softmax
    vector<Type> expv(n_omega);
    for (int q = 0; q < n_free; q++) {
      expv(q) = exp(vraw(q));
    }
    expv(n_omega - 1) = Type(1);
    Type sumv = expv.sum();
    for (int q = 0; q < n_omega; q++) {
      omega(q) = expv(q) / sumv;
    }
  }

  // 3) Gaussian error parameters (only for dist_code == 1)
  Type log_sigma0 = Type(0);
  Type log_sigma1 = Type(0);
  if (dist_code == 1) {
    log_sigma0 = theta(idx_after_b);
    idx_after_b++;
    if (homoskedastic == 1) {
      log_sigma1 = log_sigma0;
    } else {
      log_sigma1 = theta(idx_after_b);
      idx_after_b++;
    }
  }
  Type sigma0 = exp(log_sigma0);
  Type sigma1 = exp(log_sigma1);

  // 4) Build negative log-likelihood
  Type nll = Type(0);

  for (int i = 0; i < n; i++) {
    Type yi = Y(i);
    Type wi = wt(i);
    int j = z_hat(i);  // observed proxy category

    // Linear predictor base: alpha' * x_i
    Type eta_base = Type(0);
    for (int q = 0; q < r; q++) {
      eta_base += alpha(q) * X(i, q);
    }

    // Sum over true categories l = 0, ..., K-1
    // log( sum_l omega[j,l] * f(y_i | Z=l, psi) )
    // Use log-sum-exp for numerical stability

    vector<Type> log_terms(K);
    for (int l = 0; l < K; l++) {
      Type w_jl = omega(j * K + l);
      Type eta_l = eta_base + gamma_vec(l);
      Type log_f;

      switch(dist_code) {
        case 1: { // Gaussian
          // Use component-specific sigma for heteroskedastic
          Type sig;
          if (homoskedastic == 1) {
            sig = sigma0;
          } else {
            // sigma0 for Z=0, sigma1 for Z!=0
            sig = (l == 0) ? sigma0 : sigma1;
          }
          log_f = dnorm(yi, eta_l, sig, true);
          break;
        }
        case 2: { // Poisson
          Type mu_l = exp(eta_l);
          log_f = dpois(yi, mu_l, true);
          break;
        }
        case 3: { // Binomial (Bernoulli with logit link)
          Type p_l = invlogit(eta_l);
          log_f = yi * log(p_l + Type(1e-20)) + (Type(1) - yi) * log(Type(1) - p_l + Type(1e-20));
          break;
        }
        default:
          error("Unknown dist_code %d", dist_code);
          log_f = Type(0); // unreachable
      }

      // log(w_jl * f) = log(w_jl) + log_f
      if (w_jl > Type(1e-20)) {
        log_terms(l) = log(w_jl) + log_f;
      } else {
        log_terms(l) = Type(-1000);  // effectively zero contribution
      }
    }

    // log-sum-exp
    Type max_log = log_terms(0);
    for (int l = 1; l < K; l++) {
      if (log_terms(l) > max_log) max_log = log_terms(l);
    }
    Type sum_exp = Type(0);
    for (int l = 0; l < K; l++) {
      sum_exp += exp(log_terms(l) - max_log);
    }
    nll -= wi * (max_log + log(sum_exp));
  }

  return nll;
}
