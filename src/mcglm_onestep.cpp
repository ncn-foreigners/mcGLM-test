// src/mcglm_onestep.cpp
// TMB template for one-step joint estimation with misclassified covariates
// Supports GLMs (model_type=0) and multinomial logistic (model_type=1)
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {

  // ---- Common data ----
  DATA_INTEGER(model_type);   // 0=GLM, 1=multinomial logistic
  DATA_MATRIX(X);             // correctly observed covariates incl intercept (n x r)
  DATA_IVECTOR(z_hat);        // observed proxy covariate (n), values in {0,...,K-1}
  DATA_INTEGER(weights_fixed);// 1=use omega_data, 0=estimate via softmax
  DATA_VECTOR(omega_data);    // pre-computed mixture weights (K*K)
  DATA_INTEGER(K);            // number of categories for Z (2=binary)
  DATA_VECTOR(wt);            // frequency weights (n)

  // ---- Parameters ----
  PARAMETER_VECTOR(theta);

  int n = X.rows();
  int r = X.cols();
  int s = K - 1;

  // ---- Mixture weights ----
  int n_omega = K * K;
  vector<Type> omega(n_omega);

  if (model_type == 0) {
    // ===================================================================
    //  GLM model (Gaussian / Poisson / Binomial)
    // ===================================================================
    DATA_VECTOR(Y);              // response (n)
    DATA_INTEGER(dist_code);     // 1=Gaussian, 2=Poisson, 3=Binomial
    DATA_INTEGER(homoskedastic); // for Gaussian: 1=common sigma

    int d = s + r;

    // 1) Unpack regression coefficients
    vector<Type> b = theta.segment(0, d);
    vector<Type> gamma_vec(K);
    gamma_vec(0) = Type(0);
    for (int k = 1; k < K; k++) {
      gamma_vec(k) = b(k - 1);
    }
    vector<Type> alpha = b.segment(s, r);

    // 2) Mixture weights
    int idx_after_b = d;
    if (weights_fixed == 1) {
      for (int q = 0; q < n_omega; q++) omega(q) = omega_data(q);
    } else {
      int n_free = n_omega - 1;
      vector<Type> vraw = theta.segment(idx_after_b, n_free);
      idx_after_b += n_free;
      vector<Type> expv(n_omega);
      for (int q = 0; q < n_free; q++) expv(q) = exp(vraw(q));
      expv(n_omega - 1) = Type(1);
      Type sumv = expv.sum();
      for (int q = 0; q < n_omega; q++) omega(q) = expv(q) / sumv;
    }

    // 3) Gaussian error parameters
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
      int j = z_hat(i);

      Type eta_base = Type(0);
      for (int q = 0; q < r; q++) eta_base += alpha(q) * X(i, q);

      vector<Type> log_terms(K);
      for (int l = 0; l < K; l++) {
        Type w_jl = omega(j * K + l);
        Type eta_l = eta_base + gamma_vec(l);
        Type log_f;

        switch(dist_code) {
          case 1: {
            Type sig = (homoskedastic == 1) ? sigma0 : ((l == 0) ? sigma0 : sigma1);
            log_f = dnorm(yi, eta_l, sig, true);
            break;
          }
          case 2: {
            Type mu_l = exp(eta_l);
            log_f = dpois(yi, mu_l, true);
            break;
          }
          case 3: {
            Type p_l = invlogit(eta_l);
            log_f = yi * log(p_l + Type(1e-20)) + (Type(1) - yi) * log(Type(1) - p_l + Type(1e-20));
            break;
          }
          default:
            error("Unknown dist_code %d", dist_code);
            log_f = Type(0);
        }

        if (w_jl > Type(1e-20)) {
          log_terms(l) = log(w_jl) + log_f;
        } else {
          log_terms(l) = Type(-1000);
        }
      }

      Type max_log = log_terms(0);
      for (int l = 1; l < K; l++) {
        if (log_terms(l) > max_log) max_log = log_terms(l);
      }
      Type sum_exp = Type(0);
      for (int l = 0; l < K; l++) sum_exp += exp(log_terms(l) - max_log);
      nll -= wi * (max_log + log(sum_exp));
    }

    return nll;

  } else {
    // ===================================================================
    //  Multinomial logistic model (model_type == 1)
    // ===================================================================
    DATA_IVECTOR(Y_int);   // categorical response (n), values in {0,...,J-1}
    DATA_INTEGER(J);       // number of response categories

    // Parameter layout: (J-1) blocks, each of size (s + r)
    // Block j (j=0,...,J-2 for response categories 1,...,J-1):
    //   gamma_{j+1,1}, ..., gamma_{j+1,K-1}, alpha_{j+1,0}, ..., alpha_{j+1,r-1}
    int block_size = s + r;
    int d = (J - 1) * block_size;

    // Unpack parameters into gamma[j][k] and alpha[j][q] matrices
    // gamma[j][k] for j=1,...,J-1, k=0,...,K-1 (gamma[j][0] = 0)
    // alpha[j][q] for j=1,...,J-1, q=0,...,r-1
    // Response category 0 is baseline (all zeros)
    matrix<Type> gamma_mat(J, K);   // gamma_mat(j, k)
    matrix<Type> alpha_mat(J, r);   // alpha_mat(j, q)

    // Baseline response (j=0): all zeros
    for (int k = 0; k < K; k++) gamma_mat(0, k) = Type(0);
    for (int q = 0; q < r; q++) alpha_mat(0, q) = Type(0);

    // Fill from theta
    for (int jj = 0; jj < J - 1; jj++) {
      int offset = jj * block_size;
      // gamma_{jj+1, 0} = 0 (baseline Z category)
      gamma_mat(jj + 1, 0) = Type(0);
      for (int k = 1; k < K; k++) {
        gamma_mat(jj + 1, k) = theta(offset + k - 1);
      }
      for (int q = 0; q < r; q++) {
        alpha_mat(jj + 1, q) = theta(offset + s + q);
      }
    }

    // Mixture weights
    int idx_after_b = d;
    if (weights_fixed == 1) {
      for (int q = 0; q < n_omega; q++) omega(q) = omega_data(q);
    } else {
      int n_free = n_omega - 1;
      vector<Type> vraw = theta.segment(idx_after_b, n_free);
      idx_after_b += n_free;
      vector<Type> expv(n_omega);
      for (int q = 0; q < n_free; q++) expv(q) = exp(vraw(q));
      expv(n_omega - 1) = Type(1);
      Type sumv = expv.sum();
      for (int q = 0; q < n_omega; q++) omega(q) = expv(q) / sumv;
    }

    // Build negative log-likelihood
    Type nll = Type(0);

    for (int i = 0; i < n; i++) {
      int yi = Y_int(i);
      Type wi = wt(i);
      int zh = z_hat(i);

      // Compute alpha_j' * x_i for each response category j
      vector<Type> alpha_x(J);
      alpha_x(0) = Type(0);  // baseline response
      for (int jj = 1; jj < J; jj++) {
        Type ax = Type(0);
        for (int q = 0; q < r; q++) ax += alpha_mat(jj, q) * X(i, q);
        alpha_x(jj) = ax;
      }

      // For each latent Z = l, compute log P(Y = yi | Z = l, x_i)
      vector<Type> log_mixture_terms(K);
      for (int l = 0; l < K; l++) {
        Type w_jl = omega(zh * K + l);

        // Compute eta_{j,l} for all response categories j
        // eta_{j,l} = gamma_{j,l} + alpha_j' x_i
        // Then P(Y=j | Z=l) = exp(eta_{j,l}) / sum_{j'} exp(eta_{j',l})
        // Use log-sum-exp for the denominator

        // Find max eta for numerical stability
        Type max_eta = Type(0);  // eta for baseline (j=0) is 0
        for (int jj = 1; jj < J; jj++) {
          Type eta_jl = gamma_mat(jj, l) + alpha_x(jj);
          if (eta_jl > max_eta) max_eta = eta_jl;
        }

        // log(sum_j exp(eta_{j,l}))
        Type sum_exp_eta = exp(Type(0) - max_eta);  // j=0: eta=0
        for (int jj = 1; jj < J; jj++) {
          Type eta_jl = gamma_mat(jj, l) + alpha_x(jj);
          sum_exp_eta += exp(eta_jl - max_eta);
        }
        Type log_denom = max_eta + log(sum_exp_eta);

        // eta for the observed response category yi
        Type eta_yi_l;
        if (yi == 0) {
          eta_yi_l = Type(0);
        } else {
          eta_yi_l = gamma_mat(yi, l) + alpha_x(yi);
        }

        Type log_prob = eta_yi_l - log_denom;

        if (w_jl > Type(1e-20)) {
          log_mixture_terms(l) = log(w_jl) + log_prob;
        } else {
          log_mixture_terms(l) = Type(-1000);
        }
      }

      // log-sum-exp over latent Z categories
      Type max_log = log_mixture_terms(0);
      for (int l = 1; l < K; l++) {
        if (log_mixture_terms(l) > max_log) max_log = log_mixture_terms(l);
      }
      Type sum_exp = Type(0);
      for (int l = 0; l < K; l++) {
        sum_exp += exp(log_mixture_terms(l) - max_log);
      }
      nll -= wi * (max_log + log(sum_exp));
    }

    return nll;
  }
}
