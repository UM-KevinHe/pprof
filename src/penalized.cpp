// Penalized logistic regression with unpenalized provider effects: lasso and
// group-lasso paths by coordinate descent on standardized (and, per group,
// orthonormalized) covariates, with one Newton step on the provider effects
// per cycle. Ported from grplasso 1.0 (pp.lasso, grp.lasso) with three
// changes: binomial trials N, a cyclic (deterministic) update order instead
// of a random permutation, and private loop variables in the OpenMP loop.
//
// Objective, with n rows and penalty multipliers m_j (lasso) or m_g (groups):
//   -(1/n) loglik(gamma, beta) + lambda * sum_j m_j |beta_j|
//   -(1/n) loglik(gamma, beta) + lambda * sum_g m_g ||beta_g||_2
// The covariate matrix handed in is already standardized / orthonormalized,
// so that every penalized column (lasso) or group (group lasso) satisfies
// Z'Z / n = I; the R side maps the estimates back to the original scale.
#define STRICT_R_HEADERS
#define ARMA_DONT_USE_OPENMP
#include <RcppArmadillo.h>
#include <cmath>
#include "pprof_omp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec rep_vec(const arma::vec &x, const arma::vec &each); // fixed_effect.cpp

namespace {

double soft(double z, double l) {
  if (z > l) return z - l;
  if (z < -l) return z + l;
  return 0.0;
}

// Binomial log-likelihood without the constant, in a form stable for large |eta|.
double loglik_eta(const vec &N, const vec &Y, const vec &eta) {
  return sum(Y % eta - N % (clamp(eta, 0.0, datum::inf) + log1p(exp(-abs(eta)))));
}

vec expit(const vec &eta) {
  return 1.0 / (1.0 + exp(-eta));
}

// One Newton step of every provider effect given eta (which contains the
// effects), optionally with a backtracking line search; gamma and eta are
// updated in place and gamma is kept within median +/- bound.
void gamma_step(const vec &N, const vec &Y, const vec &n_prov, const uvec &first,
                const uvec &last, vec &gamma, vec &eta, double bound, bool backtrack,
                int threads) {
  const int m = gamma.n_elem;
  vec p = expit(eta);
  vec resid = Y - N % p;
  vec pq = N % p % (1 - p);
  vec score(m), d_gamma(m);
  omp_set_num_threads(threads);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < m; i++) {
    double s = sum(resid.subvec(first(i), last(i)));
    double info = std::max(1e-20, sum(pq.subvec(first(i), last(i))));
    score(i) = s;
    d_gamma(i) = s / info;
  }
  double u = 1.0;
  if (backtrack) {
    const double s_armijo = 0.01, t_shrink = 0.8;
    double ll = loglik_eta(N, Y, eta);
    double slope = dot(score, d_gamma);
    vec eta_try = eta + rep_vec(u * d_gamma, n_prov);
    while (loglik_eta(N, Y, eta_try) - ll < s_armijo * u * slope && u > 1e-10) {
      u *= t_shrink;
      eta_try = eta + rep_vec(u * d_gamma, n_prov);
    }
  }
  vec gamma_new = gamma + u * d_gamma;
  gamma_new = clamp(gamma_new, median(gamma_new) - bound, median(gamma_new) + bound);
  eta += rep_vec(gamma_new - gamma, n_prov);
  gamma = gamma_new;
}

// -2 * loglik without the constant (the deviance up to the saturated model).
double deviance_eta(const vec &N, const vec &Y, const vec &eta) {
  return -2.0 * loglik_eta(N, Y, eta);
}

uvec provider_first(const vec &n_prov) {
  uvec first(n_prov.n_elem);
  unsigned int ind = 0;
  for (unsigned int i = 0; i < n_prov.n_elem; i++) {
    first(i) = ind;
    ind += n_prov(i);
  }
  return first;
}

} // namespace

// Lasso path. `penalized` marks the columns that carry a penalty (0/1);
// `multiplier` is the per-column penalty multiplier. With mm = true the
// working weights are the constant 0.25 * max(N) (majorization), otherwise
// the IRLS weights N p (1 - p).
// [[Rcpp::export]]
List lasso_path(const arma::vec &N, const arma::vec &Y, const arma::mat &Z,
                const arma::vec &n_prov, arma::vec gamma, arma::vec beta,
                const arma::uvec &penalized, const arma::vec &multiplier,
                const arma::vec &lambda, bool mm, bool backtrack, double bound, double tol,
                int max_iter, bool active_set, int threads) {
  const int n = Z.n_rows, p = Z.n_cols, m = n_prov.n_elem, n_lambda = lambda.n_elem;
  uvec first = provider_first(n_prov);
  uvec last = first + conv_to<uvec>::from(n_prov) - 1;
  const double v = 0.25 * N.max();
  vec zz(p);
  for (int j = 0; j < p; j++) zz(j) = dot(Z.col(j), Z.col(j));

  mat beta_path(p, n_lambda), gamma_path(m, n_lambda);
  vec dev(n_lambda), iters(n_lambda);
  vec eta = rep_vec(gamma, n_prov) + Z * beta;
  uvec active(p);
  for (int j = 0; j < p; j++) active(j) = (!active_set || penalized(j) == 0 || beta(j) != 0) ? 1 : 0;

  for (int l = 0; l < n_lambda; l++) {
    R_CheckUserInterrupt();
    const double lam = lambda(l);
    int iter = 0;
    while (true) {
      // inner cycles over the active columns
      while (iter < max_iter) {
        iter++;
        gamma_step(N, Y, n_prov, first, last, gamma, eta, bound, backtrack, threads);
        vec pr = expit(eta);
        vec w = mm ? vec(n, fill::value(v)) : vec(clamp(N % pr % (1 - pr), 1e-10, datum::inf));
        vec r = (Y - N % pr) / w;
        double max_change = 0.0;
        for (int j = 0; j < p; j++) {
          if (!active(j)) continue;
          double wz2 = mm ? v * zz(j) : dot(w % Z.col(j), Z.col(j));
          double wzr = mm ? v * dot(Z.col(j), r) : dot(w % Z.col(j), r);
          double b_new = wzr / wz2 + beta(j);
          if (penalized(j)) b_new = soft(b_new, n * lam * multiplier(j) / wz2);
          double change = b_new - beta(j);
          if (change != 0.0) {
            r -= change * Z.col(j);
            eta += change * Z.col(j);
            beta(j) = b_new;
            max_change = std::max(max_change, std::fabs(change));
          }
        }
        if (max_change < tol) break;
      }
      if (!active_set) break;
      // KKT check of the inactive columns at the current fit
      vec pr = expit(eta);
      vec w = mm ? vec(n, fill::value(v)) : vec(clamp(N % pr % (1 - pr), 1e-10, datum::inf));
      vec r = (Y - N % pr) / w;
      int added = 0;
      for (int j = 0; j < p; j++) {
        if (active(j)) continue;
        double wz2 = mm ? v * zz(j) : dot(w % Z.col(j), Z.col(j));
        double wzr = mm ? v * dot(Z.col(j), r) : dot(w % Z.col(j), r);
        if (soft(wzr / wz2, n * lam * multiplier(j) / wz2) != 0.0) {
          active(j) = 1;
          added++;
        }
      }
      if (added == 0 || iter >= max_iter) break;
    }
    beta_path.col(l) = beta;
    gamma_path.col(l) = gamma;
    dev(l) = deviance_eta(N, Y, eta);
    iters(l) = iter;
  }
  return List::create(_["beta"] = beta_path, _["gamma"] = gamma_path, _["deviance"] = dev,
                      _["iter"] = iters);
}

// Group-lasso path by majorization (working weight 0.25 * max(N)). Columns of
// a penalized group are contiguous, `group_first`/`group_last` give their
// 0-based range, `unpenalized` lists the columns without a penalty, and the
// penalized groups are orthonormalized (Z_g' Z_g / n = I).
// [[Rcpp::export]]
List group_lasso_path(const arma::vec &N, const arma::vec &Y, const arma::mat &Z,
                      const arma::vec &n_prov, arma::vec gamma, arma::vec beta,
                      const arma::uvec &unpenalized, const arma::uvec &group_first,
                      const arma::uvec &group_last, const arma::vec &multiplier,
                      const arma::vec &lambda, bool backtrack, double bound, double tol,
                      int max_iter, bool active_set, int threads) {
  const int n = Z.n_rows, p = Z.n_cols, m = n_prov.n_elem, n_lambda = lambda.n_elem;
  const int n_group = group_first.n_elem;
  uvec first = provider_first(n_prov);
  uvec last = first + conv_to<uvec>::from(n_prov) - 1;
  const double v = 0.25 * N.max();
  vec zz(p);
  for (int j = 0; j < p; j++) zz(j) = dot(Z.col(j), Z.col(j));

  mat beta_path(p, n_lambda), gamma_path(m, n_lambda);
  vec dev(n_lambda), iters(n_lambda);
  vec eta = rep_vec(gamma, n_prov) + Z * beta;
  uvec active(n_group);
  for (int g = 0; g < n_group; g++) {
    active(g) = (!active_set || beta(group_first(g)) != 0) ? 1 : 0;
  }
  // soft-thresholded norm of a group's working solution
  auto group_update = [&](int g, const vec &r, double lam, bool apply, vec &r_mut) -> double {
    int k = group_last(g) - group_first(g) + 1;
    vec u(k);
    for (int j = 0; j < k; j++) {
      int col = group_first(g) + j;
      u(j) = dot(Z.col(col), r) / n + (apply ? beta(col) : 0.0);
    }
    double norm_u = norm(u, 2);
    double len = soft(norm_u, lam * multiplier(g) / v);
    if (!apply) return len;
    double max_change = 0.0;
    if (len != 0.0 || beta(group_first(g)) != 0.0) {
      for (int j = 0; j < k; j++) {
        int col = group_first(g) + j;
        double b_new = len == 0.0 ? 0.0 : len * u(j) / norm_u;
        double change = b_new - beta(col);
        if (change != 0.0) {
          r_mut -= change * Z.col(col);
          eta += change * Z.col(col);
          beta(col) = b_new;
          max_change = std::max(max_change, std::fabs(change));
        }
      }
    }
    return max_change;
  };

  for (int l = 0; l < n_lambda; l++) {
    R_CheckUserInterrupt();
    const double lam = lambda(l);
    int iter = 0;
    while (true) {
      while (iter < max_iter) {
        iter++;
        gamma_step(N, Y, n_prov, first, last, gamma, eta, bound, backtrack, threads);
        vec pr = expit(eta);
        vec r = (Y - N % pr) / v;
        double max_change = 0.0;
        for (unsigned int u = 0; u < unpenalized.n_elem; u++) {
          int j = unpenalized(u);
          double change = dot(Z.col(j), r) / zz(j);
          if (change != 0.0) {
            r -= change * Z.col(j);
            eta += change * Z.col(j);
            beta(j) += change;
            max_change = std::max(max_change, std::fabs(change));
          }
        }
        for (int g = 0; g < n_group; g++) {
          if (!active(g)) continue;
          max_change = std::max(max_change, group_update(g, r, lam, true, r));
        }
        if (max_change < tol) break;
      }
      if (!active_set) break;
      vec pr = expit(eta);
      vec r = (Y - N % pr) / v;
      int added = 0;
      for (int g = 0; g < n_group; g++) {
        if (active(g)) continue;
        if (group_update(g, r, lam, false, r) != 0.0) {
          active(g) = 1;
          added++;
        }
      }
      if (added == 0 || iter >= max_iter) break;
    }
    beta_path.col(l) = beta;
    gamma_path.col(l) = gamma;
    dev(l) = deviance_eta(N, Y, eta);
    iters(l) = iter;
  }
  return List::create(_["beta"] = beta_path, _["gamma"] = gamma_path, _["deviance"] = dev,
                      _["iter"] = iters);
}
