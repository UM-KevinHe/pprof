// Penalized survival regression for variable selection with the provider
// effects kept: (1) the Cox model stratified by provider (the provider
// baselines absorb the provider effects) with a lasso or group-lasso penalty
// on the covariates, and (2) the discrete-time hazard model
//   logit h(k | i, j) = alpha_k + gamma_i + Z_j' beta
// with unpenalized period effects alpha and provider effects gamma and a
// lasso penalty on beta. Ported from grplasso 1.0 (Strat.cox, pp.DiscSurv)
// with these changes: counting-process rows (start, stop], case weights and
// offsets in the Cox model; a cyclic (deterministic) update order; the
// discrete-time provider effects identified by a zero median (grplasso fixed
// the first provider at zero); convergence judged on every parameter block.
//
// Objectives, with n rows (Cox: the weighted partial likelihood):
//   -(1/n) loglik(beta) + lambda * sum_g m_g ||beta_g||_2        (Cox)
//   -(1/n) loglik(alpha, gamma, beta) + lambda * sum_j m_j |beta_j| (discrete)
// The covariate matrix is standardized / orthonormalized by the R side.
#define STRICT_R_HEADERS
#define ARMA_DONT_USE_OPENMP
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec rep_vec(const arma::vec &x, const arma::vec &each); // fixed_effect.cpp

namespace {

double soft_threshold(double z, double l) {
  if (z > l) return z - l;
  if (z < -l) return z + l;
  return 0.0;
}

// ---------------------------------------------------------------------------
// Cox model stratified by provider.
// Rows are grouped by stratum (row_first/row_last, 0-based, inclusive). Each
// stratum has a grid of its distinct event times (grid_first offset into the
// concatenated grid, grid_count points, d_grid the weighted events at each);
// k_stop(i) and k_start(i) count the stratum's grid points <= stop_i and
// <= start_i (0 = none), so that row i is at risk at grid point k exactly
// when k_start(i) < k <= k_stop(i).
// ---------------------------------------------------------------------------
struct CoxData {
  const vec &event, &w;
  const uvec &k_stop, &k_start, &row_first, &row_last, &grid_first, &grid_count;
  const vec &d_grid;
};

// Working residuals r_i = w_i (delta_i - exp(eta_i) H_i), with H_i the
// stratum's Breslow cumulative hazard over (start_i, stop_i], and the
// deviance -2 sum w_i delta_i (eta_i - log RS(stop_i)).
double cox_working(const CoxData &d, const vec &eta, vec &r) {
  const int n_strata = d.row_first.n_elem;
  const double shift = eta.max();
  vec whaz = d.w % exp(eta - shift);
  double deviance = 0.0;
  for (int s = 0; s < n_strata; s++) {
    const int k_n = d.grid_count(s);
    const unsigned int g0 = d.grid_first(s);
    if (k_n == 0) {
      for (unsigned int i = d.row_first(s); i <= d.row_last(s); i++) r(i) = 0.0;
      continue;
    }
    vec s_stop(k_n + 2, fill::zeros), s_start(k_n + 2, fill::zeros);
    for (unsigned int i = d.row_first(s); i <= d.row_last(s); i++) {
      s_stop(d.k_stop(i)) += whaz(i);
      s_start(d.k_start(i)) += whaz(i);
    }
    for (int k = k_n - 1; k >= 1; k--) {
      s_stop(k) += s_stop(k + 1);
      s_start(k) += s_start(k + 1);
    }
    vec cumhaz(k_n + 1, fill::zeros), log_rs(k_n + 1, fill::zeros);
    for (int k = 1; k <= k_n; k++) {
      double rs = s_stop(k) - s_start(k);
      cumhaz(k) = cumhaz(k - 1) + d.d_grid(g0 + k - 1) / rs;
      log_rs(k) = std::log(rs);
    }
    for (unsigned int i = d.row_first(s); i <= d.row_last(s); i++) {
      double h = cumhaz(d.k_stop(i)) - cumhaz(d.k_start(i));
      r(i) = d.w(i) * (d.event(i) - std::exp(eta(i) - shift) * h);
      if (d.event(i) > 0) {
        deviance -= 2.0 * d.w(i) * d.event(i) * (eta(i) - shift - log_rs(d.k_stop(i)));
      }
    }
  }
  return deviance;
}

// ---------------------------------------------------------------------------
// Discrete-time hazard model. Row j is at risk in periods 1..t_j and has the
// event in period t_j when delta_j = 1; eta_j = gamma_i + Z_j' beta + offset.
// ---------------------------------------------------------------------------
inline double expit(double x) { return 1.0 / (1.0 + std::exp(-x)); }

// Per row: the sum of the hazards over its periods and of their variances.
void discrete_sums(const uvec &t, const vec &alpha, const vec &eta, vec &p_sum, vec &pq_sum) {
  const int n = t.n_elem;
  for (int j = 0; j < n; j++) {
    double ps = 0.0, pqs = 0.0;
    for (unsigned int k = 0; k < t(j); k++) {
      double p = expit(alpha(k) + eta(j));
      ps += p;
      pqs += p * (1.0 - p);
    }
    p_sum(j) = ps;
    pq_sum(j) = pqs;
  }
}

double discrete_loglik(const uvec &t, const vec &event, const vec &alpha, const vec &eta) {
  const int n = t.n_elem;
  double ll = 0.0;
  for (int j = 0; j < n; j++) {
    for (unsigned int k = 0; k < t(j); k++) {
      double x = alpha(k) + eta(j);
      ll -= std::max(x, 0.0) + std::log1p(std::exp(-std::fabs(x)));
      if (k == t(j) - 1 && event(j) > 0) ll += x;
    }
  }
  return ll;
}

// Newton step on the period effects, with optional backtracking; alpha is
// kept within median +/- bound.
void alpha_step(const uvec &t, const vec &event, const vec &d_period, vec &alpha,
                  const vec &eta, double bound, bool backtrack) {
  const int n = t.n_elem, n_alpha = alpha.n_elem;
  vec score = d_period, info(n_alpha, fill::zeros);
  for (int j = 0; j < n; j++) {
    for (unsigned int k = 0; k < t(j); k++) {
      double p = expit(alpha(k) + eta(j));
      score(k) -= p;
      info(k) += p * (1.0 - p);
    }
  }
  vec d_alpha = score / clamp(info, 1e-10, datum::inf);
  double u = 1.0;
  if (backtrack) {
    double ll = discrete_loglik(t, event, alpha, eta);
    double slope = dot(score, d_alpha);
    while (discrete_loglik(t, event, alpha + u * d_alpha, eta) - ll < 0.01 * u * slope &&
           u > 1e-10) {
      u *= 0.8;
    }
  }
  vec alpha_new = alpha + u * d_alpha;
  alpha = clamp(alpha_new, median(alpha_new) - bound, median(alpha_new) + bound);
}

// Newton step on the provider effects given the per-row hazard sums, with
// optional backtracking; the effects are then centred at a zero median (the
// shift moves into alpha) and kept within +/- bound.
void gamma_step(const uvec &t, const vec &event, const vec &n_prov, const uvec &first,
                  const uvec &last, vec &gamma, vec &alpha, vec &eta, double bound,
                  bool backtrack) {
  const int m = gamma.n_elem, n = t.n_elem;
  vec p_sum(n), pq_sum(n);
  discrete_sums(t, alpha, eta, p_sum, pq_sum);
  vec score(m), d_gamma(m);
  for (int i = 0; i < m; i++) {
    double s = 0.0, info = 0.0;
    for (unsigned int j = first(i); j <= last(i); j++) {
      s += event(j) - p_sum(j);
      info += pq_sum(j);
    }
    score(i) = s;
    d_gamma(i) = s / std::max(info, 1e-10);
  }
  double u = 1.0;
  if (backtrack) {
    double ll = discrete_loglik(t, event, alpha, eta);
    double slope = dot(score, d_gamma);
    vec eta_try = eta + rep_vec(u * d_gamma, n_prov);
    while (discrete_loglik(t, event, alpha, eta_try) - ll < 0.01 * u * slope && u > 1e-10) {
      u *= 0.8;
      eta_try = eta + rep_vec(u * d_gamma, n_prov);
    }
  }
  vec gamma_new = gamma + u * d_gamma;
  double med = median(gamma_new);
  gamma_new -= med;
  alpha += med;
  gamma_new = clamp(gamma_new, -bound, bound);
  eta += rep_vec(gamma_new - gamma, n_prov);
  gamma = gamma_new;
}

uvec block_first(const vec &sizes) {
  uvec first(sizes.n_elem);
  unsigned int ind = 0;
  for (unsigned int i = 0; i < sizes.n_elem; i++) {
    first(i) = ind;
    ind += sizes(i);
  }
  return first;
}

} // namespace

// Working residuals of the stratified Cox model at a linear predictor (for
// the largest penalty of a path) and its deviance.
// [[Rcpp::export]]
List cox_residuals(const arma::vec &event, const arma::vec &w, const arma::uvec &k_stop,
                   const arma::uvec &k_start, const arma::uvec &row_first,
                   const arma::uvec &row_last, const arma::uvec &grid_first,
                   const arma::uvec &grid_count, const arma::vec &d_grid, const arma::vec &eta) {
  CoxData d{event, w, k_stop, k_start, row_first, row_last, grid_first, grid_count, d_grid};
  vec r(event.n_elem, fill::zeros);
  double dev = cox_working(d, eta, r);
  return List::create(_["residual"] = r, _["deviance"] = dev);
}

// Deviance of the stratified Cox model at several linear predictors (columns
// of eta), for cross-validation.
// [[Rcpp::export]]
arma::vec cox_deviance(const arma::vec &event, const arma::vec &w, const arma::uvec &k_stop,
                       const arma::uvec &k_start, const arma::uvec &row_first,
                       const arma::uvec &row_last, const arma::uvec &grid_first,
                       const arma::uvec &grid_count, const arma::vec &d_grid,
                       const arma::mat &eta) {
  CoxData d{event, w, k_stop, k_start, row_first, row_last, grid_first, grid_count, d_grid};
  vec r(event.n_elem, fill::zeros), dev(eta.n_cols);
  for (unsigned int l = 0; l < eta.n_cols; l++) dev(l) = cox_working(d, eta.col(l), r);
  return dev;
}

// Group-lasso path of the stratified Cox model (the lasso is the case of
// singleton groups). Columns of a penalized group are contiguous
// (group_first/group_last, 0-based), `unpenalized` lists the columns without
// a penalty; the groups are orthonormalized (Z_g' Z_g / n = I). The working
// weight is 1 (grplasso's majorization).
// [[Rcpp::export]]
List cox_lasso_path(const arma::vec &event, const arma::vec &w, const arma::vec &offset,
                    const arma::mat &Z, const arma::uvec &k_stop, const arma::uvec &k_start,
                    const arma::uvec &row_first, const arma::uvec &row_last,
                    const arma::uvec &grid_first, const arma::uvec &grid_count,
                    const arma::vec &d_grid, arma::vec beta, const arma::uvec &unpenalized,
                    const arma::uvec &group_first, const arma::uvec &group_last,
                    const arma::vec &multiplier, const arma::vec &lambda, double tol,
                    int max_iter, bool active_set) {
  CoxData d{event, w, k_stop, k_start, row_first, row_last, grid_first, grid_count, d_grid};
  const int n = Z.n_rows, p = Z.n_cols, n_lambda = lambda.n_elem, n_group = group_first.n_elem;
  vec zz(p);
  for (int j = 0; j < p; j++) zz(j) = dot(Z.col(j), Z.col(j));
  mat beta_path(p, n_lambda);
  vec dev(n_lambda), iters(n_lambda);
  vec eta = Z * beta + offset;
  vec r(n, fill::zeros);
  uvec active(n_group);
  for (int g = 0; g < n_group; g++) {
    active(g) = (!active_set || beta(group_first(g)) != 0) ? 1 : 0;
  }
  auto group_update = [&](int g, double lam, bool apply) -> double {
    int k = group_last(g) - group_first(g) + 1;
    vec u(k);
    for (int j = 0; j < k; j++) {
      int col = group_first(g) + j;
      u(j) = dot(Z.col(col), r) / n + (apply ? beta(col) : 0.0);
    }
    double norm_u = norm(u, 2);
    double len = soft_threshold(norm_u, lam * multiplier(g));
    if (!apply) return len;
    double max_change = 0.0;
    if (len != 0.0 || beta(group_first(g)) != 0.0) {
      for (int j = 0; j < k; j++) {
        int col = group_first(g) + j;
        double b_new = len == 0.0 ? 0.0 : len * u(j) / norm_u;
        double change = b_new - beta(col);
        if (change != 0.0) {
          r -= change * Z.col(col);
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
    double deviance = 0.0;
    while (true) {
      while (iter < max_iter) {
        iter++;
        deviance = cox_working(d, eta, r);
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
          max_change = std::max(max_change, group_update(g, lam, true));
        }
        if (max_change < tol) break;
      }
      if (!active_set) break;
      deviance = cox_working(d, eta, r);
      int added = 0;
      for (int g = 0; g < n_group; g++) {
        if (active(g)) continue;
        if (group_update(g, lam, false) != 0.0) {
          active(g) = 1;
          added++;
        }
      }
      if (added == 0 || iter >= max_iter) break;
    }
    beta_path.col(l) = beta;
    dev(l) = cox_working(d, eta, r);
    iters(l) = iter;
  }
  return List::create(_["beta"] = beta_path, _["deviance"] = dev, _["iter"] = iters);
}

// Hazards of every row in every period at risk, row by row (periods 1..t_j),
// and the per-row sums of the hazards.
// [[Rcpp::export]]
List discrete_hazards(const arma::uvec &t, const arma::vec &alpha, const arma::vec &eta) {
  const int n = t.n_elem;
  vec p_sum(n), pq_sum(n);
  discrete_sums(t, alpha, eta, p_sum, pq_sum);
  vec hazard(sum(t));
  unsigned int pos = 0;
  for (int j = 0; j < n; j++) {
    for (unsigned int k = 0; k < t(j); k++) hazard(pos++) = expit(alpha(k) + eta(j));
  }
  return List::create(_["hazard"] = hazard, _["expected"] = p_sum);
}

// Lasso path of the discrete-time hazard model. t_j is the row's exit period
// (1-based), d_period the events per period, n_prov the rows per provider
// (rows grouped by provider), `penalized` marks the columns with a penalty.
// With mm = true the working weights are 0.25 t_j (majorization), otherwise
// the IRLS weights sum_k p (1 - p). A single infinite lambda fits the
// unpenalized model (alpha, gamma and the unpenalized covariates).
// [[Rcpp::export]]
List discrete_lasso_path(const arma::uvec &t, const arma::vec &event, const arma::mat &Z,
                         const arma::vec &offset, const arma::vec &n_prov,
                         const arma::vec &d_period, arma::vec alpha, arma::vec gamma,
                         arma::vec beta, const arma::uvec &penalized,
                         const arma::vec &multiplier, const arma::vec &lambda, bool mm,
                         bool backtrack, double bound, double tol, int max_iter,
                         bool active_set) {
  const int n = Z.n_rows, p = Z.n_cols, m = n_prov.n_elem, n_lambda = lambda.n_elem;
  const int n_alpha = alpha.n_elem;
  uvec first = block_first(n_prov);
  uvec last = first + conv_to<uvec>::from(n_prov) - 1;
  mat beta_path(p, n_lambda), alpha_path(n_alpha, n_lambda), gamma_path(m, n_lambda);
  vec dev(n_lambda), iters(n_lambda);
  vec eta = rep_vec(gamma, n_prov) + Z * beta + offset;
  vec p_sum(n), pq_sum(n), w(n), r(n);
  uvec active(p);
  for (int j = 0; j < p; j++) active(j) = (!active_set || penalized(j) == 0 || beta(j) != 0) ? 1 : 0;

  auto working = [&]() {
    discrete_sums(t, alpha, eta, p_sum, pq_sum);
    if (mm) {
      for (int j = 0; j < n; j++) w(j) = 0.25 * t(j);
    } else {
      w = clamp(pq_sum, 1e-10, datum::inf);
    }
    r = (event - p_sum) / w;
  };

  for (int l = 0; l < n_lambda; l++) {
    R_CheckUserInterrupt();
    const double lam = lambda(l);
    int iter = 0;
    while (true) {
      while (iter < max_iter) {
        iter++;
        // the change of a block is measured over the whole iteration: the
        // alpha step and the re-centring of gamma move along the flat
        // direction (alpha + c, gamma - c), which is no change at all
        vec alpha_start = alpha, gamma_start = gamma;
        alpha_step(t, event, d_period, alpha, eta, bound, backtrack);
        gamma_step(t, event, n_prov, first, last, gamma, alpha, eta, bound, backtrack);
        double max_change = std::max(abs(alpha - alpha_start).max(),
                                     abs(gamma - gamma_start).max());
        working();
        for (int j = 0; j < p; j++) {
          if (!active(j)) continue;
          double wz2 = dot(w % Z.col(j), Z.col(j));
          double wzr = dot(w % Z.col(j), r);
          double b_new = wzr / wz2 + beta(j);
          if (penalized(j)) b_new = soft_threshold(b_new, n * lam * multiplier(j) / wz2);
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
      working();
      int added = 0;
      for (int j = 0; j < p; j++) {
        if (active(j)) continue;
        double wz2 = dot(w % Z.col(j), Z.col(j));
        double wzr = dot(w % Z.col(j), r);
        if (soft_threshold(wzr / wz2, n * lam * multiplier(j) / wz2) != 0.0) {
          active(j) = 1;
          added++;
        }
      }
      if (added == 0 || iter >= max_iter) break;
    }
    beta_path.col(l) = beta;
    alpha_path.col(l) = alpha;
    gamma_path.col(l) = gamma;
    dev(l) = -2.0 * discrete_loglik(t, event, alpha, eta);
    iters(l) = iter;
  }
  return List::create(_["beta"] = beta_path, _["alpha"] = alpha_path, _["gamma"] = gamma_path,
                      _["deviance"] = dev, _["iter"] = iters);
}
