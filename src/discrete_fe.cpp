// Discrete-time hazard model with fixed provider effects, unpenalized:
//   logit h(k | i, j) = alpha_k + gamma_i + Z_j' beta,   k = 1, ..., t_j,
// fitted on one row per subject. The person-periods of a subject are looped
// over and never stored. The periods are the distinct event times (a period
// without events has a zero hazard and adds nothing to the likelihood), t_j is
// the number of them at or before the subject's time, and alpha_1 = 0: the
// provider effects carry the first period's level.
//
// The Newton step uses the block structure of the information, as SerBIN does
// for the logistic model: the provider block is diagonal and
// theta = (alpha_2, ..., alpha_K, beta) is small, so the step is solved
// through the Schur complement of the provider block.
//
// Armadillo's internal OpenMP is switched off (sums do not depend on the
// number of cores); the explicit loop of the direct standardization is
// parallel over providers.
#define STRICT_R_HEADERS
#define ARMA_DONT_USE_OPENMP
#include <RcppArmadillo.h>
#include <cmath>
#include "pprof_omp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace {

inline double dfe_expit(double x) { return 1.0 / (1.0 + std::exp(-x)); }

// Provider (0-based) of every row, rows grouped by provider.
uvec dfe_rows_provider(const vec &n_prov) {
  uvec out((uword) accu(n_prov));
  uword pos = 0;
  for (uword i = 0; i < n_prov.n_elem; i++) {
    for (uword r = 0; r < (uword) n_prov(i); r++) out(pos++) = i;
  }
  return out;
}

double dfe_loglik(const uvec &t, const vec &event, const vec &alpha, const vec &eta) {
  double ll = 0.0;
  for (uword j = 0; j < t.n_elem; j++) {
    for (uword k = 0; k < t(j); k++) {
      const double x = alpha(k) + eta(j);
      ll -= std::max(x, 0.0) + std::log1p(std::exp(-std::fabs(x)));
      if (k + 1 == t(j) && event(j) > 0) ll += x;
    }
  }
  return ll;
}

// Score and information of (gamma, theta), theta = (alpha_2..alpha_K, beta).
// With clip = true the hazards are kept within [1e-10, 1 - 1e-10] in the
// information (as the variance routines of the logistic model do).
struct DfeBlocks {
  vec u_gamma, i_gamma, u_theta, H, W;
  mat i_gamma_theta, i_theta;
};

DfeBlocks dfe_blocks(const uvec &t, const vec &event, const mat &Z, const uvec &prov, uword m,
                     const vec &alpha, const vec &eta, bool clip) {
  const uword n = Z.n_rows, p = Z.n_cols, K = alpha.n_elem, q = K - 1 + p;
  const mat Zt = Z.t();
  DfeBlocks b;
  b.H.zeros(n);
  b.W.zeros(n);
  vec u_alpha(K, fill::zeros), a_alpha(K, fill::zeros);
  mat b_alpha(p, K, fill::zeros);  // column k: sum of w_jk Z_j over the subjects at risk
  mat c(m, K, fill::zeros);        // sum of w_jk by provider and period
  mat i_gamma_beta(p, m, fill::zeros);
  b.u_gamma.zeros(m);
  b.i_gamma.zeros(m);
  for (uword j = 0; j < n; j++) {
    const uword i = prov(j);
    double hs = 0.0, ws = 0.0;
    const double *zj = Zt.colptr(j);
    for (uword k = 0; k < t(j); k++) {
      const double h = dfe_expit(alpha(k) + eta(j));
      hs += h;
      u_alpha(k) -= h;
      const double hc = clip ? std::min(std::max(h, 1e-10), 1.0 - 1e-10) : h;
      const double w = hc * (1.0 - hc);
      ws += w;
      a_alpha(k) += w;
      c(i, k) += w;
      double *bk = b_alpha.colptr(k);
      for (uword l = 0; l < p; l++) bk[l] += w * zj[l];
    }
    if (event(j) > 0 && t(j) > 0) u_alpha(t(j) - 1) += 1.0;
    b.H(j) = hs;
    b.W(j) = ws;
    b.u_gamma(i) += event(j) - hs;
    b.i_gamma(i) += ws;
    double *gi = i_gamma_beta.colptr(i);
    for (uword l = 0; l < p; l++) gi[l] += ws * zj[l];
  }
  const vec resid = event - b.H;
  b.i_gamma_theta.zeros(m, q);
  b.u_theta.zeros(q);
  b.i_theta.zeros(q, q);
  if (K > 1) {
    b.i_gamma_theta.cols(0, K - 2) = c.cols(1, K - 1);
    b.u_theta.subvec(0, K - 2) = u_alpha.subvec(1, K - 1);
    b.i_theta.submat(0, 0, K - 2, K - 2) = diagmat(a_alpha.subvec(1, K - 1));
  }
  if (p > 0) {
    b.i_gamma_theta.cols(K - 1, q - 1) = i_gamma_beta.t();
    b.u_theta.subvec(K - 1, q - 1) = Z.t() * resid;
    b.i_theta.submat(K - 1, K - 1, q - 1, q - 1) = Z.t() * (Z.each_col() % b.W);
    if (K > 1) {
      mat ab = b_alpha.cols(1, K - 1);  // p x (K - 1)
      b.i_theta.submat(0, K - 1, K - 2, q - 1) = ab.t();
      b.i_theta.submat(K - 1, 0, q - 1, K - 2) = ab;
    }
  }
  return b;
}

// Per subject: the sum of the hazards over its periods and of h (1 - h).
void dfe_sums(const uvec &t, const vec &alpha, const vec &eta, vec &H, vec &W) {
  const uword n = t.n_elem;
  H.set_size(n);
  W.set_size(n);
  for (uword j = 0; j < n; j++) {
    double hs = 0.0, ws = 0.0;
    for (uword k = 0; k < t(j); k++) {
      const double h = dfe_expit(alpha(k) + eta(j));
      hs += h;
      ws += h * (1.0 - h);
    }
    H(j) = hs;
    W(j) = ws;
  }
}

}  // namespace

// Maximum likelihood fit. t: periods at risk per subject (0 = none); rows
// grouped by provider (n_prov rows each); alpha has one value per period with
// alpha(0) held at 0.
// [[Rcpp::export]]
List discrete_fe_fit(const arma::uvec &t, const arma::vec &event, const arma::mat &Z,
                     const arma::vec &n_prov, arma::vec alpha, arma::vec gamma, arma::vec beta,
                     double tol, int max_iter, double bound, bool backtrack, bool message) {
  const uword m = n_prov.n_elem, p = Z.n_cols, K = alpha.n_elem;
  const uvec prov = dfe_rows_provider(n_prov);
  alpha(0) = 0.0;
  vec eta = gamma.elem(prov) + Z * beta;
  int iter = 0;
  double crit = 100.0;
  while (iter < max_iter && crit >= tol) {
    R_CheckUserInterrupt();
    iter++;
    DfeBlocks b = dfe_blocks(t, event, Z, prov, m, alpha, eta, false);
    vec ig_inv = 1.0 / clamp(b.i_gamma, 1e-20, datum::inf);
    mat j1 = b.i_gamma_theta.each_col() % ig_inv;               // m x q
    mat schur = b.i_theta - b.i_gamma_theta.t() * j1;            // q x q
    vec rhs = b.u_theta - j1.t() * b.u_gamma;
    vec d_theta = schur.n_rows > 0 ? vec(solve(schur, rhs, solve_opts::likely_sympd)) : vec();
    vec d_gamma = ig_inv % (b.u_gamma - b.i_gamma_theta * d_theta);
    vec d_alpha(K, fill::zeros);
    if (K > 1) d_alpha.subvec(1, K - 1) = d_theta.subvec(0, K - 2);
    vec d_beta = p > 0 ? vec(d_theta.subvec(K - 1, K - 2 + p)) : vec();

    double v = 1.0;
    if (backtrack) {
      const double ll = dfe_loglik(t, event, alpha, eta);
      const double slope = dot(b.u_gamma, d_gamma) + dot(b.u_theta, d_theta);
      auto ll_at = [&](double s) {
        vec eta_s = (gamma + s * d_gamma).eval().elem(prov) +
          (p > 0 ? vec(Z * (beta + s * d_beta)) : vec(Z.n_rows, fill::zeros));
        return dfe_loglik(t, event, alpha + s * d_alpha, eta_s);
      };
      while (ll_at(v) - ll < 0.01 * v * slope && v > 1e-10) v *= 0.6;
    }
    gamma += v * d_gamma;
    const double med = median(gamma);
    gamma = clamp(gamma, med - bound, med + bound);
    alpha += v * d_alpha;
    if (p > 0) beta += v * d_beta;
    eta = gamma.elem(prov) + Z * beta;
    crit = d_theta.n_elem > 0 ? v * norm(d_theta, "inf") : 0.0;
    if (message) Rcout << "iteration " << iter << ": max change " << crit << std::endl;
  }
  return List::create(_["gamma"] = gamma, _["alpha"] = alpha, _["beta"] = beta,
                      _["iter"] = iter, _["loglik"] = dfe_loglik(t, event, alpha, eta));
}

// The information blocks at the estimate (hazards clipped as in the logistic
// variance routines), and the person-period means of the design per provider
// (the case-mix point of the Wald variance): x_bar = (share of the provider's
// person-periods in periods 2..K, person-period mean of Z).
// [[Rcpp::export]]
List discrete_fe_info(const arma::uvec &t, const arma::vec &event, const arma::mat &Z,
                      const arma::vec &n_prov, const arma::vec &alpha, const arma::vec &gamma,
                      const arma::vec &beta) {
  const uword m = n_prov.n_elem, p = Z.n_cols, K = alpha.n_elem, q = K - 1 + p;
  const uvec prov = dfe_rows_provider(n_prov);
  const vec eta = gamma.elem(prov) + Z * beta;
  DfeBlocks b = dfe_blocks(t, event, Z, prov, m, alpha, eta, true);
  mat xbar(m, q, fill::zeros);
  vec periods(m, fill::zeros);
  for (uword j = 0; j < Z.n_rows; j++) {
    const uword i = prov(j);
    periods(i) += t(j);
    for (uword k = 1; k < t(j); k++) xbar(i, k - 1) += 1.0;
    for (uword l = 0; l < p; l++) xbar(i, K - 1 + l) += t(j) * Z(j, l);
  }
  for (uword i = 0; i < m; i++) {
    if (periods(i) > 0) xbar.row(i) /= periods(i);
  }
  return List::create(_["i_gamma"] = b.i_gamma, _["i_gamma_theta"] = b.i_gamma_theta,
                      _["i_theta"] = b.i_theta, _["xbar"] = xbar, _["person_periods"] = periods);
}

// Patient-clustered meat of theta with the provider effects profiled out:
// sum over clusters of u_c u_c', u_c = sum over the cluster's subjects of
// s_j - J_i (delta_j - H_j), s_j the subject's score for theta and J_i the
// provider's row of I_gamma_theta / I_gamma. cluster: codes of the rows,
// each cluster's rows contiguous.
// [[Rcpp::export]]
arma::mat discrete_fe_theta_meat(const arma::uvec &t, const arma::vec &event,
                                 const arma::mat &Z, const arma::vec &n_prov,
                                 const arma::uvec &cluster, const arma::vec &alpha,
                                 const arma::vec &gamma, const arma::vec &beta,
                                 const arma::mat &J) {
  const uword n = Z.n_rows, p = Z.n_cols, K = alpha.n_elem, q = K - 1 + p;
  const uvec prov = dfe_rows_provider(n_prov);
  const vec eta = gamma.elem(prov) + Z * beta;
  mat meat(q, q, fill::zeros);
  vec u(q, fill::zeros);
  for (uword j = 0; j < n; j++) {
    if (j > 0 && cluster(j) != cluster(j - 1)) {
      meat += u * u.t();
      u.zeros();
    }
    double hs = 0.0;
    for (uword k = 0; k < t(j); k++) {
      const double h = dfe_expit(alpha(k) + eta(j));
      hs += h;
      if (k > 0) u(k - 1) -= h;
    }
    if (event(j) > 0 && t(j) > 1) u(t(j) - 2) += 1.0;
    const double r = event(j) - hs;
    for (uword l = 0; l < p; l++) u(K - 1 + l) += Z(j, l) * r;
    u -= J.row(prov(j)).t() * r;
  }
  if (n > 0) meat += u * u.t();
  return meat;
}

// Per subject: sum of the hazards over its periods and of h (1 - h), at the
// linear predictors eta (provider effect included).
// [[Rcpp::export]]
List discrete_subject_sums(const arma::uvec &t, const arma::vec &alpha, const arma::vec &eta) {
  vec H, W;
  dfe_sums(t, alpha, eta, H, W);
  return List::create(_["H"] = H, _["W"] = W);
}

// Direct standardization: for every provider effect g, the sums of
// expit(g + alpha_k + zb_j) and of its h (1 - h) over every person-period of
// the population, with expit(g + e) = 1 / (1 + exp(-g) exp(-e)) so that the
// exponentials are taken once per provider, period and subject.
// [[Rcpp::export]]
arma::mat discrete_direct(const arma::vec &gamma, const arma::uvec &t, const arma::vec &alpha,
                          const arma::vec &zb, int threads = 1) {
  const int m = gamma.n_elem;
  const uword n = zb.n_elem, K = alpha.n_elem;
  mat out(m, 2, fill::zeros);
  vec a = exp(-alpha), c = exp(-zb);
  omp_set_num_threads(std::max(1, threads));
#pragma omp parallel for schedule(static)
  for (int i = 0; i < m; i++) {
    const double g = gamma(i);
    long double s1 = 0.0L, s2 = 0.0L;
    const double b = std::exp(-g);
    const bool direct = !std::isfinite(b);
    for (uword j = 0; j < n; j++) {
      for (uword k = 0; k < t(j) && k < K; k++) {
        double h;
        if (direct) {
          h = dfe_expit(g + alpha(k) + zb(j));
        } else {
          const double e = b * a(k) * c(j);
          h = std::isfinite(e) ? 1.0 / (1.0 + e) : dfe_expit(g + alpha(k) + zb(j));
        }
        s1 += h;
        s2 += h * (1.0 - h);
      }
    }
    out(i, 0) = (double) s1;
    out(i, 1) = (double) s2;
  }
  return out;
}
