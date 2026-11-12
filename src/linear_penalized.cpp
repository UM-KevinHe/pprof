// Lasso and group-lasso paths of the linear model with fixed provider effects,
//   y_ij = gamma_i + Z_ij' beta + e_ij.
// The provider effects are profiled out exactly: for any beta the least-squares
// gamma_i is the provider mean of y - Z beta, so the objective
//   (1 / 2n) sum (y - gamma - Z beta)^2 + lambda * sum_g m_g ||beta_g||_2
// minimized over gamma is the same objective in beta on the outcome and
// covariates centred within providers. The path works on the Gram matrix
// G = Z'Z / n and c = Z'y / n of the centred data (covariance updates), so an
// update costs O(p) whatever the number of rows.
//
// Groups are given by column ranges; a lasso is the case of groups of one
// column. A group's update is the proximal step of the quadratic majorized by
// L_g = the largest eigenvalue of G_gg (exact for one column); unpenalized
// columns take exact coordinate steps.
#define STRICT_R_HEADERS
#define ARMA_DONT_USE_OPENMP
#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// G, c: Gram matrix and cross-products of the centred data divided by n; yy:
// the centred outcome's sum of squares; beta: start (unpenalized columns at
// their fit, penalized at zero); unpenalized: 0-based columns;
// group_first / group_last: 0-based column ranges of the penalized groups.
// Returns beta per lambda, the residual sum of squares and the iterations.
// [[Rcpp::export]]
List linear_lasso_path(const arma::mat &G, const arma::vec &c, double yy, double n,
                       arma::vec beta, const arma::uvec &unpenalized,
                       const arma::uvec &group_first, const arma::uvec &group_last,
                       const arma::vec &multiplier, const arma::vec &lambda, double tol,
                       int max_iter, bool active_set) {
  const uword p = G.n_cols, n_group = group_first.n_elem, n_lambda = lambda.n_elem;
  vec lip(n_group);
  for (uword g = 0; g < n_group; g++) {
    mat block = G.submat(group_first(g), group_first(g), group_last(g), group_last(g));
    lip(g) = group_first(g) == group_last(g) ? block(0, 0) : eig_sym(block).max();
    if (lip(g) <= 0) lip(g) = 1e-12;
  }
  vec gb = G * beta;  // G beta, kept up to date
  mat beta_path(p, n_lambda);
  vec rss(n_lambda), iters(n_lambda);
  uvec active(n_group);
  for (uword g = 0; g < n_group; g++) {
    active(g) = (!active_set || any(beta.subvec(group_first(g), group_last(g)) != 0)) ? 1 : 0;
  }

  auto move = [&](uword first, uword last, const vec &delta) {
    beta.subvec(first, last) += delta;
    gb += G.cols(first, last) * delta;
  };
  // the proximal step of group g; with apply = false only whether it would
  // leave zero (the optimality condition at zero)
  auto group_step = [&](uword g, double lam, bool apply) -> double {
    const uword f = group_first(g), l = group_last(g);
    vec grad = c.subvec(f, l) - gb.subvec(f, l);
    if (!apply) return norm(grad, 2) > lam * multiplier(g) ? 1.0 : 0.0;
    vec z = beta.subvec(f, l) + grad / lip(g);
    double nz = norm(z, 2);
    double shrink = nz > 0 ? std::max(0.0, 1.0 - lam * multiplier(g) / (lip(g) * nz)) : 0.0;
    vec delta = shrink * z - beta.subvec(f, l);
    double change = abs(delta).max();
    if (change > 0) move(f, l, delta);
    return change;
  };

  for (uword li = 0; li < n_lambda; li++) {
    R_CheckUserInterrupt();
    const double lam = lambda(li);
    int iter = 0;
    while (true) {
      while (iter < max_iter) {
        iter++;
        double max_change = 0.0;
        for (uword u = 0; u < unpenalized.n_elem; u++) {
          const uword j = unpenalized(u);
          if (G(j, j) <= 0) continue;
          vec delta(1);
          delta(0) = (c(j) - gb(j)) / G(j, j);
          if (delta(0) != 0) {
            move(j, j, delta);
            max_change = std::max(max_change, std::fabs(delta(0)));
          }
        }
        for (uword g = 0; g < n_group; g++) {
          if (!active(g)) continue;
          max_change = std::max(max_change, group_step(g, lam, true));
        }
        if (max_change < tol) break;
      }
      if (!active_set) break;
      int added = 0;
      for (uword g = 0; g < n_group; g++) {
        if (active(g)) continue;
        if (group_step(g, lam, false) > 0) {
          active(g) = 1;
          added++;
        }
      }
      if (added == 0 || iter >= max_iter) break;
    }
    beta_path.col(li) = beta;
    rss(li) = yy - 2.0 * n * dot(c, beta) + n * dot(beta, gb);
    iters(li) = iter;
  }
  return List::create(_["beta"] = beta_path, _["deviance"] = rss, _["iter"] = iters);
}
