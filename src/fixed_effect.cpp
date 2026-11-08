// Unified C++ engine for the fixed provider effects logistic model.
//
// Every routine takes the number of trials N (binomial outcome). With N = 1
// each expression reduces exactly (multiplication by 1.0 is exact) to the
// binary code previously copied into the SFR / PPPW / SRR measure folders.
//
// Armadillo's internal OpenMP is switched off so that sums inside Armadillo
// expressions do not depend on the number of cores; the explicit
// "#pragma omp" loops below remain parallel.
#define STRICT_R_HEADERS
#define ARMA_DONT_USE_OPENMP
#include <RcppArmadillo.h>
#include <cmath>
#include <iomanip>
#include <string>
#include "pprof_omp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Repeat x[i] each[i] times (provider-level value -> observation level).
arma::vec rep_vec(const arma::vec &x, const arma::vec &each) {
  arma::vec x_rep(sum(each));
  int ind = 0, m = x.n_elem;
  for (int i = 0; i < m; i++) {
    x_rep.subvec(ind, ind + each(i) - 1) = x(i) * ones(each(i));
    ind += each(i);
  }
  return x_rep;
}

void ind2uppsub(unsigned int index, unsigned int dim, unsigned int &row, unsigned int &col) {
  row = 0, col = dim - 1;
  unsigned int n = dim * (dim - 1) / 2 - (dim - row) * (dim - row - 1) / 2 + col;
  while (index > n) {
    ++row;
    n = dim * (dim - 1) / 2 - (dim - row) * (dim - row - 1) / 2 + col;
  }
  while (index < n) {
    --col;
    --n;
  }
}

arma::mat info_beta_omp(const arma::mat &Z, const arma::vec &pq, const int &threads) {
  omp_set_num_threads(threads);
  unsigned int p = Z.n_cols;
  unsigned int loops = p * (1 + p) / 2;
  arma::mat output(p, p);
#pragma omp parallel for schedule(static)
  for (unsigned int i = 0; i < loops; i++) {
    unsigned int r, c;
    ind2uppsub(i, p, r, c);
    output(r, c) = dot(Z.col(r), Z.col(c) % pq);
    output(c, r) = output(r, c);
  }
  return output;
}

// Binomial log-likelihood without the log(choose(N, Y)) constant.
double Loglkd(const arma::vec &N, const arma::vec &Y, const arma::vec &Z_beta,
              const arma::vec &gamma_obs) {
  return sum((gamma_obs + Z_beta) % Y - N % log(1 + exp(gamma_obs + Z_beta)));
}

double logist(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

// Convergence criterion for the stopping rules ported from pprof.
// "beta" is the rule the measure code has always used.
double stop_crit(const std::string &stop, double crit_beta, double d_loglkd,
                 double loglkd, double loglkd_init) {
  double crit_relch = std::abs(d_loglkd / (d_loglkd + loglkd));
  double crit_ratch = std::abs(d_loglkd / (d_loglkd + loglkd - loglkd_init));
  if (stop == "relch") return crit_relch;
  if (stop == "ratch") return crit_ratch;
  if (stop == "all") return std::max(crit_beta, std::max(crit_relch, crit_ratch));
  if (stop == "or") return std::min(crit_beta, std::min(crit_relch, crit_ratch));
  return crit_beta;
}

void print_crit(const std::string &stop, int iter, double crit, bool newline) {
  if (stop == "beta") {
    Rcout << "Iter " << iter << ": Inf norm of running diff in est reg parm is "
          << std::scientific << std::setprecision(3) << crit << ";";
  } else if (stop == "relch") {
    Rcout << "Iter " << iter << ": Relative change in est log likelihood is "
          << std::scientific << std::setprecision(3) << crit << ";";
  } else if (stop == "ratch") {
    Rcout << "Iter " << iter << ": Adjusted relative change in est log likelihood is "
          << std::scientific << std::setprecision(3) << crit << ";";
  } else if (stop == "all") {
    Rcout << "Iter " << iter << ": Maximum criterion across all checks is "
          << std::scientific << std::setprecision(3) << crit << ";";
  } else {
    Rcout << "Iter " << iter << ": Minimum criterion across all checks is "
          << std::scientific << std::setprecision(3) << crit << ";";
  }
  if (newline) Rcout << std::endl;
}

// ---------------------------------------------------------------------------
// Model fitting
// ---------------------------------------------------------------------------

// Block ascent Newton (BAN).
// [[Rcpp::export]]
List logis_fe_prov(arma::vec &N, arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec gamma,
                   arma::vec beta, bool backtrack = true, int max_iter = 10000, double bound = 10.0,
                   double tol = 1e-5, bool message = false, std::string stop = "beta") {

  int iter = 0, n = Z.n_rows, m = n_prov.n_elem, ind;
  double v;
  arma::vec gamma_obs(n);
  double crit = 100.0;
  bool need_lkd = (stop != "beta");
  double loglkd_init = 0.0, loglkd_old = 0.0;
  if (need_lkd) {
    gamma_obs = rep_vec(gamma, n_prov);
    loglkd_init = Loglkd(N, Y, Z * beta, gamma_obs);
  }

  if (message == true) {
    Rcout << "Implementing BAN algorithm (Rcpp) for fixed provider effects model ..." << std::endl;
  }

  double s = 0.01, t = 0.6; // only used for "backtrack = true"
  double lambda, d_loglkd, loglkd;
  arma::vec gamma_obs_tmp(n), gamma_tmp(m), beta_tmp(Z.n_cols);

  while (iter < max_iter) {
    if (crit < tol) {
      break;
    }
    iter++;
    // provider effect update
    gamma_obs = rep_vec(gamma, n_prov);
    arma::vec Z_beta = Z * beta;
    if (need_lkd) loglkd_old = Loglkd(N, Y, Z_beta, gamma_obs);
    arma::vec p = 1 / (1 + exp(-gamma_obs - Z_beta));
    arma::vec Yp = Y - N % p, pq = N % p % (1 - p);
    if (any(pq == 0)) {
      pq.replace(0, 1e-20);
    }
    ind = 0;
    arma::vec score_gamma(m), d_gamma(m);
    for (int i = 0; i < m; i++) {
      score_gamma(i) = sum(Yp(span(ind, ind + n_prov(i) - 1)));
      d_gamma(i) = score_gamma(i) / sum(pq(span(ind, ind + n_prov(i) - 1)));
      ind += n_prov(i);
    }
    v = 1.0; // initialize step size
    if (backtrack == true) {
      loglkd = Loglkd(N, Y, Z_beta, gamma_obs);
      gamma_tmp = gamma + v * d_gamma;
      gamma_obs_tmp = rep_vec(gamma_tmp, n_prov);
      d_loglkd = Loglkd(N, Y, Z_beta, gamma_obs_tmp) - loglkd;
      lambda = dot(score_gamma, d_gamma);
      while (d_loglkd < s * v * lambda) {
        v = t * v;
        gamma_tmp = gamma + v * d_gamma;
        gamma_obs_tmp = rep_vec(gamma_tmp, n_prov);
        d_loglkd = Loglkd(N, Y, Z_beta, gamma_obs_tmp) - loglkd;
      }
    }
    gamma += v * d_gamma;
    gamma = clamp(gamma, median(gamma) - bound, median(gamma) + bound);
    gamma_obs = rep_vec(gamma, n_prov);

    // regression parameter update
    p = 1 / (1 + exp(-gamma_obs - Z_beta));
    pq = N % p % (1 - p);
    arma::vec score_beta = Z.t() * (Y - N % p);
    arma::mat info_beta = Z.t() * (Z.each_col() % pq);
    arma::vec d_beta = solve(info_beta, score_beta, solve_opts::fast + solve_opts::likely_sympd);
    v = 1.0; // initialize step size
    if (backtrack == true) {
      loglkd = Loglkd(N, Y, Z_beta, gamma_obs);
      beta_tmp = beta + v * d_beta;
      d_loglkd = Loglkd(N, Y, Z * beta_tmp, gamma_obs) - loglkd;
      lambda = dot(score_beta, d_beta);
      while (d_loglkd < s * v * lambda) {
        v = t * v;
        beta_tmp = beta + v * d_beta;
        d_loglkd = Loglkd(N, Y, Z * beta_tmp, gamma_obs) - loglkd;
      }
    }
    beta += v * d_beta;
    double crit_beta = norm(v * d_beta, "inf");
    if (need_lkd) {
      double d_lkd = Loglkd(N, Y, Z * beta, gamma_obs) - loglkd_old;
      crit = stop_crit(stop, crit_beta, d_lkd, loglkd_old, loglkd_init);
    } else {
      crit = crit_beta;
    }
    if (message == true) {
      print_crit(stop, iter, crit, true);
    }
  }
  if (message == true) {
    Rcout << "BAN algorithm (Rcpp) converged after " << iter << " iterations!" << std::endl;
  }
  List ret = List::create(_["gamma"] = gamma, _["beta"] = beta, _["iter"] = iter);
  return ret;
}

// Serial blockwise inversion Newton (SerBIN).
// [[Rcpp::export]]
List logis_BIN_fe_prov(arma::vec &N, arma::vec &Y, arma::mat &Z, arma::vec &n_prov,
                       arma::vec gamma, arma::vec beta, int threads = 1, double tol = 1e-8,
                       int max_iter = 10000, double bound = 10.0, bool message = true,
                       bool backtrack = true, std::string stop = "beta") {

  int iter = 0, n = Z.n_rows, m = n_prov.n_elem, ind;
  double v;
  arma::vec gamma_obs(n);
  double crit = 100.0;
  bool need_lkd = (stop != "beta");
  double loglkd_init = 0.0, loglkd_cur = 0.0;
  if (need_lkd) {
    gamma_obs = rep_vec(gamma, n_prov);
    loglkd_init = Loglkd(N, Y, Z * beta, gamma_obs);
  }
  if (message == true) {
    Rcout << "Implementing SerBIN algorithm (Rcpp) for fixed provider effects model ..." << std::endl;
  }

  double s = 0.01, t = 0.6; // only used for "backtrack = true"
  double lambda, d_loglkd, loglkd;
  arma::vec gamma_obs_tmp(n), gamma_tmp(m), beta_tmp(Z.n_cols);

  while (iter <= max_iter) {
    if (crit < tol) {
      break;
    }
    iter++;
    gamma_obs = rep_vec(gamma, n_prov);
    arma::vec Z_beta = Z * beta;
    if (need_lkd) loglkd_cur = Loglkd(N, Y, Z_beta, gamma_obs);
    arma::vec p = 1 / (1 + exp(-gamma_obs - Z_beta));
    arma::vec Yp = Y - N % p, pq = N % p % (1 - p);
    if (any(pq == 0)) {
      pq.replace(0, 1e-20);
    }
    arma::vec score_gamma(m), info_gamma_inv(m);
    arma::mat info_betagamma(Z.n_cols, m);
    ind = 0;
    for (int i = 0; i < m; i++) {
      score_gamma(i) = sum(Yp(span(ind, ind + n_prov(i) - 1)));
      info_gamma_inv(i) = 1 / sum(pq(span(ind, ind + n_prov(i) - 1)));
      info_betagamma.col(i) = sum(Z.rows(ind, ind + n_prov(i) - 1).each_col() % pq.subvec(ind, ind + n_prov(i) - 1)).t();
      ind += n_prov(i);
    }
    arma::vec score_beta = Z.t() * Yp;
    arma::mat info_beta(Z.n_cols, Z.n_cols);
    if (threads > 1) {
      info_beta = info_beta_omp(Z, pq, threads);
    } else {
      info_beta = Z.t() * (Z.each_col() % pq);
    }
    arma::mat mat_tmp1 = info_betagamma.each_row() % info_gamma_inv.t();                 // p*m matrix (J_1)
    arma::mat schur_system = info_beta - mat_tmp1 * info_betagamma.t();                  // p*p matrix (S)
    arma::mat mat_tmp2 = arma::solve(schur_system, mat_tmp1, solve_opts::likely_sympd);  // p*m matrix (J_2)
    arma::vec d_gamma = info_gamma_inv % score_gamma + mat_tmp2.t() * (mat_tmp1 * score_gamma - score_beta);
    arma::vec schur_system_solve = arma::solve(schur_system, score_beta, solve_opts::likely_sympd);
    arma::vec d_beta = schur_system_solve - mat_tmp2 * score_gamma;

    v = 1.0; // initialize step size
    if (backtrack == true) {
      loglkd = Loglkd(N, Y, Z_beta, gamma_obs);
      gamma_tmp = gamma + v * d_gamma;
      gamma_obs_tmp = rep_vec(gamma_tmp, n_prov);
      arma::vec Z_beta_tmp = Z * (beta + v * d_beta);
      lambda = dot(score_gamma, d_gamma) + dot(score_beta, d_beta);
      d_loglkd = Loglkd(N, Y, Z_beta_tmp, gamma_obs_tmp) - loglkd;
      while (d_loglkd < s * v * lambda) {
        v = t * v;
        gamma_tmp = gamma + v * d_gamma;
        gamma_obs_tmp = rep_vec(gamma_tmp, n_prov);
        Z_beta_tmp = Z * (beta + v * d_beta);
        d_loglkd = Loglkd(N, Y, Z_beta_tmp, gamma_obs_tmp) - loglkd;
      }
    }
    gamma += v * d_gamma;
    gamma = clamp(gamma, median(gamma) - bound, median(gamma) + bound);
    beta += v * d_beta;
    double crit_beta = norm(v * d_beta, "inf");
    if (need_lkd) {
      double d_lkd = Loglkd(N, Y, Z * beta, rep_vec(gamma, n_prov)) - loglkd_cur;
      crit = stop_crit(stop, crit_beta, d_lkd, loglkd_cur, loglkd_init);
    } else {
      crit = crit_beta;
    }

    if (message == true) {
      print_crit(stop, iter, crit, false);
    }
  }
  if (message == true) {
    Rcout << "serBIN (Rcpp) algorithm converged after " << iter << " iterations!" << std::endl;
  }
  List ret = List::create(_["gamma"] = gamma, _["beta"] = beta, _["iter"] = iter);
  return ret;
}

// Model-based (inverse information) variance of beta and gamma.
// [[Rcpp::export]]
List logis_fe_var(arma::vec &N, arma::mat &Z, arma::vec &n_prov, arma::vec &gamma, arma::vec &beta) {
  arma::vec gamma_obs = rep_vec(gamma, n_prov);
  int m = n_prov.n_elem;
  arma::vec p = 1 / (1 + exp(-gamma_obs - Z * beta));
  p = clamp(p, 1e-10, 1 - 1e-10);
  arma::vec pq = N % p % (1 - p);
  int ind = 0;
  arma::vec info_gamma_inv(m);
  arma::mat info_betagamma(Z.n_cols, m);
  for (int i = 0; i < m; i++) {
    info_gamma_inv(i) = 1 / sum(pq.subvec(ind, ind + n_prov(i) - 1));
    info_betagamma.col(i) = sum(Z.rows(ind, ind + n_prov(i) - 1).each_col() % pq.subvec(ind, ind + n_prov(i) - 1)).t();
    ind += n_prov(i);
  }
  arma::mat info_beta = Z.t() * (Z.each_col() % pq);
  arma::mat info_beta_inv = inv_sympd(info_beta - (info_betagamma.each_row() % info_gamma_inv.t()) * info_betagamma.t());
  arma::vec var_gamma(m);
  for (int i = 0; i < m; i++) {
    arma::vec mat_tmp = info_gamma_inv(i) * info_betagamma.col(i); // J_1^T
    var_gamma(i) = info_gamma_inv(i) + as_scalar(mat_tmp.t() * info_beta_inv * mat_tmp);
  }
  List ret;
  ret["var.beta"] = info_beta_inv;
  ret["se.beta"] = sqrt(info_beta_inv.diag());
  ret["var.gamma"] = var_gamma;
  ret["se.gamma"] = sqrt(var_gamma);
  return ret;
}

// ---------------------------------------------------------------------------
// Provider-effect (gamma) tests
// ---------------------------------------------------------------------------

// Model-based Wald standard errors of selected provider effects.
// parm_indices are 1-based.
// [[Rcpp::export]]
arma::vec wald_gamma(arma::vec &N, arma::mat &Z, arma::vec &p, arma::vec &n_prov,
                     arma::vec &parm_indices, int threads = 1) {
  int m = n_prov.n_elem, ind = 0;
  arma::vec pq = N % p % (1 - p);
  if (any(pq == 0)) {
    pq.replace(0, 1e-20);
  }

  arma::vec info_gamma_inv(m);
  arma::mat info_betagamma(Z.n_cols, m);
  arma::mat info_beta(Z.n_cols, Z.n_cols);

  for (int i = 0; i < m; i++) {
    info_gamma_inv(i) = 1 / sum(pq(span(ind, ind + n_prov(i) - 1)));
    info_betagamma.col(i) = sum(Z.rows(ind, ind + n_prov(i) - 1).each_col() % pq.subvec(ind, ind + n_prov(i) - 1)).t();
    ind += n_prov(i);
  }
  if (threads > 1) {
    info_beta = info_beta_omp(Z, pq, threads);
  } else {
    info_beta = Z.t() * (Z.each_col() % pq);
  }

  arma::mat J1 = info_betagamma.each_row() % info_gamma_inv.t(); // p*m matrix (J_1)
  arma::mat schur_system = info_beta - J1 * info_betagamma.t();  // p*p matrix (S)

  arma::uvec col_indices = arma::conv_to<arma::uvec>::from(parm_indices) - 1;
  arma::mat selected_J1 = J1.cols(col_indices);
  arma::vec selected_info_gamma_inv = info_gamma_inv(col_indices);

  int n_parm = col_indices.n_elem;
  arma::vec se_gamma(n_parm);
  for (int i = 0; i < n_parm; i++) {
    arma::vec J1_i = selected_J1.col(i);
    arma::vec J2_i = arma::solve(schur_system, J1_i, solve_opts::likely_sympd);
    se_gamma(i) = sqrt(selected_info_gamma_inv(i) + dot(J1_i, J2_i));
  }
  return se_gamma;
}

// Bread of the patient-clustered sandwich for selected provider effects:
// returns sqrt(I_gamma^-1). The R side multiplies se^4 by the meat.
// [[Rcpp::export]]
arma::vec robust_wald_gamma(arma::vec &N, arma::mat &Z, arma::vec &p, arma::vec &n_prov,
                            arma::vec &parm_indices) {
  int m = n_prov.n_elem, ind = 0;
  arma::vec pq = N % p % (1 - p);
  if (any(pq == 0)) {
    pq.replace(0, 1e-20);
  }

  arma::vec info_gamma_inv(m);
  for (int i = 0; i < m; i++) {
    info_gamma_inv(i) = 1 / sum(pq(span(ind, ind + n_prov(i) - 1)));
    ind += n_prov(i);
  }

  arma::uvec col_indices = arma::conv_to<arma::uvec>::from(parm_indices) - 1;
  arma::vec selected_info_gamma_inv = info_gamma_inv(col_indices);

  int n_parm = col_indices.n_elem;
  arma::vec se_gamma(n_parm);
  for (int i = 0; i < n_parm; i++) {
    se_gamma(i) = sqrt(selected_info_gamma_inv(i));
  }
  return se_gamma;
}

// ---------------------------------------------------------------------------
// Covariate (beta) tests
// ---------------------------------------------------------------------------

// Model-based Wald test for selected covariates; indices are 1-based.
// [[Rcpp::export]]
List wald_covar(arma::vec &N, arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec &gamma,
                arma::vec &beta, arma::uvec &indices, double null, double alpha) {

  indices -= 1; // switch to C indexing
  arma::vec gamma_obs = rep_vec(gamma, n_prov);
  int m = n_prov.n_elem;
  arma::vec p = 1 / (1 + exp(-gamma_obs - Z * beta));
  p = clamp(p, 1e-10, 1 - 1e-10);
  int ind = 0;
  arma::vec info_gamma_inv(m);
  arma::mat info_betagamma(Z.n_cols, m);
  for (int i = 0; i < m; i++) {
    info_gamma_inv(i) = 1 / dot(N.subvec(ind, ind + n_prov(i) - 1) % p.subvec(ind, ind + n_prov(i) - 1), 1 - p.subvec(ind, ind + n_prov(i) - 1));
    info_betagamma.col(i) =
      sum(Z.rows(ind, ind + n_prov(i) - 1).each_col() % (N.subvec(ind, ind + n_prov(i) - 1) % p.subvec(ind, ind + n_prov(i) - 1) % (1 - p.subvec(ind, ind + n_prov(i) - 1)))).t();
    ind += n_prov(i);
  }
  arma::mat info_beta = Z.t() * (Z.each_col() % (N % p % (1 - p)));
  arma::mat info_beta_inv = inv_sympd(info_beta - (info_betagamma.each_row() % info_gamma_inv.t()) * info_betagamma.t());
  arma::vec se_beta = sqrt(info_beta_inv.diag());
  arma::vec stat = (beta(indices) - null) / se_beta(indices);
  arma::vec p_lower = normcdf(stat);
  List ret;
  ret["stat"] = stat;
  ret["p"] = 2 * min(p_lower, 1 - p_lower);
  ret["se.beta"] = 1 * se_beta(indices);
  ret["beta.lower"] = beta(indices) - R::qnorm5(1 - alpha / 2, 0, 1, true, false) * se_beta(indices);
  ret["beta.upper"] = beta(indices) + R::qnorm5(1 - alpha / 2, 0, 1, true, false) * se_beta(indices);
  return ret;
}

// Patient-clustered sandwich (robust Wald) test for selected covariates.
// A0: provider-level sum of squared patient score sums (gamma block);
// A1_t: provider x covariate cross terms; A2: covariate block of the meat.
// [[Rcpp::export]]
List robust_wald_covar(arma::vec &N, arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec &gamma,
                       arma::vec &beta, arma::uvec &indices, double null_value, double alpha,
                       arma::vec &A0, arma::mat &A1_t, arma::mat &A2) {

  indices -= 1;
  arma::vec gamma_obs = rep_vec(gamma, n_prov);
  int m = n_prov.n_elem;
  arma::vec p = 1 / (1 + exp(-gamma_obs - Z * beta));
  p = clamp(p, 1e-10, 1 - 1e-10);
  int ind = 0;
  arma::vec info_gamma_inv(m);
  arma::mat info_betagamma(Z.n_cols, m);
  for (int i = 0; i < m; i++) {
    info_gamma_inv(i) = 1 / dot(N.subvec(ind, ind + n_prov(i) - 1) % p.subvec(ind, ind + n_prov(i) - 1), 1 - p.subvec(ind, ind + n_prov(i) - 1));
    info_betagamma.col(i) =
      sum(Z.rows(ind, ind + n_prov(i) - 1).each_col() % (N.subvec(ind, ind + n_prov(i) - 1) % p.subvec(ind, ind + n_prov(i) - 1) % (1 - p.subvec(ind, ind + n_prov(i) - 1)))).t();
    ind += n_prov(i);
  }
  arma::mat info_beta = Z.t() * (Z.each_col() % (N % p % (1 - p)));
  arma::mat mat_tmp1 = info_betagamma.each_row() % info_gamma_inv.t();                 // p*m matrix (J_1)
  arma::mat schur_system = info_beta - mat_tmp1 * info_betagamma.t();                  // p*p matrix (S)
  arma::mat mat_tmp2 = arma::solve(schur_system, mat_tmp1, solve_opts::likely_sympd);  // p*m matrix (J_2; -B2)
  arma::mat schur_system_inv;                                                          // p*p matrix (S^-1; B3)
  bool success = arma::inv(schur_system_inv, schur_system);
  if (!success) {
    double eps = 1e-8 * arma::trace(schur_system) / schur_system.n_rows;
    schur_system_inv = arma::pinv(schur_system + eps * arma::eye(schur_system.n_rows, schur_system.n_cols));
  }
  arma::mat mat_tmp2_t = mat_tmp2.t();
  arma::mat var_beta = mat_tmp2 * (mat_tmp2_t.each_col() % A0) -  // B2 A0 B2
    2 * mat_tmp2 * A1_t * schur_system_inv +                        // 2 * (B2 A1^T B3)
    schur_system_inv * A2 * schur_system_inv;                       // B3 A2 B3

  arma::vec se_beta = sqrt(var_beta.diag());
  arma::vec stat = (beta(indices) - null_value) / se_beta(indices);
  arma::vec p_lower = normcdf(stat);

  List ret;
  ret["stat"] = stat;
  ret["p"] = 2 * min(p_lower, 1 - p_lower);
  ret["se.beta"] = 1 * se_beta(indices);
  ret["beta.lower"] = beta(indices) - R::qnorm5(1 - alpha / 2, 0, 1, true, false) * se_beta(indices);
  ret["beta.upper"] = beta(indices) + R::qnorm5(1 - alpha / 2, 0, 1, true, false) * se_beta(indices);
  return ret;
}

// Sum of outer products of the rows of X (covariate block of the meat).
// [[Rcpp::export]]
arma::mat compute_A2(const arma::mat &X_matrix) {
  arma::mat A2 = arma::zeros<arma::mat>(X_matrix.n_cols, X_matrix.n_cols);
  for (unsigned int i = 0; i < X_matrix.n_rows; ++i) {
    arma::rowvec Xi = X_matrix.row(i);
    A2 += Xi.t() * Xi;
  }
  return A2;
}

// ---------------------------------------------------------------------------
// Standardized measures
// ---------------------------------------------------------------------------

// For every provider effect gamma_j, sums over a standard population:
//   column 1: sum_i w_i * expit(gamma_j + Z_beta_i)
//   column 2: sum_i w_i * expit(.) * (1 - expit(.))   (delta-method SE numerator)
// w = 1 gives the binary direct standardization used by SFR / PPPW;
// w = number of trials gives the binomial (AOH) version.
//
// Two exact rearrangements keep this O(m x distinct(Z_beta)) and free of
// exp() in the inner loop:
//   * rows sharing a linear predictor are merged, their weights summed
//     (categorical covariates leave few distinct values);
//   * expit(g + e) = 1 / (1 + t) with t = exp(-g) * exp(-e), and
//     expit * (1 - expit) = t / (1 + t)^2, so each exponential is taken once
//     per provider and once per distinct row.
// Only the rounding differs from summing logist(g + e) row by row.
// [[Rcpp::export]]
arma::mat computeDirectExp(const arma::vec &gamma_prov, const arma::vec &Z_beta,
                           const arma::vec &w, const int threads = 1) {
  const int m = gamma_prov.n_elem;
  const arma::uword n = Z_beta.n_elem;
  arma::mat results(m, 2, arma::fill::zeros);
  if (m == 0 || n == 0) return results;

  // distinct linear predictors and their total weight
  std::vector<double> eta, wt;
  if (Z_beta.has_nan()) { // no ordering possible: keep the rows as they are
    eta.assign(Z_beta.begin(), Z_beta.end());
    wt.assign(w.begin(), w.end());
  } else {
    arma::uvec ord = arma::stable_sort_index(Z_beta); // ties keep row order
    double cur = Z_beta(ord(0));
    long double acc = 0.0L;
    for (arma::uword k = 0; k < n; ++k) {
      const double e = Z_beta(ord(k));
      if (e != cur) {
        eta.push_back(cur);
        wt.push_back((double) acc);
        cur = e;
        acc = 0.0L;
      }
      acc += w(ord(k));
    }
    eta.push_back(cur);
    wt.push_back((double) acc);
  }
  const std::size_t u = eta.size();
  std::vector<double> a(u);
  double a_max = 0.0;
  bool a_nan = false;
  for (std::size_t i = 0; i < u; ++i) {
    a[i] = std::exp(-eta[i]);
    if (std::isnan(a[i])) a_nan = true;
    else if (a[i] > a_max) a_max = a[i];
  }
  if (a_nan) a_max = arma::datum::nan; // sends every provider to the direct loop

  // long double accumulators: for a rate close to 1, 1 - rate would otherwise
  // be lost to the rounding of a large sum
  omp_set_num_threads(threads);
#pragma omp parallel for schedule(static)
  for (int j = 0; j < m; ++j) {
    const double g = gamma_prov[j];
    const double b = std::exp(-g);
    long double sum1 = 0.0L;
    long double sum2 = 0.0L;
    if (std::isfinite(b * a_max)) {
      for (std::size_t i = 0; i < u; ++i) {
        const double t = b * a[i];
        const double d = 1.0 / (1.0 + t);
        sum1 += wt[i] * d;
        sum2 += wt[i] * (t * d * d);
      }
    } else { // extreme or non-finite values: evaluate directly
      for (std::size_t i = 0; i < u; ++i) {
        const double r = logist(g + eta[i]);
        sum1 += wt[i] * r;
        sum2 += wt[i] * ((1.0 - r) * r);
      }
    }
    results(j, 0) = (double) sum1;
    results(j, 1) = (double) sum2;
  }
  return results;
}

// ---------------------------------------------------------------------------
// Patient-clustered sums for the robust and stabilized tests
// ---------------------------------------------------------------------------
// cluster: 1-based (provider, patient) code of each row, numbered in order of
// first appearance; prov_row: 1-based provider of each row; cluster_prov:
// 1-based provider of each cluster. Sums use long double accumulators in row
// order, i.e. exactly what base R's sum() does on the same groups.

// Per provider: U = sum of r over its rows, and
// meat = sum over its patients of (patient sum of r)^2.
// [[Rcpp::export]]
List cluster_sums(const arma::vec &r, const arma::uvec &cluster, const arma::uvec &prov_row,
                  const arma::uvec &cluster_prov, const int m) {
  const arma::uword n = r.n_elem, K = cluster_prov.n_elem;
  std::vector<long double> csum(K, 0.0L), U(m, 0.0L), meat(m, 0.0L);
  for (arma::uword i = 0; i < n; ++i) {
    csum[cluster(i) - 1] += r(i);
    U[prov_row(i) - 1] += r(i);
  }
  for (arma::uword c = 0; c < K; ++c) {
    const double s = (double) csum[c];
    meat[cluster_prov(c) - 1] += s * s;
  }
  arma::vec U_out(m), meat_out(m);
  for (int j = 0; j < m; ++j) {
    U_out(j) = (double) U[j];
    meat_out(j) = (double) meat[j];
  }
  return List::create(Named("U") = U_out, Named("meat") = meat_out);
}

// Stabilized variance of the provider score (Feiyang's stable_wald_tmp):
// residuals standardized by their own sd, e = r / sqrt(v), give one pooled
// matrix over the positions of a patient's rows (1st, 2nd, ... in data
// order), R[k, l] = sum_patients e_k e_l / #{patients with > max(k, l) rows};
// each patient's variance is then a' R a with a = sqrt(v), summed per
// provider. Rows with v = 0 drop out.
// [[Rcpp::export]]
arma::vec stabilized_meat(const arma::vec &r, const arma::vec &v, const arma::uvec &cluster,
                          const arma::uvec &cluster_prov, const int m) {
  const arma::uword n = r.n_elem, K = cluster_prov.n_elem;

  // rows of every patient, in data order (counting sort on the cluster code)
  arma::uvec start(K + 1, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i) start(cluster(i))++;
  for (arma::uword c = 0; c < K; ++c) start(c + 1) += start(c);
  arma::uvec rows(n), next = start.head(K);
  for (arma::uword i = 0; i < n; ++i) rows(next(cluster(i) - 1)++) = i;

  arma::uword dim = 0;
  for (arma::uword c = 0; c < K; ++c) dim = std::max(dim, start(c + 1) - start(c));

  arma::vec e(n), a(n);
  for (arma::uword i = 0; i < n; ++i) {
    if (v(i) > 0) {
      a(i) = std::sqrt(v(i));
      e(i) = r(i) / a(i);
    } else {
      a(i) = 0.0;
      e(i) = 0.0;
    }
  }

  // pooled matrix, entry (k, l) averaged over patients with > max(k, l) rows
  arma::mat S(dim, dim, arma::fill::zeros);
  arma::vec n_ge(dim, arma::fill::zeros);
  for (arma::uword c = 0; c < K; ++c) {
    const arma::uword b = start(c), s = start(c + 1) - b;
    for (arma::uword k = 0; k < s; ++k) {
      const double ek = e(rows(b + k));
      for (arma::uword l = 0; l < s; ++l) S(k, l) += ek * e(rows(b + l));
      n_ge(k) += 1.0;
    }
  }
  for (arma::uword k = 0; k < dim; ++k) {
    for (arma::uword l = 0; l < dim; ++l) {
      const double d = n_ge(std::max(k, l));
      S(k, l) = d > 0 ? S(k, l) / d : 0.0;
    }
  }

  arma::vec meat(m, arma::fill::zeros);
  for (arma::uword c = 0; c < K; ++c) {
    const arma::uword b = start(c), s = start(c + 1) - b;
    double W = 0.0;
    for (arma::uword k = 0; k < s; ++k) {
      const double ak = a(rows(b + k));
      if (ak == 0.0) continue;
      for (arma::uword l = 0; l < s; ++l) W += ak * S(k, l) * a(rows(b + l));
    }
    meat(cluster_prov(c) - 1) += W;
  }
  return meat;
}

// ---------------------------------------------------------------------------
// Restricted score test of H0: gamma_i = gamma_null for the requested providers
// (1-based `parm_indices`): the score of gamma_i at the null divided by the
// square root of its efficient information with every other parameter at its
// full-model estimate. The `Modified_score` of pprof 1.0.3 (score_modified =
// FALSE) extended to binomial trials. Only provider i's rows differ between the
// full and the null information, so the p x p system is updated per provider
// instead of rebuilt: O(n p^2) in total.
// [[Rcpp::export]]
arma::vec restricted_score(const arma::vec &N, const arma::vec &Y, const arma::mat &Z,
                           const arma::vec &n_prov, const arma::vec &gamma, const arma::vec &beta,
                           double gamma_null, const arma::uvec &parm_indices) {
  const int m = n_prov.n_elem;
  arma::vec gamma_obs = rep_vec(gamma, n_prov);
  arma::vec Zb = Z * beta;
  arma::vec p = 1 / (1 + exp(-gamma_obs - Zb));
  arma::vec pq = N % p % (1 - p);
  pq.replace(0, 1e-20);
  arma::vec p0 = 1 / (1 + exp(-gamma_null - Zb));
  arma::vec pq0 = N % p0 % (1 - p0);
  arma::vec r0 = Y - N % p0;

  arma::vec info_gamma_inv(m);
  arma::mat info_bg(Z.n_cols, m);
  arma::uvec start(m + 1);
  int ind = 0;
  for (int i = 0; i < m; i++) {
    start(i) = ind;
    int last = ind + n_prov(i) - 1;
    info_gamma_inv(i) = 1 / sum(pq.subvec(ind, last));
    info_bg.col(i) = Z.rows(ind, last).t() * pq.subvec(ind, last);
    ind += n_prov(i);
  }
  start(m) = ind;
  arma::mat info_beta_full = Z.t() * (Z.each_col() % pq);
  arma::mat m_full = (info_bg.each_row() % info_gamma_inv.t()) * info_bg.t();

  arma::vec z(m);
  z.fill(arma::datum::nan);
  for (arma::uword k = 0; k < parm_indices.n_elem; k++) {
    int i = parm_indices(k) - 1;
    int a = start(i), b = start(i + 1) - 1;
    arma::mat Zi = Z.rows(a, b);
    arma::vec d = pq0.subvec(a, b) - pq.subvec(a, b);
    arma::mat S = info_beta_full + Zi.t() * (Zi.each_col() % d) - m_full +
      info_bg.col(i) * info_bg.col(i).t() * info_gamma_inv(i);
    arma::mat B11;
    if (!arma::inv_sympd(B11, S)) B11 = arma::pinv(S);
    double info_alpha = sum(pq0.subvec(a, b));
    arma::vec info_ba = Zi.t() * pq0.subvec(a, b);
    double v0 = info_alpha - arma::as_scalar(info_ba.t() * B11 * info_ba);
    z(i) = sum(r0.subvec(a, b)) / sqrt(v0);
  }
  return z;
}
