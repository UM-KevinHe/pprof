//#define ARMA_NO_DEBUG
#define STRICT_R_HEADERS // needed on Windows, not on macOS
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <omp.h>
#include <iostream>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
using namespace RcppParallel;
using namespace Rcpp;
using namespace std;
using namespace arma;

//` @importFrom RcppParallel RcppParallelLibs

double Loglkd2(const arma::vec &Y, const arma::vec &Z_beta, const arma::vec &gamma_obs) {
  return arma::accu((gamma_obs + Z_beta) % Y - arma::log(1 + arma::exp(gamma_obs + Z_beta)));
}

arma::vec rep2(arma::vec &x, arma::vec &each) {
  arma::vec x_rep(sum(each));
  int ind = 0, m = x.n_elem;
  for (int i = 0; i < m; i++) {
    x_rep.subvec(ind,ind+each(i)-1) = x(i) * ones(each(i));
    ind += each(i);
  }
  return x_rep;
}

// [[Rcpp::export]]
List logis_BIN_fe_prov1(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec gamma, arma::vec beta,
                       int parallel=1, int threads=1, double tol=1e-8, int max_iter=10000,
                       double bound=10.0, bool message = true, bool backtrack = false) {

  int iter = 0, n = Z.n_rows, m = n_prov.n_elem, ind;
  double v;
  arma::vec gamma_obs(n);
  double crit = 100.0;
  if (message == true) {
    cout << "Implementing SerBIN algorithm (Rcpp) for fixed provider effects model ..." << endl;
  }

  double s = 0.01, t = 0.6; //only used for "backtrack = true"
  double lambda, d_loglkd, loglkd;
  arma::vec gamma_obs_tmp(n), gamma_tmp(m), beta_tmp(Z.n_cols);

  // while (iter <= max_iter) {
    // if (crit < tol) {
    //   break;
    // }
    iter++;
    gamma_obs = rep2(gamma, n_prov);
    arma::vec Z_beta = Z * beta;
    arma::vec p = 1 / (1 + exp(-gamma_obs-Z_beta));
    arma::vec Yp = Y - p, pq = p % (1-p);
    arma::vec score_gamma(m), info_gamma_inv(m);
    arma::mat info_betagamma(Z.n_cols,m);
    ind = 0;
    for (int i = 0; i < m; i++) {
      score_gamma(i) = sum(Yp(span(ind,ind+n_prov(i)-1)));
      info_gamma_inv(i) = 1 / sum(pq(span(ind,ind+n_prov(i)-1)));
      info_betagamma.col(i) =
        sum(Z.rows(ind,ind+n_prov(i)-1).each_col()%(p.subvec(ind,ind+n_prov(i)-1)%(1-p.subvec(ind,ind+n_prov(i)-1)))).t();
      ind += n_prov(i);
    }

    // List ret = List::create(_["score_gamma"]=score_gamma,
    //                         _["info_gamma_inv"]=info_gamma_inv,
    //                         _["info_betagamma"]=info_betagamma,
    //                         _["ind"]=ind);
    // return ret;

    arma::vec score_beta = Z.t() * Yp;
    arma::mat info_beta(Z.n_cols, Z.n_cols);
    // if (parallel==1) { // parallel
    //   info_beta = info_beta_omp(Z, pq, threads); // omp
    //   // info_beta = info_beta_tbb(Z, pq); // tbb
    // } else if (parallel==0) { // serial
      info_beta = Z.t() * (Z.each_col()%pq);
    // }
  // //
  List ret = List::create(_["score_beta"]=score_beta,
                          _["info_beta"]=info_beta);
  return ret;

    // arma::mat mat_tmp1 = trans(info_betagamma.each_row()%info_gamma_inv.t());
    // arma::mat schur_inv = inv_sympd(info_beta-mat_tmp1.t()*info_betagamma.t());
    //
    // arma::mat mat_tmp2 = mat_tmp1*schur_inv;
    // arma::vec d_gamma = info_gamma_inv%score_gamma + mat_tmp2*(mat_tmp1.t()*score_gamma-score_beta);
    // arma::vec d_beta = schur_inv*score_beta - mat_tmp2.t()*score_gamma;
  //
  //   v = 1.0; // initialize step size
  //   if (backtrack == true){
  //     loglkd = Loglkd2(Y, Z * beta, rep2(gamma, n_prov));
  //     gamma_tmp = gamma + v * d_gamma;
  //     gamma_obs_tmp = rep2(gamma_tmp, n_prov);
  //     arma::vec Z_beta_tmp = Z * (beta+v*d_beta);
  //     lambda = dot(score_gamma, d_gamma) + dot(score_beta, d_beta);
  //
  //     while (d_loglkd < s*v*lambda) {
  //       v = t*v;
  //       gamma_tmp = gamma + v * d_gamma;
  //       gamma_obs_tmp = rep2(gamma_tmp, n_prov);
  //       Z_beta_tmp = Z * (beta+v*d_beta);
  //       d_loglkd = Loglkd2(Y, Z_beta_tmp, gamma_obs_tmp) - loglkd;
  //     }
  //   }
  //
  //
  //   gamma += v * d_gamma;
  //   gamma = clamp(gamma, median(gamma)-bound, median(gamma)+bound);
  //   beta += v * d_beta;
  //   crit = norm(v*d_beta, "inf");
  //
  //   if (message == true) {
  //     cout << "Iter " << iter << ": Inf norm of running diff in est reg parm is " << scientific << setprecision(3) << crit << ";";
  //   }
  // // }
  // if (message == true) {
  //   cout << "serBIN (Rcpp) algorithm converged after " << iter << " iterations!" << endl;
  // }
  // List ret = List::create(_["gamma"]=gamma, _["beta"]=beta);
  // return ret;
}




