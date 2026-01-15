//#define ARMA_NO_DEBUG
#define STRICT_R_HEADERS // needed on Windows, not on macOS
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#define ARMA_DONT_USE_OPENMP
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include "header.h"
#include "myomp.h"
#include <iostream>
#include <chrono>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
using namespace RcppParallel;
using namespace Rcpp;
using namespace std;
using namespace arma;

//` @importFrom RcppParallel RcppParallelLibs


// [[Rcpp::export]]
List logis_LR_prov(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec gamma, arma::vec beta,
                   int n_obs, int m, const int &threads, double tol = 1e-5, int max_iter = 100,
                   double bound = 10.0, bool message = true, //bool backtrack = false,
                   const std::string stop = "beta") {

  arma::ivec indices(m + 1);
  int ind = 0;
  for (int i = 0; i < m; i++) {
    indices(i) = ind;
    ind += n_prov(i);
  }
  indices(m) = ind;

  arma::vec gamma_obs = rep(gamma, n_prov);
  double loglkd_init = Loglkd(Y, Z*beta, gamma_obs);

  int chunk_id;

  int beta_size = Z.n_cols;
  arma::vec score_beta(beta_size);
  arma::vec score_gamma(m);
  arma::vec info_gamma_inv(m);
  arma::mat J(beta_size, m);
  arma::mat J2(beta_size, m);  //I^*(\beta)^{-1} * J; where "I^*(\beta)^{-1}" is "schur_inv".
  arma::mat schur(beta_size, beta_size);  //I^*(\beta)
  arma::mat schur_inv(beta_size, beta_size);   //I^*(\beta)^{-1}
  arma::vec p(n_obs);
  arma::vec Yp(n_obs);
  arma::vec pq(n_obs);
  arma::vec H(beta_size);

  arma::vec old_gamma = gamma;


  omp_set_num_threads(threads);
  int iter = 0;
  double crit = 100.0, d_loglkd = 0.0, old_loglkd = loglkd_init, loglkd = loglkd_init;

  #pragma omp parallel shared(indices, iter, crit)
  {
    arma::mat local_schur(beta_size, beta_size, arma::fill::zeros);
    arma::vec local_H(beta_size, arma::fill::zeros);
    while (iter < max_iter && crit > tol) {
        #pragma omp single
        {
            chunk_id = 0;
            schur.zeros();
        }
        #pragma omp barrier

        while (true){
            int my_chunk;
            #pragma omp critical (get_chunk_id_1)
            {
                if (chunk_id < m) {
                    my_chunk = chunk_id;
                    chunk_id++;
                } else {
                    my_chunk = -1;
                }
            }
            if (my_chunk == -1) {
                break;
            }

            // Get the chunk of Z and Y
            // these "subviews" provide a window into the original matrices without copying data.
            const arma::mat &TDS_Z = Z.rows(indices(my_chunk), indices(my_chunk + 1) - 1);
            const arma::vec &TDS_Y = Y.subvec(indices(my_chunk), indices(my_chunk + 1) - 1);

            arma::vec TDS_gamma_obs = rep(gamma(my_chunk), n_prov(my_chunk)); // get the chunk of gamma
            arma::vec TDS_Z_beta = TDS_Z * beta;
            arma::vec TDS_p = 1 / (1 + exp(-TDS_gamma_obs-TDS_Z_beta));
            arma::vec TDS_Yp = TDS_Y - TDS_p;
            arma::vec TDS_pq = TDS_p % (1 - TDS_p);
            if (arma::any(TDS_pq == 0)) {
                TDS_pq.replace(0, 1e-10);
            }

            p.subvec(indices(my_chunk), indices(my_chunk + 1) - 1) = TDS_p;
            Yp.subvec(indices(my_chunk), indices(my_chunk + 1) - 1) = TDS_Yp;
            pq.subvec(indices(my_chunk), indices(my_chunk + 1) - 1) = TDS_pq;

            // 1. compute I(gamma, gamma) & I(beta, gamma) & I(beta, beta)
            double TDS_info_gamma_inv = 1 / sum(TDS_pq);
            info_gamma_inv(my_chunk) = TDS_info_gamma_inv;

            arma::vec TDS_info_betagamma = sum(TDS_Z.each_col() % TDS_pq).t();  //p*1 ("B_i" in notes)

            arma::vec TDS_J = TDS_info_betagamma * TDS_info_gamma_inv; //p*1
            J.col(my_chunk) = TDS_J;

            arma::mat TDS_info_beta = TDS_Z.t() * (TDS_Z.each_col() % TDS_pq);  //p*p ("C_i" in notes)
            arma::mat TDS_schur = TDS_info_beta - TDS_J * TDS_info_betagamma.t(); //p*p (S_i)
            local_schur += TDS_schur;
        }
        #pragma omp critical (update_schur)
        {
            schur += local_schur;  //sum over all threads to get the global schur matrix
        }
        local_schur.zeros();  // Reset local_schur for the next iteration
        #pragma omp barrier

        // compute I(\theta)^-1
        #pragma omp single
        {
            schur_inv = inv_sympd(schur);  //S^-1
            chunk_id = 0;
            H.zeros();
        }
        #pragma omp barrier

        while (true){
            int my_chunk;
            #pragma omp critical (get_chunk_id_2)
            {
                if (chunk_id < m) {
                    my_chunk = chunk_id;
                    chunk_id++;
                } else {
                    my_chunk = -1;
                }
            }
            if (my_chunk == -1) {
                break;
            }

            arma::vec TDS_J2 = schur_inv * J.col(my_chunk); //p*1 (J2_i)
            J2.col(my_chunk) = TDS_J2;

            const arma::vec &TDS_Yp = Yp.subvec(indices(my_chunk), indices(my_chunk + 1) - 1); //n_i*1
            const arma::mat &TDS_Z = Z.rows(indices(my_chunk), indices(my_chunk + 1) - 1);  //n_i*p

            arma::vec TDS_score_beta = TDS_Z.t() * TDS_Yp; //p*1

            double TDS_score_gamma = sum(TDS_Yp);
            score_gamma(my_chunk) = TDS_score_gamma;

            arma::vec TDS_G = J.col(my_chunk) * TDS_score_gamma; //p*1
            local_H += (TDS_G - TDS_score_beta); //H_i
        }
        #pragma omp critical (update_H)
        {
            H += local_H; //sum over all threads to get the global H
        }
        local_H.zeros();
        #pragma omp barrier

        #pragma omp single
        {
            chunk_id = 0;
        }
        #pragma omp barrier

        // 1. update gamma
        while (true){
            int my_chunk;
            #pragma omp critical (get_chunk_id_3)
            {
                if (chunk_id < m) {
                    my_chunk = chunk_id;
                    chunk_id++;
                } else {
                    my_chunk = -1;
                }
            }
            if (my_chunk == -1) {
                break;
            }

            double d_gamma_TDS = info_gamma_inv(my_chunk) * score_gamma(my_chunk) + arma::as_scalar(J2.col(my_chunk).t() * H);
            gamma(my_chunk) += d_gamma_TDS;
        }
        #pragma omp barrier


        // 2. update beta
        #pragma omp single
        {
          gamma = clamp(gamma, median(gamma) - bound, median(gamma) + bound);
          arma::vec d_beta = - schur_inv * H;
          beta += d_beta; //update beta

          if (stop == "beta"){
            crit = norm(d_beta, "inf");
            if (message == true) {
                Rcout << "Iter " << iter << ": Inf norm of running diff in est reg parm is " << setprecision(3) << scientific << crit << ";" << endl;
            }
          } else {
            gamma_obs = rep(gamma, n_prov);
            loglkd = Loglkd(Y, Z*beta, gamma_obs);  //current log-likelihood
            d_loglkd = loglkd - old_loglkd;
            old_loglkd = loglkd;
            if (stop == "relch") {
                crit = abs(d_loglkd/loglkd);
                if (message == true) {
                    Rcout << "Iter " << iter << ": Relative change in est log likelihood is " << setprecision(3) << scientific << crit << ";" << endl;
                }
            } else if (stop == "ratch") {
                crit = abs(d_loglkd/(loglkd-loglkd_init));
                if (message == true) {
                    Rcout << "Iter " << iter << ": Adjusted relative change in est log likelihood is " << setprecision(3) << scientific << crit << ";" << endl;
                }
            } else if (stop == "all") {
                arma::vec crits(3);
                crits(0) = norm(d_beta, "inf");
                crits(1) = abs(d_loglkd/loglkd);
                crits(2) = crit = abs(d_loglkd/(loglkd-loglkd_init));
                crit = crits.max();
                if (message == true) {
                    Rcout << "Iter " << iter << ": Maximum criterion across all checks is " << setprecision(3) << scientific << crit << ";" << endl;
                }
            } else if (stop == "or") {
                arma::vec crits(3);
                crits(0) = norm(d_beta, "inf");
                crits(1) = abs(d_loglkd/loglkd);
                crits(2) = crit = abs(d_loglkd/(loglkd-loglkd_init));
                crit = crits.min();
                if (message == true) {
                    Rcout << "Iter " << iter << ": Minimum criterion across all checks is " << setprecision(3) << scientific << crit << ";" << endl;
                }
            } else {
                Rcpp::stop("Argument 'stop' NOT as required!");
            }
          }

          chunk_id = 0;
          iter++;
        }
        #pragma omp barrier

    }
  }

  if (message == true) {
    std::cout << "Algorithm with " << threads <<  " cores converged after " << iter << " iterations." << endl;
  }
  List ret = List::create(_["gamma"]=gamma, _["beta"]=beta);
  return ret;
}


// [[Rcpp::export]]
List logis_firth_prov(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec gamma, arma::vec beta,
                      int n_obs, int m, const int &threads, double tol = 1e-5, int max_iter = 100,
                      double bound = 10.0, bool message = true) {
  if (message == true) {
    cout << "Implementing firth-corrected fixed provider effects model (Rcpp) ..." << endl;
  }

  // make "provider" as the basic unit of parallelization (so ideally maximum number of threads = number of providers)
  arma::ivec indices(m + 1);
  int ind = 0;
  for (int i = 0; i < m; i++) {
    indices(i) = ind;
    ind += n_prov(i);
  }
  indices(m) = ind;

  int chunk_id; //each thread will process a chunk of data

  int beta_size = Z.n_cols;
  arma::vec score_beta(beta_size);
  arma::vec score_gamma(m);
  arma::vec info_gamma_inv(m);
  arma::mat J(beta_size, m);
  arma::mat J2(beta_size, m);  //I^*(\beta)^{-1} * J; where "I^*(\beta)^{-1}" is "schur_inv".
  arma::mat schur(beta_size, beta_size);  //I^*(\beta)
  arma::mat schur_inv(beta_size, beta_size);   //I^*(\beta)^{-1}
  arma::vec p(n_obs);
  arma::vec Yp(n_obs);
  arma::vec pq(n_obs);
  arma::vec H(beta_size);

  arma::vec old_gamma = gamma;


  omp_set_num_threads(threads);
  int iter = 0;
  double crit = 100.0;

  #pragma omp parallel shared(indices, iter, crit)
  {
    arma::mat local_schur(beta_size, beta_size, arma::fill::zeros);
    arma::vec local_H(beta_size, arma::fill::zeros);
    while (iter < max_iter && crit > tol) {
        #pragma omp single
        {
            chunk_id = 0;
            schur.zeros();
        }
        #pragma omp barrier

        while (true){
            int my_chunk;
            #pragma omp critical (get_chunk_id_1)
            {
                if (chunk_id < m) {
                    my_chunk = chunk_id;
                    chunk_id++;
                } else {
                    my_chunk = -1;
                }
            }
            if (my_chunk == -1) {
                break;
            }

            const arma::mat &TDS_Z = Z.rows(indices(my_chunk), indices(my_chunk + 1) - 1);
            const arma::vec &TDS_Y = Y.subvec(indices(my_chunk), indices(my_chunk + 1) - 1);

            arma::vec TDS_gamma_obs = rep(gamma(my_chunk), n_prov(my_chunk)); // get the chunk of gamma
            arma::vec TDS_Z_beta = TDS_Z * beta;
            arma::vec TDS_p = 1 / (1 + exp(-TDS_gamma_obs-TDS_Z_beta));
            arma::vec TDS_Yp = TDS_Y - TDS_p;
            arma::vec TDS_pq = TDS_p % (1 - TDS_p);
            if (arma::any(TDS_pq == 0)) {
                TDS_pq.replace(0, 1e-10);
            }

            p.subvec(indices(my_chunk), indices(my_chunk + 1) - 1) = TDS_p;
            Yp.subvec(indices(my_chunk), indices(my_chunk + 1) - 1) = TDS_Yp;
            pq.subvec(indices(my_chunk), indices(my_chunk + 1) - 1) = TDS_pq;

            // 1. compute I(gamma, gamma) & I(beta, gamma) & I(beta, beta)
            double TDS_info_gamma_inv = 1 / sum(TDS_pq);
            info_gamma_inv(my_chunk) = TDS_info_gamma_inv;

            arma::vec TDS_info_betagamma = sum(TDS_Z.each_col() % TDS_pq).t();  //p*1 ("B_i" in notes)

            arma::vec TDS_J = TDS_info_betagamma * TDS_info_gamma_inv; //p*1
            J.col(my_chunk) = TDS_J;

            arma::mat TDS_info_beta = TDS_Z.t() * (TDS_Z.each_col() % TDS_pq);  //p*p ("C_i" in notes)
            arma::mat TDS_schur = TDS_info_beta - TDS_J * TDS_info_betagamma.t(); //p*p (S_i)
            local_schur += TDS_schur;
        }
        #pragma omp critical (update_schur)
        {
            schur += local_schur;  //sum over all threads to get the global schur matrix
        }
        local_schur.zeros();  // Reset local_schur for the next iteration
        #pragma omp barrier

        // compute I(\theta)^-1
        #pragma omp single
        {
            schur_inv = inv_sympd(schur);  //S^-1
            chunk_id = 0;
            H.zeros();
        }
        #pragma omp barrier

        while (true){
            int my_chunk;
            #pragma omp critical (get_chunk_id_2)
            {
                if (chunk_id < m) {
                    my_chunk = chunk_id;
                    chunk_id++;
                } else {
                    my_chunk = -1;
                }
            }
            if (my_chunk == -1) {
                break;
            }

            arma::vec TDS_J2 = schur_inv * J.col(my_chunk); //p*1 (J2_i)
            J2.col(my_chunk) = TDS_J2;

            const arma::vec &TDS_Yp = Yp.subvec(indices(my_chunk), indices(my_chunk + 1) - 1); //n_i*1
            const arma::mat &TDS_Z = Z.rows(indices(my_chunk), indices(my_chunk + 1) - 1);  //n_i*p
            arma::vec TDS_YpA;

            const arma::vec &TDS_p = p.subvec(indices(my_chunk), indices(my_chunk + 1) - 1);
            const arma::vec &TDS_pq = pq.subvec(indices(my_chunk), indices(my_chunk + 1) - 1);

            double TDS_prod = arma::as_scalar(J.col(my_chunk).t() * TDS_J2); // 1*1
            double TDS_diag_prod = info_gamma_inv(my_chunk) + TDS_prod;
            arma::vec TDS_c1 = rep(TDS_diag_prod, n_prov(my_chunk)); // n_i*1

            arma::vec TDS_c2 = -TDS_Z * TDS_J2; // n_i*1
            int n_i = n_prov(my_chunk);
            arma::vec TDS_c3(n_i);
            for (int i = 0; i < n_i; i++) {
                TDS_c3(i) = arma::as_scalar(TDS_Z.row(i) * schur_inv * TDS_Z.row(i).t());
            }
            TDS_YpA = TDS_Yp + TDS_pq % (TDS_c1 + TDS_c2 + TDS_c2 + TDS_c3) % (0.5 - TDS_p); // n_i*1
            arma::vec TDS_score_beta = TDS_Z.t() * TDS_YpA; //p*1

            double TDS_score_gamma = sum(TDS_YpA);
            score_gamma(my_chunk) = TDS_score_gamma;

            arma::vec TDS_G = J.col(my_chunk) * TDS_score_gamma; //p*1
            local_H += (TDS_G - TDS_score_beta); //H_i
        }
        #pragma omp critical (update_H)
        {
            H += local_H; //sum over all threads to get the global H
        }
        local_H.zeros();
        #pragma omp barrier

        #pragma omp single
        {
          arma::vec d_beta = - schur_inv * H;
          beta += d_beta; //update beta
          crit = norm(d_beta, "inf");
          chunk_id = 0;
          iter++;
        }
        #pragma omp barrier

        // update gamma
        while (true){
            int my_chunk;
            #pragma omp critical (get_chunk_id_3)
            {
                if (chunk_id < m) {
                    my_chunk = chunk_id;
                    chunk_id++;
                } else {
                    my_chunk = -1;
                }
            }
            if (my_chunk == -1) {
                break;
            }

            double d_gamma_TDS = info_gamma_inv(my_chunk) * score_gamma(my_chunk) + arma::as_scalar(J2.col(my_chunk).t() * H);
            gamma(my_chunk) += d_gamma_TDS;
        }
        #pragma omp barrier

        #pragma omp single
        {
          gamma = clamp(gamma, median(gamma) - bound, median(gamma) + bound);
        }
    }
  }

  if (message == true) {
    std::cout << "Algorithm with " << threads <<  " cores converged after " << iter << " iterations." << endl;
  }
  List ret = List::create(_["gamma"]=gamma, _["beta"]=beta);
  return ret;
}

