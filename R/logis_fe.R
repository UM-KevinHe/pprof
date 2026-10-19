#' Fit the fixed provider effects logistic model
#'
#' Fits \eqn{logit(p_{ij}) = \gamma_i + Z_{ij}^T\beta} with a binary or
#' binomial outcome (\eqn{Y_{ij} \sim Bin(N_{ij}, p_{ij})}) by the Serial
#' blockwise inversion Newton (SerBIN) or block ascent Newton (BAN)
#' algorithm.
#'
#' @param data.prep an object from [fe_data_prep()].
#' @param algorithm `"SerBIN"` (default; Wu et al. 2022) or `"BAN"` (He et al. 2013).
#' @param max.iter maximum number of iterations. Default 10,000.
#' @param tol convergence tolerance. Default 1e-5. The KECC measures use
#'   1e-5 (SFR, SRR), 1e-10 (AOH) and 1e-15 (PPPW).
#' @param bound provider effects are kept within `median(gamma) +/- bound`
#'   during the iterations. Default 10.
#' @param backtrack a logical; use backtracking line search. Default `TRUE`.
#' @param Rcpp kept for backward compatibility; the fit always uses the
#'   compiled code.
#' @param threads number of OpenMP threads for the information matrix
#'   (SerBIN only). Default 1.
#' @param AUC a logical; also report the area under the ROC curve (binary
#'   outcomes only). Default `FALSE`.
#' @param message a logical; print the iteration history. Default `FALSE`.
#' @param stop stopping rule: `"beta"` (default; infinity norm of the change
#'   in beta), `"relch"` (relative change in log-likelihood), `"ratch"`
#'   (change relative to the total change since the start), `"all"` (all
#'   three below `tol`) or `"or"` (any of the three below `tol`).
#'
#' @details SerBIN updates provider effects and covariate coefficients jointly
#' using the block structure of the information matrix; BAN alternates between
#' the two blocks. SerBIN is recommended; BAN is a fallback when the SerBIN
#' information matrix is numerically singular.
#'
#' @return An object of class `logis_fe`:
#' \item{beta}{covariate coefficients (matrix with one column).}
#' \item{gamma}{provider effects (matrix with one column, rows named by provider).}
#' \item{linear_pred}{\eqn{Z\beta} for every included row.}
#' \item{pred}{fitted probabilities \eqn{p_{ij}}.}
#' \item{neg2Loglkd, AIC, BIC}{fit statistics (binomial likelihood without the
#'   constant \eqn{\log\binom{N}{Y}}).}
#' \item{obs, prov}{outcome and provider of every included row.}
#' \item{df.prov}{provider-level observed events and \eqn{\hat\gamma}.}
#' \item{char_list}{column names from `data.prep`.}
#' \item{data_include}{the rows used for fitting.}
#' \item{control}{the settings used for the fit.}
#' \item{AUC}{if requested.}
#'
#' @references
#' He K, Kalbfleisch J, Li Y, et al. (2013) Evaluating hospital readmission
#' rates in dialysis facilities; adjusting for hospital effects. \emph{Lifetime
#' Data Analysis}, 19: 490-512.
#'
#' Wu W, Yang Y, Kang J, He K. (2022) Improving large-scale estimation and
#' inference for profiling health care providers. \emph{Statistics in
#' Medicine}, 41(15): 2840-2853.
#'
#' @export
logis_fe <- function(data.prep, algorithm = "SerBIN", max.iter = 10000, tol = 1e-5, bound = 10,
                     backtrack = TRUE, Rcpp = TRUE, threads = 1, AUC = FALSE, message = FALSE,
                     stop = "beta") {
  if (missing(data.prep)) stop("Argument 'data.prep' is required!", call. = FALSE)
  if (!inherits(data.prep, "data_prep")) stop("Object 'data.prep' should be generated from 'fe_data_prep' function!", call. = FALSE)
  if (!is.logical(backtrack)) stop("Argument 'backtrack' NOT as required!", call. = FALSE)
  if (!stop %in% c("beta", "relch", "ratch", "all", "or")) stop("Argument 'stop' NOT as required!", call. = FALSE)
  threads <- max(1L, as.integer(threads))

  data <- data.prep$data
  N.char <- data.prep$char_list$N.char
  Y.char <- data.prep$char_list$Y.char
  prov.char <- data.prep$char_list$prov.char
  Z.char <- data.prep$char_list$Z.char

  # only providers with "included == 1" (provider size cutoff) enter the fit
  data <- data[data$included == 1, ]
  n.prov <- sapply(split(data[, Y.char], data[, prov.char]), length) # provider-specific number of rows
  Z <- as.matrix(data[, Z.char])
  N <- trials_of(data, N.char)
  Y_bar <- sum(data[, Y.char]) / sum(N)
  gamma.prov <- rep(log(Y_bar / (1 - Y_bar)), length(n.prov))
  beta <- rep(0, NCOL(Z))

  if (algorithm == "SerBIN") {
    ls <- logis_BIN_fe_prov(as.matrix(N), as.matrix(data[, Y.char]), Z, n.prov, gamma.prov, beta,
                            threads, tol, max.iter, bound, message, backtrack, stop)
  } else if (algorithm == "BAN") {
    ls <- logis_fe_prov(as.matrix(N), as.matrix(data[, Y.char]), Z, n.prov, gamma.prov, beta,
                        backtrack, max.iter, bound, tol, message, stop)
  } else {
    stop("Argument 'algorithm' NOT as required!", call. = FALSE)
  }
  gamma.prov <- as.numeric(ls$gamma)
  beta <- as.numeric(ls$beta)

  gamma.obs <- rep(gamma.prov, n.prov)
  neg2Loglkd <- -2 * sum((gamma.obs + Z %*% beta) * data[, Y.char] - N * log(1 + exp(gamma.obs + Z %*% beta)))
  AIC <- neg2Loglkd + 2 * (length(gamma.prov) + length(beta))
  BIC <- neg2Loglkd + log(nrow(data)) * (length(gamma.prov) + length(beta))

  df.prov <- data.frame(Obs_provider = sapply(split(data[, Y.char], data[, prov.char]), sum),
                        gamma_est = gamma.prov) # original gamma-hat, for internal use
  linear_pred <- Z %*% beta
  pred <- as.numeric(stats::plogis(gamma.obs + linear_pred))

  beta <- matrix(beta)
  gamma.prov <- matrix(gamma.prov)
  dimnames(beta) <- list(Z.char, "beta")
  dimnames(gamma.prov) <- list(names(n.prov), "gamma")

  return_ls <- structure(list(beta = beta,
                              gamma = gamma.prov,       # provider effects
                              linear_pred = linear_pred, # linear predictor
                              pred = pred,              # predicted probability
                              neg2Loglkd = neg2Loglkd,
                              AIC = AIC,
                              BIC = BIC,
                              obs = data[, Y.char],     # row-level outcome
                              prov = data[, prov.char]),
                         class = "logis_fe")
  if (AUC) {
    return_ls$AUC <- if (all(N == 1)) auc_binary(data[, Y.char], pred) else {
      warning("AUC is only defined for binary outcomes (N = 1); returning NA.", call. = FALSE)
      NA_real_
    }
  }
  return_ls$df.prov <- df.prov
  return_ls$char_list <- data.prep$char_list
  return_ls$data_include <- data
  return_ls$control <- list(algorithm = algorithm, max.iter = max.iter, tol = tol, bound = bound,
                            backtrack = backtrack, threads = threads, stop = stop, iter = ls$iter)
  return_ls
}

# Empirical AUC (Mann-Whitney), with the direction chosen as pROC does by
# default: cases are expected to have higher predictions unless their median
# prediction is lower than the controls'.
auc_binary <- function(y, pred) {
  cases <- pred[y == 1]
  controls <- pred[y == 0]
  n1 <- length(cases)
  n0 <- length(controls)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(c(cases, controls))
  auc <- (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  if (stats::median(controls) > stats::median(cases)) auc <- 1 - auc
  auc
}

#' Extract coefficients from a `logis_fe` fit
#'
#' @param object a model fitted by [logis_fe()].
#' @param ... unused.
#' @return A list with `gamma` (provider effects) and `beta` (covariate coefficients).
#' @export
coef.logis_fe <- function(object, ...) {
  check_fit(object)
  list(gamma = object$gamma, beta = object$beta)
}

#' Model-based variance of the covariate coefficients
#'
#' Inverse of the observed information for \eqn{\beta}, after profiling out
#' the provider effects. The model-based variances of the provider effects
#' are returned as attribute `"var.gamma"`.
#'
#' @param object a model fitted by [logis_fe()].
#' @param ... unused.
#' @return A covariance matrix for `beta`.
#' @export
vcov.logis_fe <- function(object, ...) {
  check_fit(object)
  data <- object$data_include
  cl <- object$char_list
  n.prov <- sapply(split(data[, cl$Y.char], data[, cl$prov.char]), length)
  N <- trials_of(data, cl$N.char)
  v <- logis_fe_var(N, as.matrix(data[, cl$Z.char]), n.prov,
                    object$df.prov$gamma_est, as.numeric(object$beta))
  out <- v$var.beta
  dimnames(out) <- list(cl$Z.char, cl$Z.char)
  attr(out, "var.gamma") <- stats::setNames(as.numeric(v$var.gamma), names(n.prov))
  out
}
