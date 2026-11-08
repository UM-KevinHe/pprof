#' Fit the logistic model with fixed provider effects
#'
#' Fits \eqn{\mathrm{logit}\, p_{ij} = \gamma_i + Z_{ij}^\top\beta} with a
#' binary or binomial outcome (\eqn{Y_{ij} \sim \mathrm{Bin}(N_{ij}, p_{ij})})
#' by the serial blockwise inversion Newton (SerBIN) or block ascent Newton
#' (BAN) algorithm, on the providers that [pp_data()] marked as included,
#' optionally with Firth's penalty.
#'
#' @param data an object from [pp_data()] with a binary or binomial outcome.
#' @param algorithm `"SerBIN"` (default; Wu et al. 2022) or `"BAN"` (He et
#'   al. 2013).
#' @param penalty `"none"` (default; maximum likelihood) or `"firth"`
#'   (penalized likelihood with Jeffreys' prior, Firth 1993), which keeps the
#'   effects of providers with no or all events finite.
#' @param max_iter maximum number of iterations. Default 10,000.
#' @param tol convergence tolerance. Default 1e-5.
#' @param bound provider effects are kept within `median(gamma) +/- bound`
#'   during the iterations. Default 10.
#' @param backtrack use a backtracking line search (not used with the Firth
#'   penalty). Default `TRUE`.
#' @param threads number of OpenMP threads for the information matrix
#'   (SerBIN and Firth). Results do not depend on it.
#' @param stop stopping rule: `"beta"` (default; infinity norm of the change
#'   in beta), `"relch"` (relative change in the log-likelihood), `"ratch"`
#'   (change relative to the total change since the start), `"all"` (all
#'   three below `tol`) or `"or"` (any of the three below `tol`). With the
#'   Firth penalty the rules read the penalized log-likelihood.
#' @param auc also compute the area under the ROC curve (binary outcomes).
#' @param message print the iteration history.
#'
#' @details SerBIN updates provider effects and covariate coefficients
#' jointly using the block structure of the information matrix; BAN
#' alternates between the two blocks. Providers with no or all events are
#' fitted (their effects sit at the bound); the tests treat them as the data
#' allow.
#'
#' With `penalty = "firth"` the penalized log-likelihood \eqn{\ell(\gamma,
#' \beta) + \frac12 \log\det I(\gamma, \beta)} is maximized by a joint
#' Newton step through the same block structure; the modified score adds
#' \eqn{h_{ij}(\frac12 - p_{ij})} to each row's residual, with \eqn{h_{ij}}
#' the row's leverage \eqn{N_{ij} p_{ij}(1 - p_{ij})\, x_{ij}^\top I^{-1}
#' x_{ij}}. The estimates stay finite when a provider has no or all events.
#' `loglik`, `aic` and `bic` report the unpenalized log-likelihood at the
#' penalized estimate; the penalized value is `penalized_loglik`.
#'
#' @return An object of class `logistic_fe`:
#' \item{coefficients}{a list with `gamma` (provider effects, named by
#'   provider) and `beta` (covariate coefficients, named by covariate).}
#' \item{fitted}{fitted probabilities of the included rows.}
#' \item{linear_predictor}{\eqn{Z\hat\beta} of the included rows.}
#' \item{loglik, aic, bic}{binomial log-likelihood (without the constant
#'   \eqn{\log\binom{N}{Y}}) and the information criteria.}
#' \item{penalized_loglik}{the penalized log-likelihood (Firth only).}
#' \item{auc}{the AUC, when requested (`NA` for binomial outcomes).}
#' \item{data}{the included rows, as sorted by [pp_data()].}
#' \item{vars, sizes, cutoff}{from [pp_data()] (`sizes` restricted to the
#'   included providers).}
#' \item{rows}{rows per provider, in the order of `coefficients$gamma`.}
#' \item{control}{the settings and the number of iterations.}
#' \item{converged}{whether the stopping rule was met before `max_iter`.}
#'
#' @references
#' Firth D (1993). Bias reduction of maximum likelihood estimates.
#' *Biometrika*, 80(1), 27-38.
#'
#' He K, Kalbfleisch JD, Li Y, Li Y (2013). Evaluating hospital readmission
#' rates in dialysis facilities; adjusting for hospital effects. *Lifetime
#' Data Analysis*, 19, 490-512.
#'
#' Wu W, Yang Y, Kang J, He K (2022). Improving large-scale estimation and
#' inference for profiling health care providers. *Statistics in Medicine*,
#' 41(15), 2840-2853.
#'
#' @seealso [test_providers()], [test_covariates()], [standardize()],
#'   [confint.logistic_fe()]
#' @examples
#' d <- sim_kecc(n_providers = 20, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_logistic_fe(pd)
#' fit
#' coef(fit)$beta
#' fit_firth <- fit_logistic_fe(pd, penalty = "firth")
#' @export
fit_logistic_fe <- function(data, algorithm = c("SerBIN", "BAN"), penalty = c("none", "firth"),
                            max_iter = 10000, tol = 1e-5, bound = 10, backtrack = TRUE,
                            threads = 1, stop = c("beta", "relch", "ratch", "all", "or"),
                            auc = FALSE, message = FALSE) {
  if (!inherits(data, "pp_data")) stop("`data` must come from pp_data().", call. = FALSE)
  if (outcome_of(data) == "continuous") {
    stop("The logistic model needs a binary or binomial outcome; see fit_linear_fe().",
         call. = FALSE)
  }
  algorithm <- match.arg(algorithm)
  penalty <- match.arg(penalty)
  stop <- match.arg(stop)
  if (!is.logical(backtrack) || length(backtrack) != 1) {
    stop("`backtrack` must be TRUE or FALSE.", call. = FALSE)
  }
  threads <- max(1L, as.integer(threads))
  v <- data$vars
  d <- data$data[data$data$included, , drop = FALSE]
  if (nrow(d) == 0) stop("No provider reaches the cutoff.", call. = FALSE)
  rownames(d) <- NULL
  prov <- d[[v$provider]]
  levels <- unique(prov)
  n_prov <- as.numeric(table(factor(prov, levels = levels)))
  y <- as.numeric(d[[v$y]])
  n <- if (is.null(v$n)) rep(1, nrow(d)) else as.numeric(d[[v$n]])
  z <- as.matrix(d[, v$x, drop = FALSE])
  storage.mode(z) <- "double"

  y_bar <- sum(y) / sum(n)
  gamma0 <- rep(log(y_bar / (1 - y_bar)), length(levels))
  beta0 <- rep(0, ncol(z))
  res <- if (penalty == "firth") {
    logis_firth_prov(n, y, z, n_prov, gamma0, beta0, threads, tol, max_iter, bound, message, stop)
  } else if (algorithm == "SerBIN") {
    logis_BIN_fe_prov(n, y, z, n_prov, gamma0, beta0, threads, tol, max_iter, bound, message,
                      backtrack, stop)
  } else {
    logis_fe_prov(n, y, z, n_prov, gamma0, beta0, backtrack, max_iter, bound, tol, message, stop)
  }
  gamma <- stats::setNames(as.numeric(res$gamma), levels)
  beta <- stats::setNames(as.numeric(res$beta), v$x)

  lp <- as.numeric(z %*% beta)
  eta <- rep(gamma, n_prov) + lp
  fitted <- stats::plogis(eta)
  loglik <- sum(eta * y - n * (pmax(eta, 0) + log1p(exp(-abs(eta)))))
  k <- length(gamma) + length(beta)

  out <- structure(list(coefficients = list(gamma = gamma, beta = beta),
                        fitted = fitted, linear_predictor = lp,
                        loglik = loglik, aic = -2 * loglik + 2 * k,
                        bic = -2 * loglik + log(nrow(d)) * k,
                        penalized_loglik = if (penalty == "firth") res$loglik_pen else NULL,
                        data = d, vars = v,
                        sizes = data$sizes[data$sizes$included, , drop = FALSE],
                        cutoff = data$cutoff, rows = n_prov,
                        control = list(algorithm = algorithm, penalty = penalty,
                                       max_iter = max_iter, tol = tol, bound = bound,
                                       backtrack = backtrack, threads = threads, stop = stop,
                                       iterations = res$iter),
                        converged = res$iter < max_iter),
                   class = "logistic_fe")
  rownames(out$sizes) <- NULL
  if (auc) {
    out$auc <- if (all(n == 1)) {
      auc_binary(y, fitted)
    } else {
      warning("The AUC is only defined for binary outcomes; returning NA.", call. = FALSE)
      NA_real_
    }
  }
  if (!out$converged) {
    warning("The fit stopped at max_iter = ", max_iter, " without meeting the stopping rule.",
            call. = FALSE)
  }
  out
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

#' @export
print.logistic_fe <- function(x, ...) {
  cat("Logistic model with fixed provider effects (",
      if (identical(x$control$penalty, "firth")) "Firth penalty" else x$control$algorithm,
      ")\n", sep = "")
  cat("  ", nrow(x$data), " rows, ", length(x$coefficients$gamma), " providers, ",
      length(x$coefficients$beta), " covariates; ",
      if (is.null(x$vars$n)) "binary outcome" else "binomial outcome", "\n", sep = "")
  cat("  log-likelihood ", format(x$loglik, digits = 6), "; iterations ", x$control$iterations,
      if (!x$converged) " (not converged)", "\n", sep = "")
  cat("  provider effects: median ", format(stats::median(x$coefficients$gamma), digits = 4),
      ", range ", paste(format(range(x$coefficients$gamma), digits = 4), collapse = " to "),
      "\n", sep = "")
  invisible(x)
}

#' @export
coef.logistic_fe <- function(object, ...) {
  object$coefficients
}

#' Model-based variance of a logistic fixed-effect fit
#'
#' The inverse of the observed information for \eqn{\beta} after profiling
#' out the provider effects. The model-based variances of the provider
#' effects are returned as the attribute `"var_gamma"`; see
#' [gamma_variance()] for the choice of variance used by the Wald test.
#'
#' @param object a model fitted by [fit_logistic_fe()].
#' @param ... unused.
#' @return A covariance matrix for `beta` with attribute `"var_gamma"`.
#' @export
vcov.logistic_fe <- function(object, ...) {
  m <- fit_matrices(object)
  v <- logis_fe_var(m$n, m$z, m$n_prov, object$coefficients$gamma, object$coefficients$beta)
  out <- v$var.beta
  dimnames(out) <- list(object$vars$x, object$vars$x)
  attr(out, "var_gamma") <- stats::setNames(as.numeric(v$var.gamma), m$levels)
  out
}

#' Variance of the provider effects for a Wald test
#'
#' Two model-based variances of a fixed provider effect \eqn{\hat\gamma_i}
#' are available. The `"case_mix"` variance is that of \eqn{\hat\gamma_i +
#' \bar x_i^\top\hat\beta}, the provider's outcome (log-odds or mean) at its
#' own average case mix \eqn{\bar x_i}; it does not depend on where the
#' covariates are centred and is the default of [test_providers()]. The
#' `"full"` variance is the \eqn{(\gamma_i,\gamma_i)} element of the inverse
#' information, \eqn{1/I_{ii} + \bar z_i^\top \mathrm{Var}(\hat\beta)\bar
#' z_i}; it is the variance of the outcome at \eqn{Z = 0} and changes with
#' the origin of the covariates, so it is only appropriate when the null
#' value is an external constant stated in the same coding.
#'
#' For the logistic model \eqn{\bar z_i} is the \eqn{N p (1 - p)}-weighted
#' mean covariate and the case-mix variance is \eqn{1/I_{ii} + (\bar z_i -
#' \bar x_i)^\top \mathrm{Var}(\hat\beta) (\bar z_i - \bar x_i)}; for the
#' linear model \eqn{\bar z_i = \bar x_i} and the case-mix variance is
#' exactly \eqn{\hat\sigma^2 / n_i}.
#'
#' @param fit a model fitted by [fit_logistic_fe()] or [fit_linear_fe()].
#' @param type `"case_mix"` (default) or `"full"`.
#' @param ... unused.
#' @return A named numeric vector, one variance per provider.
#' @export
gamma_variance <- function(fit, type = c("case_mix", "full"), ...) {
  UseMethod("gamma_variance")
}

#' @rdname gamma_variance
#' @export
gamma_variance.logistic_fe <- function(fit, type = c("case_mix", "full"), ...) {
  type <- match.arg(type)
  m <- fit_matrices(fit)
  v <- logis_fe_var(m$n, m$z, m$n_prov, fit$coefficients$gamma, fit$coefficients$beta)
  if (type == "full") return(stats::setNames(as.numeric(v$var.gamma), m$levels))
  p <- fit$fitted
  pq <- m$n * p * (1 - p)
  code <- provider_codes(m$prov, m$levels)
  i_gamma <- as.numeric(rowsum(pq, code, reorder = TRUE))
  j <- rowsum(pq * m$z, code, reorder = TRUE) / i_gamma
  xbar <- rowsum(m$n * m$z, code, reorder = TRUE) / as.numeric(rowsum(m$n, code, reorder = TRUE))
  d <- j - xbar
  stats::setNames(1 / i_gamma + rowSums((d %*% v$var.beta) * d), m$levels)
}

# The pieces of a fit every test reads: outcome, trials, covariate matrix,
# provider and patient of the included rows, provider levels and rows per
# provider (in the order of the effects).
fit_matrices <- function(fit) {
  d <- fit$data
  v <- fit$vars
  z <- as.matrix(d[, v$x, drop = FALSE])
  storage.mode(z) <- "double"
  list(y = as.numeric(d[[v$y]]),
       n = if (is.null(v$n)) rep(1, nrow(d)) else as.numeric(d[[v$n]]),
       z = z, prov = d[[v$provider]],
       patient = if (is.null(v$patient)) NULL else d[[v$patient]],
       levels = fit$sizes$provider, n_prov = fit$rows)
}

check_fit <- function(fit, class = "logistic_fe") {
  if (!inherits(fit, class)) stop("`fit` must be a ", class, " object.", call. = FALSE)
  invisible(TRUE)
}
