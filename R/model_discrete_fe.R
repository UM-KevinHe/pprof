#' Fit the discrete-time hazard model with fixed provider effects
#'
#' Fits \eqn{\mathrm{logit}\, h_{ijk} = \alpha_k + \gamma_i + Z_{ij}^\top\beta},
#' the hazard of the event in period \eqn{k} for subject \eqn{j} of provider
#' \eqn{i}, by maximum likelihood on one row per subject. The person-periods
#' are looped over inside the computation and never stored.
#'
#' @param data an object from [pp_data()] with a survival outcome: `time` is
#'   the subject's last period (the event period, or the last period observed
#'   without the event) and `event` the event indicator. One row per subject,
#'   without `start`, `weights` or `offset`.
#' @param max_iter maximum number of Newton iterations. Default 10,000.
#' @param tol convergence tolerance on the period effects and covariate
#'   coefficients. Default 1e-5.
#' @param bound provider effects are kept within `median(gamma) +/- bound`.
#'   Default 10.
#' @param backtrack use a backtracking line search. Default `TRUE`.
#' @param message print the iteration history.
#'
#' @details The periods are the distinct values of `time` among the events: a
#' period in which no subject has the event has a hazard of zero and adds
#' nothing to the likelihood. A subject is at risk in every period up to its
#' own time and has the event, if any, in its last one. The first period's
#' effect is fixed at zero, so \eqn{\gamma_i} is the log-odds of the hazard in
#' the first period at \eqn{Z = 0}; tests and standardized measures do not
#' depend on this choice.
#'
#' The likelihood is that of the logistic model on the person-period data,
#' and every result (coefficients, variances, tests, standardized measures)
#' equals what [fit_logistic_fe()] gives on the expanded data with the period
#' effects as covariates. The Newton step solves through the Schur complement
#' of the diagonal provider block of the information, as SerBIN does.
#'
#' @return An object of class `discrete_fe`:
#' \item{coefficients}{a list with `gamma` (provider effects), `alpha`
#'   (period effects, named by period, the first fixed at zero) and `beta`
#'   (covariate coefficients).}
#' \item{linear_predictor}{\eqn{Z\hat\beta} of the included subjects.}
#' \item{loglik, aic, bic}{log-likelihood and information criteria (the BIC
#'   counts person-periods).}
#' \item{periods}{the number of periods each subject is at risk.}
#' \item{time_grid}{the event periods.}
#' \item{person_periods}{the total number of person-periods at risk.}
#' \item{data, vars, sizes, cutoff, rows, control, converged}{as for
#'   [fit_logistic_fe()]; `rows` counts subjects.}
#'
#' @seealso [test.discrete_fe()], [standardize.discrete_fe()],
#'   [confint.discrete_fe()], [summary.discrete_fe()], and
#'   `var_select(model = "discrete")` for selecting the covariates.
#' @examples
#' s <- sim_kecc(n_providers = 20, provider_size = c(20, 40), rows_per_patient = c(1, 1),
#'               type = "survival", seed = 1)
#' s$month <- pmin(ceiling(s$stop * 12), 12)
#' s$event[s$stop * 12 > 12] <- 0
#' pd <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider", patient = "patient",
#'               time = "month", event = "event", message = FALSE)
#' fit <- fit_discrete_fe(pd)
#' fit
#' coef(fit)$alpha
#' @export
fit_discrete_fe <- function(data, max_iter = 10000, tol = 1e-5, bound = 10, backtrack = TRUE,
                            message = FALSE) {
  check_model_data(data, "survival", "The discrete-time model", "fit_logistic_fe()")
  v <- data$vars
  if (!is.null(v$start)) {
    stop("The discrete-time model takes one row per subject, without `start`.", call. = FALSE)
  }
  if (!is.logical(backtrack) || length(backtrack) != 1) {
    stop("`backtrack` must be TRUE or FALSE.", call. = FALSE)
  }
  d <- data$data[data$data$included, , drop = FALSE]
  if (nrow(d) == 0) stop("No provider reaches the cutoff.", call. = FALSE)
  rownames(d) <- NULL
  prov <- d[[v$provider]]
  levels <- unique(prov)
  n_prov <- as.numeric(table(factor(prov, levels = levels)))
  time <- as.numeric(d[[v$time]])
  event <- as.numeric(d[[v$event]])
  if (!all(event %in% c(0, 1))) stop("`event` must be 0 or 1.", call. = FALSE)
  grid <- sort(unique(time[event > 0]))
  if (length(grid) == 0) stop("There are no events.", call. = FALSE)
  t_index <- findInterval(time, grid)
  z <- as.matrix(d[, v$x, drop = FALSE])
  storage.mode(z) <- "double"
  n_period <- length(grid)
  if (n_period == 1 && ncol(z) == 0) {
    stop("With a single event period the model needs at least one covariate.", call. = FALSE)
  }

  at_risk <- rev(cumsum(rev(tabulate(t_index[t_index > 0], n_period))))
  events <- tabulate(t_index[event > 0], n_period)
  crude <- stats::qlogis(pmin(pmax(events / at_risk, 1e-6), 1 - 1e-6))
  gamma0 <- rep(crude[1], length(levels))
  alpha0 <- crude - crude[1]
  res <- discrete_fe_fit(t_index, event, z, n_prov, alpha0, gamma0, rep(0, ncol(z)), tol,
                         max(1L, as.integer(max_iter)), bound, backtrack, message)
  gamma <- stats::setNames(as.numeric(res$gamma), levels)
  alpha <- stats::setNames(as.numeric(res$alpha), format(grid, trim = TRUE))
  beta <- stats::setNames(as.numeric(res$beta), v$x)
  person_periods <- sum(t_index)
  k <- length(gamma) + n_period - 1 + length(beta)

  out <- structure(list(coefficients = list(gamma = gamma, alpha = alpha, beta = beta),
                        linear_predictor = as.numeric(z %*% beta),
                        loglik = res$loglik, aic = -2 * res$loglik + 2 * k,
                        bic = -2 * res$loglik + log(person_periods) * k,
                        periods = as.integer(t_index), time_grid = grid,
                        person_periods = person_periods, data = d, vars = v,
                        sizes = data$sizes[data$sizes$included, , drop = FALSE],
                        cutoff = data$cutoff, rows = n_prov,
                        control = list(max_iter = max_iter, tol = tol, bound = bound,
                                       backtrack = backtrack, iterations = res$iter),
                        converged = res$iter < max_iter),
                   class = "discrete_fe")
  rownames(out$sizes) <- NULL
  if (!out$converged) {
    warning("The fit stopped at max_iter = ", max_iter, " without converging.", call. = FALSE)
  }
  out
}

#' @export
print.discrete_fe <- function(x, ...) {
  cat("Discrete-time hazard model with fixed provider effects\n")
  cat("  ", nrow(x$data), " subjects, ", x$person_periods, " person-periods, ",
      length(x$coefficients$gamma), " providers, ", length(x$time_grid), " event periods, ",
      length(x$coefficients$beta), " covariates\n", sep = "")
  cat("  log-likelihood ", format(x$loglik, digits = 6), "; iterations ", x$control$iterations,
      if (!x$converged) " (not converged)", "\n", sep = "")
  cat("  provider effects: median ", format(stats::median(x$coefficients$gamma), digits = 4),
      ", range ", paste(format(range(x$coefficients$gamma), digits = 4), collapse = " to "),
      "\n", sep = "")
  invisible(x)
}

#' @export
coef.discrete_fe <- function(object, ...) {
  object$coefficients
}

#' Model-based variance of a discrete-time fixed-effect fit
#'
#' The inverse of the observed information for the covariate coefficients
#' after profiling out the provider effects (the period effects are estimated
#' jointly with them). The model-based variances of the provider effects are
#' returned as the attribute `"var_gamma"`, those of the period effects
#' (periods 2 onwards) as `"var_alpha"`.
#'
#' @param object a model fitted by [fit_discrete_fe()].
#' @param ... unused.
#' @return A covariance matrix for `beta` with attributes `"var_gamma"` and
#'   `"var_alpha"`.
#' @export
vcov.discrete_fe <- function(object, ...) {
  inf <- discrete_information(object)
  p <- length(object$coefficients$beta)
  k <- length(object$time_grid) - 1
  out <- inf$theta_cov[k + seq_len(p), k + seq_len(p), drop = FALSE]
  dimnames(out) <- list(object$vars$x, object$vars$x)
  attr(out, "var_gamma") <- gamma_variance(object, "full")
  attr(out, "var_alpha") <- stats::setNames(diag(inf$theta_cov)[seq_len(k)],
                                            names(object$coefficients$alpha)[-1])
  out
}

#' @rdname gamma_variance
#' @export
gamma_variance.discrete_fe <- function(fit, type = c("case_mix", "full"), ...) {
  type <- match.arg(type)
  inf <- discrete_information(fit)
  jbar <- inf$i_gamma_theta / inf$i_gamma
  d <- if (type == "full") jbar else jbar - inf$xbar
  stats::setNames(1 / inf$i_gamma + rowSums((d %*% inf$theta_cov) * d), fit$sizes$provider)
}

# Information blocks at the estimate and the covariance of theta =
# (alpha_2..alpha_K, beta) with the provider effects profiled out.
discrete_information <- function(fit) {
  m <- discrete_matrices(fit)
  cf <- fit$coefficients
  inf <- discrete_fe_info(m$t, m$event, m$z, m$n_prov, cf$alpha, cf$gamma, cf$beta)
  i_gamma <- as.numeric(inf$i_gamma)
  schur <- inf$i_theta - crossprod(inf$i_gamma_theta, inf$i_gamma_theta / i_gamma)
  list(i_gamma = i_gamma, i_gamma_theta = inf$i_gamma_theta, xbar = inf$xbar,
       theta_cov = solve(schur), person_periods = as.numeric(inf$person_periods))
}

# The pieces of a discrete-time fit the tests read: periods at risk, events,
# covariates, provider and patient of the included subjects, provider levels
# and subjects per provider (in the order of the effects).
discrete_matrices <- function(fit) {
  d <- fit$data
  v <- fit$vars
  z <- as.matrix(d[, v$x, drop = FALSE])
  storage.mode(z) <- "double"
  list(t = fit$periods, event = as.numeric(d[[v$event]]), z = z, prov = d[[v$provider]],
       patient = if (is.null(v$patient)) NULL else d[[v$patient]],
       levels = fit$sizes$provider, n_prov = fit$rows)
}

# Per subject, the expected events and their variance at provider effects
# `gamma_row` (one per subject).
discrete_sums_at <- function(fit, m, gamma_row) {
  discrete_subject_sums(m$t, fit$coefficients$alpha, gamma_row + fit$linear_predictor)
}

# Cluster codes for the clustered variances: the patient within provider, or
# the subject when there is no patient column.
discrete_clusters <- function(m) {
  if (is.null(m$patient)) {
    list(code = seq_along(m$prov), first = rep(TRUE, length(m$prov)))
  } else {
    cluster_codes(m$prov, m$patient)
  }
}
