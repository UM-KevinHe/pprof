#' Cox model with provider strata or a pooled baseline
#'
#' Fits the Cox proportional hazards model of the survival measures (SMR,
#' SHR, STrR) through [survival::coxph()]. By default the provider is a
#' stratification variable, so that each provider has its own baseline
#' hazard and the covariate effects are estimated within providers (He and
#' Schaubel 2014); the national baseline is then the Breslow estimator over
#' every row with the fitted linear predictor as an offset. This is the
#' two-stage model of the KECC measures: stage one estimates the covariate
#' effects in the stratified model, stage two estimates the optional
#' covariates `x2` with the stage-one linear predictor as offset and gives
#' the baseline. With `strata = "none"` one pooled model is fitted, without
#' provider effects.
#'
#' @param data an object from [pp_data()] with a survival outcome (`time`,
#'   `event`, optional `start`, `weights` and `offset`).
#' @param strata `"provider"` (default) or `"none"`.
#' @param x2 optional names of covariate columns fitted in the second stage
#'   with the first-stage linear predictor as offset (calendar-year
#'   indicators in the KECC measures, say); they must be in `data`'s
#'   covariates `x` or be further numeric columns of the original data frame
#'   named in `x`. `NULL` (default) for no second-stage covariates.
#' @param ties `"breslow"` (default, the KECC convention) or `"efron"`.
#' @param message print the fit. Default `FALSE`.
#'
#' @details Rows with `included = FALSE` are left out; set `cutoff = 0` in
#' [pp_data()] to fit every provider, as the measures do, and choose the
#' providers to report in [test()].
#'
#' When `data` has a patient column it is passed as `cluster` to
#' [survival::coxph()], so that the robust (patient-clustered) covariance is
#' available from [vcov.cox()] and [summary()].
#'
#' The expected number of events of every row is
#' \eqn{E_i = \exp(\eta_i)\,[\Lambda_0(t_i) - \Lambda_0(s_i)]} with
#' \eqn{\eta_i} the full linear predictor (covariates of both stages and the
#' offset) and \eqn{\Lambda_0} the Breslow estimator of the national
#' cumulative baseline hazard at \eqn{\eta = 0}, weighted by the case weights;
#' the row's martingale residual is \eqn{\delta_i - E_i}. For `ties =
#' "breslow"` these are the residuals [survival::residuals.coxph()] returns
#' for the second-stage model, and `expected` is what the KECC measures
#' compute as `event - residual`.
#'
#' @return An object of class `cox`: a list with
#' \item{coefficients}{`beta` (first stage) and `beta2` (second stage, or
#'   `NULL`).}
#' \item{stage1, stage2}{the [survival::coxph()] fits (`stage2` is `NULL`
#'   without `x2`).}
#' \item{linear_predictor}{the full linear predictor of every included row,
#'   including the offset.}
#' \item{expected, residuals}{expected events and martingale residuals per
#'   row.}
#' \item{baseline}{a data frame with the distinct event times and the
#'   national cumulative baseline hazard.}
#' \item{data, vars, sizes, rows, control, converged}{as in
#'   [fit_logistic_fe()].}
#'
#' @references
#' He K, Schaubel DE (2014). Standardized mortality ratio for evaluating
#' center-specific mortality: assessment and alternative. *Statistics in
#' Biosciences*, 6(2), 1-19.
#'
#' Liu D, Schaubel DE, Kalbfleisch JD (2012). Computationally efficient
#' marginal models for clustered recurrent event data. *Biometrics*, 68(2),
#' 637-647.
#' @seealso [standardize.cox()], [test.cox()],
#'   [summary.cox()]
#' @examples
#' s <- sim_kecc(n_providers = 20, type = "survival", seed = 1)
#' ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", time = "stop", start = "start",
#'               event = "event", cutoff = 0, message = FALSE)
#' fit <- fit_cox(ps)
#' fit
#' coef(fit)$beta
#' @importFrom survival coxph Surv strata
#' @export
fit_cox <- function(data, strata = c("provider", "none"), x2 = NULL,
                    ties = c("breslow", "efron"), message = FALSE) {
  if (!inherits(data, "pp_data")) stop("`data` must come from pp_data().", call. = FALSE)
  if (outcome_of(data) != "survival") {
    stop("fit_cox() needs a survival outcome: give `time` and `event` to pp_data().",
         call. = FALSE)
  }
  strata <- match.arg(strata)
  ties <- match.arg(ties)
  v <- data$vars
  if (!is.null(x2)) {
    if (!is.character(x2) || !all(x2 %in% v$x)) {
      stop("`x2` must name covariates of pp_data()'s `x`.", call. = FALSE)
    }
    if (strata == "none") stop("`x2` is only used with strata = \"provider\".", call. = FALSE)
  }
  x1 <- setdiff(v$x, x2)
  d <- data$data[data$data$included, , drop = FALSE]
  if (nrow(d) == 0) stop("No provider reaches the cutoff.", call. = FALSE)
  rownames(d) <- NULL
  n <- nrow(d)
  sizes <- data$sizes[data$sizes$included, , drop = FALSE]
  rownames(sizes) <- NULL
  levels <- sizes$provider
  code <- provider_codes(d[[v$provider]], levels)
  start <- if (is.null(v$start)) rep(0, n) else as.numeric(d[[v$start]])
  stop_t <- as.numeric(d[[v$time]])
  event <- as.numeric(d[[v$event]])
  w <- if (is.null(v$weights)) rep(1, n) else as.numeric(d[[v$weights]])
  offset <- if (is.null(v$offset)) rep(0, n) else as.numeric(d[[v$offset]])

  dd <- d[, v$x, drop = FALSE]
  dd$.start <- start
  dd$.time <- stop_t
  dd$.event <- event
  dd$.provider <- factor(code, levels = seq_along(levels))
  dd$.offset <- offset
  if (!is.null(v$weights)) dd$.weights <- w
  if (!is.null(v$patient)) dd$.cluster <- d[[v$patient]]
  response <- "Surv(.start, .time, .event)"

  rhs1 <- c(if (strata == "provider") "strata(.provider)", backtick(x1), "offset(.offset)")
  s1 <- cox_stage(rhs1, response, dd, ties)
  beta <- stats::setNames(as.numeric(stats::coef(s1$fit)), x1)
  z1 <- as.matrix(dd[, x1, drop = FALSE])
  storage.mode(z1) <- "double"
  eta <- as.numeric(z1 %*% beta) + offset
  s2 <- NULL
  beta2 <- NULL
  if (!is.null(x2)) {
    dd$.eta1 <- eta
    s2 <- cox_stage(c(backtick(x2), "offset(.eta1)"), response, dd, ties)
    beta2 <- stats::setNames(as.numeric(stats::coef(s2$fit)), x2)
    z2 <- as.matrix(dd[, x2, drop = FALSE])
    storage.mode(z2) <- "double"
    eta <- eta + as.numeric(z2 %*% beta2)
  }
  ex <- cox_expectations(eta, start, stop_t, event, w, code, length(levels), direct = FALSE)

  out <- structure(list(coefficients = list(beta = beta, beta2 = beta2),
                        stage1 = s1$fit, stage2 = if (is.null(s2)) NULL else s2$fit,
                        linear_predictor = eta, expected = ex$rows,
                        residuals = event - ex$rows, baseline = ex$baseline,
                        data = d, vars = v, sizes = sizes,
                        rows = as.integer(tabulate(code, length(levels))),
                        control = list(strata = strata, x2 = x2, ties = ties,
                                       robust = !is.null(v$patient),
                                       iterations = c(s1$fit$iter,
                                                      if (!is.null(s2)) s2$fit$iter)),
                        converged = s1$converged && (is.null(s2) || s2$converged)),
                   class = "cox")
  if (message) print(out)
  out
}

backtick <- function(x) if (length(x)) paste0("`", x, "`") else NULL

# One coxph() stage: the formula is built from `rhs` and `response` and
# evaluated with the data `dd`, whose columns .weights and .cluster (if
# present) are the case weights and the cluster. A fit that runs out of
# iterations is recorded, other warnings pass through.
cox_stage <- function(rhs, response, dd, ties) {
  f <- stats::reformulate(rhs, response = response)
  environment(f) <- environment()
  args <- list(quote(survival::coxph), formula = f, data = quote(dd), ties = ties,
               model = FALSE, x = FALSE, y = FALSE)
  if (".weights" %in% names(dd)) args$weights <- as.name(".weights")
  if (".cluster" %in% names(dd)) args$cluster <- as.name(".cluster")
  converged <- TRUE
  fit <- withCallingHandlers(
    eval(as.call(args)),
    warning = function(w) {
      if (grepl("Ran out of iterations", conditionMessage(w))) {
        converged <<- FALSE
        invokeRestart("muffleWarning")
      }
    }
  )
  list(fit = fit, converged = converged)
}

#' @export
print.cox <- function(x, ...) {
  cat("Cox model", if (x$control$strata == "provider") " stratified by provider (two-stage)"
      else " with a pooled baseline (no provider effects)", "\n", sep = "")
  cat("  ", nrow(x$data), " rows, ", nrow(x$sizes), " providers, ",
      length(x$coefficients$beta), " covariate", if (length(x$coefficients$beta) != 1) "s",
      if (!is.null(x$coefficients$beta2)) {
        paste0(" + ", length(x$coefficients$beta2), " in the second stage")
      },
      "; ", format(sum(x$data[[x$vars$event]])), " events; ties ", x$control$ties,
      if (!x$converged) " (not converged)", "\n", sep = "")
  cat("  expected events sum to ", format(sum(x$expected), digits = 6), "\n", sep = "")
  invisible(x)
}

#' @export
coef.cox <- function(object, ...) {
  object$coefficients
}

#' Covariance of the covariate effects of a Cox fit
#'
#' @param object a model fitted by [fit_cox()].
#' @param variance `"naive"` (default; the inverse information) or
#'   `"robust"` (the patient-clustered sandwich, available when [pp_data()]
#'   has a patient column).
#' @param stage `1` (default) or `2`, the stage whose coefficients are
#'   wanted.
#' @param ... unused.
#' @return A covariance matrix.
#' @export
vcov.cox <- function(object, variance = c("naive", "robust"), stage = 1, ...) {
  variance <- match.arg(variance)
  fit <- cox_stage_fit(object, stage)
  if (variance == "robust") {
    if (!object$control$robust) {
      stop("The robust covariance needs a patient column in pp_data().", call. = FALSE)
    }
    v <- fit$var
  } else {
    v <- if (is.null(fit$naive.var)) fit$var else fit$naive.var
  }
  nm <- if (stage == 1) names(object$coefficients$beta) else names(object$coefficients$beta2)
  dimnames(v) <- list(nm, nm)
  v
}

cox_stage_fit <- function(object, stage) {
  if (!inherits(object, "cox")) stop("`fit` must be a cox object.", call. = FALSE)
  if (!stage %in% c(1, 2)) stop("`stage` must be 1 or 2.", call. = FALSE)
  if (stage == 2 && is.null(object$stage2)) {
    stop("The fit has no second stage (`x2` was NULL).", call. = FALSE)
  }
  if (stage == 1) object$stage1 else object$stage2
}
