#' Tests of the covariate coefficients of a Cox fit
#'
#' Wald tests of the log hazard ratios of a stage of [fit_cox()].
#'
#' @param object a model fitted by [fit_cox()].
#' @param variance `"naive"` (default; the inverse information) or
#'   `"robust"` (the patient-clustered sandwich; needs a patient column).
#' @param covariates optional subset of covariate names.
#' @param level confidence level for the Wald intervals.
#' @param null the null value of each coefficient. Default 0.
#' @param stage `1` (default) or `2`.
#' @param ... unused.
#'
#' @return A data frame with one row per covariate: `covariate`, `estimate`
#'   (the log hazard ratio), `se`, `statistic` (z), `p`, `lower` and `upper`.
#' @examples
#' s <- sim_kecc(n_providers = 20, type = "survival", seed = 1)
#' ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", time = "stop", start = "start",
#'               event = "event", cutoff = 0, message = FALSE)
#' summary(fit_cox(ps), variance = "robust")
#' @export
summary.cox <- function(object, variance = c("naive", "robust"), covariates = NULL,
                                level = 0.95, null = 0, stage = 1, ...) {
  check_fit(object, "cox")
  variance <- match.arg(variance)
  check_level(level)
  beta <- if (stage == 1) object$coefficients$beta else object$coefficients$beta2
  v <- vcov(object, variance = variance, stage = stage)
  x <- names(beta)
  ind <- if (is.null(covariates)) seq_along(x) else match(covariates, x)
  if (anyNA(ind)) {
    stop("Unknown covariate(s): ", paste(covariates[is.na(ind)], collapse = ", "), call. = FALSE)
  }
  se <- sqrt(diag(v))[ind]
  z <- (beta[ind] - null) / se
  q <- stats::qnorm(1 - (1 - level) / 2)
  pp_summary(data.frame(covariate = x[ind], estimate = unname(beta[ind]), se = unname(se),
             statistic = unname(z), p = unname(2 * stats::pnorm(-abs(z))),
             lower = unname(beta[ind] - q * se), upper = unname(beta[ind] + q * se),
             stringsAsFactors = FALSE))
}
