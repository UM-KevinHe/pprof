#' Tests of the fixed effects of a random-effect fit
#'
#' Wald tests of the fixed-effect coefficients (including the intercept)
#' with the covariance matrix of the lme4 object, referred to the standard
#' normal.
#'
#' @param object a model from [fit_logistic_re()], [fit_logistic_cre()],
#'   [fit_linear_re()] or [fit_linear_cre()].
#' @param covariates optional subset of term names (as in
#'   `names(coef(object)$beta)`).
#' @param level confidence level for the intervals.
#' @param null the null value of each coefficient.
#' @param ... unused.
#'
#' @return A data frame with one row per term: `covariate`, `estimate`,
#'   `se`, `statistic`, `p`, `lower`, `upper`.
#' @examples
#' \donttest{
#' d <- sim_kecc(n_providers = 20, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' summary(fit_logistic_re(pd))
#' }
#' @name summary_mixed
NULL

#' @rdname summary_mixed
#' @export
summary.logistic_re <- function(object, covariates = NULL, level = 0.95, null = 0, ...) {
  summary_mixed(object, covariates, level, null)
}

#' @rdname summary_mixed
#' @export
summary.linear_re <- function(object, covariates = NULL, level = 0.95, null = 0, ...) {
  summary_mixed(object, covariates, level, null)
}

summary_mixed <- function(object, covariates, level, null) {
  check_mixed(object)
  check_level(level)
  beta <- object$coefficients$beta
  terms <- names(beta)
  ind <- if (is.null(covariates)) seq_along(terms) else match(covariates, terms)
  if (anyNA(ind)) {
    stop("Unknown term(s): ", paste(covariates[is.na(ind)], collapse = ", "), call. = FALSE)
  }
  se <- sqrt(diag(object$variance$beta))[ind]
  stat <- (beta[ind] - null) / se
  decision <- z_decision(stat, level, "two.sided")
  limits <- wald_limits(beta[ind], se, level, "two.sided")
  pp_summary(data.frame(covariate = terms[ind], estimate = unname(beta[ind]), se = unname(se),
             statistic = unname(stat), p = decision$p, lower = limits$lower,
             upper = limits$upper, stringsAsFactors = FALSE))
}
