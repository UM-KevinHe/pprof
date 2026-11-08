#' Tests of the fixed effects of a random-effect fit
#'
#' Wald tests of the fixed-effect coefficients (including the intercept)
#' with the covariance matrix of the lme4 fit, referred to the standard
#' normal.
#'
#' @param fit a model from [fit_logistic_re()], [fit_logistic_cre()],
#'   [fit_linear_re()] or [fit_linear_cre()].
#' @param covariates optional subset of term names (as in
#'   `names(coef(fit)$beta)`).
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
#' test_covariates(fit_logistic_re(pd))
#' }
#' @name test_covariates_mixed
NULL

#' @rdname test_covariates_mixed
#' @export
test_covariates.logistic_re <- function(fit, covariates = NULL, level = 0.95, null = 0, ...) {
  test_covariates_mixed(fit, covariates, level, null)
}

#' @rdname test_covariates_mixed
#' @export
test_covariates.linear_re <- function(fit, covariates = NULL, level = 0.95, null = 0, ...) {
  test_covariates_mixed(fit, covariates, level, null)
}

test_covariates_mixed <- function(fit, covariates, level, null) {
  check_mixed(fit)
  check_level(level)
  beta <- fit$coefficients$beta
  terms <- names(beta)
  ind <- if (is.null(covariates)) seq_along(terms) else match(covariates, terms)
  if (anyNA(ind)) {
    stop("Unknown term(s): ", paste(covariates[is.na(ind)], collapse = ", "), call. = FALSE)
  }
  se <- sqrt(diag(fit$variance$beta))[ind]
  stat <- (beta[ind] - null) / se
  decision <- z_decision(stat, level, "two.sided")
  limits <- wald_limits(beta[ind], se, level, "two.sided")
  data.frame(covariate = terms[ind], estimate = unname(beta[ind]), se = unname(se),
             statistic = unname(stat), p = decision$p, lower = limits$lower,
             upper = limits$upper, stringsAsFactors = FALSE)
}
