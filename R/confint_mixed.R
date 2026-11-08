#' Confidence intervals for the provider effects of a random-effect fit
#'
#' Wald intervals \eqn{\hat\alpha_i \pm q \sqrt{v_i}} from the predicted
#' random effects and their conditional variances.
#'
#' @param object a model from [fit_logistic_re()], [fit_logistic_cre()],
#'   [fit_linear_re()] or [fit_linear_cre()].
#' @param parm optional subset of provider identifiers.
#' @param level confidence level. Default 0.95.
#' @param alternative `"two.sided"` (default), `"greater"` (lower limit
#'   only) or `"less"` (upper limit only).
#' @param ... unused.
#'
#' @return A data frame with `provider`, `estimate`, `lower`, `upper`.
#' @examples
#' \donttest{
#' d <- sim_kecc(n_providers = 10, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' confint(fit_logistic_re(pd))
#' }
#' @name confint_mixed
NULL

#' @rdname confint_mixed
#' @export
confint.logistic_re <- function(object, parm = NULL, level = 0.95,
                                alternative = c("two.sided", "greater", "less"), ...) {
  confint_mixed(object, parm, level, match.arg(alternative))
}

#' @rdname confint_mixed
#' @export
confint.linear_re <- function(object, parm = NULL, level = 0.95,
                              alternative = c("two.sided", "greater", "less"), ...) {
  confint_mixed(object, parm, level, match.arg(alternative))
}

confint_mixed <- function(object, parm, level, alternative) {
  res <- test(object, providers = parm, level = level, alternative = alternative)
  data.frame(provider = res$provider, estimate = res$estimate, lower = res$lower,
             upper = res$upper, stringsAsFactors = FALSE)
}
