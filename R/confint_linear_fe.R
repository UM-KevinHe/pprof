#' Confidence intervals for the provider effects of a linear fixed-effect fit
#'
#' Wald intervals for \eqn{\gamma_i} with the model-based standard error
#' (t reference with \eqn{n - m - p} degrees of freedom) or the
#' patient-clustered standard error (normal reference).
#'
#' @param object a model fitted by [fit_linear_fe()].
#' @param parm optional subset of provider identifiers.
#' @param level confidence level. Default 0.95.
#' @param type `"wald"` (default) or `"robust_wald"`.
#' @param wald_variance for `type = "wald"`: `"case_mix"` (default) or
#'   `"full"`, see [gamma_variance()].
#' @param alternative `"two.sided"` (default), `"greater"` (lower limit
#'   only) or `"less"` (upper limit only).
#' @param ... unused.
#'
#' @return A data frame with `provider`, `estimate`, `lower`, `upper`.
#' @examples
#' d <- sim_kecc(n_providers = 10, type = "continuous", seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' confint(fit_linear_fe(pd))
#' @export
confint.linear_fe <- function(object, parm = NULL, level = 0.95, type = c("wald", "robust_wald"),
                              wald_variance = c("case_mix", "full"),
                              alternative = c("two.sided", "greater", "less"), ...) {
  check_fit(object, "linear_fe")
  type <- match.arg(type)
  wald_variance <- match.arg(wald_variance)
  alternative <- match.arg(alternative)
  res <- test_providers(object, variance = if (type == "wald") "naive" else "robust",
                        providers = parm, level = level, alternative = alternative,
                        wald_variance = wald_variance)
  data.frame(provider = res$provider, estimate = res$estimate, lower = res$lower,
             upper = res$upper, stringsAsFactors = FALSE)
}
