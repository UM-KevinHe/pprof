#' Tests of the provider effects of a random-effect fit
#'
#' Tests \eqn{H_0: \alpha_i = \alpha_0} for every provider (or a subset) with
#' the statistic \eqn{(\hat\alpha_i - \alpha_0) / \sqrt{v_i}}, where
#' \eqn{\hat\alpha_i} is the predicted random effect and \eqn{v_i} its
#' conditional variance given the data (`post_var`), referred to the
#' standard normal.
#'
#' @param fit a model from [fit_logistic_re()], [fit_logistic_cre()],
#'   [fit_linear_re()] or [fit_linear_cre()].
#' @param null the null value \eqn{\alpha_0}: `"median"` (default) or `"mean"` (weighted by
#'   rows) of the predicted effects, or a number (0 is the population
#'   intercept).
#' @param providers optional subset of provider identifiers to report.
#' @param level confidence level; the significance level is `1 - level`.
#' @param alternative `"two.sided"` (default), `"greater"` or `"less"`.
#' @param null_model optional null model (a `pp_null` object or a function of
#'   the result returning one), applied through [calibrate()].
#' @param ... unused.
#'
#' @return A data frame of class `pp_test` with the columns described in
#'   [test.logistic_fe()]; `estimate` is \eqn{\hat\alpha_i}.
#' @examples
#' \donttest{
#' d <- sim_kecc(n_providers = 20, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' test(fit_logistic_re(pd))
#' }
#' @name test_mixed
NULL

#' @rdname test_mixed
#' @export
test.logistic_re <- function(fit, null = "median", providers = NULL, level = 0.95,
                                       alternative = c("two.sided", "greater", "less"),
                                       null_model = NULL, ...) {
  test_mixed(fit, null, providers, level, match.arg(alternative), null_model)
}

#' @rdname test_mixed
#' @export
test.linear_re <- function(fit, null = "median", providers = NULL, level = 0.95,
                                     alternative = c("two.sided", "greater", "less"),
                                       null_model = NULL, ...) {
  test_mixed(fit, null, providers, level, match.arg(alternative), null_model)
}

test_mixed <- function(fit, null, providers, level, alternative, null_model) {
  check_mixed(fit)
  check_level(level)
  alpha <- fit$coefficients$alpha
  levels <- fit$sizes$provider
  a0 <- null_value(null, alpha, sizes = fit$rows)
  idx <- provider_positions(providers, levels)
  se <- sqrt(unname(fit$post_var[idx]))
  stat <- (alpha[idx] - a0) / se
  stat[!is.finite(stat)] <- NA_real_
  decision <- z_decision(stat, level, alternative)
  limits <- wald_limits(alpha[idx], se, level, alternative)
  out <- pp_test_frame(levels[idx], alpha[idx], se, a0, stat, decision, limits$lower,
                       limits$upper,
                test = list(statistic = "wald", variance = "conditional", method = NA,
                            null = null, null_value = a0, level = level,
                            alternative = alternative),
                sizes = data.frame(provider = levels[idx], rows = fit$rows[idx],
                                   patients = fit$sizes$patients[idx], stringsAsFactors = FALSE))
  if (!is.null(null_model)) out <- calibrate(out, null_model)
  out
}
