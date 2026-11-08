#' Tests of the provider effects of a linear fixed-effect fit
#'
#' Wald tests of \eqn{H_0: \gamma_i = \gamma_0} for every provider (or a
#' subset) with the model-based or the patient-clustered standard error.
#'
#' @param fit a model fitted by [fit_linear_fe()].
#' @param variance `"naive"` (default; model-based, see [gamma_variance()])
#'   or `"robust"` (patient-clustered sandwich with \eqn{\hat\beta} held
#'   fixed; needs a patient column in [pp_data()]).
#' @param null the null value \eqn{\gamma_0}: `"median"` (default) or
#'   `"mean"` (weighted by rows) of the fitted effects, or a number.
#' @param providers optional subset of provider identifiers to report.
#' @param level confidence level; the significance level is `1 - level`.
#' @param alternative `"two.sided"` (default), `"greater"` or `"less"`.
#' @param wald_variance for `variance = "naive"`: `"case_mix"` (default) or
#'   `"full"`, see [gamma_variance()].
#' @param ... unused.
#'
#' @details The naive statistic \eqn{(\hat\gamma_i - \gamma_0) / se_i} is
#' referred to the t distribution with \eqn{n - m - p} degrees of freedom,
#' the robust statistic to the standard normal. The robust standard error is
#' \eqn{\sqrt{\sum_{k \in i} (\sum_{j \in k} e_{ij})^2} / n_i} over the
#' provider's patients \eqn{k}, with \eqn{e_{ij}} the residuals: the
#' sandwich variance of \eqn{\bar Y_i - \bar Z_i^\top\hat\beta}.
#'
#' @return A data frame of class `pp_test` with the columns described in
#'   [test_providers.logistic_fe()].
#' @examples
#' d <- sim_kecc(n_providers = 20, type = "continuous", seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_linear_fe(pd)
#' test_providers(fit, variance = "robust")
#' @export
test_providers.linear_fe <- function(fit, variance = c("naive", "robust"), null = "median",
                                     providers = NULL, level = 0.95,
                                     alternative = c("two.sided", "greater", "less"),
                                     wald_variance = c("case_mix", "full"), ...) {
  check_fit(fit, "linear_fe")
  variance <- match.arg(variance)
  alternative <- match.arg(alternative)
  wald_variance <- match.arg(wald_variance)
  check_level(level)
  if (variance == "robust" && is.null(fit$vars$patient)) {
    stop("The robust variance needs a patient column in pp_data().", call. = FALSE)
  }
  m <- fit_matrices(fit)
  gamma <- fit$coefficients$gamma
  g0 <- null_value(null, gamma, sizes = m$n_prov)
  idx <- provider_positions(providers, m$levels)

  if (variance == "naive") {
    se_all <- sqrt(gamma_variance(fit, wald_variance))
    df <- fit$df_residual
  } else {
    meat <- robust_gamma_meat(fit$residuals, m$prov, m$patient, m$levels)
    se_all <- sqrt(meat) / m$n_prov
    df <- Inf
  }
  se <- unname(se_all[idx])
  stat <- (gamma[idx] - g0) / se
  stat[!is.finite(stat)] <- NA_real_
  decision <- z_decision(stat, level, alternative, df)
  limits <- wald_limits(gamma[idx], se, level, alternative, df)
  pp_test_frame(m$levels[idx], gamma[idx], se, g0, stat, decision, limits$lower, limits$upper,
                test = list(statistic = "wald", variance = variance,
                            method = if (variance == "naive") wald_variance else NA,
                            null = null, null_value = g0, level = level,
                            alternative = alternative, df = df),
                sizes = data.frame(provider = m$levels[idx], rows = m$n_prov[idx],
                                   patients = fit$sizes$patients[idx], stringsAsFactors = FALSE))
}
