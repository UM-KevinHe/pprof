#' Standardized differences of a linear fixed-effect fit
#'
#' Indirect or direct standardization of every provider's mean outcome, as
#' a difference from the expected mean.
#'
#' @param fit a model fitted by [fit_linear_fe()].
#' @param type `"indirect"` (default) or `"direct"`.
#' @param null for indirect standardization, the norm of the provider effect:
#'   `"median"` (default), `"mean"` (weighted by rows) or a number.
#' @param providers optional subset of provider identifiers.
#' @param gamma_se optional standard errors of the provider effects, named by
#'   provider (for instance the `se` column of a [test_providers()] result);
#'   with them `se_difference` is returned.
#' @param level if given, confidence limits of the difference are added,
#'   mapped from Wald limits of \eqn{\gamma_i}.
#' @param ci the standard error of the limits: `"wald"` (default; model-based,
#'   with `wald_variance`) or `"robust_wald"`, see [confint.linear_fe()].
#' @param wald_variance `"case_mix"` (default) or `"full"`, see
#'   [gamma_variance()].
#' @param alternative `"two.sided"` (default), `"greater"` or `"less"`.
#' @param ... unused.
#'
#' @details With \eqn{\mu_{ij}(\gamma) = \gamma + Z_{ij}^\top\hat\beta}:
#' * indirect: `observed` \eqn{O_i = \sum_{j \in i} Y_{ij}}, `expected`
#'   \eqn{E_i = \sum_{j \in i} \mu_{ij}(\gamma_0)}, difference
#'   \eqn{(O_i - E_i) / n_i = \hat\gamma_i - \gamma_0};
#' * direct: `expected` \eqn{E_i = \sum_k \mu_k(\hat\gamma_i)} over every
#'   fitted row, `observed` \eqn{O = \sum_k Y_k}, difference
#'   \eqn{(E_i - O) / n}, the deviation of \eqn{\hat\gamma_i} from the
#'   row-weighted mean of the effects.
#'
#' Both differences are linear in \eqn{\hat\gamma_i}, so their standard
#' error with \eqn{\hat\beta} held fixed is \eqn{se(\hat\gamma_i)} and the
#' limits are those of \eqn{\gamma_i} shifted.
#'
#' @return A data frame with one row per provider: `provider`, `observed`,
#'   `expected`, `difference`, and when requested `se_difference`, `lower`,
#'   `upper`. Attributes `type`, `null_value`, `level`.
#' @examples
#' d <- sim_kecc(n_providers = 20, type = "continuous", seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_linear_fe(pd)
#' head(standardize(fit, level = 0.95))
#' @export
standardize.linear_fe <- function(fit, type = c("indirect", "direct"), null = "median",
                                  providers = NULL, gamma_se = NULL, level = NULL,
                                  ci = c("wald", "robust_wald"),
                                  wald_variance = c("case_mix", "full"),
                                  alternative = c("two.sided", "greater", "less"), ...) {
  check_fit(fit, "linear_fe")
  type <- match.arg(type)
  ci <- match.arg(ci)
  wald_variance <- match.arg(wald_variance)
  alternative <- match.arg(alternative)
  m <- fit_matrices(fit)
  gamma <- fit$coefficients$gamma
  idx <- provider_positions(providers, m$levels)
  code <- provider_codes(m$prov, m$levels)
  lp <- fit$linear_predictor

  if (type == "indirect") {
    g0 <- null_value(null, gamma, sizes = m$n_prov)
    observed <- as.numeric(rowsum(m$y, code, reorder = TRUE))
    expected <- as.numeric(rowsum(g0 + lp, code, reorder = TRUE))
    difference <- (observed - expected) / m$n_prov
  } else {
    g0 <- NA_real_
    observed <- rep(sum(m$y), length(gamma))
    expected <- unname(gamma) * length(m$y) + sum(lp)
    difference <- (expected - observed) / length(m$y)
  }
  out <- data.frame(provider = m$levels, observed = observed, expected = expected,
                    difference = difference, stringsAsFactors = FALSE)

  if (!is.null(gamma_se)) {
    if (!is.null(names(gamma_se))) {
      gamma_se <- gamma_se[m$levels]
    } else if (length(gamma_se) != length(m$levels)) {
      stop("`gamma_se` must be named by provider or have one value per fitted provider.",
           call. = FALSE)
    }
    out$se_difference <- as.numeric(gamma_se)
  }
  if (!is.null(level)) {
    limits <- stats::confint(fit, level = level, type = ci, wald_variance = wald_variance,
                             alternative = alternative)
    out$lower <- difference - (unname(gamma) - limits$lower)
    out$upper <- difference + (limits$upper - unname(gamma))
  }
  out <- out[idx, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "type") <- type
  attr(out, "null_value") <- g0
  attr(out, "level") <- level
  out
}
