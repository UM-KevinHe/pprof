#' Standardized ratios and rates of a logistic fixed-effect fit
#'
#' Indirect or direct standardization of every provider's outcome, with the
#' observed and expected counts, the delta-method numerator of the direct
#' standard error, and optional confidence intervals mapped from the limits
#' of the provider effects.
#'
#' @param fit a model fitted by [fit_logistic_fe()].
#' @param type `"indirect"` (default) or `"direct"`.
#' @param null for indirect standardization, the norm of the provider effect:
#'   `"median"` (default), `"mean"` (weighted by rows) or a number.
#' @param n_pop denominator of the rates: the number of trials in the
#'   standard population. Defaults to the number of trials among the fitted
#'   rows.
#' @param providers optional subset of provider identifiers.
#' @param gamma_se optional standard errors of the provider effects, named by
#'   provider (for instance the `se` column of a Wald [test_providers()]
#'   result); with them the direct rate's standard error
#'   \eqn{\mathrm{numer\_se}/n_{pop} \times se(\hat\gamma_i)} is returned.
#' @param level if given, confidence limits of the ratio and the rate are
#'   added, mapped from the limits of \eqn{\gamma_i} obtained with `ci`.
#' @param ci the test inverted for the limits, see [confint.logistic_fe()].
#' @param alternative `"two.sided"` (default), `"greater"` or `"less"`.
#' @param threads OpenMP threads for direct standardization.
#' @param ... unused.
#'
#' @details With \eqn{p_{ij}(\gamma) = \mathrm{expit}(\gamma + Z_{ij}^\top
#' \hat\beta)} and \eqn{N_{ij}} trials per row:
#' * indirect: \eqn{E_i = \sum_{j \in i} N_{ij}\, p_{ij}(\gamma_0)}, ratio
#'   \eqn{O_i / E_i};
#' * direct: \eqn{E_i = \sum_k N_k\, p_k(\hat\gamma_i)} over every fitted row,
#'   ratio \eqn{E_i / \sum O}.
#'
#' In both cases rate = ratio \eqn{\times \sum O / n_{pop}}, truncated to
#' \[0, 1\]. The direct `numer_se` is \eqn{\sum_k N_k p_k(1 - p_k)} at
#' \eqn{\hat\gamma_i}; the direct rate's standard error with \eqn{\hat\beta}
#' held fixed is `numer_se / n_pop * se(gamma_i)`. For an arbitrary standard
#' population use [direct_expectation()].
#'
#' @return A data frame with one row per provider: `provider`, `observed`,
#'   `expected`, `ratio`, `rate`, `variance_expected` (indirect:
#'   \eqn{\sum N p_0(1 - p_0)}; direct: `numer_se`), and when requested
#'   `se_rate`, `lower_ratio`, `upper_ratio`, `lower_rate`, `upper_rate`.
#'   Attributes `type`, `null_value`, `n_pop` and `population_rate`.
#' @references He K (2019). Indirect and direct standardization for
#'   evaluating transplant centers. *Journal of Hospital Administration*,
#'   8(1), 9-14.
#' @examples
#' d <- sim_kecc(n_providers = 20, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_logistic_fe(pd)
#' head(standardize(fit, type = "direct"))
#' @export
standardize.logistic_fe <- function(fit, type = c("indirect", "direct"), null = "median",
                                    n_pop = NULL, providers = NULL, gamma_se = NULL, level = NULL,
                                    ci = c("exact", "score", "wald", "robust_wald"),
                                    alternative = c("two.sided", "greater", "less"), threads = 1,
                                    ...) {
  check_fit(fit)
  type <- match.arg(type)
  ci <- match.arg(ci)
  alternative <- match.arg(alternative)
  m <- fit_matrices(fit)
  gamma <- fit$coefficients$gamma
  idx <- provider_positions(providers, m$levels)
  code <- provider_codes(m$prov, m$levels)
  lp <- fit$linear_predictor
  total_obs <- sum(m$y)
  total_trials <- sum(m$n)
  if (is.null(n_pop)) n_pop <- total_trials
  population_rate <- total_obs / n_pop
  observed <- as.numeric(rowsum(m$y, code, reorder = TRUE))
  clip <- function(x) pmax(pmin(x, 1), 0)

  if (type == "indirect") {
    g0 <- null_value(null, gamma, sizes = m$n_prov)
    p0 <- stats::plogis(g0 + lp)
    expected <- as.numeric(rowsum(m$n * p0, code, reorder = TRUE))
    variance_expected <- as.numeric(rowsum(m$n * p0 * (1 - p0), code, reorder = TRUE))
    ratio <- observed / expected
  } else {
    g0 <- NA_real_
    res <- computeDirectExp(unname(gamma), lp, m$n, max(1L, as.integer(threads)))
    expected <- res[, 1]
    variance_expected <- res[, 2]
    ratio <- expected / total_obs
  }
  out <- data.frame(provider = m$levels, observed = observed, expected = expected, ratio = ratio,
                    rate = clip(ratio * population_rate), variance_expected = variance_expected,
                    stringsAsFactors = FALSE)

  if (!is.null(gamma_se)) {
    if (type == "indirect") {
      stop("`gamma_se` gives the standard error of the direct rate; use type = \"direct\".",
           call. = FALSE)
    }
    if (!is.null(names(gamma_se))) {
      gamma_se <- gamma_se[m$levels]
    } else if (length(gamma_se) != length(m$levels)) {
      stop("`gamma_se` must be named by provider or have one value per fitted provider.",
           call. = FALSE)
    }
    out$se_rate <- variance_expected / n_pop * as.numeric(gamma_se)
  }

  if (!is.null(level)) {
    limits <- gamma_limits(fit, m, ci, seq_along(m$levels), level, alternative)
    trials_of_provider <- as.numeric(rowsum(m$n, code, reorder = TRUE))
    if (type == "indirect") {
      rows_of <- split(seq_along(code), factor(code, levels = seq_along(m$levels)))
      expected_at <- function(g) {
        mapply(function(r, gi) {
          if (!is.finite(gi)) return(if (gi > 0) sum(m$n[r]) else 0)
          sum(m$n[r] * stats::plogis(gi + lp[r]))
        }, rows_of, g)
      }
      lower <- expected_at(limits$lower) / expected
      upper <- expected_at(limits$upper) / expected
      lower[fit$sizes$no_events] <- 0
      upper[fit$sizes$all_events] <- (trials_of_provider / expected)[fit$sizes$all_events]
      if (alternative == "less") lower <- rep(0, length(lower))
      if (alternative == "greater") upper <- trials_of_provider / expected
    } else {
      direct_at <- function(g) {
        val <- rep(NA_real_, length(g))
        fin <- is.finite(g)
        if (any(fin)) val[fin] <- computeDirectExp(g[fin], lp, m$n, 1L)[, 1]
        val[!fin & g > 0] <- total_trials
        val[!fin & g < 0] <- 0
        val
      }
      lower <- direct_at(limits$lower) / total_obs
      upper <- direct_at(limits$upper) / total_obs
      lower[fit$sizes$no_events] <- 0
      upper[fit$sizes$all_events] <- total_trials / total_obs
      if (alternative == "less") lower <- rep(0, length(lower))
      if (alternative == "greater") upper <- rep(total_trials / total_obs, length(upper))
    }
    out$lower_ratio <- lower
    out$upper_ratio <- upper
    out$lower_rate <- clip(lower * population_rate)
    out$upper_rate <- clip(upper * population_rate)
  }

  out <- out[idx, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "type") <- type
  attr(out, "null_value") <- g0
  attr(out, "n_pop") <- n_pop
  attr(out, "population_rate") <- population_rate
  attr(out, "level") <- level
  out
}

#' Direct standardization over an arbitrary standard population
#'
#' For each provider effect \eqn{\gamma_i}, the expected number of events
#' \eqn{E_i = \sum_k N_k\, \mathrm{expit}(\gamma_i + \eta_k)} over the rows of
#' a standard population with linear predictors \eqn{\eta_k} and trials
#' \eqn{N_k}, and \eqn{V_i = \sum_k N_k\, p_k (1 - p_k)}. The direct
#' standardized rate is \eqn{E_i / \sum_k N_k}; its delta-method standard
#' error with \eqn{\hat\beta} held fixed is \eqn{V_i / \sum_k N_k \times
#' se(\hat\gamma_i)}.
#'
#' @param gamma provider effects.
#' @param linear_predictor \eqn{Z\hat\beta} of the standard population's rows.
#' @param n number of trials per row of the standard population (default 1).
#' @param threads number of OpenMP threads.
#' @return A data frame with `expected` and `variance`, one row per element of
#'   `gamma`.
#' @examples
#' direct_expectation(c(-1, 0, 1), linear_predictor = rnorm(100))
#' @export
direct_expectation <- function(gamma, linear_predictor, n = 1, threads = 1) {
  linear_predictor <- as.numeric(linear_predictor)
  if (length(n) == 1) n <- rep(n, length(linear_predictor))
  if (length(n) != length(linear_predictor)) {
    stop("`n` must have one value per row of the standard population.", call. = FALSE)
  }
  res <- computeDirectExp(as.numeric(gamma), linear_predictor, as.numeric(n),
                          max(1L, as.integer(threads)))
  data.frame(expected = res[, 1], variance = res[, 2])
}
