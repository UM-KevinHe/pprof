#' Confidence intervals for the provider effects of a logistic fixed-effect fit
#'
#' Intervals for \eqn{\gamma_i} obtained by inverting a test: the exact
#' Poisson-binomial test, the score test with \eqn{\hat\beta} held fixed
#' (naive, robust or stabilized variance), or a Wald test with the naive or
#' patient-clustered standard error. Intervals for standardized measures are
#' produced by [standardize()], which maps these limits through the
#' standardization.
#'
#' @param object a model fitted by [fit_logistic_fe()].
#' @param parm optional subset of provider identifiers.
#' @param level confidence level. Default 0.95.
#' @param type `"exact"` (default), `"score"`, `"robust_score"`,
#'   `"stabilized_score"`, `"wald"` or `"robust_wald"`.
#' @param wald_variance for `type = "wald"`: `"case_mix"` (default) or
#'   `"full"`, see [gamma_variance()].
#' @param alternative `"two.sided"` (default), `"greater"` (lower limit
#'   only) or `"less"` (upper limit only).
#' @param ... unused.
#'
#' @details Providers with no events have no finite lower limit and providers
#' with all events no finite upper limit; for `"exact"` and `"score"` their
#' finite limit is one-sided at level `1 - level` (the convention of pprof
#' 1.0), while `"robust_score"` and `"stabilized_score"` (and the limits that
#' [test()] attaches to its score statistics) use the two-sided threshold
#' throughout, so that a limit excludes the null value exactly when the
#' provider is flagged. The exact inversion works on the individual Bernoulli
#' trials and is slow for providers with very many trials.
#'
#' @return A data frame with `provider`, `estimate`, `lower`, `upper`.
#' @examples
#' d <- sim_kecc(n_providers = 10, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_logistic_fe(pd)
#' confint(fit, type = "score")
#' @export
confint.logistic_fe <- function(object, parm = NULL, level = 0.95,
                                type = c("exact", "score", "robust_score", "stabilized_score",
                                         "wald", "robust_wald"),
                                wald_variance = c("case_mix", "full"),
                                alternative = c("two.sided", "greater", "less"), ...) {
  check_fit(object)
  type <- match.arg(type)
  wald_variance <- match.arg(wald_variance)
  alternative <- match.arg(alternative)
  m <- fit_matrices(object)
  idx <- provider_positions(parm, m$levels)
  limits <- if (type %in% c("robust_score", "stabilized_score")) {
    r <- test(object, statistic = "score", variance = sub("_score$", "", type),
              providers = m$levels[idx], level = level, alternative = alternative)
    list(lower = r$lower, upper = r$upper)
  } else {
    gamma_limits(object, m, type, idx, level, alternative, wald_variance)
  }
  data.frame(provider = m$levels[idx], estimate = unname(object$coefficients$gamma[idx]),
             lower = limits$lower, upper = limits$upper, stringsAsFactors = FALSE)
}

# Confidence limits of the provider effects at positions `idx`.
gamma_limits <- function(fit, m, type, idx, level, alternative, wald_variance = "case_mix") {
  gamma <- fit$coefficients$gamma
  alpha <- 1 - level
  if (type %in% c("wald", "robust_wald")) {
    se <- if (type == "wald") {
      sqrt(gamma_variance(fit, wald_variance))[idx]
    } else {
      test(fit, statistic = "wald", variance = "robust", providers = m$levels[idx])$se
    }
    q <- normal_quantile(level, alternative)
    lower <- if (alternative == "less") rep(-Inf, length(idx)) else gamma[idx] - q * se
    upper <- if (alternative == "greater") rep(Inf, length(idx)) else gamma[idx] + q * se
    return(list(lower = unname(lower), upper = unname(upper)))
  }

  code <- provider_codes(m$prov, m$levels)
  rows_of <- split(seq_along(code), factor(code, levels = seq_along(m$levels)))
  no_events <- fit$sizes$no_events
  all_events <- fit$sizes$all_events
  lp <- fit$linear_predictor
  q2 <- stats::qnorm(alpha / 2, lower.tail = FALSE)
  q1 <- stats::qnorm(alpha, lower.tail = FALSE)

  limits_one <- function(i) {
    r <- rows_of[[i]]
    zb <- lp[r]
    nt <- m$n[r]
    obs <- sum(m$y[r])
    g_hat <- gamma[i]
    if (type == "exact") {
      trials <- function(g) {
        p <- stats::plogis(g + zb)
        if (all(nt == 1)) p else rep(p, nt)
      }
      upper_fun <- function(g) {
        p <- trials(g)
        if (alternative == "two.sided") {
          poibin::ppoibin(obs - 1, p) + 0.5 * poibin::dpoibin(obs, p) - alpha / 2
        } else {
          poibin::ppoibin(obs, p) - alpha
        }
      }
      lower_fun <- function(g) {
        p <- trials(g)
        if (alternative == "two.sided") {
          1 - poibin::ppoibin(obs, p) + 0.5 * poibin::dpoibin(obs, p) - alpha / 2
        } else {
          1 - poibin::ppoibin(obs - 1, p) - alpha
        }
      }
      none_fun <- function(g) exp(sum(nt * stats::plogis(-(g + zb), log.p = TRUE))) / 2 - alpha
      all_fun <- function(g) exp(sum(nt * stats::plogis(g + zb, log.p = TRUE))) / 2 - alpha
    } else {
      z_of <- function(g) {
        p <- stats::plogis(g + zb)
        (obs - sum(nt * p)) / sqrt(sum(nt * p * (1 - p)))
      }
      q <- if (alternative == "two.sided") q2 else q1
      upper_fun <- function(g) z_of(g) + q
      lower_fun <- function(g) z_of(g) - q
      none_fun <- function(g) {
        p <- stats::plogis(g + zb)
        q1 - sum(nt * p) / sqrt(sum(nt * p * (1 - p)))
      }
      all_fun <- function(g) {
        p <- stats::plogis(g + zb)
        sum(nt * (1 - p)) / sqrt(max(sum(nt * p * (1 - p)), 1e-20)) - q1
      }
    }
    span <- 10 + max(abs(zb))
    if (no_events[i]) return(c(-Inf, find_root(none_fun, -span, span)))
    if (all_events[i]) return(c(find_root(all_fun, -span, span), Inf))
    lo <- if (alternative == "less") -Inf else find_root(lower_fun, g_hat - span, g_hat)
    up <- if (alternative == "greater") Inf else find_root(upper_fun, g_hat, g_hat + span)
    c(lo, up)
  }
  res <- vapply(idx, limits_one, numeric(2))
  list(lower = res[1, ], upper = res[2, ])
}

# Root of a monotone function on [a, b], widening the interval if needed.
find_root <- function(f, a, b, max_widen = 5) {
  width <- b - a
  for (k in 0:max_widen) {
    fa <- f(a)
    fb <- f(b)
    if (is.finite(fa) && is.finite(fb) && sign(fa) != sign(fb)) {
      return(stats::uniroot(f, c(a, b), tol = 1e-8)$root)
    }
    a <- a - width
    b <- b + width
  }
  NA_real_
}
