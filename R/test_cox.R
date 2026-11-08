#' Tests of the standardized ratios of a Cox fit
#'
#' Tests every provider's indirect standardized ratio \eqn{O_j / E_j} (see
#' [standardize.cox()]) against 1, treating the observed count \eqn{O_j} as
#' Poisson with mean \eqn{E_j}, and flags the providers with more or fewer
#' events than expected. These are the tests of the survival measures (SMR,
#' SHR, STrR, SWR).
#'
#' @param fit a model fitted by [fit_cox()].
#' @param statistic `"midp"` (default; the two-sided mid-p Poisson test
#'   turned into a z statistic, the SMR convention), `"exact"` (the two-sided
#'   Poisson test with exact or Byar limits, the SWR convention; theoretical
#'   null only) or `"wald"` (\eqn{\log(O_j / E_j)} over its standard error
#'   \eqn{\sqrt{\phi / O_j}}, the SHR convention).
#' @param dispersion for `statistic = "wald"`: `"none"` (default; \eqn{\phi =
#'   1}) or `"deviance"` (\eqn{\phi} = the Poisson deviance of the
#'   patient-level counts about \eqn{(O_j / E_j) e_k} over its degrees of
#'   freedom, as `PROC GENMOD`'s `DSCALE`; the units are the patients when
#'   [pp_data()] has a patient column, the rows otherwise).
#' @param providers optional subset of provider identifiers to report.
#' @param level confidence level; the significance level is `1 - level`.
#' @param alternative `"two.sided"` (default), `"greater"` or `"less"`.
#' @param min_expected providers with fewer than `min_expected` expected
#'   events are dropped from the result (the measures use 3, so that an
#'   empirical null fitted on the result does not see them). Default 0.
#' @param null_model optional null model: a `pp_null` object or a function
#'   of the result returning one, for instance
#'   `function(t) null_model("empirical", t, size = t$person_time, grouping = "quantile")`
#'   (the survival measures group providers by quartiles of person-time);
#'   the result is then passed through [calibrate()]. Not available for the
#'   exact test.
#' @param ... unused.
#'
#' @details The mid-p statistic is
#' \eqn{z = \pm\Phi^{-1}(p)} with \eqn{p = \max(10^{-6}, \min(P_{low},
#' P_{high}) / 2)}, \eqn{P_{low} = 2 F(O_j) - f(O_j)}, \eqn{P_{high} = 2 [1 -
#' F(O_j - 1)] - f(O_j)}, \eqn{F} and \eqn{f} the Poisson distribution and
#' mass functions with mean \eqn{E_j}, the sign negative when the lower tail
#' is the smaller. Its confidence limits invert the (calibrated) test: the
#' Poisson means at which the provider's count is just significant, divided
#' by \eqn{E_j}, so that a limit excludes 1 exactly when the provider is
#' flagged.
#'
#' The exact test's p-value is \eqn{2 P(X \ge O_j)} for a ratio above 1 and
#' \eqn{2 P(X \le O_j)} otherwise, capped at 0.999; its limits are the exact
#' chi-square limits when \eqn{E_j < 100} and Byar's approximation otherwise.
#'
#' The Wald limits are \eqn{\exp[\log(O_j / E_j) \pm c\, se]}; a provider
#' with no events has no Wald statistic.
#'
#' @return A data frame of class `pp_test` with one row per tested provider:
#'   `provider`, `estimate` (the ratio \eqn{O_j / E_j}), `se` (Wald only;
#'   the standard error of the log ratio), `null_value` (1), `statistic`, `p`,
#'   `flag` (`1` more events than expected, `-1` fewer, `0` not different, `NA`
#'   not testable), `lower`, `upper` (limits for the ratio), `observed`,
#'   `expected`, `person_time`. Attributes `test` (the settings, including
#'   `phi`) and `sizes`.
#' @references
#' He K, Schaubel DE (2015). Standardized mortality ratio for evaluating
#' center-specific mortality: assessment and alternative. *Statistics in
#' Biosciences*, 7, 296-321.
#'
#' Breslow NE, Day NE (1987). *Statistical Methods in Cancer Research, Volume
#' II: The Design and Analysis of Cohort Studies*. IARC.
#' @examples
#' s <- sim_kecc(n_providers = 30, type = "survival", seed = 1)
#' ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", time = "stop", start = "start",
#'               event = "event", cutoff = 0, message = FALSE)
#' fit <- fit_cox(ps)
#' test(fit)
#' test(fit, statistic = "wald", dispersion = "deviance")
#' @export
test.cox <- function(fit, statistic = c("midp", "exact", "wald"),
                               dispersion = c("none", "deviance"), providers = NULL,
                               level = 0.95, alternative = c("two.sided", "greater", "less"),
                               min_expected = 0, null_model = NULL, ...) {
  check_fit(fit, "cox")
  statistic <- match.arg(statistic)
  dispersion <- match.arg(dispersion)
  alternative <- match.arg(alternative)
  check_level(level)
  if (statistic == "exact" && !is.null(null_model)) {
    stop("The exact Poisson test uses the theoretical null; calibrate the mid-p statistic ",
         "(statistic = \"midp\") instead.", call. = FALSE)
  }
  if (!is.numeric(min_expected) || length(min_expected) != 1 || min_expected < 0) {
    stop("`min_expected` must be a single non-negative number.", call. = FALSE)
  }
  sm <- standardize(fit, type = "indirect", providers = providers)
  keep <- sm$expected >= min_expected
  n_dropped <- sum(!keep)
  sm <- sm[keep, , drop = FALSE]
  obs <- sm$observed
  expected <- sm$expected
  k <- nrow(sm)
  q <- normal_quantile(level, alternative)
  se <- rep(NA_real_, k)
  phi <- NA_real_

  if (statistic == "midp") {
    z <- midp_z(obs, expected)
    decision <- z_decision(z, level, alternative)
    lim <- midp_limits(obs, expected, rep(0, k), rep(1, k), q, alternative)
  } else if (statistic == "exact") {
    ex <- exact_poisson(obs, expected, level, alternative)
    z <- ex$z
    decision <- list(flag = ex$flag, p = ex$p)
    lim <- ex[c("lower", "upper")]
  } else {
    phi <- if (dispersion == "deviance") deviance_dispersion(fit, sm$provider) else 1
    log_ratio <- ifelse(obs > 0, log(obs / expected), NA_real_)
    se <- ifelse(obs > 0, sqrt(phi / obs), NA_real_)
    z <- log_ratio / se
    decision <- z_decision(z, level, alternative)
    wl <- wald_limits(log_ratio, se, level, alternative)
    lim <- list(lower = exp(wl$lower), upper = exp(wl$upper))
  }
  settings <- list(model = "cox", statistic = statistic, variance = "poisson",
                   method = NA_character_, dispersion = if (statistic == "wald") dispersion else NA,
                   phi = phi, null = "ratio", null_value = 1, level = level,
                   alternative = alternative, min_expected = min_expected, n_dropped = n_dropped,
                   limits = switch(statistic, midp = "midp", exact = "exact", wald = "log_wald"))
  sizes <- fit$sizes[match(sm$provider, fit$sizes$provider), , drop = FALSE]
  rownames(sizes) <- NULL
  out <- pp_test_frame(sm$provider, obs / expected, se, 1, z, decision, lim$lower, lim$upper,
                       settings, sizes)
  out$observed <- obs
  out$expected <- expected
  out$person_time <- sm$person_time
  if (!is.null(null_model)) out <- calibrate(out, null_model)
  out
}

# Mid-p z statistic of Poisson counts against their expected means (the SMR
# convention; the floor 1e-6 caps |z| at 4.75).
midp_z <- function(observed, expected) {
  p_low <- 2 * stats::ppois(observed, expected) - stats::dpois(observed, expected)
  p_high <- 2 * (1 - stats::ppois(observed - 1, expected)) - stats::dpois(observed, expected)
  p <- pmax(1e-6, pmin(p_low, p_high) / 2)
  ifelse(p_low <= p_high, stats::qnorm(p), -stats::qnorm(p))
}

# Limits for the ratio O/E from inverting the mid-p test calibrated by a null
# N(mean, sd^2) with critical value q: the Poisson means at which z(O, mean)
# equals mean + q sd (lower) and mean - q sd (upper), over E. z(O, t) falls
# with t, so each equation has at most one root; 0 (lower) and Inf (upper)
# where it has none.
midp_limits <- function(observed, expected, mean, sd, q, alternative) {
  k <- length(observed)
  lower <- upper <- rep(NA_real_, k)
  for (i in seq_len(k)) {
    o <- observed[i]
    e <- expected[i]
    if (!is.finite(o) || !is.finite(e) || e <= 0 || !is.finite(mean[i]) || !is.finite(sd[i])) {
      next
    }
    lower[i] <- if (alternative == "less") 0 else midp_root(o, e, mean[i] + q * sd[i])
    upper[i] <- if (alternative == "greater") Inf else midp_root(o, e, mean[i] - q * sd[i])
  }
  list(lower = lower / expected, upper = upper / expected)
}

# The Poisson mean t at which the decreasing function midp_z(o, t) equals
# `target`: 0 when it is at or below the target already at the smallest mean,
# Inf when it never falls to it.
midp_root <- function(o, e, target) {
  lo <- 1e-10 * max(e, 1)
  if (midp_z(o, lo) <= target) return(0)
  hi <- 10 * (o + e + 10)
  while (midp_z(o, hi) > target) {
    if (hi > 1e12) return(Inf)
    hi <- hi * 10
  }
  stats::uniroot(function(t) midp_z(o, t) - target, c(lo, hi),
                 tol = 1e-10 * max(e, 1))$root
}

# The exact Poisson test of counts against expected means: p-values (doubled
# tail, capped at 0.999), the normal quantile of the two-sided p signed by
# the direction as the statistic, flags, and exact (E < 100) or Byar limits
# for the ratio.
exact_poisson <- function(observed, expected, level, alternative, normal_threshold = 100) {
  alpha <- 1 - level
  high <- observed / expected > 1
  p_high <- 1 - stats::ppois(observed - 1, expected)
  p_low <- stats::ppois(observed, expected)
  p_two <- pmin(0.999, 2 * ifelse(high, p_high, p_low))
  z <- sign(observed - expected) * stats::qnorm(1 - p_two / 2)
  if (alternative == "two.sided") {
    p <- p_two
    flag <- ifelse(p < alpha, sign(observed - expected), 0)
    a_low <- a_high <- alpha / 2
  } else if (alternative == "greater") {
    p <- p_high
    flag <- ifelse(p < alpha, 1, 0)
    a_low <- alpha
    a_high <- 0
  } else {
    p <- p_low
    flag <- ifelse(p < alpha, -1, 0)
    a_low <- 0
    a_high <- alpha
  }
  flag <- as.integer(flag)
  approx <- expected >= normal_threshold
  lower <- ifelse(observed == 0, 0, stats::qchisq(a_low, 2 * observed) / 2)
  upper <- stats::qchisq(1 - a_high, 2 * (observed + 1)) / 2
  if (any(approx)) {
    o <- observed[approx]
    z_low <- stats::qnorm(1 - a_low)
    z_high <- stats::qnorm(1 - a_high)
    lower[approx] <- ifelse(o == 0, 0, o * (1 - 1 / (9 * o) - z_low / (3 * sqrt(o)))^3)
    upper[approx] <- (o + 1) * (1 - 1 / (9 * (o + 1)) + z_high / (3 * sqrt(o + 1)))^3
  }
  if (alternative == "less") lower[] <- 0
  if (alternative == "greater") upper[] <- Inf
  bad <- !is.finite(observed) | !is.finite(expected) | expected <= 0
  z[bad] <- NA_real_
  p[bad] <- NA_real_
  flag[bad] <- NA_integer_
  list(z = z, p = p, flag = flag, lower = lower / expected, upper = upper / expected)
}

# Deviance dispersion of the tested providers' ratios: the Poisson deviance
# of the unit-level counts (patients, or rows) about ratio_j * e_k, over the
# number of units minus the number of providers.
deviance_dispersion <- function(fit, providers) {
  m <- cox_matrices(fit)
  unit <- if (is.null(m$patient)) seq_along(m$code) else cluster_codes(m$code, m$patient)$code
  o_unit <- as.numeric(rowsum(m$w * m$event, unit, reorder = TRUE))
  e_unit <- as.numeric(rowsum(m$w * fit$expected, unit, reorder = TRUE))
  prov_unit <- m$code[!duplicated(unit)][order(unique(unit))]
  tested <- match(providers, m$levels)
  keep <- prov_unit %in% tested
  o_prov <- as.numeric(rowsum(o_unit[keep], prov_unit[keep], reorder = TRUE))
  e_prov <- as.numeric(rowsum(e_unit[keep], prov_unit[keep], reorder = TRUE))
  ratio <- stats::setNames(o_prov / e_prov, sort(unique(prov_unit[keep])))
  mu <- ratio[as.character(prov_unit[keep])] * e_unit[keep]
  o <- o_unit[keep]
  dev <- 2 * (ifelse(o > 0, o * log(o / mu), 0) - (o - mu))
  df <- sum(keep) - length(tested)
  if (df <= 0) stop("The deviance dispersion needs more units than providers.", call. = FALSE)
  sum(dev) / df
}

#' Confidence limits of the standardized ratios of a Cox fit
#'
#' The limits of [test.cox()] for the ratio \eqn{O_j / E_j}.
#'
#' @param object a model fitted by [fit_cox()].
#' @param parm optional subset of provider identifiers.
#' @param level confidence level. Default 0.95.
#' @param statistic,dispersion,alternative,min_expected,null_model as in
#'   [test.cox()].
#' @param ... unused.
#' @return A data frame with `provider`, `estimate` (the ratio), `lower`,
#'   `upper`.
#' @examples
#' s <- sim_kecc(n_providers = 20, type = "survival", seed = 1)
#' ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", time = "stop", start = "start",
#'               event = "event", cutoff = 0, message = FALSE)
#' confint(fit_cox(ps))
#' @export
confint.cox <- function(object, parm = NULL, level = 0.95, statistic = c("midp", "exact", "wald"),
                        dispersion = c("none", "deviance"),
                        alternative = c("two.sided", "greater", "less"), min_expected = 0,
                        null_model = NULL, ...) {
  statistic <- match.arg(statistic)
  dispersion <- match.arg(dispersion)
  alternative <- match.arg(alternative)
  tst <- test(object, statistic = statistic, dispersion = dispersion, providers = parm,
                        level = level, alternative = alternative, min_expected = min_expected,
                        null_model = null_model)
  data.frame(provider = tst$provider, estimate = tst$estimate, lower = tst$lower,
             upper = tst$upper, stringsAsFactors = FALSE)
}
