#' Confidence intervals for provider effects and standardized measures
#'
#' Ported from `pprof` and extended to binomial outcomes. Intervals for the
#' provider effects are obtained by inverting a test; intervals for the
#' standardized ratios / rates map the provider-effect limits through the
#' standardization.
#'
#' @param object a model fitted by [logis_fe()].
#' @param parm optional subset of providers.
#' @param level confidence level. Default 0.95.
#' @param test `"exact"` (default; inverts the exact Poisson-binomial test),
#'   `"score"` (inverts the modified score test), `"wald"` (model-based
#'   standard error) or `"robust_wald"` (patient-clustered standard error).
#' @param option `"SM"` (default) for standardized measures or `"gamma"` for
#'   provider effects.
#' @param stdz for `option = "SM"`: `"indirect"` (default), `"direct"`, or both.
#' @param null the population norm for indirect standardization.
#' @param measure `"rate"`, `"ratio"` or both (default).
#' @param alternative `"two.sided"` (default), `"greater"` (lower limit only)
#'   or `"less"` (upper limit only). Provider-effect intervals are two-sided.
#' @param N_pop denominator of the rates, as in [SR_output()].
#' @param ... unused.
#'
#' @details Providers with no events have no finite lower limit for the
#' provider effect and providers with all events have no finite upper limit;
#' as in `pprof`, their finite limit is one-sided at level `1 - level`. The
#' exact test works on the individual Bernoulli trials, so it can be slow for
#' providers with a very large number of trials; the score or Wald versions
#' are then preferable.
#'
#' @return For `option = "gamma"`, a data frame with `gamma`, `gamma.lower`,
#'   `gamma.upper`. For `option = "SM"`, a list with the requested
#'   `CI.indirect_ratio`, `CI.indirect_rate`, `CI.direct_ratio`,
#'   `CI.direct_rate` data frames.
#'
#' @export
confint.logis_fe <- function(object, parm, level = 0.95, test = "exact", option = "SM",
                             stdz = "indirect", null = "median", measure = c("rate", "ratio"),
                             alternative = "two.sided", N_pop, ...) {
  check_fit(object)
  if (!option %in% c("gamma", "SM")) stop("Argument 'option' should be 'gamma' or 'SM'", call. = FALSE)
  if (!test %in% c("exact", "score", "wald", "robust_wald")) stop("Argument 'test' NOT as required!", call. = FALSE)
  if (!"indirect" %in% stdz & !"direct" %in% stdz) stop("Argument 'stdz' NOT as required!", call. = FALSE)
  alternative <- match_alternative(alternative)
  if (option == "gamma" && alternative != "two.sided")
    stop("Provider effect (option = 'gamma') only supports two-sided confidence intervals.", call. = FALSE)
  alpha <- 1 - level

  cl <- object$char_list
  data <- object$data_include
  N <- trials_of(data, cl$N.char)
  prov <- data[, cl$prov.char]
  n.prov <- sapply(split(data[, cl$Y.char], prov), length)
  prov_levels <- names(n.prov)
  indices <- provider_indices(if (missing(parm)) NULL else parm, prov, prov_levels)
  selected <- prov_levels[indices]
  gamma <- object$df.prov$gamma_est
  Z_beta <- as.numeric(object$linear_pred)

  ci_gamma <- gamma_limits(object, test, indices, alpha, alternative, data, N, prov, prov_levels, Z_beta)
  rownames(ci_gamma) <- selected
  if (option == "gamma") {
    attr(ci_gamma, "description") <- "Provider Effects"
    return(ci_gamma)
  }

  if (missing(N_pop)) N_pop <- sum(N)
  population_rate <- sum(object$obs) / N_pop
  no_events <- as.logical(tapply(data$no.events, factor(as.character(prov), levels = prov_levels), max))[indices]
  all_events <- as.logical(tapply(data$all.events, factor(as.character(prov), levels = prov_levels), max))[indices]
  return_ls <- list()

  make_output <- function(est, lo, up, label) {
    out <- list()
    type <- switch(alternative, greater = "upper one-sided", less = "lower one-sided", "two-sided")
    ratio <- data.frame(est, lo, up, row.names = selected)
    colnames(ratio) <- c(paste0(label, "_ratio"), "CI_ratio.lower", "CI_ratio.upper")
    attr(ratio, "confidence_level") <- paste(level * 100, "%")
    attr(ratio, "type") <- type
    if ("ratio" %in% measure) out[[paste0("CI.", label, "_ratio")]] <- ratio
    if ("rate" %in% measure) {
      clip <- function(x) pmax(pmin(x * population_rate, 1), 0)
      rate <- data.frame(clip(est), clip(lo), clip(up), row.names = selected)
      colnames(rate) <- c(paste0(label, "_rate"), "CI_rate.lower", "CI_rate.upper")
      attr(rate, "confidence_level") <- paste(level * 100, "%")
      attr(rate, "type") <- type
      attr(rate, "population_rate") <- population_rate
      out[[paste0("CI.", label, "_rate")]] <- rate
    }
    out
  }

  if ("indirect" %in% stdz) {
    gamma.null <- gamma_null_value(null, object$gamma)
    grp <- factor(as.character(prov), levels = prov_levels)
    Exp <- sapply(split(N * stats::plogis(gamma.null + Z_beta), grp), sum)[indices]
    Obs <- object$df.prov$Obs_provider[indices]
    rows_of <- split(seq_along(prov), grp)[indices]
    expected_at <- function(g) mapply(function(r, gi) {
      if (!is.finite(gi)) return(if (gi > 0) sum(N[r]) else 0)
      sum(N[r] * stats::plogis(gi + Z_beta[r]))
    }, rows_of, g)
    lo <- expected_at(ci_gamma$gamma.lower) / Exp
    up <- expected_at(ci_gamma$gamma.upper) / Exp
    lo[no_events] <- 0
    up[all_events] <- (sapply(rows_of, function(r) sum(N[r])) / Exp)[all_events]
    if (alternative == "less") lo <- rep(0, length(lo))
    if (alternative == "greater") up <- sapply(rows_of, function(r) sum(N[r])) / Exp
    return_ls <- c(return_ls, make_output(Obs / Exp, lo, up, "indirect"))
  }

  if ("direct" %in% stdz) {
    Obs_all <- sum(object$obs)
    direct_at <- function(g) {
      out <- rep(NA_real_, length(g))
      fin <- is.finite(g)
      out[fin] <- computeDirectExp(g[fin], Z_beta, N, 1L)[, 1]
      out[!fin & g > 0] <- sum(N)
      out[!fin & g < 0] <- 0
      out
    }
    est <- computeDirectExp(gamma[indices], Z_beta, N, 1L)[, 1] / Obs_all
    lo <- direct_at(ci_gamma$gamma.lower) / Obs_all
    up <- direct_at(ci_gamma$gamma.upper) / Obs_all
    lo[no_events] <- 0
    up[all_events] <- sum(N) / Obs_all
    if (alternative == "less") lo <- rep(0, length(lo))
    if (alternative == "greater") up <- rep(sum(N) / Obs_all, length(up))
    return_ls <- c(return_ls, make_output(est, lo, up, "direct"))
  }
  return_ls
}

# Confidence limits for the provider effects of the providers `indices`.
gamma_limits <- function(object, test, indices, alpha, alternative, data, N, prov, prov_levels, Z_beta) {
  cl <- object$char_list
  gamma <- object$df.prov$gamma_est
  if (test %in% c("wald", "robust_wald")) {
    if (test == "wald") {
      se <- sqrt(attr(vcov.logis_fe(object), "var.gamma"))[indices]
    } else {
      se <- test_fe(object, test = "robust_wald")$se.gamma[indices]
    }
    q <- if (alternative == "two.sided") stats::qnorm(1 - alpha / 2) else stats::qnorm(1 - alpha)
    lo <- gamma[indices] - q * se
    up <- gamma[indices] + q * se
    if (alternative == "less") lo <- rep(-Inf, length(lo))
    if (alternative == "greater") up <- rep(Inf, length(up))
    return(data.frame(gamma = gamma[indices], gamma.lower = lo, gamma.upper = up))
  }

  grp <- factor(as.character(prov), levels = prov_levels)
  rows_of <- split(seq_along(prov), grp)
  Y <- data[, cl$Y.char]
  q2 <- stats::qnorm(alpha / 2, lower.tail = FALSE)
  q1 <- stats::qnorm(alpha, lower.tail = FALSE)
  no_events <- tapply(data$no.events, grp, max) == 1
  all_events <- tapply(data$all.events, grp, max) == 1

  limits_one <- function(i) {
    r <- rows_of[[i]]
    zb <- Z_beta[r]
    nt <- N[r]
    obs <- sum(Y[r])
    g_hat <- gamma[i]
    if (test == "exact") {
      trials <- function(G) {
        p <- stats::plogis(G + zb)
        if (all(nt == 1)) p else rep(p, nt)
      }
      upper_fun <- function(G) {
        p <- trials(G)
        if (alternative == "two.sided") poibin::ppoibin(obs - 1, p) + 0.5 * poibin::dpoibin(obs, p) - alpha / 2
        else poibin::ppoibin(obs, p) - alpha
      }
      lower_fun <- function(G) {
        p <- trials(G)
        if (alternative == "two.sided") 1 - poibin::ppoibin(obs, p) + 0.5 * poibin::dpoibin(obs, p) - alpha / 2
        else 1 - poibin::ppoibin(obs - 1, p) - alpha
      }
      # one-sided limits at level alpha for providers at the boundary
      none_fun <- function(G) exp(sum(nt * stats::plogis(-(G + zb), log.p = TRUE))) / 2 - alpha
      all_fun <- function(G) exp(sum(nt * stats::plogis(G + zb, log.p = TRUE))) / 2 - alpha
    } else {
      z_of <- function(G) {
        p <- stats::plogis(G + zb)
        (obs - sum(nt * p)) / sqrt(sum(nt * p * (1 - p)))
      }
      q <- if (alternative == "two.sided") q2 else q1
      upper_fun <- function(G) z_of(G) + q
      lower_fun <- function(G) z_of(G) - q
      none_fun <- function(G) {
        p <- stats::plogis(G + zb)
        q1 - sum(nt * p) / sqrt(sum(nt * p * (1 - p)))
      }
      all_fun <- function(G) {
        p <- stats::plogis(G + zb)
        sum(nt * (1 - p)) / sqrt(max(sum(nt * p * (1 - p)), 1e-20)) - q1
      }
    }
    span <- 10 + max(abs(zb))
    if (no_events[i]) {
      return(c(g_hat, -Inf, find_root(none_fun, -span, span)))
    }
    if (all_events[i]) {
      return(c(g_hat, find_root(all_fun, -span, span), Inf))
    }
    lo <- if (alternative == "less") -Inf else find_root(lower_fun, g_hat - span, g_hat)
    up <- if (alternative == "greater") Inf else find_root(upper_fun, g_hat, g_hat + span)
    c(g_hat, lo, up)
  }
  res <- t(vapply(indices, limits_one, numeric(3)))
  data.frame(gamma = res[, 1], gamma.lower = res[, 2], gamma.upper = res[, 3])
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
