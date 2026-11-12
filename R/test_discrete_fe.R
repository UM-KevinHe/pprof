#' Tests of the provider effects of a discrete-time fixed-effect fit
#'
#' Tests \eqn{H_0: \gamma_i = \gamma_0} for every provider (or a subset) with
#' the score or the Wald statistic and a model-based or clustered variance,
#' and flags the providers above or below the null value.
#'
#' @param fit a model fitted by [fit_discrete_fe()].
#' @param statistic `"score"` (default) or `"wald"`.
#' @param variance `"naive"` (default; model-based) or `"robust"` (sandwich
#'   clustered by patient, or by subject when [pp_data()] has no patient
#'   column).
#' @param null the null value \eqn{\gamma_0}: `"median"` (default) or
#'   `"mean"` (weighted by person-periods) of the fitted effects, or a number.
#' @param providers optional subset of provider identifiers to report.
#' @param level confidence level; the significance level is `1 - level`.
#' @param alternative `"two.sided"` (default), `"greater"` or `"less"`.
#' @param wald_variance for `statistic = "wald"` with the naive variance:
#'   `"case_mix"` (default) or `"full"`, see [gamma_variance()].
#' @param null_model optional null model, as in [test.logistic_fe()]; the
#'   result is then passed through [calibrate()].
#' @param interval for the score statistic, compute the confidence limits of
#'   \eqn{\gamma_i} by inverting the test (default `TRUE`).
#' @param ... unused.
#'
#' @details With \eqn{H_j(g) = \sum_{k \le t_j} \mathrm{expit}(g + \hat\alpha_k
#' + Z_j^\top\hat\beta)} the expected events of subject \eqn{j} over its
#' periods at risk and \eqn{W_j(g)} the sum of \eqn{h(1 - h)}:
#' * score: \eqn{(O_i - \sum_{j \in i} H_j(\gamma_0)) / \sqrt{V_i}} with the
#'   naive \eqn{V_i = \sum_{j \in i} W_j(\gamma_0)} or the clustered
#'   \eqn{V_i = \sum_c (\sum_{j \in c} \delta_j - H_j(\gamma_0))^2};
#' * Wald: \eqn{(\hat\gamma_i - \gamma_0)/se}, with the variance of
#'   [gamma_variance()] or the clustered standard error \eqn{\sqrt{V_i} /
#'   \sum_{j \in i} W_j(\hat\gamma_i)}.
#'
#' These are the statistics of [test.logistic_fe()] on the person-period
#' data. The exact and likelihood-ratio statistics are not offered: a
#' subject leaves at its event, so a provider's event count is not a sum of
#' independent trials. Under the model the person-periods of a subject are
#' uncorrelated, so there is no stabilized variance.
#'
#' The limits of the score statistic invert the test as in
#' [test.logistic_fe()]; a provider with no events has no finite lower limit.
#'
#' @return A data frame of class `pp_test`, with the columns of
#'   [test.logistic_fe()].
#' @examples
#' s <- sim_kecc(n_providers = 20, provider_size = c(20, 40), rows_per_patient = c(1, 1),
#'               type = "survival", seed = 1)
#' s$month <- pmin(ceiling(s$stop * 12), 12)
#' s$event[s$stop * 12 > 12] <- 0
#' pd <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider", patient = "patient",
#'               time = "month", event = "event", message = FALSE)
#' fit <- fit_discrete_fe(pd)
#' test(fit)
#' test(fit, statistic = "wald", variance = "robust")
#' @export
test.discrete_fe <- function(fit, statistic = c("score", "wald"),
                             variance = c("naive", "robust"), null = "median", providers = NULL,
                             level = 0.95, alternative = c("two.sided", "greater", "less"),
                             wald_variance = c("case_mix", "full"), null_model = NULL,
                             interval = TRUE, ...) {
  check_fit(fit, "discrete_fe")
  if (length(statistic) == 1 && statistic %in% c("exact", "lr")) {
    stop("The discrete-time model offers the score and Wald statistics: a subject leaves at its ",
         "event, so a provider's event count is not a sum of independent trials.", call. = FALSE)
  }
  if (length(variance) == 1 && variance == "stabilized") {
    stop("The person-periods of a subject are uncorrelated under the model; use variance = ",
         "\"naive\" or \"robust\".", call. = FALSE)
  }
  statistic <- match.arg(statistic)
  variance <- match.arg(variance)
  alternative <- match.arg(alternative)
  wald_variance <- match.arg(wald_variance)
  check_level(level)

  m <- discrete_matrices(fit)
  gamma <- fit$coefficients$gamma
  code <- provider_codes(m$prov, m$levels)
  person_periods <- as.numeric(rowsum(m$t, code, reorder = TRUE))
  g0 <- null_value(null, gamma, sizes = person_periods)
  idx <- provider_positions(providers, m$levels)
  sel <- m$levels[idx]
  q <- normal_quantile(level, alternative)
  se <- rep(NA_real_, length(idx))
  lower <- upper <- rep(NA_real_, length(idx))
  extra <- inversion <- NULL
  cl <- if (variance == "robust") discrete_clusters(m) else NULL
  clustered_meat <- function(r) {
    cluster_sums(r, cl$code, code, code[cl$first], length(m$levels))$meat
  }

  if (statistic == "score") {
    s0 <- discrete_sums_at(fit, m, rep(g0, length(m$t)))
    r0 <- m$event - s0$H
    obs_all <- as.numeric(rowsum(m$event, code, reorder = TRUE))
    exp_all <- as.numeric(rowsum(s0$H, code, reorder = TRUE))
    v_all <- if (variance == "naive") {
      as.numeric(rowsum(s0$W, code, reorder = TRUE))
    } else {
      as.numeric(clustered_meat(r0))
    }
    z <- ((obs_all - exp_all) / sqrt(v_all))[idx]
    z[!is.finite(z)] <- NA_real_
    extra <- list(observed = obs_all[idx], expected = exp_all[idx], variance = v_all[idx])
    inversion <- list(fit = fit, variance = variance, pooled = NULL)
    if (isTRUE(interval)) {
      tg <- inversion_targets(rep(0, length(idx)), rep(1, length(idx)), q, alternative)
      lim <- score_limits(fit, variance, idx, tg$z_lower, tg$z_upper)
      lower <- lim$lower
      upper <- lim$upper
    }
  } else {
    se_all <- if (variance == "naive") {
      sqrt(gamma_variance(fit, wald_variance))
    } else {
      s <- discrete_sums_at(fit, m, gamma[code])
      i_gamma <- as.numeric(rowsum(s$W, code, reorder = TRUE))
      sqrt(as.numeric(clustered_meat(m$event - s$H))) / i_gamma
    }
    se <- unname(se_all[idx])
    z <- (gamma[idx] - g0) / se
    z[!is.finite(z)] <- NA_real_
    lower <- if (alternative == "less") rep(-Inf, length(idx)) else gamma[idx] - q * se
    upper <- if (alternative == "greater") rep(Inf, length(idx)) else gamma[idx] + q * se
  }
  decision <- z_decision(z, level, alternative)

  out <- data.frame(provider = sel, estimate = unname(gamma[idx]), se = se, null_value = g0,
                    statistic = unname(z), p = decision$p, flag = decision$flag,
                    lower = unname(lower), upper = unname(upper), stringsAsFactors = FALSE)
  if (!is.null(extra)) {
    out$observed <- extra$observed
    out$expected <- extra$expected
    out$variance <- extra$variance
  }
  attr(out, "test") <- list(statistic = statistic, variance = variance,
                            method = if (statistic == "score") {
                              "modified"
                            } else if (variance == "naive") {
                              wald_variance
                            } else {
                              NA
                            },
                            null = null, null_value = g0, level = level, alternative = alternative)
  attr(out, "sizes") <- data.frame(provider = sel, rows = m$n_prov[idx],
                                   patients = fit$sizes$patients[idx], stringsAsFactors = FALSE)
  if (!is.null(inversion)) attr(out, "inversion") <- inversion
  class(out) <- c("pp_test", "data.frame")
  if (!is.null(null_model)) out <- calibrate(out, null_model)
  out
}

# Score-test limits of the discrete-time model (see score_limits()): z_i(g)
# with the provider's expected events and variance re-evaluated at g.
score_limits_discrete <- function(fit, variance, idx, z_lower, z_upper) {
  m <- discrete_matrices(fit)
  gamma <- fit$coefficients$gamma
  alpha <- fit$coefficients$alpha
  lp <- fit$linear_predictor
  code <- provider_codes(m$prov, m$levels)
  rows_of <- split(seq_along(code), factor(code, levels = seq_along(m$levels)))
  cluster <- if (variance == "naive") NULL else discrete_clusters(m)$code
  lower <- upper <- rep(NA_real_, length(idx))
  for (k in seq_along(idx)) {
    i <- idx[k]
    r <- rows_of[[i]]
    tt <- m$t[r]
    zb <- lp[r]
    ev <- m$event[r]
    obs <- sum(ev)
    cl <- if (is.null(cluster)) NULL else match(cluster[r], unique(cluster[r]))
    z_of <- function(g) {
      s <- discrete_subject_sums(tt, alpha, g + zb)
      u <- obs - sum(s$H)
      v <- if (is.null(cl)) sum(s$W) else sum(rowsum(ev - s$H, cl, reorder = FALSE)^2)
      if (!is.finite(v) || v <= 0) return(sign(u) * 1e6)
      u / sqrt(v)
    }
    lower[k] <- monotone_root(z_of, gamma[i], z_lower[k])
    upper[k] <- monotone_root(z_of, gamma[i], z_upper[k])
  }
  list(lower = unname(lower), upper = unname(upper))
}

#' Confidence intervals for the provider effects of a discrete-time fit
#'
#' Intervals for \eqn{\gamma_i} by inverting the score test (naive or
#' clustered variance) or from a Wald test, as attached by
#' [test.discrete_fe()].
#'
#' @param object a model fitted by [fit_discrete_fe()].
#' @param parm optional subset of provider identifiers.
#' @param level confidence level. Default 0.95.
#' @param type `"score"` (default), `"robust_score"`, `"wald"` or
#'   `"robust_wald"`.
#' @param wald_variance for `type = "wald"`: `"case_mix"` (default) or
#'   `"full"`, see [gamma_variance()].
#' @param alternative `"two.sided"` (default), `"greater"` (lower limit
#'   only) or `"less"` (upper limit only).
#' @param ... unused.
#' @return A data frame with `provider`, `estimate`, `lower`, `upper`.
#' @export
confint.discrete_fe <- function(object, parm = NULL, level = 0.95,
                                type = c("score", "robust_score", "wald", "robust_wald"),
                                wald_variance = c("case_mix", "full"),
                                alternative = c("two.sided", "greater", "less"), ...) {
  check_fit(object, "discrete_fe")
  type <- match.arg(type)
  wald_variance <- match.arg(wald_variance)
  alternative <- match.arg(alternative)
  m <- discrete_matrices(object)
  idx <- provider_positions(parm, m$levels)
  r <- test(object, statistic = if (grepl("score", type)) "score" else "wald",
            variance = if (grepl("^robust", type)) "robust" else "naive",
            providers = m$levels[idx], level = level, alternative = alternative,
            wald_variance = wald_variance)
  data.frame(provider = m$levels[idx], estimate = unname(object$coefficients$gamma[idx]),
             lower = r$lower, upper = r$upper, stringsAsFactors = FALSE)
}

#' Standardized ratios and rates of a discrete-time fixed-effect fit
#'
#' Indirect or direct standardization of every provider's events over the
#' person-periods at risk, with the observed and expected events and
#' optional confidence intervals mapped from the limits of the provider
#' effects.
#'
#' @param fit a model fitted by [fit_discrete_fe()].
#' @param type `"indirect"` (default) or `"direct"`.
#' @param null for indirect standardization, the norm of the provider effect:
#'   `"median"` (default), `"mean"` (weighted by person-periods) or a number.
#' @param n_pop denominator of the rates. Defaults to the number of
#'   person-periods at risk.
#' @param providers optional subset of provider identifiers.
#' @param gamma_se optional standard errors of the provider effects, named by
#'   provider; with them the direct rate's standard error is returned.
#' @param level if given, confidence limits of the ratio and the rate are
#'   added, mapped from the limits of \eqn{\gamma_i} obtained with `ci`.
#' @param ci the interval used for the limits, see [confint.discrete_fe()].
#' @param alternative `"two.sided"` (default), `"greater"` or `"less"`.
#' @param threads OpenMP threads for direct standardization.
#' @param ... unused.
#'
#' @details With \eqn{H_j(g)} the expected events of subject \eqn{j} over its
#' periods at risk under the provider effect \eqn{g}:
#' * indirect: \eqn{E_i = \sum_{j \in i} H_j(\gamma_0)}, ratio \eqn{O_i /
#'   E_i};
#' * direct: \eqn{E_i = \sum_j H_j(\hat\gamma_i)} over every subject, ratio
#'   \eqn{E_i / \sum O}.
#'
#' In both cases rate = ratio \eqn{\times \sum O / n_{pop}}, the events per
#' person-period. These are the measures of [standardize.logistic_fe()] on
#' the person-period data.
#'
#' @return A data frame with the columns of [standardize.logistic_fe()].
#' @export
standardize.discrete_fe <- function(fit, type = c("indirect", "direct"), null = "median",
                                    n_pop = NULL, providers = NULL, gamma_se = NULL, level = NULL,
                                    ci = c("score", "robust_score", "wald", "robust_wald"),
                                    alternative = c("two.sided", "greater", "less"), threads = 1,
                                    ...) {
  check_fit(fit, "discrete_fe")
  type <- match.arg(type)
  ci <- match.arg(ci)
  alternative <- match.arg(alternative)
  m <- discrete_matrices(fit)
  gamma <- fit$coefficients$gamma
  alpha <- fit$coefficients$alpha
  lp <- fit$linear_predictor
  idx <- provider_positions(providers, m$levels)
  code <- provider_codes(m$prov, m$levels)
  person_periods <- as.numeric(rowsum(m$t, code, reorder = TRUE))
  total_obs <- sum(m$event)
  total_pp <- sum(m$t)
  if (is.null(n_pop)) n_pop <- total_pp
  population_rate <- total_obs / n_pop
  observed <- as.numeric(rowsum(m$event, code, reorder = TRUE))
  clip <- function(x) pmax(pmin(x, 1), 0)
  direct_at <- function(g) {
    val <- rep(NA_real_, length(g))
    fin <- is.finite(g)
    if (any(fin)) {
      val[fin] <- discrete_direct(g[fin], m$t, alpha, lp, max(1L, as.integer(threads)))[, 1]
    }
    val[!fin & g > 0] <- total_pp
    val[!fin & g < 0] <- 0
    val
  }

  if (type == "indirect") {
    g0 <- null_value(null, gamma, sizes = person_periods)
    s0 <- discrete_sums_at(fit, m, rep(g0, length(m$t)))
    expected <- as.numeric(rowsum(s0$H, code, reorder = TRUE))
    variance_expected <- as.numeric(rowsum(s0$W, code, reorder = TRUE))
    ratio <- observed / expected
  } else {
    g0 <- NA_real_
    res <- discrete_direct(unname(gamma), m$t, alpha, lp, max(1L, as.integer(threads)))
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
    limits <- stats::confint(fit, level = level, type = ci, alternative = alternative)
    if (type == "indirect") {
      rows_of <- split(seq_along(code), factor(code, levels = seq_along(m$levels)))
      expected_at <- function(g) {
        mapply(function(r, gi) {
          if (!is.finite(gi)) return(if (gi > 0) sum(m$t[r]) else 0)
          sum(discrete_subject_sums(m$t[r], alpha, gi + lp[r])$H)
        }, rows_of, g)
      }
      lower <- expected_at(limits$lower) / expected
      upper <- expected_at(limits$upper) / expected
    } else {
      lower <- direct_at(limits$lower) / total_obs
      upper <- direct_at(limits$upper) / total_obs
    }
    lower[fit$sizes$no_events] <- 0
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
  as_pp_measure(out)
}

#' Tests of the covariates of a discrete-time fixed-effect fit
#'
#' Wald tests of the covariate coefficients with the model-based variance or
#' the sandwich clustered by patient (by subject without a patient column),
#' the provider effects profiled out and the period effects estimated
#' jointly.
#'
#' @param object a model fitted by [fit_discrete_fe()].
#' @param test `"wald"` (default) or `"robust_wald"`.
#' @param covariates optional subset of covariate names.
#' @param level confidence level. Default 0.95.
#' @param null the null value of the coefficients. Default 0.
#' @param ... unused.
#' @return A data frame of class `pp_summary` with `covariate`, `estimate`,
#'   `se`, `statistic`, `p`, `lower`, `upper`.
#' @export
summary.discrete_fe <- function(object, test = c("wald", "robust_wald"), covariates = NULL,
                                level = 0.95, null = 0, ...) {
  check_fit(object, "discrete_fe")
  test <- match.arg(test)
  check_level(level)
  x <- object$vars$x
  ind <- if (is.null(covariates)) seq_along(x) else match(covariates, x)
  if (anyNA(ind)) {
    stop("Unknown covariate(s): ", paste(covariates[is.na(ind)], collapse = ", "), call. = FALSE)
  }
  inf <- discrete_information(object)
  k <- length(object$time_grid) - 1
  pos <- k + ind
  cov_theta <- if (test == "wald") {
    inf$theta_cov
  } else {
    m <- discrete_matrices(object)
    cf <- object$coefficients
    cl <- discrete_clusters(m)
    meat <- discrete_fe_theta_meat(m$t, m$event, m$z, m$n_prov, cl$code, cf$alpha, cf$gamma,
                                   cf$beta, inf$i_gamma_theta / inf$i_gamma)
    inf$theta_cov %*% meat %*% inf$theta_cov
  }
  beta <- object$coefficients$beta[ind]
  se <- sqrt(diag(cov_theta)[pos])
  stat <- (beta - null) / se
  q <- stats::qnorm(1 - (1 - level) / 2)
  pp_summary(data.frame(covariate = x[ind], estimate = unname(beta), se = se,
                        statistic = unname(stat), p = 2 * stats::pnorm(-abs(unname(stat))),
                        lower = unname(beta - q * se), upper = unname(beta + q * se),
                        stringsAsFactors = FALSE))
}
