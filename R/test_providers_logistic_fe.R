#' Tests of the provider effects of a logistic fixed-effect fit
#'
#' Tests \eqn{H_0: \gamma_i = \gamma_0} for every provider (or a subset) with
#' a chosen statistic and variance, and flags the providers above or below
#' the null value.
#'
#' @param fit a model fitted by [fit_logistic_fe()].
#' @param statistic `"exact"` (default), `"score"`, `"wald"` or `"lr"` (the
#'   likelihood-ratio statistic of \eqn{\gamma_i} with every other parameter
#'   held at its full-model estimate; the penalized likelihood of a Firth
#'   fit, as in ParaFLR).
#' @param variance `"naive"` (default; model-based), `"robust"`
#'   (patient-clustered sandwich) or `"stabilized"` (sandwich with one
#'   within-patient correlation matrix pooled over all patients). The
#'   clustered variances need a patient column in [pp_data()] and are not
#'   available for the exact test.
#' @param null the null value \eqn{\gamma_0}: `"median"` (default) or
#'   `"mean"` (weighted by rows) of the fitted effects, or a number.
#' @param providers optional subset of provider identifiers to report. Every
#'   provider still enters the stabilized pooling.
#' @param level confidence level; the significance level is `1 - level`.
#' @param alternative `"two.sided"` (default), `"greater"` or `"less"`.
#' @param exact_method `"poibin"` (default; the Poisson-binomial distribution
#'   of the provider's event count, each row contributing `n` Bernoulli
#'   trials) or `"bootstrap"` (the same by simulation).
#' @param score_type `"modified"` (default; \eqn{\hat\beta} held at the
#'   full-model estimate) or `"restricted"` (the other parameters re-solved
#'   under the null; naive variance only).
#' @param wald_variance for `statistic = "wald"` with the naive variance:
#'   `"case_mix"` (default) or `"full"`, see [gamma_variance()].
#' @param n_resample number of draws for `exact_method = "bootstrap"`.
#' @param null_model optional null model: a `pp_null` object or a function
#'   of the result returning one, for instance
#'   `function(t) null_empirical(t, n_groups = 4)`; the result is then passed
#'   through [calibrate()].
#' @param ... unused.
#'
#' @details Row residuals \eqn{r_j = Y_j - N_j p_j} and variances
#' \eqn{v_j = N_j p_j(1 - p_j)} are evaluated at the fitted \eqn{\hat p_j}
#' for the Wald statistics and at \eqn{p_j = \mathrm{expit}(\gamma_0 +
#' Z_j^\top\hat\beta)} for the score statistics. With
#' \eqn{U_i = \sum_{j \in i} r_j} and \eqn{I_i = \sum_{j \in i} v_j}:
#' * naive score: \eqn{U_i / \sqrt{I_i}};
#' * robust score: \eqn{U_i / \sqrt{V_i}}, \eqn{V_i = \sum_{k \in i}
#'   (\sum_{j \in k} r_j)^2} over the provider's patients \eqn{k};
#' * stabilized score: \eqn{V_i} replaced by \eqn{\sum_{k \in i} a_k^\top
#'   \hat R a_k} with \eqn{a_k} the \eqn{\sqrt{v_j}} of patient \eqn{k}'s rows
#'   and \eqn{\hat R} pooled over all patients by row position;
#' * naive Wald: \eqn{(\hat\gamma_i - \gamma_0)/\sqrt{\mathrm{Var}}} with the
#'   variance of [gamma_variance()];
#' * robust and stabilized Wald: standard error \eqn{\sqrt{V_i} / I_i};
#' * likelihood ratio: \eqn{-2[\ell^*(\hat\gamma_{(i)}, \hat\beta) -
#'   \ell^*(\hat\gamma, \hat\beta)]} with \eqn{\hat\gamma_{(i)}} the fitted
#'   effects except \eqn{\gamma_i = \gamma_0}, and \eqn{\ell^*} the
#'   (penalized, for a Firth fit) log-likelihood; the other parameters are
#'   not re-estimated under the null. The statistic reported is the signed
#'   square root, so that two-sided p-values equal the chi-square tail.
#'
#' The stabilized statistics take a patient's rows in their order in the
#' data, so sort by time within patient before [pp_data()]. A patient with
#' more than 5,000 rows stops the stabilized statistics.
#'
#' @return A data frame of class `pp_test` with one row per tested provider:
#'   `provider`, `estimate` (\eqn{\hat\gamma_i}), `se` (Wald statistics only),
#'   `null_value`, `statistic` (a z statistic; for the exact test the normal
#'   quantile of the mid-p, for `"lr"` the signed root of the chi-square
#'   statistic), `p`, `flag` (`1` above the null, `-1` below,
#'   `0` not different, `NA` not testable) and, for the Wald statistics,
#'   `lower` and `upper` limits for \eqn{\gamma_i}. Attributes `test` (the
#'   settings) and `sizes` (rows and patients per provider).
#'
#' @references
#' Wu W, Yang Y, Kang J, He K (2022). Improving large-scale estimation and
#' inference for profiling health care providers. *Statistics in Medicine*,
#' 41(15), 2840-2853.
#' @examples
#' d <- sim_kecc(n_providers = 20, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_logistic_fe(pd)
#' test_providers(fit, statistic = "score", variance = "robust")
#' @export
test_providers.logistic_fe <- function(fit, statistic = c("exact", "score", "wald", "lr"),
                                       variance = c("naive", "robust", "stabilized"),
                                       null = "median", providers = NULL, level = 0.95,
                                       alternative = c("two.sided", "greater", "less"),
                                       exact_method = c("poibin", "bootstrap"),
                                       score_type = c("modified", "restricted"),
                                       wald_variance = c("case_mix", "full"), n_resample = 10000,
                                       null_model = NULL, ...) {
  check_fit(fit)
  statistic <- match.arg(statistic)
  variance <- match.arg(variance)
  alternative <- match.arg(alternative)
  exact_method <- match.arg(exact_method)
  score_type <- match.arg(score_type)
  wald_variance <- match.arg(wald_variance)
  if (statistic == "exact" && variance != "naive") {
    stop("The exact test assumes independent trials; use variance = \"naive\".", call. = FALSE)
  }
  if (variance != "naive" && is.null(fit$vars$patient)) {
    stop("The clustered variances need a patient column in pp_data().", call. = FALSE)
  }
  if (score_type == "restricted" && variance != "naive") {
    stop("The restricted score test is only available with the naive variance.", call. = FALSE)
  }
  if (statistic == "lr" && variance != "naive") {
    stop("The likelihood-ratio statistic has no clustered variance; use variance = \"naive\".",
         call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("`level` must be a number between 0 and 1.", call. = FALSE)
  }

  m <- fit_matrices(fit)
  gamma <- fit$coefficients$gamma
  beta <- fit$coefficients$beta
  g0 <- null_value(null, gamma, sizes = m$n_prov)
  idx <- provider_positions(providers, m$levels)
  sel <- m$levels[idx]
  code <- provider_codes(m$prov, m$levels)
  lp <- fit$linear_predictor
  se <- rep(NA_real_, length(idx))
  z <- rep(NA_real_, length(idx))
  lower <- upper <- rep(NA_real_, length(idx))

  if (statistic == "exact") {
    p0 <- stats::plogis(g0 + lp)
    rows_of <- split(seq_along(code), factor(code, levels = seq_along(m$levels)))[idx]
    res <- vapply(rows_of, function(r) {
      obs <- sum(m$y[r])
      if (exact_method == "poibin") {
        pr <- if (all(m$n[r] == 1)) p0[r] else rep(p0[r], m$n[r])
        exact_poibin(obs, pr, level, alternative)
      } else {
        exact_bootstrap(obs, p0[r], m$n[r], n_resample, level, alternative)
      }
    }, numeric(3))
    z <- res[3, ]
    decision <- list(flag = as.integer(res[1, ]), p = res[2, ])
  } else if (statistic == "score") {
    p0 <- stats::plogis(g0 + lp)
    r0 <- m$y - m$n * p0
    if (variance == "naive") {
      if (score_type == "modified") {
        u <- as.numeric(rowsum(r0, code, reorder = TRUE))
        v0 <- as.numeric(rowsum(m$n * p0 * (1 - p0), code, reorder = TRUE))
        z <- (u / sqrt(v0))[idx]
      } else {
        z <- restricted_score(m$n, m$y, m$z, m$n_prov, gamma, beta, g0, idx)[idx]
      }
    } else {
      cl <- cluster_codes(m$prov, m$patient)
      prov_cl <- code[cl$first]
      s <- cluster_sums(r0, cl$code, code, prov_cl, length(m$levels))
      meat <- if (variance == "robust") {
        s$meat
      } else {
        stabilized_variance(r0, m$n * p0 * (1 - p0), cl, prov_cl, length(m$levels))
      }
      z <- (s$U / sqrt(meat))[idx]
    }
    z[!is.finite(z)] <- NA_real_
    decision <- z_decision(z, level, alternative)
  } else if (statistic == "lr") {
    firth <- identical(fit$control$penalty, "firth")
    lr <- penalized_lr(m$n, m$y, m$z, m$n_prov, gamma, beta, g0, idx, firth)[idx]
    z <- sign(gamma[idx] - g0) * sqrt(pmax(lr, 0))
    z[!is.finite(z)] <- NA_real_
    decision <- z_decision(z, level, alternative)
  } else {
    if (variance == "naive") {
      se_all <- sqrt(gamma_variance(fit, wald_variance))
    } else {
      p <- fit$fitted
      r <- m$y - m$n * p
      i_gamma <- as.numeric(rowsum(m$n * p * (1 - p), code, reorder = TRUE))
      meat <- if (variance == "robust") {
        robust_gamma_meat(r, m$prov, m$patient, m$levels)
      } else {
        cl <- cluster_codes(m$prov, m$patient)
        stabilized_variance(r, m$n * p * (1 - p), cl, code[cl$first], length(m$levels))
      }
      se_all <- sqrt(meat) / i_gamma
    }
    se <- unname(se_all[idx])
    z <- (gamma[idx] - g0) / se
    z[!is.finite(z)] <- NA_real_
    decision <- z_decision(z, level, alternative)
    q <- normal_quantile(level, alternative)
    lower <- if (alternative == "less") rep(-Inf, length(idx)) else gamma[idx] - q * se
    upper <- if (alternative == "greater") rep(Inf, length(idx)) else gamma[idx] + q * se
  }

  out <- data.frame(provider = sel, estimate = unname(gamma[idx]), se = se, null_value = g0,
                    statistic = unname(z), p = decision$p, flag = decision$flag,
                    lower = unname(lower), upper = unname(upper), stringsAsFactors = FALSE)
  attr(out, "test") <- list(statistic = statistic, variance = variance,
                            method = switch(statistic, exact = exact_method, score = score_type,
                                            wald = if (variance == "naive") wald_variance else NA,
                                            lr = if (identical(fit$control$penalty, "firth")) {
                                              "penalized"
                                            } else {
                                              "likelihood"
                                            }),
                            null = null, null_value = g0, level = level, alternative = alternative)
  attr(out, "sizes") <- data.frame(provider = sel, rows = m$n_prov[idx],
                                   patients = fit$sizes$patients[idx], stringsAsFactors = FALSE)
  class(out) <- c("pp_test", "data.frame")
  if (!is.null(null_model)) out <- calibrate(out, null_model)
  out
}

# Pooled stabilized variance per provider; stops when a patient has too many rows.
stabilized_variance <- function(residual, v, cl, prov_cl, m) {
  rows <- max(tabulate(cl$code))
  if (rows > 5000) {
    stop("The stabilized variance pools a matrix over the rows of a patient; the largest patient ",
         "has ", rows, " rows.", call. = FALSE)
  }
  stabilized_meat(residual, v, cl$code, prov_cl, m)
}

# Exact test for one provider from the Poisson-binomial distribution of its
# event count; `probs` holds one success probability per Bernoulli trial.
# Returns c(flag, p-value, z-score).
exact_poibin <- function(obs, probs, level, alternative) {
  alpha <- 1 - level
  if (alternative == "two.sided") {
    p_upper <- 1 - poibin::ppoibin(obs, probs) + 0.5 * poibin::dpoibin(obs, probs)
    c(flag_two_sided(p_upper, alpha), 2 * min(p_upper, 1 - p_upper),
      stats::qnorm(p_upper, lower.tail = FALSE))
  } else if (alternative == "greater") {
    p_upper <- 1 - poibin::ppoibin(obs - 1, probs)
    c(if (p_upper < alpha) 1L else 0L, p_upper, stats::qnorm(p_upper, lower.tail = FALSE))
  } else {
    p_lower <- poibin::ppoibin(obs, probs)
    c(if (p_lower < alpha) -1L else 0L, p_lower, stats::qnorm(p_lower, lower.tail = TRUE))
  }
}

# The same test by simulation.
exact_bootstrap <- function(obs, probs, trials, n_resample, level, alternative) {
  alpha <- 1 - level
  size <- if (all(trials == 1)) 1 else rep(trials, times = n_resample)
  draws <- stats::rbinom(length(probs) * n_resample, size = size,
                         prob = rep(probs, times = n_resample))
  sums <- colSums(matrix(draws, ncol = n_resample))
  if (alternative == "two.sided") {
    p_upper <- (sum(sums > obs) + 0.5 * sum(sums == obs)) / n_resample
    c(flag_two_sided(p_upper, alpha), 2 * min(p_upper, 1 - p_upper),
      stats::qnorm(p_upper, lower.tail = FALSE))
  } else if (alternative == "greater") {
    p_upper <- sum(sums >= obs) / n_resample
    c(if (p_upper < alpha) 1L else 0L, p_upper, stats::qnorm(p_upper, lower.tail = FALSE))
  } else {
    p_lower <- sum(sums <= obs) / n_resample
    c(if (p_lower < alpha) -1L else 0L, p_lower, stats::qnorm(p_lower, lower.tail = TRUE))
  }
}

#' @export
print.pp_test <- function(x, digits = 4, ...) {
  s <- attr(x, "test")
  cat("Provider tests: ", s$statistic, " statistic, ", s$variance, " variance",
      if (!is.na(s$method)) paste0(" (", s$method, ")"), "; null ",
      if (is.character(s$null)) s$null else "value", " = ", format(s$null_value, digits = digits),
      "; ", s$alternative, ", level ", s$level, "\n", sep = "")
  f <- x$flag
  cat("  ", nrow(x), " providers: ", sum(f == 1, na.rm = TRUE), " above, ",
      sum(f == -1, na.rm = TRUE), " below, ", sum(f == 0, na.rm = TRUE), " not different, ",
      sum(is.na(f)), " not tested\n", sep = "")
  print.data.frame(utils::head(as.data.frame(x), 6), digits = digits, row.names = FALSE)
  if (nrow(x) > 6) cat("  ... ", nrow(x) - 6, " more rows\n", sep = "")
  invisible(x)
}
