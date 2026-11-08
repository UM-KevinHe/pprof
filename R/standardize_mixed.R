#' Standardized measures of a random-effect fit
#'
#' Indirect or direct standardization of every provider's outcome from a
#' random-effect model: ratios and rates for the logistic models,
#' differences for the linear models.
#'
#' @param fit a model from [fit_logistic_re()], [fit_logistic_cre()],
#'   [fit_linear_re()] or [fit_linear_cre()].
#' @param type `"indirect"` (default) or `"direct"`.
#' @param numerator for indirect standardization, the provider's
#'   `"observed"` outcome (default) or its `"predicted"` outcome under the
#'   fitted model with the provider's random effect, \eqn{\sum_{j \in i}
#'   N_{ij}\, \hat p_{ij}} (a shrunken numerator, as in predicted-to-expected
#'   ratios).
#' @param null for indirect standardization, the norm of the random effect:
#'   `"median"` (default) or `"mean"` (weighted by rows) of the predicted
#'   effects, or a number (0 is the population intercept).
#' @param n_pop denominator of the rates (logistic models); defaults to the
#'   number of trials among the fitted rows.
#' @param providers optional subset of provider identifiers.
#' @param alpha_se optional standard errors of the random effects, named by
#'   provider; with them the direct rate's standard error (logistic) or the
#'   difference's standard error (linear) is returned.
#' @param level if given, confidence limits are added, mapped from Wald
#'   limits of \eqn{\alpha_i} (see [confint.logistic_re()]); for the
#'   logistic models only with `numerator = "predicted"`, whose limits are
#'   the numerator evaluated at the limits of \eqn{\alpha_i}.
#' @param alternative `"two.sided"` (default), `"greater"` or `"less"`.
#' @param threads OpenMP threads for direct standardization.
#' @param ... unused.
#'
#' @details With \eqn{\eta_{ij} = \hat\mu + Z_{ij}^\top\hat\beta} the fixed
#' part of the linear predictor and \eqn{p_{ij}(a) = \mathrm{expit}(a +
#' \eta_{ij})}:
#' * logistic, indirect: \eqn{E_i = \sum_{j \in i} N_{ij}\, p_{ij}(\alpha_0)},
#'   ratio = numerator / \eqn{E_i};
#' * logistic, direct: \eqn{E_i = \sum_k N_k\, p_k(\hat\alpha_i)} over every
#'   fitted row, ratio \eqn{E_i / \sum O};
#' * linear, indirect: \eqn{E_i = \sum_{j \in i} (\alpha_0 + \eta_{ij})},
#'   difference = (numerator \eqn{- E_i) / n_i};
#' * linear, direct: \eqn{E_i = \sum_k (\hat\alpha_i + \eta_k)}, difference
#'   \eqn{(E_i - \sum O) / n}.
#'
#' Rates are ratio \eqn{\times \sum O / n_{pop}}, truncated to \[0, 1\].
#'
#' @return A data frame with one row per provider. Logistic models:
#'   `provider`, `observed` (the numerator), `expected`, `ratio`, `rate`,
#'   `variance_expected` and, when requested, `se_rate`, `lower_ratio`,
#'   `upper_ratio`, `lower_rate`, `upper_rate`. Linear models: `provider`,
#'   `observed`, `expected`, `difference` and, when requested,
#'   `se_difference`, `lower`, `upper`. Attributes `type`, `numerator`,
#'   `null_value`, `level` (and `n_pop`, `population_rate` for the logistic
#'   models).
#' @examples
#' \donttest{
#' d <- sim_kecc(n_providers = 20, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_logistic_re(pd)
#' head(standardize(fit, numerator = "predicted", level = 0.95))
#' }
#' @name standardize_mixed
NULL

#' @rdname standardize_mixed
#' @export
standardize.logistic_re <- function(fit, type = c("indirect", "direct"),
                                    numerator = c("observed", "predicted"), null = "median",
                                    n_pop = NULL, providers = NULL, alpha_se = NULL,
                                    level = NULL,
                                    alternative = c("two.sided", "greater", "less"),
                                    threads = 1, ...) {
  check_mixed(fit)
  type <- match.arg(type)
  numerator <- match.arg(numerator)
  alternative <- match.arg(alternative)
  m <- fit_matrices(fit)
  alpha <- fit$coefficients$alpha
  idx <- provider_positions(providers, m$levels)
  code <- provider_codes(m$prov, m$levels)
  lp <- fit$linear_predictor
  total_obs <- sum(m$y)
  total_trials <- sum(m$n)
  if (is.null(n_pop)) n_pop <- total_trials
  population_rate <- total_obs / n_pop
  clip <- function(x) pmax(pmin(x, 1), 0)
  predicted <- as.numeric(rowsum(m$n * fit$fitted, code, reorder = TRUE))
  threads <- max(1L, as.integer(threads))

  if (type == "indirect") {
    a0 <- null_value(null, alpha, sizes = m$n_prov)
    p0 <- stats::plogis(a0 + lp)
    expected <- as.numeric(rowsum(m$n * p0, code, reorder = TRUE))
    variance_expected <- as.numeric(rowsum(m$n * p0 * (1 - p0), code, reorder = TRUE))
    observed <- if (numerator == "observed") {
      as.numeric(rowsum(m$y, code, reorder = TRUE))
    } else {
      predicted
    }
    ratio <- observed / expected
  } else {
    a0 <- NA_real_
    res <- computeDirectExp(unname(alpha), lp, m$n, threads)
    expected <- res[, 1]
    variance_expected <- res[, 2]
    observed <- rep(total_obs, length(alpha))
    ratio <- expected / total_obs
  }
  out <- data.frame(provider = m$levels, observed = observed, expected = expected, ratio = ratio,
                    rate = clip(ratio * population_rate), variance_expected = variance_expected,
                    stringsAsFactors = FALSE)

  if (!is.null(alpha_se)) {
    if (type == "indirect") {
      stop("`alpha_se` gives the standard error of the direct rate; use type = \"direct\".",
           call. = FALSE)
    }
    out$se_rate <- variance_expected / n_pop * effect_se_vector(alpha_se, m$levels)
  }

  if (!is.null(level)) {
    if (type == "indirect" && numerator != "predicted") {
      stop("Limits of the indirect ratio need numerator = \"predicted\".", call. = FALSE)
    }
    limits <- wald_limits(alpha, sqrt(fit$post_var), level, alternative)
    if (type == "indirect") {
      rows_of <- split(seq_along(code), factor(code, levels = seq_along(m$levels)))
      numerator_at <- function(a) {
        mapply(function(r, ai) {
          if (!is.finite(ai)) return(if (ai > 0) sum(m$n[r]) else 0)
          sum(m$n[r] * stats::plogis(ai + lp[r]))
        }, rows_of, a)
      }
      lower <- numerator_at(limits$lower) / expected
      upper <- numerator_at(limits$upper) / expected
    } else {
      direct_at <- function(a) {
        val <- rep(NA_real_, length(a))
        fin <- is.finite(a)
        if (any(fin)) val[fin] <- computeDirectExp(a[fin], lp, m$n, threads)[, 1]
        val[!fin & a > 0] <- total_trials
        val[!fin & a < 0] <- 0
        val
      }
      lower <- direct_at(limits$lower) / total_obs
      upper <- direct_at(limits$upper) / total_obs
    }
    out$lower_ratio <- lower
    out$upper_ratio <- upper
    out$lower_rate <- clip(lower * population_rate)
    out$upper_rate <- clip(upper * population_rate)
  }

  out <- out[idx, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "type") <- type
  attr(out, "numerator") <- if (type == "indirect") numerator else "predicted"
  attr(out, "null_value") <- a0
  attr(out, "n_pop") <- n_pop
  attr(out, "population_rate") <- population_rate
  attr(out, "level") <- level
  out
}

#' @rdname standardize_mixed
#' @export
standardize.linear_re <- function(fit, type = c("indirect", "direct"),
                                  numerator = c("observed", "predicted"), null = "median",
                                  providers = NULL, alpha_se = NULL, level = NULL,
                                  alternative = c("two.sided", "greater", "less"), ...) {
  check_mixed(fit)
  type <- match.arg(type)
  numerator <- match.arg(numerator)
  alternative <- match.arg(alternative)
  m <- fit_matrices(fit)
  alpha <- fit$coefficients$alpha
  idx <- provider_positions(providers, m$levels)
  code <- provider_codes(m$prov, m$levels)
  lp <- fit$linear_predictor
  value <- if (numerator == "observed") m$y else fit$fitted

  if (type == "indirect") {
    a0 <- null_value(null, alpha, sizes = m$n_prov)
    observed <- as.numeric(rowsum(value, code, reorder = TRUE))
    expected <- as.numeric(rowsum(a0 + lp, code, reorder = TRUE))
    difference <- (observed - expected) / m$n_prov
  } else {
    a0 <- NA_real_
    observed <- rep(sum(value), length(alpha))
    expected <- unname(alpha) * length(lp) + sum(lp)
    difference <- (expected - observed) / length(lp)
  }
  out <- data.frame(provider = m$levels, observed = observed, expected = expected,
                    difference = difference, stringsAsFactors = FALSE)
  if (!is.null(alpha_se)) out$se_difference <- effect_se_vector(alpha_se, m$levels)
  if (!is.null(level)) {
    limits <- wald_limits(alpha, sqrt(fit$post_var), level, alternative)
    out$lower <- difference - (unname(alpha) - limits$lower)
    out$upper <- difference + (limits$upper - unname(alpha))
  }
  out <- out[idx, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "type") <- type
  attr(out, "numerator") <- numerator
  attr(out, "null_value") <- a0
  attr(out, "level") <- level
  out
}

# Standard errors of the provider effects in the order of `levels`, from a
# named vector or one value per provider.
effect_se_vector <- function(se, levels) {
  if (!is.null(names(se))) {
    se <- se[levels]
  } else if (length(se) != length(levels)) {
    stop("The standard errors must be named by provider or have one value per fitted provider.",
         call. = FALSE)
  }
  as.numeric(se)
}
