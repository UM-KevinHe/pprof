#' Test a standardized measure against a reference value
#'
#' Measure-level z-test shared by the KECC measures. With a transformation
#' \eqn{g} (identity or logit), the statistic is
#' \deqn{z_i = \frac{g(M_i) - g(M_{ref})}{g'(M_i)\, SE(M_i)}.}
#' Under the null, \eqn{z_i} is taken to follow \eqn{N(location, scale^2)}:
#' `location = 0, scale = 1` is the theoretical null, a fixed inflation such
#' as `scale = 1.81` is the PPPW convention, and empirical-null estimates
#' computed by the measure code can be passed in the same way.
#'
#' @param measure standardized measure \eqn{M_i} (e.g. a direct standardized
#'   rate), one value per provider.
#' @param se standard error of `measure`.
#' @param ref reference value \eqn{M_{ref}} (scalar or one value per provider),
#'   e.g. the mean of the measure or the population rate.
#' @param transform `"identity"` (default) or `"logit"`.
#' @param location,scale location and scale of the null distribution of z
#'   (scalars or one value per provider). Defaults 0 and 1.
#' @param level confidence level. Default 0.95.
#' @param crit critical value. Defaults to the normal quantile for `level`;
#'   the measures use the rounded value 1.96.
#' @param ci.location location used to centre the confidence interval.
#'   Defaults to `location`; PPPW centres its interval at 0 even when an
#'   empirical-null location is estimated.
#' @param alternative `"two.sided"` (default), `"greater"` or `"less"`.
#'
#' @return A data frame with one row per provider:
#' \item{est, se.t, ref.t}{the measure, its SE and the reference on the
#'   transformed scale.}
#' \item{z}{the test statistic.}
#' \item{prob}{\eqn{P(Z \le z)} under the null distribution.}
#' \item{p}{p-value.}
#' \item{flag}{1 = significantly higher than the reference, -1 = lower,
#'   0 = not significant; `NA` when z is missing.}
#' \item{lower.t, upper.t}{confidence limits on the transformed scale.}
#' \item{lower, upper}{confidence limits on the measure scale.}
#'
#' @details The interval inverts the calibrated test:
#' \eqn{g(M_i) - (ci.location \pm crit \cdot scale) \cdot SE_t}, mapped back
#' with \eqn{g^{-1}}. One-sided alternatives give one-sided intervals.
#'
#' @export
test_SM <- function(measure, se, ref, transform = "identity", location = 0, scale = 1,
                    level = 0.95, crit = NULL, ci.location = location, alternative = "two.sided") {
  if (!transform %in% c("identity", "logit")) stop("Argument 'transform' should be 'identity' or 'logit'.", call. = FALSE)
  alternative <- match_alternative(alternative)
  if (is.null(crit)) {
    crit <- if (alternative == "two.sided") stats::qnorm(1 - (1 - level) / 2) else stats::qnorm(level)
  }

  if (transform == "identity") {
    est <- measure
    se.t <- se
    ref.t <- ref
    back <- function(x) x
  } else {
    est <- log(measure / (1 - measure))
    se.t <- se / (measure * (1 - measure))
    ref.t <- log(ref / (1 - ref))
    back <- function(x) ifelse(x > 700, 1, exp(x) / (1 + exp(x))) # guard exp() overflow
  }
  z <- (est - ref.t) / se.t
  prob <- stats::pnorm(z, mean = location, sd = scale)
  zc <- (z - location) / scale

  if (alternative == "two.sided") {
    p <- 2 * pmin(prob, 1 - prob)
    flag <- ifelse(zc > crit, 1, ifelse(zc < -crit, -1, 0))
    lower.t <- est - (ci.location + crit * scale) * se.t
    upper.t <- est - (ci.location - crit * scale) * se.t
  } else if (alternative == "greater") {
    p <- 1 - prob
    flag <- ifelse(zc > crit, 1, 0)
    lower.t <- est - (ci.location + crit * scale) * se.t
    upper.t <- rep(Inf, length(est))
  } else {
    p <- prob
    flag <- ifelse(zc < -crit, -1, 0)
    lower.t <- rep(-Inf, length(est))
    upper.t <- est - (ci.location - crit * scale) * se.t
  }
  n <- max(length(est), length(ref.t))
  data.frame(est = rep_len(est, n), se.t = rep_len(se.t, n), ref.t = rep_len(ref.t, n),
             z = z, prob = prob, p = p, flag = flag,
             lower.t = lower.t, upper.t = upper.t,
             lower = back(lower.t), upper = back(upper.t))
}
