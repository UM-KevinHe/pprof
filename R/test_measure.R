#' Test standardized measures against a reference value
#'
#' Measure-level z test. With a transformation \eqn{g} (identity, logit or
#' log) the statistic is
#' \deqn{z_i = \frac{g(M_i) - g(M_{ref})}{g'(M_i)\, SE_i}.}
#' Under the null, \eqn{z_i} follows \eqn{N(\mathrm{location},
#' \mathrm{scale}^2)}: `location = 0, scale = 1` is the theoretical null, a
#' fixed inflation such as `scale = 1.81` is a known-scale null, and the
#' location and scale of an empirical null can be passed per provider.
#'
#' @param estimate the standardized measures \eqn{M_i}, one per provider.
#' @param se their standard errors.
#' @param reference the reference value \eqn{M_{ref}} (one number or one
#'   value per provider), for instance the population rate.
#' @param transform `"identity"` (default), `"logit"` or `"log"`.
#' @param location,scale location and scale of the null distribution of
#'   \eqn{z} (one number or one value per provider).
#' @param level confidence level. Default 0.95.
#' @param critical critical value; defaults to the normal quantile for
#'   `level`.
#' @param ci_location location used to centre the confidence interval;
#'   defaults to `location`.
#' @param alternative `"two.sided"` (default), `"greater"` or `"less"`.
#' @param providers optional provider identifiers for the output.
#'
#' @details The interval inverts the calibrated test:
#' \eqn{g(M_i) - (\mathrm{ci\_location} \pm \mathrm{critical} \cdot
#' \mathrm{scale})\, SE_i^{(g)}}, mapped back with \eqn{g^{-1}}. One-sided
#' alternatives give one-sided intervals.
#'
#' @return A data frame with one row per provider: `provider` (if given),
#'   `estimate`, `se`, `reference`, `estimate_t`, `se_t`, `reference_t` (on
#'   the transformed scale), `z`, `p`, `flag` (`1` above the reference, `-1`
#'   below, `0` not different, `NA` when `z` is missing), `lower`, `upper`
#'   (measure scale) and `lower_t`, `upper_t`.
#' @examples
#' test_measure(estimate = c(0.12, 0.30, 0.22), se = c(0.03, 0.04, 0.05),
#'              reference = 0.2, transform = "logit")
#' @export
test_measure <- function(estimate, se, reference, transform = c("identity", "logit", "log"),
                         location = 0, scale = 1, level = 0.95, critical = NULL,
                         ci_location = location, alternative = c("two.sided", "greater", "less"),
                         providers = NULL) {
  transform <- match.arg(transform)
  alternative <- match.arg(alternative)
  n <- length(estimate)
  if (length(se) != n) stop("`se` must have one value per estimate.", call. = FALSE)
  for (arg in c("reference", "location", "scale", "ci_location")) {
    len <- length(get(arg))
    if (!len %in% c(1, n)) {
      stop("`", arg, "` must be one number or one value per estimate.", call. = FALSE)
    }
  }
  if (is.null(critical)) critical <- normal_quantile(level, alternative)

  if (transform == "identity") {
    est_t <- estimate
    se_t <- se
    ref_t <- rep_len(reference, n)
    back <- function(x) x
  } else if (transform == "logit") {
    est_t <- stats::qlogis(estimate)
    se_t <- se / (estimate * (1 - estimate))
    ref_t <- stats::qlogis(rep_len(reference, n))
    back <- stats::plogis
  } else {
    est_t <- log(estimate)
    se_t <- se / estimate
    ref_t <- log(rep_len(reference, n))
    back <- exp
  }
  z <- (est_t - ref_t) / se_t
  location <- rep_len(location, n)
  scale <- rep_len(scale, n)
  ci_location <- rep_len(ci_location, n)
  prob <- stats::pnorm(z, mean = location, sd = scale)
  zc <- (z - location) / scale

  if (alternative == "two.sided") {
    p <- 2 * pmin(prob, 1 - prob)
    flag <- ifelse(zc > critical, 1L, ifelse(zc < -critical, -1L, 0L))
    lower_t <- est_t - (ci_location + critical * scale) * se_t
    upper_t <- est_t - (ci_location - critical * scale) * se_t
  } else if (alternative == "greater") {
    p <- 1 - prob
    flag <- ifelse(zc > critical, 1L, 0L)
    lower_t <- est_t - (ci_location + critical * scale) * se_t
    upper_t <- rep(Inf, n)
  } else {
    p <- prob
    flag <- ifelse(zc < -critical, -1L, 0L)
    lower_t <- rep(-Inf, n)
    upper_t <- est_t - (ci_location - critical * scale) * se_t
  }
  flag[!is.finite(z)] <- NA_integer_
  out <- data.frame(estimate = estimate, se = se, reference = rep_len(reference, n),
                    estimate_t = est_t, se_t = se_t, reference_t = ref_t, z = z, p = p,
                    flag = flag, lower = back(lower_t), upper = back(upper_t),
                    lower_t = lower_t, upper_t = upper_t, stringsAsFactors = FALSE)
  if (!is.null(providers)) out <- cbind(provider = as.character(providers), out,
                                        stringsAsFactors = FALSE)
  attr(out, "test") <- list(transform = transform, level = level, critical = critical,
                            alternative = alternative)
  out
}
