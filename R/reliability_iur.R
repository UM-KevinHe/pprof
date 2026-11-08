#' Inter-unit reliability of a provider measure
#'
#' The fraction of the variation of a provider-level measure that reflects
#' true differences between providers rather than within-provider sampling
#' noise, \eqn{\mathrm{IUR} = \sigma^2_b / (\sigma^2_b + \sigma^2_w / n')},
#' from an ANOVA-type decomposition of the measure's variance (He et al.
#' 2019). Three estimators, chosen by `method`:
#' * `"direct"` uses provider-level estimates and their standard errors
#'   (for instance a standardized rate and the `se_rate` of [standardize()]);
#' * `"bootstrap"` resamples the rows of every provider with replacement
#'   and takes the variance of the resampled measures as the within-provider
#'   variance;
#' * `"split_half"` splits every provider's rows into two random halves
#'   and converts the agreement of the two halves' measures across providers
#'   into a reliability by the Spearman-Brown formula \eqn{2r / (1 + r)}; six
#'   agreement coefficients are reported (weighted kappa, Kendall's tau and
#'   Spearman's rho of the measures cut into ordered categories, and Pearson's
#'   r, Kendall's tau and Spearman's rho of the measures themselves).
#'
#' @param method `"direct"`, `"bootstrap"` or `"split_half"`.
#' @param estimate,se provider-level estimates and standard errors (direct
#'   method).
#' @param size number of observations per provider.
#' @param obs,expected observed and expected outcome per row (patient;
#'   bootstrap and split-half methods).
#' @param provider provider identifier per row.
#' @param measure a function `(obs, expected, provider)` returning one value
#'   per provider in sorted provider order; by default the ratio of sums
#'   \eqn{\sum O / \sum E}.
#' @param n_boot number of bootstrap replicates. Default 100.
#' @param n_iter number of random splits. Default 10.
#' @param category_probs quantile probabilities that cut the measures into
#'   ordered categories for the categorical coefficients.
#' @param seed optional seed; the session's random-number state is restored
#'   afterwards.
#'
#' @details With sizes \eqn{n_i}, within variances \eqn{v_i} of the
#' provider-level measure and measures \eqn{T_i}: the pooled within variance
#' \eqn{s^2_{w} = \sum (n_i - 1) v_i / \sum (n_i - 1)}, the effective size
#' \eqn{n' = (N - \sum n_i^2 / N) / (m - 1)}, the total variance
#' \eqn{s^2_t = \sum n_i (T_i - \bar T)^2 / (n'(m - 1))} with \eqn{\bar T}
#' the size-weighted mean, the between variance \eqn{s^2_b = s^2_t - s^2_w},
#' and \eqn{\mathrm{IUR} = s^2_b / s^2_t}. The reported `s2_within` is
#' \eqn{n' s^2_w}, the record-level within variance, so that the reliability
#' of a provider of size \eqn{n} is \eqn{s^2_b / (s^2_b + s2\_within / n)}
#' (see [iur_curve()]).
#'
#' @return An object of class `pp_iur`: a list with `method`, `iur` (for the
#'   split-half method a named vector of the six coefficients' reliabilities),
#'   and for the direct and bootstrap methods `s2_between`, `s2_within`,
#'   `n_prime`, `n_providers`, `providers`, `sizes`, `measure`,
#'   `iur_providers` (each provider's reliability at its own size) and, for
#'   the bootstrap, the `replicates` matrix; for the split-half method
#'   `iterations` (one row per split) and `excluded` (providers with fewer
#'   than two rows).
#' @references He K, Kalbfleisch JD, Yang Y, Li Y (2019). Inter-unit
#'   reliability for nonlinear models. *Statistics in Medicine*, 38(5),
#'   844-854.
#' @seealso [iur_curve()], [iur_stratified()]
#' @examples
#' d <- sim_kecc(n_providers = 30, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_logistic_fe(pd)
#' expected <- stats::median(coef(fit)$gamma) + fit$linear_predictor
#' expected <- stats::plogis(expected)
#' boot <- iur("bootstrap", obs = fit$data$y, expected = expected,
#'             provider = fit$data$provider, n_boot = 50, seed = 1)
#' boot
#' iur_curve(boot)
#' @name iur
NULL

#' @rdname iur
#' @export
iur <- function(method = c("direct", "bootstrap", "split_half"), estimate = NULL, se = NULL,
                size = NULL, obs = NULL, expected = NULL, provider = NULL, n_boot = 100,
                n_iter = 10, category_probs = c(0.1, 0.3, 0.7, 0.9), measure = NULL,
                seed = NULL) {
  method <- match.arg(method)
  if (method == "direct") {
    if (is.null(estimate) || is.null(se) || is.null(size)) {
      stop("method = \"direct\" needs `estimate`, `se` and `size`.", call. = FALSE)
    }
    return(iur_direct(estimate, se, size))
  }
  if (is.null(obs) || is.null(expected) || is.null(provider)) {
    stop("method = \"", method, "\" needs `obs`, `expected` and `provider`.", call. = FALSE)
  }
  if (method == "bootstrap") {
    return(iur_bootstrap(obs, expected, provider, n_boot = n_boot, measure = measure, seed = seed))
  }
  iur_split_half(obs, expected, provider, n_iter = n_iter, category_probs = category_probs,
                 measure = measure, seed = seed)
}

iur_direct <- function(estimate, se, size) {
  if (length(estimate) != length(se) || length(estimate) != length(size)) {
    stop("`estimate`, `se` and `size` must have one value per provider.", call. = FALSE)
  }
  dec <- iur_decomposition(as.numeric(size), as.numeric(se)^2, as.numeric(estimate))
  structure(c(list(method = "direct", providers = names(estimate), sizes = as.numeric(size),
                   measure = as.numeric(estimate), within_variance = as.numeric(se)^2), dec),
            class = "pp_iur")
}

iur_bootstrap <- function(obs, expected, provider, n_boot = 100, measure = NULL, seed = NULL) {
  d <- iur_rows(obs, expected, provider)
  fn <- if (is.null(measure)) ratio_measure else measure
  point <- fn(d$obs, d$expected, d$provider)
  if (length(point) != length(d$levels)) {
    stop("`measure` must return one value per provider.", call. = FALSE)
  }
  replicates <- with_seed(seed, {
    r <- matrix(NA_real_, length(d$levels), n_boot)
    for (b in seq_len(n_boot)) {
      idx <- stratified_resample(d$first, d$sizes)
      r[, b] <- fn(d$obs[idx], d$expected[idx], d$provider[idx])
    }
    r
  })
  dimnames(replicates) <- list(d$levels, NULL)
  within <- apply(replicates, 1, stats::var)
  dec <- iur_decomposition(d$sizes, within, point)
  structure(c(list(method = "bootstrap", providers = d$levels, sizes = d$sizes,
                   measure = stats::setNames(point, d$levels), within_variance = within,
                   replicates = replicates), dec),
            class = "pp_iur")
}

iur_split_half <- function(obs, expected, provider, n_iter = 10,
                           category_probs = c(0.1, 0.3, 0.7, 0.9), measure = NULL, seed = NULL) {
  d <- iur_rows(obs, expected, provider)
  fn <- if (is.null(measure)) ratio_measure else measure
  splittable <- d$sizes >= 2
  excluded <- d$levels[!splittable]
  if (length(excluded)) {
    warning(length(excluded), " of ", length(d$levels), " providers have fewer than two rows ",
            "and cannot be split; they are left out.", call. = FALSE)
    keep <- splittable[d$code]
    d <- iur_rows(d$obs[keep], d$expected[keep], d$provider[keep])
  }
  metrics <- c("kappa", "kendall_cat", "spearman_cat", "pearson", "kendall", "spearman")
  iterations <- with_seed(seed, {
    it <- matrix(NA_real_, n_iter, length(metrics), dimnames = list(NULL, metrics))
    for (k in seq_len(n_iter)) {
      half1 <- split_half_sample(d$first, d$sizes)
      half2 <- setdiff(seq_along(d$obs), half1)
      m1 <- fn(d$obs[half1], d$expected[half1], d$provider[half1])
      m2 <- fn(d$obs[half2], d$expected[half2], d$provider[half2])
      it[k, ] <- split_half_agreement(m1, m2, category_probs)
    }
    it
  })
  iur <- 2 * iterations / (1 + iterations)
  iur[is.nan(iur)] <- 0
  structure(list(method = "split_half", iur = colMeans(iur), iterations = iur,
                 correlations = iterations, n_providers = length(d$levels),
                 providers = d$levels, excluded = excluded),
            class = "pp_iur")
}

# Rows sorted by provider with the provider levels, sizes and first rows.
iur_rows <- function(obs, expected, provider) {
  if (length(obs) != length(expected) || length(obs) != length(provider)) {
    stop("`obs`, `expected` and `provider` must have the same length.", call. = FALSE)
  }
  provider <- as.character(provider)
  ord <- order(provider, method = "radix")
  provider <- provider[ord]
  levels <- unique(provider)
  code <- match(provider, levels)
  sizes <- as.numeric(tabulate(code, length(levels)))
  first <- c(1L, cumsum(sizes)[-length(sizes)] + 1L)
  list(obs = as.numeric(obs)[ord], expected = as.numeric(expected)[ord], provider = provider,
       levels = levels, code = code, sizes = sizes, first = as.integer(first))
}

#' @rdname iur
#' @export
ratio_measure <- function(obs, expected, provider) {
  provider <- as.character(provider)
  levels <- sort(unique(provider), method = "radix")
  code <- match(provider, levels)
  o <- as.numeric(rowsum(as.numeric(obs), code, reorder = TRUE))
  e <- as.numeric(rowsum(as.numeric(expected), code, reorder = TRUE))
  stats::setNames(ifelse(e != 0, o / e, NA_real_), levels)
}

stratified_resample <- function(first, sizes) {
  unlist(lapply(seq_along(sizes), function(i) {
    first[i] - 1L + sample.int(sizes[i], sizes[i], replace = TRUE)
  }), use.names = FALSE)
}

split_half_sample <- function(first, sizes) {
  unlist(lapply(seq_along(sizes), function(i) {
    first[i] - 1L + sort(sample.int(sizes[i], sizes[i] %/% 2))
  }), use.names = FALSE)
}

# The six agreement coefficients of two halves' measures.
split_half_agreement <- function(m1, m2, probs) {
  c1 <- categorise(m1, probs)
  c2 <- categorise(m2, probs)
  c(kappa = if (identical(c1, c2)) 1 else weighted_kappa(c1, c2),
    kendall_cat = stats::cor(c1, c2, method = "kendall"),
    spearman_cat = stats::cor(c1, c2, method = "spearman"),
    pearson = if (stats::sd(m1) == 0 || stats::sd(m2) == 0) 0 else stats::cor(m1, m2),
    kendall = stats::cor(m1, m2, method = "kendall"),
    spearman = stats::cor(m1, m2, method = "spearman"))
}

categorise <- function(values, probs) {
  breaks <- stats::quantile(values, probs, type = 7, names = FALSE)
  findInterval(values, breaks, left.open = TRUE) + 1L
}

# Cohen's kappa with linear weights over the sorted union of the labels.
weighted_kappa <- function(a, b) {
  labels <- sort(unique(c(a, b)))
  k <- length(labels)
  confusion <- table(factor(a, levels = labels), factor(b, levels = labels))
  w <- abs(outer(seq_len(k), seq_len(k), "-"))
  expected <- outer(rowSums(confusion), colSums(confusion)) / sum(confusion)
  1 - sum(w * confusion) / sum(w * expected)
}

# ANOVA-type decomposition shared by the direct and bootstrap estimators.
iur_decomposition <- function(sizes, within_variance, measure) {
  m <- length(sizes)
  total <- sum(sizes)
  weights <- pmax(sizes - 1, 0)
  s2_w <- if (sum(weights) > 0) sum(weights * within_variance) / sum(weights) else 0
  n_prime <- (total - sum(sizes^2) / total) / (m - 1)
  t_bar <- sum(sizes * measure) / total
  s2_t <- sum(sizes * (measure - t_bar)^2) / (n_prime * (m - 1))
  s2_b <- s2_t - s2_w
  scale <- if (m) max(abs(measure)) else 0
  if (!(s2_t > (64 * .Machine$double.eps * scale)^2)) {
    warning("The measures do not vary between providers; the reliability is undefined.",
            call. = FALSE)
    iur <- NA_real_
    iur_providers <- rep(NA_real_, m)
  } else {
    iur <- s2_b / s2_t
    denom <- s2_b + s2_w * n_prime / sizes
    iur_providers <- ifelse(denom != 0, s2_b / denom, 0)
  }
  list(iur = iur, n_providers = m, iur_providers = iur_providers, s2_between = s2_b,
       s2_within = s2_w * n_prime, n_prime = n_prime, s2_total = s2_t)
}

#' Reliability at representative provider sizes
#'
#' Evaluates the fitted reliability curve \eqn{s^2_b / (s^2_b + s2\_within /
#' n)} at the smallest size, the mean size of each quantile group of sizes
#' and the largest size.
#'
#' @param x a `pp_iur` object from [iur()].
#' @param sizes sizes to group; default the providers' sizes.
#' @param n_quantiles number of quantile groups. Default 10.
#' @return A data frame with `size_group` (`"min"`, `"quantile 1"`, ...,
#'   `"max"`), `size` and `iur`.
#' @examples
#' r <- iur("direct", estimate = c(0.2, 0.3, 0.25, 0.4), se = c(0.05, 0.04, 0.06, 0.03),
#'                 size = c(40, 60, 30, 90))
#' iur_curve(r, n_quantiles = 2)
#' @export
iur_curve <- function(x, sizes = NULL, n_quantiles = 10) {
  check_iur_components(x)
  if (is.null(sizes)) sizes <- x$sizes
  breaks <- unique(stats::quantile(sizes, seq(0, 1, length.out = n_quantiles + 1), type = 7,
                                   names = FALSE))
  inner <- breaks[-c(1, length(breaks))]
  labels <- findInterval(sizes, inner, left.open = TRUE) + 1L
  groups <- sort(unique(labels))
  means <- vapply(groups, function(g) mean(sizes[labels == g]), numeric(1))
  at <- c(min(sizes), means, max(sizes))
  denom <- x$s2_between + x$s2_within / at
  data.frame(size_group = c("min", paste("quantile", seq_along(groups)), "max"), size = at,
             iur = ifelse(denom != 0, x$s2_between / denom, 0), stringsAsFactors = FALSE)
}

#' Reliability within strata of providers
#'
#' Repeats the variance decomposition of [iur()]
#' within strata of providers.
#'
#' @param x a `pp_iur` object from [iur()].
#' @param by one stratification value per provider; default the sizes.
#' @param cuts right-closed cut points of `by`; default its deciles.
#' @return A data frame with one row for all providers and one per stratum:
#'   `stratum`, `n_providers`, `iur`.
#' @examples
#' r <- iur("direct", estimate = c(0.2, 0.3, 0.25, 0.4, 0.1, 0.35),
#'                 se = c(0.05, 0.04, 0.06, 0.03, 0.07, 0.03),
#'                 size = c(40, 60, 30, 90, 20, 80))
#' iur_stratified(r, cuts = 50)
#' @export
iur_stratified <- function(x, by = NULL, cuts = NULL) {
  check_iur_components(x)
  if (is.null(by)) by <- x$sizes
  if (length(by) != length(x$sizes)) stop("`by` must have one value per provider.", call. = FALSE)
  if (is.null(cuts)) cuts <- stats::quantile(by, seq(0.1, 0.9, by = 0.1), type = 7, names = FALSE)
  cuts <- sort(as.numeric(cuts))
  bounds <- c(-Inf, cuts, Inf)
  rows <- list(data.frame(stratum = "all", n_providers = x$n_providers, iur = x$iur,
                          stringsAsFactors = FALSE))
  for (i in seq_len(length(bounds) - 1)) {
    sel <- by > bounds[i] & by <= bounds[i + 1]
    if (!any(sel)) next
    lab <- if (i == 1) {
      paste0(format(min(by[sel])), " <= size <= ", format(cuts[1]))
    } else if (i == length(bounds) - 1) {
      paste0(format(cuts[length(cuts)]), " < size <= ", format(max(by[sel])))
    } else {
      paste0(format(cuts[i - 1]), " < size <= ", format(cuts[i]))
    }
    dec <- iur_decomposition(x$sizes[sel], x$within_variance[sel], x$measure[sel])
    rows[[length(rows) + 1]] <- data.frame(stratum = lab, n_providers = sum(sel), iur = dec$iur,
                                           stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

check_iur_components <- function(x) {
  if (!inherits(x, "pp_iur") || x$method == "split_half") {
    stop("`x` must come from iur() with the direct or bootstrap method.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @export
print.pp_iur <- function(x, digits = 4, ...) {
  if (x$method == "split_half") {
    cat("Split-half inter-unit reliability (", nrow(x$iterations), " splits, ", x$n_providers,
        " providers)\n", sep = "")
    print(round(x$iur, digits))
  } else {
    cat(if (x$method == "direct") "Direct" else "Bootstrap",
        " inter-unit reliability: ", format(x$iur, digits = digits), " (", x$n_providers,
        " providers)\n", sep = "")
    cat("  between variance ", format(x$s2_between, digits = digits), ", within variance ",
        format(x$s2_within, digits = digits), ", effective size ",
        format(x$n_prime, digits = digits), "\n", sep = "")
  }
  invisible(x)
}
