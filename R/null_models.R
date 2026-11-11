#' Null distributions for provider statistics
#'
#' A null model gives, for every provider, the mean and standard deviation of
#' the normal distribution its z statistic is referred to:
#' \eqn{z_{adj} = (z - \mathrm{null\_mean}) / \mathrm{null\_sd}}. Three kinds
#' exist: the theoretical null \eqn{N(0, 1)}, a fixed null with values given by
#' the caller (for instance a known scale inflation), and the empirical null
#' estimated robustly from the statistics themselves, overall or within
#' groups of providers of similar size. Any of them can be handed to
#' [calibrate()], [test()] or [test_measure()].
#'
#' @param kind `"empirical"` (default), `"theoretical"` or `"fixed"`.
#' @param mean,sd mean and standard deviation of a fixed null.
#' @param z for the empirical null, the provider statistics: a numeric
#'   vector, optionally named by provider, or a test result from [test()] or
#'   [test_measure()], whose statistic, provider identifiers and (for
#'   [test()]) sizes are used.
#' @param size one size per provider used to form the groups (patients, say);
#'   taken from a test result when `z` is one. Ignored when `groups` is given.
#' @param groups explicit group label per provider; missing labels get no
#'   parameters.
#' @param n_groups number of size groups. Default 4; 1 fits one overall null.
#' @param grouping `"rank"` (default; the providers sorted by size are dealt
#'   into `n_groups` consecutive blocks of `ceiling(n / n_groups)`, the last
#'   block holding the remainder, as the KECC production code for the PPPW and
#'   SFR measures and pprof_py do), `"quantile"` (groups by the empirical
#'   quantiles of `size`, a size equal to a break going to the lower group, as
#'   the production SMR and SHR code and `EmpiNull` do) or `"cut"` (the sorted
#'   providers cut by [cut()] into `n_groups` intervals of equal length, the
#'   `Group_EN` of the KECC measures paper).
#' @param estimator robust location-scale estimator of each group, all fitted
#'   with [MASS::rlm()] on an intercept-only model: `"huber"` (default;
#'   Huber's psi with the MAD scale, the production choice), `"bisquare"`
#'   (Tukey's biweight with the MAD scale) or `"mm"` (MM-estimation: an
#'   S-estimate of scale followed by bisquare steps).
#' @param tuning tuning constant of the psi function; `MASS` defaults 1.345
#'   (Huber) and 4.685 (bisquare).
#' @param maxit,acc iteration limit and convergence tolerance of [MASS::rlm()]
#'   (defaults 20 and 1e-4; `EmpiNull`'s robust calibration uses 1000 and
#'   1e-8 with the bisquare).
#' @param common_mean `NULL` (default; each group's own fitted location),
#'   `TRUE` (the mean of all eligible statistics) or a number: replaces the
#'   groups' locations, the scales are kept as fitted.
#' @param fit_mask logical, one per provider: the providers that enter the
#'   fits (`TRUE`); the others are still calibrated with their group's
#'   parameters.
#' @param min_group_size the smallest number of eligible statistics a group
#'   may have. Default 3.
#' @param small_group `"error"` (default) or `"theoretical"` (a group below
#'   `min_group_size` gets the theoretical null, with a warning).
#'
#' @return An object of class `pp_null`: a list with `kind`
#'   (`"theoretical"`, `"fixed"` or `"empirical"`), `mean`, `sd` and `group`
#'   (one per provider for the empirical null, named when the statistics
#'   were; single values otherwise), `diagnostics` (one row per group:
#'   `group`, `n_providers`, `n_fitted`, `min_size`, `median_size`,
#'   `max_size`, `mean_z`, `sd_z`, `null_mean`, `null_sd`, `converged`,
#'   `n_iter`, `note`) and `settings`.
#' @references
#' Efron B (2004). Large-scale simultaneous hypothesis testing: the choice
#' of a null hypothesis. *Journal of the American Statistical Association*,
#' 99(465), 96-104.
#'
#' Kalbfleisch JD, Wolfe RA (2013). On monitoring outcomes of medical
#' providers. *Statistics in Biosciences*, 5, 286-302.
#' @seealso [calibrate()], [adjust_z()]
#' @examples
#' d <- sim_kecc(n_providers = 40, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_logistic_fe(pd)
#' tst <- test(fit, statistic = "score", variance = "robust")
#' nul <- null_model("empirical", tst, n_groups = 2)
#' nul
#' head(calibrate(tst, nul))
#' @name null_models
NULL

#' @rdname null_models
#' @export
null_model <- function(kind = c("empirical", "theoretical", "fixed"), z = NULL, size = NULL,
                       groups = NULL, n_groups = 4, grouping = c("rank", "quantile", "cut"),
                       estimator = c("huber", "bisquare", "mm"), tuning = NULL, maxit = 20,
                       acc = 1e-4, common_mean = NULL, fit_mask = NULL, min_group_size = 3,
                       small_group = c("error", "theoretical"), mean = 0, sd = 1) {
  kind <- match.arg(kind)
  if (kind == "theoretical") return(null_theoretical())
  if (kind == "fixed") return(null_fixed(mean, sd))
  if (is.null(z)) stop("The empirical null needs the statistics `z`.", call. = FALSE)
  null_empirical(z, size = size, groups = groups, n_groups = n_groups, grouping = grouping,
                 estimator = estimator, tuning = tuning, maxit = maxit, acc = acc,
                 common_mean = common_mean, fit_mask = fit_mask,
                 min_group_size = min_group_size, small_group = small_group)
}

null_theoretical <- function() {
  new_pp_null("theoretical", mean = 0, sd = 1, group = NA_integer_)
}

null_fixed <- function(mean = 0, sd = 1) {
  if (!is.numeric(mean) || !is.numeric(sd) || length(mean) != 1 || length(sd) != 1 ||
        !is.finite(mean) || !is.finite(sd) || sd <= 0) {
    stop("`mean` must be a finite number and `sd` a positive finite number.", call. = FALSE)
  }
  new_pp_null("fixed", mean = mean, sd = sd, group = NA_integer_)
}

null_empirical <- function(z, size = NULL, groups = NULL, n_groups = 4,
                           grouping = c("rank", "quantile", "cut"),
                           estimator = c("huber", "bisquare", "mm"), tuning = NULL, maxit = 20,
                           acc = 1e-4, common_mean = NULL, fit_mask = NULL, min_group_size = 3,
                           small_group = c("error", "theoretical")) {
  grouping <- match.arg(grouping)
  estimator <- match.arg(estimator)
  small_group <- match.arg(small_group)
  if (!is.numeric(min_group_size) || length(min_group_size) != 1 || min_group_size < 2) {
    stop("`min_group_size` must be at least 2.", call. = FALSE)
  }
  zs <- statistics_of(z)
  stat <- zs$z
  providers <- zs$providers
  if (is.null(size) && is.null(groups)) size <- zs$size
  n <- length(stat)
  if (!is.null(groups)) {
    if (length(groups) != n) stop("`groups` must have one label per provider.", call. = FALSE)
    labels <- groups
  } else if (!is.null(size) && n_groups > 1) {
    if (length(size) != n) stop("`size` must have one value per provider.", call. = FALSE)
    labels <- assign_groups(size, n_groups, grouping)
  } else {
    labels <- rep(1L, n)
  }
  if (is.null(fit_mask)) fit_mask <- rep(TRUE, n)
  if (!is.logical(fit_mask) || length(fit_mask) != n) {
    stop("`fit_mask` must be logical with one value per provider.", call. = FALSE)
  }
  eligible <- fit_mask & is.finite(stat)
  levels <- sort(unique(labels[!is.na(labels)]))
  null_mean <- rep(NA_real_, n)
  null_sd <- rep(NA_real_, n)
  rows <- vector("list", length(levels))
  problems <- character(0)
  for (k in seq_along(levels)) {
    g <- levels[k]
    members <- !is.na(labels) & labels == g
    values <- stat[members & eligible]
    sizes_g <- if (is.null(size)) rep(NA_real_, sum(members)) else size[members]
    row <- data.frame(group = g, n_providers = sum(members), n_fitted = length(values),
                      min_size = suppressWarnings(min(sizes_g, na.rm = TRUE)),
                      median_size = stats::median(sizes_g, na.rm = TRUE),
                      max_size = suppressWarnings(max(sizes_g, na.rm = TRUE)),
                      mean_z = if (length(values)) mean(values) else NA_real_,
                      sd_z = if (length(values) > 1) stats::sd(values) else NA_real_,
                      null_mean = NA_real_, null_sd = NA_real_, converged = NA,
                      n_iter = NA_integer_, note = "", stringsAsFactors = FALSE)
    if (length(values) < min_group_size) {
      if (small_group == "error") {
        stop("Group ", g, " has ", length(values), " eligible statistics; at least ",
             min_group_size, " are needed.", call. = FALSE)
      }
      row$null_mean <- 0
      row$null_sd <- 1
      row$note <- "fewer eligible statistics than min_group_size; theoretical null used"
      problems <- c(problems, paste0("group ", g, ": ", row$note))
    } else {
      fit <- robust_location_scale(values, estimator, tuning, maxit, acc)
      if (!is.finite(fit$location) || !is.finite(fit$scale) || fit$scale <= 0) {
        stop("Group ", g, ": the fitted location or scale is not finite and positive.",
             call. = FALSE)
      }
      row$null_mean <- fit$location
      row$null_sd <- fit$scale
      row$converged <- fit$converged
      row$n_iter <- fit$n_iter
      if (!fit$converged) {
        row$note <- "did not converge"
        problems <- c(problems, paste0("group ", g, ": the fit did not converge in ", fit$n_iter,
                                       " iterations"))
      }
    }
    rows[[k]] <- row
  }
  diagnostics <- do.call(rbind, rows)
  rownames(diagnostics) <- NULL
  if (!is.null(common_mean) && !identical(common_mean, FALSE)) {
    pooled <- if (isTRUE(common_mean)) {
      if (!any(eligible)) stop("No eligible statistic for a pooled mean.", call. = FALSE)
      mean(stat[eligible])
    } else {
      as.numeric(common_mean)
    }
    diagnostics$null_mean <- pooled
  }
  for (k in seq_len(nrow(diagnostics))) {
    members <- !is.na(labels) & labels == diagnostics$group[k]
    null_mean[members] <- diagnostics$null_mean[k]
    null_sd[members] <- diagnostics$null_sd[k]
  }
  if (anyNA(labels)) {
    problems <- c(problems, paste0(sum(is.na(labels)),
                                   " provider(s) have no group and are not calibrated"))
  }
  if (length(problems)) {
    warning("Empirical null: ", paste(problems, collapse = "; "), ".", call. = FALSE)
  }
  names(null_mean) <- names(null_sd) <- names(labels) <- providers
  new_pp_null("empirical", mean = null_mean, sd = null_sd, group = labels,
              diagnostics = diagnostics,
              settings = list(n_groups = if (is.null(groups)) n_groups else NA,
                              grouping = if (is.null(groups)) grouping else "given",
                              estimator = estimator, tuning = tuning, maxit = maxit, acc = acc,
                              common_mean = common_mean, min_group_size = min_group_size),
              z = stats::setNames(stat, providers))
}

new_pp_null <- function(kind, mean, sd, group, diagnostics = NULL, settings = list(), z = NULL) {
  structure(list(kind = kind, mean = mean, sd = sd, group = group, diagnostics = diagnostics,
                 settings = settings, z = z),
            class = "pp_null")
}

# Statistic, provider identifiers and sizes of a test result or a vector.
statistics_of <- function(z) {
  if (inherits(z, "pp_test")) {
    sizes <- attr(z, "sizes")
    size <- if (is.null(sizes)) {
      NULL
    } else if (all(is.na(sizes$patients))) {
      sizes$rows
    } else {
      sizes$patients
    }
    return(list(z = z$statistic, providers = z$provider, size = size))
  }
  if (inherits(z, "pp_measure_test")) {
    return(list(z = z$z, providers = if ("provider" %in% names(z)) z$provider else NULL,
                size = NULL))
  }
  if (!is.numeric(z)) stop("`z` must be a numeric vector or a test result.", call. = FALSE)
  list(z = as.numeric(z), providers = names(z), size = NULL)
}

#' Group providers by size
#'
#' @param size one size per provider.
#' @param n_groups number of groups.
#' @param grouping `"rank"`, `"quantile"` or `"cut"`, see [null_model()].
#' @return Integer group labels `1, ..., n_groups`.
#' @examples
#' assign_groups(c(5, 50, 12, 30, 8, 100), n_groups = 3)
#' @export
assign_groups <- function(size, n_groups = 4, grouping = c("rank", "quantile", "cut")) {
  grouping <- match.arg(grouping)
  if (!is.numeric(n_groups) || length(n_groups) != 1 || n_groups < 1 ||
        n_groups != round(n_groups)) {
    stop("`n_groups` must be a positive whole number.", call. = FALSE)
  }
  size <- as.numeric(size)
  n <- length(size)
  if (n == 0) return(integer(0))
  if (n_groups == 1) return(rep(1L, n))
  if (grouping == "rank") {
    labels <- integer(n)
    block <- ceiling(n / n_groups)
    labels[order(size)] <- pmin((seq_len(n) - 1L) %/% block, n_groups - 1L) + 1L
    return(labels)
  }
  if (grouping == "cut") {
    labels <- integer(n)
    labels[order(size)] <- cut(seq_len(n), breaks = n_groups, labels = FALSE)
    return(labels)
  }
  if (!all(is.finite(size))) stop("grouping = \"quantile\" needs finite sizes.", call. = FALSE)
  probs <- seq(0, 1, length.out = n_groups + 1)[-c(1, n_groups + 1)]
  breaks <- stats::quantile(size, probs, type = 7, names = FALSE)
  if (anyDuplicated(breaks)) {
    stop("The quantile breaks are not distinct; use grouping = \"rank\" or fewer groups.",
         call. = FALSE)
  }
  findInterval(size, breaks, left.open = TRUE) + 1L
}

# Intercept-only MASS::rlm fit: location, scale, convergence.
robust_location_scale <- function(values, estimator, tuning, maxit, acc) {
  args <- list(x = matrix(1, length(values), 1), y = values,
               method = if (estimator == "mm") "MM" else "M",
               psi = if (estimator == "huber") MASS::psi.huber else MASS::psi.bisquare,
               scale.est = "MAD", maxit = maxit, acc = acc)
  if (!is.null(tuning)) {
    args[[if (estimator == "huber") "k" else "c"]] <- tuning
  }
  fit <- do.call(MASS::rlm, args)
  list(location = unname(stats::coef(fit)[1]), scale = unname(fit$s),
       converged = isTRUE(fit$converged), n_iter = length(fit$conv))
}

# Null parameters aligned with `providers` (NULL = positional).
null_parameters <- function(null, providers = NULL, n = length(providers)) {
  if (!inherits(null, "pp_null")) stop("`null` must be a pp_null object.", call. = FALSE)
  if (null$kind != "empirical") {
    return(list(mean = rep(null$mean, n), sd = rep(null$sd, n), group = rep(NA_integer_, n)))
  }
  if (!is.null(providers) && !is.null(names(null$mean))) {
    pos <- match(as.character(providers), names(null$mean))
    if (anyNA(pos)) {
      stop("The null was not fitted for provider(s): ",
           paste(providers[is.na(pos)], collapse = ", "), call. = FALSE)
    }
  } else {
    if (length(null$mean) != n) {
      stop("The null has ", length(null$mean), " providers but ", n, " statistics were given.",
           call. = FALSE)
    }
    pos <- seq_len(n)
  }
  list(mean = unname(null$mean[pos]), sd = unname(null$sd[pos]), group = unname(null$group[pos]))
}

# A null model from a pp_null object or a function of the test result.
resolve_null <- function(null_model, x) {
  if (is.character(null_model) && length(null_model) == 1) {
    null_model <- pprof::null_model(null_model, z = x)
  }
  if (is.function(null_model)) null_model <- null_model(x)
  if (!inherits(null_model, "pp_null")) {
    stop("`null_model` must be a pp_null object or a function returning one.", call. = FALSE)
  }
  null_model
}

#' @export
print.pp_null <- function(x, digits = 4, ...) {
  if (x$kind == "theoretical") {
    cat("Theoretical null: N(0, 1)\n")
  } else if (x$kind == "fixed") {
    cat("Fixed null: N(", format(x$mean, digits = digits), ", ", format(x$sd, digits = digits),
        "^2)\n", sep = "")
  } else {
    s <- x$settings
    cat("Empirical null: ", nrow(x$diagnostics), " group(s) by ", s$grouping, ", ", s$estimator,
        " estimator\n", sep = "")
    print(x$diagnostics[, c("group", "n_providers", "n_fitted", "min_size", "max_size",
                            "null_mean", "null_sd", "converged")],
          digits = digits, row.names = FALSE)
  }
  invisible(x)
}

#' Adjust statistics with a null model
#'
#' @param z statistics (a numeric vector, optionally named by provider).
#' @param null a `pp_null` object.
#' @return \eqn{(z - \mathrm{null\_mean}) / \mathrm{null\_sd}}.
#' @examples
#' adjust_z(c(a = 2, b = -1), null_model("fixed", sd = 1.5))
#' @export
adjust_z <- function(z, null) {
  p <- null_parameters(null, names(z), length(z))
  stats::setNames((as.numeric(z) - p$mean) / p$sd, names(z))
}

#' Calibrate a test result with a null model
#'
#' Refers the statistics of a [test()] or [test_measure()] result
#' to a null model: the adjusted statistic \eqn{(z - \mathrm{null\_mean}) /
#' \mathrm{null\_sd}} replaces the standard normal in the p-values and
#' flags, and the confidence limits are those of the inverted calibrated test,
#' \eqn{\hat\theta - (\mathrm{null\_mean} \pm c\, \mathrm{null\_sd})\, se},
#' so that a limit excludes the null value exactly when the provider is
#' flagged.
#'
#' @param x a test result.
#' @param null_model a `pp_null` object, or a function of `x` returning one
#'   (for instance `function(t) null_model("empirical", t, n_groups = 4)`).
#' @param level,alternative default to the settings of `x`.
#' @param critical optional critical value for the flags and limits instead
#'   of the normal quantile for `level`.
#' @param ci_location location at which the confidence limits are centred:
#'   `NULL` (default) uses the null mean, so that a limit excludes the null
#'   value exactly when the provider is flagged; a number or one value per
#'   provider centres the limits elsewhere, for instance `0` for limits of the
#'   form \eqn{\hat\theta \pm c\, \mathrm{null\_sd}\, se} that ignore the
#'   null location.
#'
#' @details For a Wald-type test the limits are
#' \eqn{\hat\theta - (\mathrm{null\_mean} \pm c\, \mathrm{null\_sd})\, se}; for a
#' modified score test of [test.logistic_fe()] they invert the calibrated
#' statistic (the set of null values at which the adjusted statistic lies
#' within \eqn{\pm c}), for a mid-p test the calibrated mid-p rule; the exact
#' Poisson test cannot be calibrated.
#' @return `x` with the columns `z_adjusted`, `null_mean`, `null_sd` and
#'   `null_group` added, `p` and `flag` recomputed from the adjusted statistic,
#'   and `lower`/`upper` recomputed where a standard error is available (for
#'   [test_measure()] results on the transformed scale, mapped back; for the
#'   log-ratio Wald and mid-p tests of [test.cox()] on the ratio
#'   scale, the mid-p limits by root finding); the settings are recorded in
#'   the `test` attribute. The exact Poisson test cannot be calibrated.
#' @examples
#' d <- sim_kecc(n_providers = 40, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' tst <- test(fit_logistic_fe(pd), statistic = "wald", variance = "robust")
#' calibrate(tst, null_model("fixed", sd = 1.5))
#' @export
calibrate <- function(x, null_model, level = NULL, alternative = NULL, critical = NULL,
                      ci_location = NULL) {
  UseMethod("calibrate")
}

#' @export
calibrate.pp_test <- function(x, null_model, level = NULL, alternative = NULL, critical = NULL,
                              ci_location = NULL) {
  null <- resolve_null(null_model, x)
  settings <- attr(x, "test")
  if (is.null(level)) level <- settings$level
  if (is.null(alternative)) alternative <- settings$alternative
  check_level(level)
  alternative <- match_alternative(alternative)
  limits_kind <- if (is.null(settings$limits)) "wald" else settings$limits
  if (limits_kind == "exact") {
    stop("The exact Poisson test uses the theoretical null; calibrate a mid-p test instead.",
         call. = FALSE)
  }
  p <- null_parameters(null, x$provider)
  z_adj <- (x$statistic - p$mean) / p$sd
  decision <- null_decision(z_adj, level, alternative, critical)
  q <- if (is.null(critical)) normal_quantile(level, alternative) else critical
  loc <- ci_centre(ci_location, p$mean, nrow(x))
  se <- x$se
  inversion <- attr(x, "inversion")
  if (limits_kind == "midp") {
    lim <- midp_limits(x$observed, x$expected, loc, p$sd, q, alternative)
    lower <- lim$lower
    upper <- lim$upper
  } else if (identical(settings$statistic, "score") && !is.null(inversion) &&
               any(is.finite(x$lower) | is.finite(x$upper))) {
    # the calibrated score test inverted: null values where the adjusted statistic is within +/- q
    idx <- provider_positions(x$provider, inversion$fit$sizes$provider)
    tg <- inversion_targets(loc, p$sd, q, alternative)
    lim <- score_limits(inversion$fit, inversion$variance, idx, tg$z_lower, tg$z_upper,
                        inversion$pooled)
    lower <- lim$lower
    upper <- lim$upper
  } else {
    centre <- if (limits_kind == "log_wald") log(x$estimate) else x$estimate
    lower <- centre - (loc + q * p$sd) * se
    upper <- centre - (loc - q * p$sd) * se
    if (alternative == "greater") upper[!is.na(se)] <- Inf
    if (alternative == "less") lower[!is.na(se)] <- -Inf
    if (limits_kind == "log_wald") {
      lower <- exp(lower)
      upper <- exp(upper)
    }
  }
  attrs <- attributes(x)
  out <- as.data.frame(x)
  out$z_adjusted <- z_adj
  out$null_mean <- p$mean
  out$null_sd <- p$sd
  out$null_group <- p$group
  out$p <- decision$p
  out$flag <- decision$flag
  out$lower <- lower
  out$upper <- upper
  fixed <- c("provider", "estimate", "se", "null_value", "statistic", "z_adjusted", "null_mean",
             "null_sd", "null_group", "p", "flag", "lower", "upper")
  out <- out[, c(fixed, setdiff(names(out), fixed))]
  settings$level <- level
  settings$alternative <- alternative
  settings$critical <- critical
  settings$ci_location <- ci_location
  settings$null_model <- describe_null(null)
  attr(out, "test") <- settings
  attr(out, "sizes") <- attrs$sizes
  if (!is.null(inversion)) attr(out, "inversion") <- inversion
  class(out) <- c("pp_test", "data.frame")
  out
}

#' @export
calibrate.pp_measure_test <- function(x, null_model, level = NULL, alternative = NULL,
                                      critical = NULL, ci_location = NULL) {
  null <- resolve_null(null_model, x)
  settings <- attr(x, "test")
  if (is.null(level)) level <- settings$level
  if (is.null(alternative)) alternative <- settings$alternative
  check_level(level)
  alternative <- match_alternative(alternative)
  providers <- if ("provider" %in% names(x)) x$provider else NULL
  p <- null_parameters(null, providers, nrow(x))
  z_adj <- (x$z - p$mean) / p$sd
  decision <- null_decision(z_adj, level, alternative, critical)
  q <- if (is.null(critical)) normal_quantile(level, alternative) else critical
  loc <- ci_centre(ci_location, p$mean, nrow(x))
  back <- switch(settings$transform, identity = function(v) v, logit = stats::plogis, log = exp)
  lower_t <- x$estimate_t - (loc + q * p$sd) * x$se_t
  upper_t <- x$estimate_t - (loc - q * p$sd) * x$se_t
  if (alternative == "greater") upper_t <- rep(Inf, nrow(x))
  if (alternative == "less") lower_t <- rep(-Inf, nrow(x))
  out <- as.data.frame(x)
  out$z_adjusted <- z_adj
  out$null_mean <- p$mean
  out$null_sd <- p$sd
  out$null_group <- p$group
  out$p <- decision$p
  out$flag <- decision$flag
  out$lower_t <- lower_t
  out$upper_t <- upper_t
  out$lower <- back(lower_t)
  out$upper <- back(upper_t)
  settings$level <- level
  settings$alternative <- alternative
  settings$critical <- critical
  settings$ci_location <- ci_location
  settings$null_model <- describe_null(null)
  attr(out, "test") <- settings
  class(out) <- c("pp_measure_test", "data.frame")
  out
}

# Centre of the calibrated confidence limits: the null mean unless given.
ci_centre <- function(ci_location, null_mean, n) {
  if (is.null(ci_location)) return(null_mean)
  if (!length(ci_location) %in% c(1L, n)) {
    stop("`ci_location` must be one number or one value per provider.", call. = FALSE)
  }
  rep_len(as.numeric(ci_location), n)
}

# p-values and flags of adjusted statistics, by level or by a critical value.
null_decision <- function(z_adj, level, alternative, critical) {
  decision <- z_decision(z_adj, level, alternative)
  if (!is.null(critical)) {
    up <- z_adj > critical
    down <- z_adj < -critical
    if (alternative == "greater") down <- rep(FALSE, length(z_adj))
    if (alternative == "less") up <- rep(FALSE, length(z_adj))
    flag <- ifelse(up, 1L, ifelse(down, -1L, 0L))
    flag[!is.finite(z_adj)] <- NA_integer_
    decision$flag <- flag
  }
  decision
}

describe_null <- function(null) {
  if (null$kind == "empirical") {
    c(list(kind = "empirical"), null$settings, list(groups = null$diagnostics))
  } else {
    list(kind = null$kind, mean = null$mean, sd = null$sd)
  }
}
