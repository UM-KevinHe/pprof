#' Plots of profiling results
#'
#' `plot()` methods for the objects of the package, drawn with ggplot2 and
#' returned as ggplot objects to be printed, modified or saved:
#' * a [test()] result: `type = "caterpillar"` (every provider's estimate
#'   and interval in order of the estimate, coloured by flag),
#'   `"funnel"` (estimates against precision, with the control limits of
#'   the same test, so that a provider lies outside the limits exactly when
#'   it is flagged: on the precision axis \eqn{1 / se^2} for a test with
#'   standard errors, with limits \eqn{\theta_0 \pm c\, se}; on the
#'   precision axis \eqn{E^2 / V} for a score test, the observed over
#'   expected events with limits \eqn{1 + (m \pm c\, s) / \sqrt{E^2 / V}}
#'   where \eqn{V} is the variance of that test (naive, robust or
#'   stabilized) and \eqn{(m, s)} its null; on the expected-count axis for
#'   the Poisson tests of [test.cox()], with the counts the test flags at
#'   every expected count; one curve per group when the test was calibrated
#'   with an empirical null) or `"bar"` (the shares of providers flagged
#'   above, below and not different, by provider-size group);
#' * a [standardize()] result: `type = "observed_expected"` (observed
#'   against expected with the diagonal, coloured by the flags of a test when
#'   one is given) or `"funnel"` (the ratio against its precision
#'   \eqn{E^2 / V}; without a test, the normal limits of the naive score
#'   statistic and no colours; with `test =` a score test, that test's own
#'   funnel, so that the limits and the colours come from the same test. A
#'   Wald-type test is refused here: draw it with `plot(test, type =
#'   "funnel")`);
#' * a [summary()] result: a forest plot of the coefficients;
#' * an empirical [null_model()]: the statistics against the theoretical and
#'   the fitted null, one panel per group;
#' * an [iur()] result: the reliability curve against provider size;
#' * a [var_select()] path: the coefficient paths; its [cv()]: the
#'   cross-validation curve.
#'
#' @param x the object.
#' @param type the display, see above.
#' @param highlight provider identifiers to label.
#' @param levels for the funnel plot of a test with standard errors, the
#'   levels of the limit curves; the first is the test's own level, the
#'   others are drawn dashed.
#' @param n_groups for `type = "bar"`, the number of provider-size groups
#'   (quartiles by default).
#' @param test for a [standardize()] result, a [test()] result whose flags
#'   colour the providers.
#' @param exponentiate for a [summary()] result, show `exp(estimate)` (odds
#'   or hazard ratios) on a log axis.
#' @param ... unused.
#'
#' @return A ggplot object.
#' @examples
#' d <- sim_kecc(n_providers = 40, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_logistic_fe(pd)
#' tst <- test(fit, statistic = "wald", variance = "robust")
#' plot(tst)
#' plot(tst, type = "funnel")
#' plot(tst, type = "bar")
#' plot(standardize(fit), test = tst)
#' plot(summary(fit))
#' plot(null_model("empirical", tst, n_groups = 2))
#' @name plots
NULL

#' @rdname plots
#' @export
plot.pp_test <- function(x, type = c("caterpillar", "funnel", "bar"), highlight = NULL,
                         levels = c(0.95, 0.998), n_groups = 4, ...) {
  type <- match.arg(type)
  switch(type,
         caterpillar = caterpillar_plot(x, highlight),
         funnel = funnel_plot_test(x, highlight, levels),
         bar = bar_plot(x, n_groups))
}

#' @rdname plots
#' @export
plot.pp_measure <- function(x, type = c("observed_expected", "funnel"), test = NULL,
                            highlight = NULL, levels = c(0.95, 0.998), ...) {
  type <- match.arg(type)
  if (attr(x, "type") != "indirect") {
    stop("The plots of a standardized measure need the indirect standardization.", call. = FALSE)
  }
  status <- if (is.null(test)) {
    factor(rep("not tested", nrow(x)), levels = status_levels)
  } else {
    if (!inherits(test, "pp_test")) stop("`test` must be a test() result.", call. = FALSE)
    status_of(test$flag[match(x$provider, test$provider)])
  }
  d <- data.frame(provider = x$provider, observed = x$observed, expected = x$expected,
                  ratio = x$observed / x$expected, status = status, stringsAsFactors = FALSE)
  if (type == "observed_expected") {
    top <- max(c(d$observed, d$expected), na.rm = TRUE)
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$expected, y = .data$observed)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "dashed") +
      ggplot2::geom_point(ggplot2::aes(colour = .data$status, shape = .data$status), size = 2) +
      status_scales() +
      ggplot2::coord_cartesian(xlim = c(0, top), ylim = c(0, top)) +
      ggplot2::labs(x = "Expected", y = "Observed", colour = NULL, shape = NULL) +
      theme_pp()
    return(label_points(p, d, highlight, "expected", "observed"))
  }
  if (!is.null(test)) {
    if (!all(c("observed", "expected", "variance") %in% names(test)) ||
          !identical(attr(test, "test")$statistic, "score")) {
      stop("The limits of a funnel plot must come from the test that colours it: pass a score ",
           "test, or draw a Wald-type test with plot(test, type = \"funnel\").", call. = FALSE)
    }
    return(funnel_plot_test(test, highlight, levels))
  }
  precision <- if ("variance_expected" %in% names(x) && any(is.finite(x$variance_expected))) {
    x$expected^2 / x$variance_expected
  } else {
    x$expected
  }
  d$precision <- precision
  grid <- exp(seq(log(max(min(precision[precision > 0]), 1e-8)), log(max(precision)),
                  length.out = 200))
  curves <- do.call(rbind, lapply(levels, function(l) {
    q <- stats::qnorm(1 - (1 - l) / 2)
    data.frame(x = grid, lower = pmax(1 - q / sqrt(grid), 0), upper = 1 + q / sqrt(grid),
               level = l)
  }))
  funnel_draw(d, "precision", "ratio", 1, curves, highlight,
              if (identical(precision, x$expected)) "Expected events" else "Precision, E\u00b2 / V",
              "Observed / expected")
}

#' @rdname plots
#' @export
plot.pp_summary <- function(x, exponentiate = FALSE, ...) {
  d <- data.frame(covariate = factor(x$covariate, levels = rev(x$covariate)),
                  estimate = x$estimate, lower = x$lower, upper = x$upper,
                  stringsAsFactors = FALSE)
  if (is.null(d$lower)) d$lower <- d$upper <- NA_real_
  null_line <- 0
  if (exponentiate) {
    d$estimate <- exp(d$estimate)
    d$lower <- exp(d$lower)
    d$upper <- exp(d$upper)
    null_line <- 1
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(y = .data$covariate, x = .data$estimate)) +
    ggplot2::geom_vline(xintercept = null_line, colour = "grey55", linetype = "dashed") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$lower, xmax = .data$upper), height = 0.2,
                            colour = "#4D4D4D", na.rm = TRUE) +
    ggplot2::geom_point(size = 2.2, colour = "#0072B2") +
    ggplot2::labs(x = if (exponentiate) "Ratio (log scale)" else "Coefficient", y = NULL) +
    theme_pp()
  if (exponentiate) p <- p + ggplot2::scale_x_log10()
  p
}

#' @rdname plots
#' @export
plot.pp_null <- function(x, ...) {
  if (x$kind != "empirical" || is.null(x$z)) {
    stop("Only an empirical null fitted from statistics can be plotted.", call. = FALSE)
  }
  d <- data.frame(z = as.numeric(x$z), group = x$group, null_mean = unname(x$mean),
                  null_sd = unname(x$sd))
  d <- d[is.finite(d$z) & !is.na(d$group), , drop = FALSE]
  d$group <- factor(d$group)
  rng <- range(c(d$z, -4, 4))
  grid <- seq(rng[1], rng[2], length.out = 200)
  params <- unique(d[, c("group", "null_mean", "null_sd")])
  curves <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
    rbind(data.frame(group = params$group[i], z = grid, density = stats::dnorm(grid),
                     curve = "theoretical N(0, 1)"),
          data.frame(group = params$group[i], z = grid,
                     density = stats::dnorm(grid, params$null_mean[i], params$null_sd[i]),
                     curve = "fitted empirical null"))
  }))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$z)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30,
                            fill = "grey80", colour = "white") +
    ggplot2::geom_line(data = curves,
                       ggplot2::aes(y = .data$density, colour = .data$curve,
                                    linetype = .data$curve), linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = c("theoretical N(0, 1)" = "grey40",
                                            "fitted empirical null" = "#D55E00")) +
    ggplot2::scale_linetype_manual(values = c("theoretical N(0, 1)" = "dashed",
                                              "fitted empirical null" = "solid")) +
    ggplot2::facet_wrap(~ group, labeller = ggplot2::labeller(group = function(g) {
      paste("group", g)
    })) +
    ggplot2::labs(x = "Statistic", y = "Density", colour = NULL, linetype = NULL) +
    theme_pp()
}

#' @rdname plots
#' @export
plot.pp_iur <- function(x, ...) {
  if (x$method == "split_half") {
    d <- data.frame(coefficient = factor(names(x$iur), levels = names(x$iur)), iur = x$iur)
    return(ggplot2::ggplot(d, ggplot2::aes(x = .data$coefficient, y = .data$iur)) +
             ggplot2::geom_col(fill = "#0072B2", width = 0.6) +
             ggplot2::coord_cartesian(ylim = c(0, 1)) +
             ggplot2::labs(x = NULL, y = "Reliability (Spearman-Brown)") +
             theme_pp())
  }
  curve <- iur_curve(x)
  sizes <- x$sizes
  grid <- seq(min(sizes), max(sizes), length.out = 200)
  smooth <- data.frame(size = grid, iur = x$s2_between / (x$s2_between + x$s2_within / grid))
  points <- curve[!curve$size_group %in% c("min", "max"), , drop = FALSE]
  ggplot2::ggplot(smooth, ggplot2::aes(x = .data$size, y = .data$iur)) +
    ggplot2::geom_line(colour = "#0072B2", linewidth = 0.9) +
    ggplot2::geom_point(data = points, colour = "#D55E00", size = 2) +
    ggplot2::geom_hline(yintercept = x$iur, colour = "grey55", linetype = "dashed") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Provider size", y = "Reliability") +
    theme_pp()
}

#' @rdname plots
#' @export
plot.pp_lasso <- function(x, ...) {
  beta <- x$beta
  d <- data.frame(lambda = rep(x$lambda, each = nrow(beta)),
                  covariate = rep(rownames(beta), ncol(beta)), beta = as.numeric(beta),
                  penalized = rep(unname(x$penalized), ncol(beta)), stringsAsFactors = FALSE)
  ends <- d[d$lambda == min(d$lambda), , drop = FALSE]
  ggplot2::ggplot(d, ggplot2::aes(x = .data$lambda, y = .data$beta, group = .data$covariate,
                                  colour = .data$covariate,
                                  linetype = ifelse(.data$penalized, "penalized", "unpenalized"))) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey75") +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_text(data = ends, ggplot2::aes(label = .data$covariate), hjust = 1.1, size = 3,
                       show.legend = FALSE) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_linetype_manual(values = c(penalized = "solid", unpenalized = "dashed")) +
    ggplot2::labs(x = expression(lambda), y = "Coefficient", colour = NULL, linetype = NULL) +
    theme_pp() +
    ggplot2::theme(legend.position = "none")
}

#' @rdname plots
#' @export
plot.pp_cv <- function(x, ...) {
  d <- data.frame(lambda = x$lambda, cve = x$cve, cvse = x$cvse)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$lambda, y = .data$cve)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$cve - .data$cvse,
                                        ymax = .data$cve + .data$cvse),
                           width = 0, colour = "grey65") +
    ggplot2::geom_point(colour = "#0072B2", size = 1.8) +
    ggplot2::geom_vline(xintercept = x$lambda_min, colour = "#D55E00", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = x$lambda_1se, colour = "#D55E00", linetype = "dotted") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = expression(lambda), y = "Cross-validated deviance") +
    theme_pp()
}

#' Flags of several tests side by side
#'
#' Collects the flags of several [test()] results (scenarios: other
#' statistics, variances, references, standardizations) for the providers
#' flagged in at least one of them, the "flags under alternatives" of a
#' measure report.
#'
#' @param ... [test()] results, named by scenario.
#' @param all_providers keep every provider rather than those flagged at
#'   least once. Default `FALSE`.
#' @return A data frame of class `pp_stability`: `provider`, one column of
#'   flags per scenario, and `n_flagged`. Its `plot()` method draws the
#'   providers against the scenarios.
#' @examples
#' d <- sim_kecc(n_providers = 40, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_logistic_fe(pd)
#' st <- flag_stability(score = test(fit, statistic = "score"),
#'                      robust = test(fit, statistic = "score", variance = "robust"),
#'                      empirical = test(fit, statistic = "score", variance = "robust",
#'                                       null_model = "empirical"))
#' st
#' plot(st)
#' @export
flag_stability <- function(..., all_providers = FALSE) {
  tests <- list(...)
  if (length(tests) < 1 || !all(vapply(tests, inherits, logical(1), "pp_test"))) {
    stop("Give one or more test() results.", call. = FALSE)
  }
  if (is.null(names(tests)) || any(names(tests) == "")) {
    names(tests) <- paste0("scenario_", seq_along(tests))
  }
  providers <- unique(unlist(lapply(tests, function(t) t$provider)))
  flags <- vapply(tests, function(t) t$flag[match(providers, t$provider)],
                  integer(length(providers)))
  flags <- matrix(flags, nrow = length(providers), dimnames = list(NULL, names(tests)))
  n_flagged <- rowSums(flags != 0, na.rm = TRUE)
  keep <- if (all_providers) rep(TRUE, length(providers)) else n_flagged > 0
  out <- data.frame(provider = providers[keep], as.data.frame(flags[keep, , drop = FALSE]),
                    n_flagged = as.integer(n_flagged[keep]), stringsAsFactors = FALSE,
                    check.names = FALSE)
  out <- out[order(-out$n_flagged, out$provider), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "scenarios") <- names(tests)
  class(out) <- c("pp_stability", "data.frame")
  out
}

#' @rdname plots
#' @export
plot.pp_stability <- function(x, ...) {
  scenarios <- attr(x, "scenarios")
  d <- data.frame(provider = rep(x$provider, length(scenarios)),
                  scenario = rep(scenarios, each = nrow(x)),
                  flag = unlist(lapply(scenarios, function(s) x[[s]]), use.names = FALSE),
                  stringsAsFactors = FALSE)
  d$status <- status_of(d$flag)
  d$provider <- factor(d$provider, levels = rev(x$provider))
  d$scenario <- factor(d$scenario, levels = scenarios)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$scenario, y = .data$provider, fill = .data$status)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::scale_fill_manual(values = status_colours, drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
    theme_pp() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

#' @export
print.pp_stability <- function(x, ...) {
  cat("Flags under ", length(attr(x, "scenarios")), " scenarios: ", nrow(x), " providers",
      " (", sum(x$n_flagged == length(attr(x, "scenarios"))), " flagged under every scenario)\n",
      sep = "")
  print.data.frame(utils::head(as.data.frame(x), 10), row.names = FALSE)
  if (nrow(x) > 10) cat("  ... ", nrow(x) - 10, " more rows\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------

utils::globalVariables("density") # ggplot2::after_stat(density) in plot.pp_null

status_levels <- c("above", "below", "not different", "not tested")
status_colours <- c(above = "#D55E00", below = "#0072B2", `not different` = "#8C8C8C",
                    `not tested` = "#D9D9D9")
status_shapes <- c(above = 17, below = 15, `not different` = 16, `not tested` = 4)

status_of <- function(flag) {
  factor(ifelse(is.na(flag), "not tested",
                ifelse(flag > 0, "above", ifelse(flag < 0, "below", "not different"))),
         levels = status_levels)
}

status_scales <- function() {
  list(ggplot2::scale_colour_manual(values = status_colours, drop = FALSE),
       ggplot2::scale_shape_manual(values = status_shapes, drop = FALSE))
}

theme_pp <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom",
                   plot.title.position = "plot")
}

as_pp_measure <- function(out) {
  class(out) <- c("pp_measure", "data.frame")
  out
}

pp_summary <- function(out) {
  class(out) <- c("pp_summary", "data.frame")
  out
}

label_points <- function(p, d, highlight, xvar, yvar) {
  if (is.null(highlight)) return(p)
  lab <- d[d$provider %in% as.character(highlight), , drop = FALSE]
  if (!nrow(lab)) return(p)
  p + ggplot2::geom_text(data = lab,
                         ggplot2::aes(x = .data[[xvar]], y = .data[[yvar]],
                                      label = .data$provider),
                         size = 3, vjust = -0.8, inherit.aes = FALSE)
}

# The y-axis label of a test result's estimate.
estimate_label <- function(x) {
  s <- attr(x, "test")
  if (identical(s$model, "cox")) return("Observed / expected")
  "Provider effect"
}

caterpillar_plot <- function(x, highlight) {
  d <- data.frame(provider = x$provider, estimate = x$estimate, lower = x$lower, upper = x$upper,
                  status = status_of(x$flag), stringsAsFactors = FALSE)
  d <- d[order(d$estimate), , drop = FALSE]
  d$rank <- seq_len(nrow(d))
  null_value <- x$null_value[1]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$rank, y = .data$estimate))
  if (any(is.finite(d$lower) | is.finite(d$upper))) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                                                 colour = .data$status),
                                    width = 0, linewidth = 0.4, alpha = 0.7, na.rm = TRUE)
  }
  p <- p +
    ggplot2::geom_hline(yintercept = null_value, colour = "grey45", linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$status, shape = .data$status), size = 1.6) +
    status_scales() +
    ggplot2::labs(x = "Providers in order of the estimate", y = estimate_label(x),
                  colour = NULL, shape = NULL) +
    theme_pp() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank())
  if (identical(attr(x, "test")$model, "cox")) p <- p + ggplot2::scale_y_log10()
  label_points(p, d, highlight, "rank", "estimate")
}

funnel_plot_test <- function(x, highlight, levels) {
  s <- attr(x, "test")
  d <- data.frame(provider = x$provider, estimate = x$estimate, status = status_of(x$flag),
                  stringsAsFactors = FALSE)
  d$group <- if ("null_group" %in% names(x)) x$null_group else NA_integer_
  calibrated <- "null_mean" %in% names(x)
  null_value <- x$null_value[1]
  if (identical(s$model, "cox")) {
    if (identical(s$statistic, "wald")) {
      return(funnel_wald(x, d, highlight, levels, "Observed events", log_scale = TRUE))
    }
    d$x <- x$expected
    grid <- exp(seq(log(max(min(x$expected), 0.05)), log(max(x$expected) * 1.05),
                    length.out = 120))
    groups <- if (calibrated) unique(x[, c("null_group", "null_mean", "null_sd")]) else
      data.frame(null_group = NA_integer_, null_mean = 0, null_sd = 1)
    q <- if (!is.null(s$critical)) s$critical else normal_quantile(s$level, s$alternative)
    curves <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
      lim <- poisson_limit_curve(grid, s$statistic, groups$null_mean[i], groups$null_sd[i], q,
                                 s$alternative, s$level)
      data.frame(x = grid, lower = lim$lower, upper = lim$upper, level = s$level,
                 group = groups$null_group[i])
    }))
    return(funnel_draw(d, "x", "estimate", null_value, curves, highlight, "Expected events",
                       "Observed / expected", log_scale = TRUE, facet = calibrated))
  }
  if (identical(s$statistic, "score") && all(c("observed", "expected", "variance") %in% names(x))) {
    return(funnel_score(x, d, highlight, levels))
  }
  funnel_wald(x, d, highlight, levels, "Precision, 1 / se\u00b2", log_scale = FALSE)
}

# Funnel of a score test: observed over expected events against the precision
# E^2 / V, with V the variance the statistic was referred to, so that the
# limits 1 + (m +/- q s) / sqrt(E^2 / V) hold exactly the providers the test
# flags (z = (O - E) / sqrt(V), referred to N(m, s^2)).
funnel_score <- function(x, d, highlight, levels) {
  s <- attr(x, "test")
  ok <- is.finite(x$variance) & x$variance > 0 & is.finite(x$expected) & x$expected > 0
  if (!any(ok)) {
    stop("The funnel plot of a score test needs finite expected counts and variances ",
         "(the restricted score test has none).", call. = FALSE)
  }
  calibrated <- "null_mean" %in% names(x)
  d$estimate <- x$observed / x$expected
  d$x <- x$expected^2 / x$variance
  d <- d[ok, , drop = FALSE]
  grid <- exp(seq(log(min(d$x)), log(max(d$x) * 1.05), length.out = 150))
  groups <- if (calibrated) unique(x[ok, c("null_group", "null_mean", "null_sd")]) else
    data.frame(null_group = NA_integer_, null_mean = 0, null_sd = 1)
  curves <- do.call(rbind, lapply(seq_along(levels), function(k) {
    l <- levels[k]
    q <- if (k == 1 && !is.null(s$critical)) s$critical else normal_quantile(l, s$alternative)
    do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
      m <- groups$null_mean[i]
      sd <- groups$null_sd[i]
      data.frame(x = grid, lower = pmax(1 + (m - q * sd) / sqrt(grid), 0),
                 upper = 1 + (m + q * sd) / sqrt(grid), level = l, group = groups$null_group[i])
    }))
  }))
  if (s$alternative == "greater") curves$lower <- NA_real_
  if (s$alternative == "less") curves$upper <- NA_real_
  funnel_draw(d, "x", "estimate", 1, curves, highlight, "Precision, E\u00b2 / V",
              "Observed / expected", facet = calibrated)
}

# Funnel of a test with standard errors: x = 1 / se^2 (or the observed
# events for the log-ratio Wald test, whose se is sqrt(phi / O)).
funnel_wald <- function(x, d, highlight, levels, x_label, log_scale) {
  s <- attr(x, "test")
  if (all(!is.finite(x$se))) {
    stop("The funnel plot needs standard errors (a Wald-type test) or expected counts.",
         call. = FALSE)
  }
  calibrated <- "null_mean" %in% names(x)
  null_value <- x$null_value[1]
  d$x <- if (log_scale) x$observed else 1 / x$se^2
  ok <- is.finite(d$x) & d$x > 0
  grid <- exp(seq(log(min(d$x[ok])), log(max(d$x[ok]) * 1.05), length.out = 150))
  se_of <- if (log_scale) {
    phi <- if (is.null(s$phi) || !is.finite(s$phi)) 1 else s$phi
    function(g) sqrt(phi / g)
  } else {
    function(g) 1 / sqrt(g)
  }
  groups <- if (calibrated) unique(x[, c("null_group", "null_mean", "null_sd")]) else
    data.frame(null_group = NA_integer_, null_mean = 0, null_sd = 1)
  centre <- if (log_scale) log(null_value) else null_value
  back <- if (log_scale) exp else function(v) v
  curves <- do.call(rbind, lapply(seq_along(levels), function(k) {
    l <- levels[k]
    q <- if (k == 1 && !is.null(s$critical)) s$critical else normal_quantile(l, s$alternative)
    do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
      m <- groups$null_mean[i]
      sd <- groups$null_sd[i]
      data.frame(x = grid, lower = back(centre + (m - q * sd) * se_of(grid)),
                 upper = back(centre + (m + q * sd) * se_of(grid)), level = l,
                 group = groups$null_group[i])
    }))
  }))
  if (s$alternative == "greater") curves$lower <- NA_real_
  if (s$alternative == "less") curves$upper <- NA_real_
  funnel_draw(d, "x", "estimate", null_value, curves, highlight, x_label, estimate_label(x),
              log_scale = log_scale, facet = calibrated)
}

# The counts a Poisson test flags at expected counts `e`, as ratio limits
# half-way between the last count not flagged and the first flagged.
poisson_limit_curve <- function(e, statistic, mean, sd, q, alternative, level) {
  lower <- upper <- rep(NA_real_, length(e))
  for (i in seq_along(e)) {
    o <- 0:ceiling(e[i] + 8 * sqrt(e[i]) + 10)
    if (statistic == "midp") {
      z <- (midp_z(o, e[i]) - mean) / sd
      flag <- z_decision(z, level, alternative)$flag
    } else {
      flag <- exact_poisson(o, rep(e[i], length(o)), level, alternative)$flag
    }
    hi <- o[flag == 1 & !is.na(flag)]
    lo <- o[flag == -1 & !is.na(flag)]
    if (length(hi)) upper[i] <- (min(hi) - 0.5) / e[i]
    if (length(lo)) lower[i] <- (max(lo) + 0.5) / e[i]
  }
  list(lower = lower, upper = upper)
}

funnel_draw <- function(d, xvar, yvar, null_value, curves, highlight, x_label, y_label,
                        log_scale = FALSE, facet = FALSE) {
  curves$level <- factor(curves$level)
  first <- levels(curves$level)[1]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[xvar]], y = .data[[yvar]])) +
    ggplot2::geom_hline(yintercept = null_value, colour = "grey45") +
    ggplot2::geom_line(data = curves, ggplot2::aes(x = .data$x, y = .data$lower,
                                                   linetype = .data$level),
                       colour = "grey35", na.rm = TRUE) +
    ggplot2::geom_line(data = curves, ggplot2::aes(x = .data$x, y = .data$upper,
                                                   linetype = .data$level),
                       colour = "grey35", na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$status, shape = .data$status), size = 1.8,
                        na.rm = TRUE) +
    status_scales() +
    ggplot2::scale_linetype_manual(values = stats::setNames(
      c("solid", rep("dashed", nlevels(curves$level) - 1)), levels(curves$level)),
      name = "Level") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = x_label, y = y_label, colour = NULL, shape = NULL) +
    theme_pp()
  if (log_scale) p <- p + ggplot2::scale_y_log10()
  if (facet && length(unique(stats::na.omit(d$group))) > 1) {
    p <- p + ggplot2::facet_wrap(~ group, labeller = ggplot2::labeller(group = function(g) {
      paste("null group", g)
    }))
  }
  if (nlevels(curves$level) == 1 && first == "0.95") p <- p + ggplot2::guides(linetype = "none")
  label_points(p, d, highlight, xvar, yvar)
}

bar_plot <- function(x, n_groups) {
  sizes <- attr(x, "sizes")
  if (is.null(sizes)) stop("The test result carries no provider sizes.", call. = FALSE)
  size <- if (all(is.na(sizes$patients))) sizes$rows else sizes$patients
  size <- size[match(x$provider, sizes$provider)]
  group <- assign_groups(size, n_groups, "quantile")
  d <- data.frame(group = c(paste0("Q", group), rep("all", length(group))),
                  status = status_of(rep(x$flag, 2)))
  d$group <- factor(d$group, levels = c(paste0("Q", seq_len(n_groups)), "all"))
  counts <- as.data.frame(table(group = d$group, status = d$status))
  counts <- counts[counts$status != "not tested" | counts$Freq > 0, , drop = FALSE]
  totals <- tapply(counts$Freq, counts$group, sum)
  counts$share <- counts$Freq / as.numeric(totals[as.character(counts$group)])
  ggplot2::ggplot(counts, ggplot2::aes(x = .data$group, y = .data$share, fill = .data$status)) +
    ggplot2::geom_col(width = 0.7, colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(.data$share > 0.02,
                                                   sprintf("%.1f%%", 100 * .data$share), "")),
                       position = ggplot2::position_stack(vjust = 0.5), size = 3) +
    ggplot2::scale_fill_manual(values = status_colours, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = function(v) paste0(100 * v, "%")) +
    ggplot2::labs(x = "Provider-size group", y = "Share of providers", fill = NULL) +
    theme_pp()
}
