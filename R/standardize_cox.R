#' Standardized ratios of a Cox fit
#'
#' Indirect or direct standardization of every provider's event count from a
#' Cox model, in the form of the survival measures (SMR, SHR).
#'
#' @param fit a model fitted by [fit_cox()].
#' @param type `"indirect"` (default) or `"direct"`.
#' @param providers optional subset of provider identifiers; the baselines
#'   always use every fitted row.
#' @param ... unused.
#'
#' @details With \eqn{w_i} the case weights, \eqn{\delta_i} the event
#' indicator, \eqn{(s_i, t_i]} the at-risk interval and \eqn{\eta_i} the
#' fitted linear predictor of row \eqn{i}:
#' * indirect: \eqn{O_j = \sum_{i \in j} w_i \delta_i} and
#'   \eqn{E_j = \sum_{i \in j} w_i\, e^{\eta_i} [\Lambda_0(t_i) -
#'   \Lambda_0(s_i)]}, the events expected for provider \eqn{j}'s rows at the
#'   national baseline \eqn{\Lambda_0} (the Breslow estimator over every row);
#'   ratio \eqn{O_j / E_j}; the expected counts add up to the observed ones;
#' * direct: \eqn{E^{(j)} = \sum_i w_i\, e^{\eta_i} [\Lambda_{0j}(t_i) -
#'   \Lambda_{0j}(s_i)]} over every row, the events expected if every patient
#'   had provider \eqn{j}'s own baseline \eqn{\Lambda_{0j}} (the Breslow
#'   estimator over provider \eqn{j}'s rows at the same coefficients); ratio
#'   \eqn{E^{(j)} / O} with \eqn{O} the total observed events. Each Breslow
#'   increment of provider \eqn{j} being one weighted event over its risk-set
#'   sum \eqn{RS_j(t)}, \eqn{E^{(j)}} equals \eqn{\sum} over provider
#'   \eqn{j}'s events of \eqn{w_i\, RS(t_i) / RS_j(t_i)}, the ratio of the
#'   population's to the provider's risk-set sums of \eqn{w e^\eta}.
#'
#' Both baselines are Breslow estimators whatever `ties` fitted the
#' coefficients. The person-time is \eqn{\sum_{i \in j} (t_i - s_i)}, the
#' grouping variable of the survival measures' empirical null.
#'
#' @return A data frame with one row per provider. Indirect: `provider`,
#'   `observed`, `expected`, `ratio`, `person_time`. Direct: `provider`,
#'   `observed` (the total \eqn{O}), `expected` (\eqn{E^{(j)}}), `ratio`,
#'   `n_pop` (the number of fitted rows). Attribute `type`.
#' @references He K, Schaubel DE (2015). Standardized mortality ratio for
#'   evaluating center-specific mortality: assessment and alternative.
#'   *Statistics in Biosciences*, 7, 296-321.
#' @examples
#' s <- sim_kecc(n_providers = 20, type = "survival", seed = 1)
#' ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", time = "stop", start = "start",
#'               event = "event", cutoff = 0, message = FALSE)
#' fit <- fit_cox(ps)
#' head(standardize(fit))
#' head(standardize(fit, type = "direct"))
#' @export
standardize.cox <- function(fit, type = c("indirect", "direct"), providers = NULL, ...) {
  check_fit(fit, "cox")
  type <- match.arg(type)
  m <- cox_matrices(fit)
  idx <- provider_positions(providers, m$levels)
  if (type == "indirect") {
    ex <- cox_expectations(fit$linear_predictor, m$start, m$stop, m$event, m$w, m$code,
                           length(m$levels), direct = FALSE, rows = fit$expected)
    out <- data.frame(provider = m$levels[idx], observed = ex$observed[idx],
                      expected = ex$expected[idx], ratio = ex$observed[idx] / ex$expected[idx],
                      person_time = ex$person_time[idx], stringsAsFactors = FALSE)
  } else {
    ex <- cox_expectations(fit$linear_predictor, m$start, m$stop, m$event, m$w, m$code,
                           length(m$levels), direct = TRUE, rows = fit$expected)
    total <- sum(ex$observed)
    out <- data.frame(provider = m$levels[idx], observed = rep(total, length(idx)),
                      expected = ex$direct[idx], ratio = ex$direct[idx] / total,
                      n_pop = rep(length(m$code), length(idx)), stringsAsFactors = FALSE)
  }
  rownames(out) <- NULL
  attr(out, "type") <- type
  as_pp_measure(out)
}

# The columns of a cox fit as vectors: at-risk interval, event, weights,
# provider codes and levels.
cox_matrices <- function(fit) {
  d <- fit$data
  v <- fit$vars
  n <- nrow(d)
  levels <- fit$sizes$provider
  list(start = if (is.null(v$start)) rep(0, n) else as.numeric(d[[v$start]]),
       stop = as.numeric(d[[v$time]]), event = as.numeric(d[[v$event]]),
       w = if (is.null(v$weights)) rep(1, n) else as.numeric(d[[v$weights]]),
       code = provider_codes(d[[v$provider]], levels), levels = levels,
       patient = if (is.null(v$patient)) NULL else d[[v$patient]])
}

# Observed and expected events per provider for the indirect (national
# baseline) and direct (provider baselines) standardized ratios, the expected
# events per row and the national Breslow baseline. `rows` supplies the row
# expectations when already computed.
cox_expectations <- function(eta, start, stop, event, w, code, m, direct = FALSE, rows = NULL) {
  shift <- max(eta)
  r <- w * exp(eta - shift)
  is_event <- event > 0
  times <- sort(unique(stop[is_event]))
  d_t <- as.numeric(rowsum(w[is_event] * event[is_event], match(stop[is_event], times),
                           reorder = TRUE))
  rs_t <- at_risk_sums(times, start, stop, r)
  cumhaz <- cumsum(d_t / rs_t)
  national <- function(t) {
    k <- findInterval(t, times)
    ifelse(k > 0, cumhaz[pmax(k, 1)], 0)
  }
  if (is.null(rows)) rows <- exp(eta - shift) * (national(stop) - national(start))
  out <- list(rows = rows,
              baseline = data.frame(time = times, hazard = cumhaz * exp(-shift)),
              observed = as.numeric(rowsum(w * event, code, reorder = TRUE)),
              expected = as.numeric(rowsum(w * rows, code, reorder = TRUE)),
              person_time = as.numeric(rowsum(stop - start, code, reorder = TRUE)))
  if (!direct) return(out)
  ev_rows <- which(is_event)
  rs_pop <- at_risk_sums(stop[ev_rows], start, stop, r)
  contribution <- numeric(length(ev_rows))
  rows_of <- split(seq_along(code), factor(code, levels = seq_len(m)))
  ev_of <- split(seq_along(ev_rows), factor(code[ev_rows], levels = seq_len(m)))
  for (j in seq_len(m)) {
    e <- ev_of[[j]]
    if (!length(e)) next
    own <- rows_of[[j]]
    rs_own <- at_risk_sums(stop[ev_rows[e]], start[own], stop[own], r[own])
    contribution[e] <- w[ev_rows[e]] * event[ev_rows[e]] * rs_pop[e] / rs_own
  }
  out$direct <- as.numeric(rowsum(contribution, factor(code[ev_rows], levels = seq_len(m)),
                                  reorder = TRUE))
  out
}

# sum of r over the rows at risk at each query time: start < q <= stop.
at_risk_sums <- function(query, start, stop, r) {
  suffix <- function(keys) {
    o <- order(keys)
    list(k = keys[o], s = c(rev(cumsum(rev(r[o]))), 0))
  }
  ex <- suffix(stop)
  en <- suffix(start)
  ex$s[findInterval(query, ex$k, left.open = TRUE) + 1] -
    en$s[findInterval(query, en$k, left.open = TRUE) + 1]
}
