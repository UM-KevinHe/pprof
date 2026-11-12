# Confidence limits of the provider effects by inverting a score test.
#
# For provider i and a null value g, the score statistic with beta held at the
# full-model estimate is
#   z_i(g) = (O_i - E_i(g)) / sqrt(V_i(g)),   E_i(g) = sum_j N_ij expit(g + Z_ij' beta),
# with V_i(g) the naive variance sum_j N_ij p_ij(g) (1 - p_ij(g)), the
# patient-clustered sandwich sum_k (sum_{j in k} (Y_ij - N_ij p_ij(g)))^2, or
# the stabilized variance sum_k a_k(g)' S a_k(g), a_k = sqrt(N p (1 - p)) over
# patient k's rows and S the correlation matrix pooled over all patients at the
# test's null value (held fixed, so that at g = gamma_0 the statistic is the
# test's own). E_i(g) increases in g, so z_i(g) decreases: the lower limit is
# the root of z_i(g) = z_lower and the upper limit the root of z_i(g) = z_upper
# (for the theoretical null z_lower = c, z_upper = -c; after calibration
# m_g +/- c s_g). A limit is infinite when z_i never reaches its target: the
# lower limit of a provider with no events, the upper limit of one with all
# events, and either limit of a patient-clustered statistic, which is bounded
# by the square root of the number of events (non-events).
score_limits <- function(fit, variance, idx, z_lower, z_upper, pooled = NULL) {
  if (inherits(fit, "discrete_fe")) {
    return(score_limits_discrete(fit, variance, idx, z_lower, z_upper))
  }
  m <- fit_matrices(fit)
  gamma <- fit$coefficients$gamma
  lp <- fit$linear_predictor
  code <- provider_codes(m$prov, m$levels)
  rows_of <- split(seq_along(code), factor(code, levels = seq_along(m$levels)))
  patient <- if (variance == "naive") NULL else cluster_codes(m$prov, m$patient)$code
  if (variance == "stabilized" && is.null(pooled)) {
    stop("The stabilized inversion needs the pooled correlation matrix.", call. = FALSE)
  }
  lower <- upper <- rep(NA_real_, length(idx))
  for (k in seq_along(idx)) {
    i <- idx[k]
    r <- rows_of[[i]]
    zb <- lp[r]
    nt <- m$n[r]
    y <- m$y[r]
    obs <- sum(y)
    pat <- if (is.null(patient)) NULL else match(patient[r], unique(patient[r]))
    pairs <- if (variance == "stabilized") stabilized_pairs(pat, pooled) else NULL
    z_of <- function(g) {
      p <- stats::plogis(g + zb)
      u <- obs - sum(nt * p)
      v <- switch(variance,
                  naive = sum(nt * p * (1 - p)),
                  robust = sum(rowsum(y - nt * p, pat, reorder = FALSE)^2),
                  stabilized = {
                    a <- sqrt(pmax(nt * p * (1 - p), 0))
                    sum(a[pairs$j] * pairs$s * a[pairs$l])
                  })
      if (!is.finite(v) || v <= 0) return(sign(u) * 1e6)
      u / sqrt(v)
    }
    lower[k] <- monotone_root(z_of, gamma[i], z_lower[k])
    upper[k] <- monotone_root(z_of, gamma[i], z_upper[k])
  }
  list(lower = unname(lower), upper = unname(upper))
}

# The row pairs (j, l) of the same patient and the pooled correlation between
# their positions, for a quadratic form sum_k a_k' S a_k evaluated as
# sum(a[j] * s * a[l]).
stabilized_pairs <- function(pat, pooled) {
  pos <- stats::ave(seq_along(pat), pat, FUN = seq_along)
  by_patient <- split(seq_along(pat), pat)
  j <- unlist(lapply(by_patient, function(r) rep(r, each = length(r))), use.names = FALSE)
  l <- unlist(lapply(by_patient, function(r) rep(r, times = length(r))), use.names = FALSE)
  list(j = j, l = l, s = pooled[cbind(pos[j], pos[l])])
}

# The value x at which a decreasing function f equals `target`, searched from
# x0 outwards (the step doubles up to a span of 60), by uniroot once bracketed;
# +Inf / -Inf when f never reaches the target on that side; an infinite target
# gives the infinite limit on the corresponding side.
monotone_root <- function(f, x0, target, step = 1, max_span = 60) {
  if (!is.finite(target)) return(if (target > 0) -Inf else Inf)
  g <- function(x) f(x) - target
  f0 <- g(x0)
  if (!is.finite(f0)) return(NA_real_)
  if (f0 == 0) return(x0)
  dir <- if (f0 > 0) 1 else -1
  b <- x0 + dir * step
  fb <- g(b)
  while (sign(fb) == sign(f0) && abs(b - x0) < max_span) {
    step <- step * 2
    b <- x0 + dir * step
    fb <- g(b)
  }
  if (sign(fb) == sign(f0)) return(dir * Inf)
  stats::uniroot(g, sort(c(x0, b)), tol = 1e-8)$root
}

# The targets of the inversion for a level, an alternative and a null with
# location `loc` and scale `sd` (vectors, one per provider): the lower limit is
# the root at loc + q sd, the upper at loc - q sd; a one-sided alternative
# leaves the other limit infinite.
inversion_targets <- function(loc, sd, q, alternative) {
  z_lower <- loc + q * sd
  z_upper <- loc - q * sd
  if (alternative == "greater") z_upper <- rep(-Inf, length(z_lower))
  if (alternative == "less") z_lower <- rep(Inf, length(z_upper))
  list(z_lower = z_lower, z_upper = z_upper)
}
