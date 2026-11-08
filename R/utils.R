# Internal helpers shared across components.

# The provider effect under the null hypothesis: the median or the
# size-weighted mean of the fitted effects, or a number.
null_value <- function(null, gamma, sizes = NULL) {
  if (identical(null, "median")) return(stats::median(gamma))
  if (identical(null, "mean")) {
    if (is.null(sizes)) return(mean(gamma))
    return(sum(gamma * sizes) / sum(sizes))
  }
  if (is.numeric(null) && length(null) == 1) return(null)
  stop("`null` must be \"median\", \"mean\" or a number.", call. = FALSE)
}

match_alternative <- function(alternative) {
  if (length(alternative) != 1 || !alternative %in% c("two.sided", "greater", "less")) {
    stop("`alternative` must be \"two.sided\", \"greater\" or \"less\".", call. = FALSE)
  }
  alternative
}

# p-value and flag of a z statistic (positive z = above the reference).
z_decision <- function(z, level, alternative) {
  alpha <- 1 - level
  p_upper <- stats::pnorm(z, lower.tail = FALSE)
  if (alternative == "two.sided") {
    flag <- ifelse(p_upper < alpha / 2, 1L, ifelse(p_upper <= 1 - alpha / 2, 0L, -1L))
    p <- 2 * pmin(p_upper, 1 - p_upper)
  } else if (alternative == "greater") {
    flag <- ifelse(p_upper < alpha, 1L, 0L)
    p <- p_upper
  } else {
    p_lower <- 1 - p_upper
    flag <- ifelse(p_lower < alpha, -1L, 0L)
    p <- p_lower
  }
  flag[!is.finite(z)] <- NA_integer_
  p[!is.finite(z)] <- NA_real_
  list(flag = flag, p = p)
}

# Integer codes 1..K for (provider, patient) clusters in order of first
# appearance, and the first row of each cluster. The pair is encoded as one
# number, exact below 2^53.
cluster_codes <- function(prov, patient) {
  pc <- match(prov, unique(prov))
  qc <- match(patient, unique(patient))
  key <- (pc - 1) * max(qc, 0L) + qc
  cl <- match(key, unique(key))
  list(code = cl, first = !duplicated(cl))
}

# Integer codes 1..m of `prov` following `levels`.
provider_codes <- function(prov, levels) {
  code <- match(as.character(prov), levels)
  if (anyNA(code)) stop("Provider identifiers do not match the fitted providers.", call. = FALSE)
  code
}

# Positions of the providers named in `providers` within `levels` (all when NULL).
provider_positions <- function(providers, levels) {
  if (is.null(providers)) return(seq_along(levels))
  pos <- match(as.character(providers), levels)
  if (anyNA(pos)) {
    stop("Unknown provider(s): ", paste(providers[is.na(pos)], collapse = ", "), call. = FALSE)
  }
  pos
}
