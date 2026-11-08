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

# Flag from an upper-tail probability (exact tests): 1 above, -1 below, 0 otherwise.
flag_two_sided <- function(p_upper, alpha) {
  ifelse(p_upper < alpha / 2, 1L, ifelse(p_upper <= 1 - alpha / 2, 0L, -1L))
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

# Provider-level "meat" of the patient-clustered sandwich for gamma: for each
# provider (in the order `levels`) the sum over its patients of (patient sum of
# residuals)^2. cluster_sums() adds in row order with long double accumulators.
robust_gamma_meat <- function(residual, prov, patient, levels) {
  cl <- cluster_codes(prov, patient)
  prov_row <- provider_codes(prov, levels)
  meat <- cluster_sums(residual, cl$code, prov_row, prov_row[cl$first], length(levels))$meat
  stats::setNames(meat, levels)
}

# The three blocks of the sandwich "meat" for beta (see robust_wald_covar):
#   A0   : per provider, sum over patients of r_k^2          (length m)
#   A1_t : per provider, sum over patients of r_k * U_k      (m x p)
#   A2   : sum over patients of U_k U_k^T                    (p x p)
# where r_k and U_k are the patient-level sums of the residual and of
# Z * residual. Providers follow `levels`.
robust_beta_meat <- function(residual, z, prov, patient, levels) {
  cl <- cluster_codes(prov, patient)
  prov_cl <- provider_codes(prov[cl$first], levels)
  r_k <- as.numeric(rowsum(residual, cl$code, reorder = FALSE))
  u_k <- rowsum(z * residual, cl$code, reorder = FALSE)
  a0 <- as.numeric(rowsum(r_k^2, prov_cl, reorder = TRUE))
  a1_t <- rowsum(r_k * u_k, prov_cl, reorder = TRUE)
  dimnames(a1_t) <- NULL
  dimnames(u_k) <- NULL
  list(A0 = a0, A1_t = a1_t, A2 = compute_A2(u_k))
}

# Normal quantile for a confidence level and alternative.
normal_quantile <- function(level, alternative) {
  if (alternative == "two.sided") stats::qnorm(1 - (1 - level) / 2) else stats::qnorm(level)
}
