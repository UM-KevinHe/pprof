# Internal helpers shared by the testing and standardization functions.

# Check that `fit` is a logis_fe object.
check_fit <- function(fit) {
  if (!inherits(fit, "logis_fe")) stop("Object fit is not of the classes 'logis_fe'!", call. = FALSE)
  invisible(TRUE)
}

# Number of trials for the rows of `data`; 1 when there is no N column
# (fits created without a number of trials).
trials_of <- function(data, N.char) {
  if (is.null(N.char) || !N.char %in% colnames(data)) return(rep(1, NROW(data)))
  as.numeric(data[, N.char])
}

# Resolve the provider effect under the null hypothesis.
gamma_null_value <- function(null, gamma) {
  if (identical(null, "median")) return(stats::median(gamma))
  if (is.numeric(null)) return(null[1])
  stop("Argument 'null' NOT as required!", call. = FALSE)
}

# Provider-level flags from the upper-tail probability `p_upper`
# (1 = significantly higher than expected, -1 = lower, 0 = as expected).
flag_two_sided <- function(p_upper, alpha) {
  ifelse(p_upper < alpha / 2, 1, ifelse(p_upper <= 1 - alpha / 2, 0, -1))
}

match_alternative <- function(alternative) {
  if (length(alternative) != 1 || !alternative %in% c("two.sided", "greater", "less"))
    stop("Argument 'alternative' should be one of 'two.sided', 'greater', or 'less'", call. = FALSE)
  alternative
}

# p-value and flag for a z statistic (upper tail = "higher than expected").
z_test_output <- function(z, alpha, alternative) {
  p_upper <- stats::pnorm(z, lower.tail = FALSE)
  if (alternative == "two.sided") {
    list(flag = flag_two_sided(p_upper, alpha), p = 2 * pmin(p_upper, 1 - p_upper))
  } else if (alternative == "greater") {
    list(flag = ifelse(p_upper < alpha, 1, 0), p = p_upper)
  } else {
    p_lower <- stats::pnorm(z, lower.tail = TRUE)
    list(flag = ifelse(p_lower < alpha, -1, 0), p = p_lower)
  }
}

# Integer codes 1..K for (provider, patient) clusters, numbered in order of
# first appearance, and the first row of each cluster. The pair is encoded as
# one number, (provider code - 1) * #patients + patient code, which is exact
# below 2^53 and avoids building a string per row.
cluster_codes <- function(prov, patient) {
  pc <- match(prov, unique(prov))
  qc <- match(patient, unique(patient))
  key <- (pc - 1) * max(qc, 0L) + qc
  cl <- match(key, unique(key))
  list(code = cl, first = !duplicated(cl))
}

# Provider-level "meat" of the patient-clustered sandwich for gamma: for each
# provider (in the order `prov_levels`), the sum over its patients of
# (patient sum of residuals)^2. cluster_sums() adds in the same order and
# precision as base-R sum(), which matches the dplyr summarise() this replaces.
robust_gamma_meat <- function(residual, prov, patient, prov_levels) {
  cl <- cluster_codes(prov, patient)
  prov_row <- provider_codes(prov, prov_levels)
  meat <- cluster_sums(residual, cl$code, prov_row, prov_row[cl$first], length(prov_levels))$meat
  stats::setNames(meat, prov_levels)
}

# The three blocks of the sandwich "meat" for beta (see robust_wald_covar):
#   A0   : per provider, sum over patients of r_k^2          (length m)
#   A1_t : per provider, sum over patients of r_k * U_k      (m x p)
#   A2   : sum over patients of U_k U_k^T                    (p x p)
# where r_k and U_k are the patient-level sums of the residual and of
# Z * residual. Providers follow `prov_levels`.
robust_beta_meat <- function(residual, Z, prov, patient, prov_levels) {
  cl <- cluster_codes(prov, patient)
  prov_cl <- provider_codes(prov[cl$first], prov_levels)
  r_k <- as.numeric(rowsum(residual, cl$code, reorder = FALSE))
  U_k <- rowsum(Z * residual, cl$code, reorder = FALSE)
  A0 <- as.numeric(rowsum(r_k^2, prov_cl, reorder = TRUE))
  A1_t <- rowsum(r_k * U_k, prov_cl, reorder = TRUE)
  dimnames(A1_t) <- NULL
  dimnames(U_k) <- NULL
  list(A0 = A0, A1_t = A1_t, A2 = compute_A2(U_k))
}

# Integer codes 1..m of `prov` following `prov_levels`; rowsum(reorder = TRUE)
# on these codes returns providers in `prov_levels` order.
provider_codes <- function(prov, prov_levels) {
  code <- match(as.character(prov), prov_levels)
  if (anyNA(code)) stop("Provider identifiers do not match the fitted providers.", call. = FALSE)
  code
}

# Resolve `parm` (subset of providers) into positions within `prov_levels`.
provider_indices <- function(parm, prov_values, prov_levels) {
  if (missing(parm) || is.null(parm)) return(seq_along(prov_levels))
  if (is.numeric(parm)) parm <- as.numeric(parm)
  if (!identical(class(parm), class(prov_values)) &&
      !(is.numeric(parm) && is.numeric(prov_values)))
    stop("Argument 'parm' includes invalid elements!", call. = FALSE)
  which(prov_levels %in% as.character(parm))
}
