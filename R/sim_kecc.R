#' Simulate data shaped like a KECC quality measure
#'
#' Patients nested in providers with repeated observations per patient, a
#' provider effect, a patient random effect and a few covariates. The outcome
#' is binary per row (a patient-month, say) or binomial out of a number of
#' trials per row (days at risk, say).
#'
#' @param n_providers number of providers.
#' @param provider_size range of the number of patients per provider; sizes
#'   are drawn uniformly between the two values.
#' @param rows_per_patient range of the number of rows per patient.
#' @param type `"binary"` (one trial per row) or `"binomial"` (`n` trials per
#'   row).
#' @param trials range of the number of trials per row for `type =
#'   "binomial"`.
#' @param n_covariates number of covariates (at least 1): `x1` is a
#'   patient-level standard normal, `x2` a patient-level 0/1 with probability
#'   0.4, and any further covariates are row-level standard normals.
#' @param beta covariate coefficients; by default `0.5, -0.4, 0.3` recycled.
#' @param gamma_sd standard deviation of the provider effects around the
#'   logit of `base_rate`.
#' @param patient_sd standard deviation of the patient random effect.
#' @param base_rate event probability of an average patient at the median
#'   provider.
#' @param seed optional random seed.
#'
#' @return A data frame with columns `provider`, `patient`, `n` (trials),
#'   `y` (events) and the covariates, one row per observation, in random
#'   order. The true parameters are kept in `attr(, "truth")`: `gamma`
#'   (named by provider), `beta`, `patient_effect` (named by patient).
#' @examples
#' d <- sim_kecc(n_providers = 5, seed = 1)
#' str(attr(d, "truth"))
#' @export
sim_kecc <- function(n_providers = 50, provider_size = c(15, 60), rows_per_patient = c(1, 6),
                     type = c("binary", "binomial"), trials = c(20, 31), n_covariates = 3,
                     beta = NULL, gamma_sd = 0.5, patient_sd = 0.5, base_rate = 0.25, seed = NULL) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)
  if (n_covariates < 1) stop("`n_covariates` must be at least 1.", call. = FALSE)
  if (is.null(beta)) beta <- rep_len(c(0.5, -0.4, 0.3), n_covariates)
  if (length(beta) != n_covariates) stop("`beta` must have `n_covariates` values.", call. = FALSE)

  providers <- sprintf("P%03d", seq_len(n_providers))
  n_pat <- sample(seq(provider_size[1], provider_size[2]), n_providers, replace = TRUE)
  patient_provider <- rep(providers, n_pat)
  patients <- sprintf("%s-%03d", patient_provider, sequence(n_pat))
  n_rows <- sample(seq(rows_per_patient[1], rows_per_patient[2]), length(patients), replace = TRUE)

  gamma <- stats::qlogis(base_rate) + stats::rnorm(n_providers, 0, gamma_sd)
  gamma <- stats::setNames(gamma, providers)
  b <- stats::setNames(stats::rnorm(length(patients), 0, patient_sd), patients)
  x_patient <- cbind(x1 = stats::rnorm(length(patients)),
                     x2 = stats::rbinom(length(patients), 1, 0.4))

  row_patient <- rep(seq_along(patients), n_rows)
  total <- length(row_patient)
  x <- matrix(NA_real_, total, n_covariates,
              dimnames = list(NULL, paste0("x", seq_len(n_covariates))))
  x[, 1] <- x_patient[row_patient, 1]
  if (n_covariates >= 2) x[, 2] <- x_patient[row_patient, 2]
  if (n_covariates >= 3) x[, 3:n_covariates] <- stats::rnorm(total * (n_covariates - 2))
  n <- if (type == "binary") {
    rep(1L, total)
  } else {
    sample(seq(trials[1], trials[2]), total, replace = TRUE)
  }
  eta <- gamma[patient_provider[row_patient]] + as.numeric(x %*% beta) + b[row_patient]
  y <- stats::rbinom(total, n, stats::plogis(eta))

  d <- data.frame(provider = patient_provider[row_patient], patient = patients[row_patient],
                  n = n, y = y, x, stringsAsFactors = FALSE)
  d <- d[sample(total), , drop = FALSE]
  rownames(d) <- NULL
  attr(d, "truth") <- list(gamma = gamma, beta = stats::setNames(beta, colnames(x)),
                           patient_effect = b)
  d
}
