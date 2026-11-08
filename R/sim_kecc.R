#' Simulate data shaped like a KECC quality measure
#'
#' Patients nested in providers with repeated observations per patient, a
#' provider effect, a patient random effect and a few covariates. The outcome
#' is binary per row (a patient-month, say), binomial out of a number of
#' trials per row (days at risk, say), continuous, or a survival time with
#' censoring, split into counting-process rows.
#'
#' @param n_providers number of providers.
#' @param provider_size range of the number of patients per provider; sizes
#'   are drawn uniformly between the two values.
#' @param rows_per_patient range of the number of rows per patient.
#' @param type `"binary"` (one trial per row), `"binomial"` (`n` trials per
#'   row), `"continuous"` or `"survival"`.
#' @param trials range of the number of trials per row for `type =
#'   "binomial"`.
#' @param base_hazard,follow_up for `type = "survival"`: the hazard of an
#'   average patient at the median provider (exponential event times), and
#'   the range of the censoring time, drawn uniformly. Each patient's
#'   follow-up is cut into `rows_per_patient` consecutive rows `(start,
#'   stop]`, the last carrying the event if it came before censoring.
#' @param n_covariates number of covariates (at least 1): `x1` is a
#'   patient-level standard normal, `x2` a patient-level 0/1 with probability
#'   0.4, and any further covariates are row-level standard normals
#'   (patient-level for a survival outcome).
#' @param beta covariate coefficients; by default `0.5, -0.4, 0.3` recycled.
#' @param gamma_sd standard deviation of the provider effects around the
#'   logit of `base_rate` (binary, binomial) or around 0 (continuous).
#' @param patient_sd standard deviation of the patient random effect.
#' @param base_rate event probability of an average patient at the median
#'   provider (binary and binomial outcomes).
#' @param sigma residual standard deviation of a continuous outcome.
#' @param seed optional random seed.
#'
#' @return A data frame with columns `provider`, `patient`, `n` (trials; 1
#'   for binary and continuous outcomes), `y` and the covariates, one row per
#'   observation, in random order; for a survival outcome the columns
#'   `start`, `stop` and `event` replace `n` and `y`. The true parameters are
#'   kept in `attr(, "truth")`: `gamma` (named by provider), `beta`,
#'   `patient_effect` (named by patient).
#' @examples
#' d <- sim_kecc(n_providers = 5, seed = 1)
#' str(attr(d, "truth"))
#' s <- sim_kecc(n_providers = 5, type = "survival", seed = 1)
#' head(s)
#' @export
sim_kecc <- function(n_providers = 50, provider_size = c(15, 60), rows_per_patient = c(1, 6),
                     type = c("binary", "binomial", "continuous", "survival"), trials = c(20, 31),
                     n_covariates = 3, beta = NULL, gamma_sd = 0.5, patient_sd = 0.5,
                     base_rate = 0.25, sigma = 1, base_hazard = 0.5, follow_up = c(0.5, 2),
                     seed = NULL) {
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

  centre <- if (type %in% c("continuous", "survival")) 0 else stats::qlogis(base_rate)
  gamma <- stats::setNames(centre + stats::rnorm(n_providers, 0, gamma_sd), providers)
  b <- stats::setNames(stats::rnorm(length(patients), 0, patient_sd), patients)
  x_patient <- cbind(x1 = stats::rnorm(length(patients)),
                     x2 = stats::rbinom(length(patients), 1, 0.4))

  row_patient <- rep(seq_along(patients), n_rows)
  total <- length(row_patient)
  x <- matrix(NA_real_, total, n_covariates,
              dimnames = list(NULL, paste0("x", seq_len(n_covariates))))
  x[, 1] <- x_patient[row_patient, 1]
  if (n_covariates >= 2) x[, 2] <- x_patient[row_patient, 2]
  if (n_covariates >= 3) {
    x[, 3:n_covariates] <- if (type == "survival") {
      matrix(stats::rnorm(length(patients) * (n_covariates - 2)), ncol = n_covariates - 2)[
        row_patient, , drop = FALSE]
    } else {
      stats::rnorm(total * (n_covariates - 2))
    }
  }
  if (type == "survival") {
    first <- !duplicated(row_patient)
    eta_p <- gamma[patient_provider] + as.numeric(x[first, , drop = FALSE] %*% beta) + b
    t_event <- stats::rexp(length(patients), base_hazard * exp(eta_p))
    t_cens <- stats::runif(length(patients), follow_up[1], follow_up[2])
    t_end <- pmin(t_event, t_cens)
    cuts <- stats::runif(total)
    stop <- numeric(total)
    start <- numeric(total)
    event <- numeric(total)
    for (k in seq_along(patients)) {
      rows <- which(row_patient == k)
      u <- sort(cuts[rows])
      s <- t_end[k] * u / u[length(u)]
      stop[rows] <- s
      start[rows] <- c(0, s[-length(s)])
      event[rows[length(rows)]] <- as.numeric(t_event[k] <= t_cens[k])
    }
    d <- data.frame(provider = patient_provider[row_patient], patient = patients[row_patient],
                    start = start, stop = stop, event = event, x, stringsAsFactors = FALSE)
    d <- d[sample(total), , drop = FALSE]
    rownames(d) <- NULL
    attr(d, "truth") <- list(gamma = gamma, beta = stats::setNames(beta, colnames(x)),
                             patient_effect = b, base_hazard = base_hazard)
    return(d)
  }
  n <- if (type == "binomial") {
    sample(seq(trials[1], trials[2]), total, replace = TRUE)
  } else {
    rep(1L, total)
  }
  eta <- gamma[patient_provider[row_patient]] + as.numeric(x %*% beta) + b[row_patient]
  y <- if (type == "continuous") {
    eta + stats::rnorm(total, 0, sigma)
  } else {
    stats::rbinom(total, n, stats::plogis(eta))
  }

  d <- data.frame(provider = patient_provider[row_patient], patient = patients[row_patient],
                  n = n, y = y, x, stringsAsFactors = FALSE)
  d <- d[sample(total), , drop = FALSE]
  rownames(d) <- NULL
  attr(d, "truth") <- list(gamma = gamma, beta = stats::setNames(beta, colnames(x)),
                           patient_effect = b)
  d
}
