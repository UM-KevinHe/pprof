#' Tests of the covariate coefficients of a logistic fixed-effect fit
#'
#' @param fit a model fitted by [fit_logistic_fe()].
#' @param test `"wald"` (default; model-based), `"robust_wald"`
#'   (patient-clustered sandwich; needs a patient column), `"lr"` (likelihood
#'   ratio; refits the model without the covariate) or `"score"` (refits the
#'   model without the covariate).
#' @param covariates optional subset of covariate names.
#' @param level confidence level for the Wald intervals.
#' @param null the null value of each coefficient; `"lr"` and `"score"` only
#'   allow 0.
#' @param ... unused.
#'
#' @return A data frame with one row per covariate: `covariate`, `estimate`,
#'   `se`, `statistic` (a z statistic for the Wald tests, a chi-square
#'   statistic with one degree of freedom for `"lr"` and `"score"`), `p`,
#'   `lower` and `upper` (Wald tests only).
#' @examples
#' d <- sim_kecc(n_providers = 20, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_logistic_fe(pd)
#' test_covariates(fit, test = "robust_wald")
#' @export
test_covariates.logistic_fe <- function(fit, test = c("wald", "robust_wald", "lr", "score"),
                                        covariates = NULL, level = 0.95, null = 0, ...) {
  check_fit(fit)
  test <- match.arg(test)
  x <- fit$vars$x
  ind <- if (is.null(covariates)) seq_along(x) else match(covariates, x)
  if (anyNA(ind)) {
    stop("Unknown covariate(s): ", paste(covariates[is.na(ind)], collapse = ", "), call. = FALSE)
  }
  m <- fit_matrices(fit)
  gamma <- fit$coefficients$gamma
  beta <- fit$coefficients$beta
  alpha <- 1 - level

  if (test %in% c("wald", "robust_wald")) {
    if (test == "wald") {
      res <- wald_covar(m$n, m$y, m$z, m$n_prov, gamma, beta, ind, null, alpha)
    } else {
      if (is.null(m$patient)) {
        stop("test = \"robust_wald\" needs a patient column in pp_data().", call. = FALSE)
      }
      residual <- m$y - m$n * fit$fitted
      meat <- robust_beta_meat(residual, m$z, m$prov, m$patient, m$levels)
      res <- robust_wald_covar(m$n, m$y, m$z, m$n_prov, gamma, beta, ind, null, alpha,
                               meat$A0, meat$A1_t, meat$A2)
    }
    return(data.frame(covariate = x[ind], estimate = unname(beta[ind]),
                      se = as.numeric(res$se.beta),
                      statistic = as.numeric(res$stat), p = as.numeric(res$p),
                      lower = as.numeric(res$beta.lower), upper = as.numeric(res$beta.upper),
                      stringsAsFactors = FALSE))
  }

  if (null != 0) stop("`null` must be 0 for the lr and score tests.", call. = FALSE)
  if (length(x) < 2) stop("The lr and score tests need at least two covariates.", call. = FALSE)
  ctrl <- fit$control
  refit_without <- function(j) {
    pd <- structure(list(data = fit$data, vars = utils::modifyList(fit$vars, list(x = x[-j])),
                         sizes = fit$sizes, cutoff = fit$cutoff), class = "pp_data")
    fit_logistic_fe(pd, algorithm = ctrl$algorithm, max_iter = ctrl$max_iter, tol = ctrl$tol,
                    bound = ctrl$bound, backtrack = ctrl$backtrack, threads = ctrl$threads,
                    stop = ctrl$stop)
  }
  if (test == "lr") {
    stat <- vapply(ind, function(j) {
      null_fit <- refit_without(j)
      -2 * (null_fit$loglik - fit$loglik)
    }, numeric(1))
  } else {
    code <- provider_codes(m$prov, m$levels)
    stat <- vapply(ind, function(j) {
      null_fit <- refit_without(j)
      z_null <- m$z[, -j, drop = FALSE]
      z_add <- m$z[, j]
      p_null <- stats::plogis(rep(null_fit$coefficients$gamma, m$n_prov) +
                                as.numeric(z_null %*% null_fit$coefficients$beta))
      p_null <- pmin(pmax(p_null, 1e-10), 1 - 1e-10)
      w <- m$n * p_null * (1 - p_null)
      info_gamma_inv <- 1 / as.numeric(rowsum(w, code, reorder = TRUE))
      info_betagamma <- rowsum(w * z_null, code, reorder = TRUE)
      info_beta <- t(z_null) %*% (w * z_null)
      mat_tmp <- info_gamma_inv * info_betagamma
      schur_inv <- solve(info_beta - t(info_betagamma) %*% mat_tmp)
      info_add_gamma <- as.numeric(rowsum(w * z_add, code, reorder = TRUE))
      info_add_beta <- t(z_add) %*% (w * z_null)
      # (gamma, gamma) block of the inverse applied to info_add_gamma, without the m x m matrix
      t1 <- mat_tmp %*% (schur_inv %*% (t(mat_tmp) %*% info_add_gamma))
      quad_gamma <- sum(info_add_gamma^2 * info_gamma_inv) + sum(info_add_gamma * t1)
      cross <- -info_add_beta %*% schur_inv %*% (t(mat_tmp) %*% info_add_gamma)
      info_add <- sum(z_add * w * z_add) - quad_gamma - 2 * as.numeric(cross) -
        as.numeric(info_add_beta %*% schur_inv %*% t(info_add_beta))
      score_add <- sum(z_add * (m$y - m$n * p_null))
      score_add^2 / info_add
    }, numeric(1))
  }
  data.frame(covariate = x[ind], estimate = unname(beta[ind]), se = NA_real_, statistic = stat,
             p = stats::pchisq(stat, 1, lower.tail = FALSE), lower = NA_real_, upper = NA_real_,
             stringsAsFactors = FALSE)
}
