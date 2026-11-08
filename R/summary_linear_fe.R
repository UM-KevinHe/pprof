#' Tests of the covariate coefficients of a linear fixed-effect fit
#'
#' @param object a model fitted by [fit_linear_fe()].
#' @param test `"wald"` (default; model-based t test with \eqn{n - m - p}
#'   degrees of freedom) or `"robust_wald"` (patient-clustered sandwich,
#'   standard normal reference; needs a patient column).
#' @param covariates optional subset of covariate names.
#' @param level confidence level for the intervals.
#' @param null the null value of each coefficient.
#' @param ... unused.
#'
#' @return A data frame with one row per covariate: `covariate`, `estimate`,
#'   `se`, `statistic`, `p`, `lower`, `upper`.
#' @examples
#' d <- sim_kecc(n_providers = 20, type = "continuous", seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' object <- fit_linear_fe(pd)
#' summary(object)
#' @export
summary.linear_fe <- function(object, test = c("wald", "robust_wald"), covariates = NULL,
                                      level = 0.95, null = 0, ...) {
  check_fit(object, "linear_fe")
  test <- match.arg(test)
  check_level(level)
  x <- object$vars$x
  ind <- if (is.null(covariates)) seq_along(x) else match(covariates, x)
  if (anyNA(ind)) {
    stop("Unknown covariate(s): ", paste(covariates[is.na(ind)], collapse = ", "), call. = FALSE)
  }
  beta <- object$coefficients$beta
  if (test == "wald") {
    v <- object$vcov_beta
    df <- object$df_residual
  } else {
    m <- fit_matrices(object)
    if (is.null(m$patient)) {
      stop("test = \"robust_wald\" needs a patient column in pp_data().", call. = FALSE)
    }
    code <- provider_codes(m$prov, m$levels)
    z_w <- m$z - provider_means(m$z, code, m$n_prov)[code, , drop = FALSE]
    cl <- cluster_codes(m$prov, m$patient)
    u <- rowsum(z_w * object$residuals, cl$code, reorder = FALSE)
    bread <- object$vcov_beta / object$sigma^2
    v <- bread %*% crossprod(u) %*% bread
    df <- Inf
  }
  se <- sqrt(diag(v))[ind]
  stat <- (beta[ind] - null) / se
  decision <- z_decision(stat, level, "two.sided", df)
  limits <- wald_limits(beta[ind], se, level, "two.sided", df)
  pp_summary(data.frame(covariate = x[ind], estimate = unname(beta[ind]), se = unname(se),
             statistic = unname(stat), p = decision$p, lower = limits$lower,
             upper = limits$upper, stringsAsFactors = FALSE))
}
