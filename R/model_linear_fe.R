#' Fit the linear model with fixed provider effects
#'
#' Fits \eqn{Y_{ij} = \gamma_i + Z_{ij}^\top\beta + \epsilon_{ij}} by least
#' squares on the providers that [pp_data()] marked as included. The
#' covariate coefficients come from the within-provider regression (every
#' variable centred at its provider mean), so no provider dummy variables
#' are formed; the provider effects are \eqn{\hat\gamma_i = \bar Y_i -
#' \bar Z_i^\top\hat\beta}.
#'
#' @param data an object from [pp_data()] with a continuous (or any numeric)
#'   outcome.
#'
#' @return An object of class `linear_fe`:
#' \item{coefficients}{a list with `gamma` (provider effects, named by
#'   provider) and `beta` (covariate coefficients, named by covariate).}
#' \item{sigma}{the residual standard deviation, \eqn{\hat\sigma^2 = RSS /
#'   (n - m - p)}.}
#' \item{fitted, linear_predictor, residuals}{per included row;
#'   `linear_predictor` is \eqn{Z\hat\beta}.}
#' \item{loglik, aic, bic}{Gaussian log-likelihood at the maximum-likelihood
#'   variance and the information criteria with \eqn{m + p + 1} parameters.}
#' \item{df_residual}{\eqn{n - m - p}.}
#' \item{vcov_beta}{\eqn{\hat\sigma^2 (\tilde Z^\top \tilde Z)^{-1}} with
#'   \eqn{\tilde Z} the within-centred covariates.}
#' \item{data, vars, sizes, cutoff, rows}{as in [fit_logistic_fe()].}
#'
#' @seealso [test_providers.linear_fe()], [test_covariates.linear_fe()],
#'   [standardize.linear_fe()], [confint.linear_fe()], [gamma_variance()]
#' @examples
#' d <- sim_kecc(n_providers = 20, type = "continuous", seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_linear_fe(pd)
#' fit
#' @export
fit_linear_fe <- function(data) {
  if (!inherits(data, "pp_data")) stop("`data` must come from pp_data().", call. = FALSE)
  v <- data$vars
  d <- data$data[data$data$included, , drop = FALSE]
  if (nrow(d) == 0) stop("No provider reaches the cutoff.", call. = FALSE)
  rownames(d) <- NULL
  prov <- d[[v$provider]]
  levels <- unique(prov)
  code <- match(prov, levels)
  n_prov <- as.numeric(table(factor(prov, levels = levels)))
  y <- as.numeric(d[[v$y]])
  z <- as.matrix(d[, v$x, drop = FALSE])
  storage.mode(z) <- "double"
  n <- length(y)
  m <- length(levels)
  p <- ncol(z)
  if (n - m - p < 1) {
    stop("Not enough rows to estimate the residual variance.", call. = FALSE)
  }

  y_bar <- as.numeric(provider_means(y, code, n_prov))
  z_bar <- provider_means(z, code, n_prov)
  y_w <- y - y_bar[code]
  z_w <- z - z_bar[code, , drop = FALSE]
  xtx <- crossprod(z_w)
  xtx_inv <- tryCatch(chol2inv(chol(xtx)), error = function(e) NULL)
  if (is.null(xtx_inv)) {
    stop("The covariates are collinear within providers; drop a covariate.", call. = FALSE)
  }
  beta <- as.numeric(xtx_inv %*% crossprod(z_w, y_w))
  beta <- stats::setNames(beta, v$x)
  gamma <- stats::setNames(y_bar - as.numeric(z_bar %*% beta), levels)
  lp <- as.numeric(z %*% beta)
  fitted <- gamma[code] + lp
  residuals <- y - unname(fitted)
  df_residual <- n - m - p
  sigma2 <- sum(residuals^2) / df_residual
  vcov_beta <- sigma2 * xtx_inv
  dimnames(vcov_beta) <- list(v$x, v$x)
  loglik <- -n / 2 * (log(2 * pi) + log(sum(residuals^2) / n) + 1)
  k <- m + p + 1

  out <- structure(list(coefficients = list(gamma = gamma, beta = beta), sigma = sqrt(sigma2),
                        fitted = unname(fitted), linear_predictor = lp, residuals = residuals,
                        loglik = loglik, aic = -2 * loglik + 2 * k, bic = -2 * loglik + log(n) * k,
                        df_residual = df_residual, vcov_beta = vcov_beta,
                        data = d, vars = v,
                        sizes = data$sizes[data$sizes$included, , drop = FALSE],
                        cutoff = data$cutoff, rows = n_prov),
                   class = "linear_fe")
  rownames(out$sizes) <- NULL
  out
}

#' @export
print.linear_fe <- function(x, ...) {
  cat("Linear model with fixed provider effects\n")
  cat("  ", nrow(x$data), " rows, ", length(x$coefficients$gamma), " providers, ",
      length(x$coefficients$beta), " covariates\n", sep = "")
  cat("  residual standard deviation ", format(x$sigma, digits = 4), " on ", x$df_residual,
      " degrees of freedom\n", sep = "")
  cat("  provider effects: median ", format(stats::median(x$coefficients$gamma), digits = 4),
      ", range ", paste(format(range(x$coefficients$gamma), digits = 4), collapse = " to "),
      "\n", sep = "")
  invisible(x)
}

#' @export
coef.linear_fe <- function(object, ...) {
  object$coefficients
}

#' Model-based variance of a linear fixed-effect fit
#'
#' @param object a model fitted by [fit_linear_fe()].
#' @param ... unused.
#' @return The covariance matrix of `beta`, with the `"full"` variances of
#'   the provider effects (see [gamma_variance()]) as the attribute
#'   `"var_gamma"`.
#' @export
vcov.linear_fe <- function(object, ...) {
  out <- object$vcov_beta
  attr(out, "var_gamma") <- gamma_variance(object, type = "full")
  out
}

#' @rdname gamma_variance
#' @export
gamma_variance.linear_fe <- function(fit, type = c("case_mix", "full"), ...) {
  type <- match.arg(type)
  m <- fit_matrices(fit)
  sigma2 <- fit$sigma^2
  if (type == "case_mix") return(stats::setNames(sigma2 / m$n_prov, m$levels))
  code <- provider_codes(m$prov, m$levels)
  z_bar <- provider_means(m$z, code, m$n_prov)
  xtx_inv <- fit$vcov_beta / sigma2
  stats::setNames(sigma2 * (1 / m$n_prov + rowSums((z_bar %*% xtx_inv) * z_bar)), m$levels)
}
