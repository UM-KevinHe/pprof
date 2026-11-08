#' Fit random-effect and correlated-random-effect provider models
#'
#' Wrappers around [lme4::glmer()] (binary and binomial outcomes) and
#' [lme4::lmer()] (continuous outcomes) that fit a random provider intercept,
#' \deqn{g(E[Y_{ij}]) = \mu + \alpha_i + Z_{ij}^\top\beta, \qquad \alpha_i
#' \sim N(0, \sigma_\alpha^2),}
#' on the providers that [pp_data()] marked as included. The correlated
#' random-effect (CRE) models add the Mundlak within-between decomposition of
#' selected covariates: for each decomposed covariate the provider mean
#' (`<name>_between`) and the deviation from it (`<name>_within`) enter the
#' model instead of the covariate itself, which lets the random intercept be
#' correlated with the covariates' provider means.
#'
#' @param data an object from [pp_data()]: a binary or binomial outcome for
#'   the logistic models, a continuous outcome for the linear models.
#' @param decompose names of the covariates to decompose (CRE models);
#'   default all covariates.
#' @param ... passed to [lme4::glmer()] or [lme4::lmer()] (for instance
#'   `nAGQ` or `control`).
#'
#' @return An object of class `logistic_re`, `logistic_cre` (which inherits
#'   from `logistic_re`), `linear_re` or `linear_cre` (inherits from
#'   `linear_re`):
#' \item{coefficients}{a list with `alpha` (the predicted random effects,
#'   named by provider) and `beta` (the fixed effects, starting with
#'   `(Intercept)`).}
#' \item{variance}{a list with `alpha` (\eqn{\hat\sigma_\alpha^2}) and
#'   `beta` (the covariance matrix of the fixed effects).}
#' \item{post_var}{the conditional variance of each \eqn{\hat\alpha_i} given
#'   the data (lme4's `postVar`), named by provider.}
#' \item{sigma}{the residual standard deviation (linear models).}
#' \item{fitted}{fitted means (probabilities for the logistic models) of the
#'   included rows, including the random effects.}
#' \item{linear_predictor}{the fixed part \eqn{\hat\mu + Z\hat\beta} of the
#'   included rows.}
#' \item{residuals}{response residuals \eqn{Y - N\hat p} or \eqn{Y - \hat
#'   \mu}.}
#' \item{loglik, aic, bic}{from the lme4 fit.}
#' \item{model}{the lme4 fit.}
#' \item{data}{the included rows, with the decomposed columns added (CRE).}
#' \item{vars}{as in [pp_data()]; for the CRE models `x` holds the names of
#'   the model terms and `decompose` the covariates decomposed.}
#' \item{sizes, cutoff, rows}{as in [fit_logistic_fe()].}
#'
#' @references
#' Bates D, Maechler M, Bolker B, Walker S (2015). Fitting linear
#' mixed-effects models using lme4. *Journal of Statistical Software*,
#' 67(1), 1-48.
#'
#' Mundlak Y (1978). On the pooling of time series and cross section data.
#' *Econometrica*, 46(1), 69-85.
#' @seealso [test.logistic_re()], [standardize.logistic_re()],
#'   [confint.logistic_re()], [summary.logistic_re()]
#' @examples
#' \donttest{
#' d <- sim_kecc(n_providers = 20, seed = 1)
#' pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
#'               patient = "patient", message = FALSE)
#' fit <- fit_logistic_re(pd)
#' fit
#' fit_cre <- fit_logistic_cre(pd, decompose = "x1")
#' }
#' @name fit_mixed
NULL

#' @rdname fit_mixed
#' @export
fit_logistic_re <- function(data, ...) {
  fit_mixed(data, family = "binomial", decompose = NULL, ...)
}

#' @rdname fit_mixed
#' @export
fit_logistic_cre <- function(data, decompose = NULL, ...) {
  if (is.null(decompose)) decompose <- data$vars$x
  fit_mixed(data, family = "binomial", decompose = decompose, ...)
}

#' @rdname fit_mixed
#' @export
fit_linear_re <- function(data, ...) {
  fit_mixed(data, family = "gaussian", decompose = NULL, ...)
}

#' @rdname fit_mixed
#' @export
fit_linear_cre <- function(data, decompose = NULL, ...) {
  if (is.null(decompose)) decompose <- data$vars$x
  fit_mixed(data, family = "gaussian", decompose = decompose, ...)
}

fit_mixed <- function(data, family, decompose, ...) {
  outcome <- if (family == "binomial") {
    check_model_data(data, c("binary", "binomial"), "The logistic model", "fit_linear_re()")
  } else {
    check_model_data(data, c("continuous", "binary", "binomial"), "The linear model")
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package lme4 is needed for the random-effect models.", call. = FALSE)
  }
  v <- data$vars
  d <- data$data[data$data$included, , drop = FALSE]
  if (nrow(d) == 0) stop("No provider reaches the cutoff.", call. = FALSE)
  rownames(d) <- NULL
  prov <- d[[v$provider]]
  levels <- unique(prov)
  code <- match(prov, levels)
  n_prov <- as.numeric(table(factor(prov, levels = levels)))

  x_terms <- v$x
  if (!is.null(decompose)) {
    bad <- setdiff(decompose, v$x)
    if (length(bad)) {
      stop("Not covariates: ", paste(bad, collapse = ", "), call. = FALSE)
    }
    for (name in decompose) {
      means <- as.numeric(provider_means(d[[name]], code, n_prov))
      d[[paste0(name, "_between")]] <- means[code]
      d[[paste0(name, "_within")]] <- d[[name]] - means[code]
    }
    x_terms <- c(paste0(decompose, "_within"), paste0(decompose, "_between"),
                 setdiff(v$x, decompose))
  }
  v$x <- x_terms
  v$decompose <- decompose

  d_fit <- d
  d_fit[[v$provider]] <- factor(prov, levels = levels)
  bt <- function(x) paste0("`", x, "`")
  lhs <- if (identical(outcome, "binomial")) {
    paste0("cbind(", bt(v$y), ", ", bt(v$n), " - ", bt(v$y), ")")
  } else {
    bt(v$y)
  }
  f <- stats::as.formula(paste(lhs, "~", paste(bt(x_terms), collapse = " + "), "+ (1 |",
                               bt(v$provider), ")"))
  model <- if (family == "binomial") {
    lme4::glmer(f, data = d_fit, family = stats::binomial(), ...)
  } else {
    lme4::lmer(f, data = d_fit, ...)
  }

  beta <- lme4::fixef(model)
  re <- lme4::ranef(model, condVar = TRUE)[[v$provider]]
  alpha <- stats::setNames(re[[1]], rownames(re))[levels]
  post_var <- stats::setNames(as.numeric(attr(re, "postVar")[1, 1, ]), rownames(re))[levels]
  vc <- lme4::VarCorr(model)
  sigma_alpha2 <- as.numeric(vc[[v$provider]][1, 1])
  vcov_beta <- as.matrix(stats::vcov(model))
  dimnames(vcov_beta) <- list(names(beta), names(beta))
  lp <- as.numeric(stats::model.matrix(model) %*% beta)
  fitted <- as.numeric(stats::fitted(model))
  y <- as.numeric(d[[v$y]])
  trials <- if (is.null(v$n)) 1 else as.numeric(d[[v$n]])

  out <- structure(list(coefficients = list(alpha = alpha, beta = beta),
                        variance = list(alpha = sigma_alpha2, beta = vcov_beta),
                        post_var = post_var,
                        sigma = if (family == "gaussian") stats::sigma(model) else NULL,
                        fitted = fitted, linear_predictor = lp,
                        residuals = y - trials * fitted,
                        loglik = as.numeric(stats::logLik(model)), aic = stats::AIC(model),
                        bic = stats::BIC(model), model = model,
                        data = d, vars = v,
                        sizes = data$sizes[data$sizes$included, , drop = FALSE],
                        cutoff = data$cutoff, rows = n_prov),
                   class = c(if (!is.null(decompose)) paste0(mixed_prefix(family), "_cre"),
                             paste0(mixed_prefix(family), "_re")))
  rownames(out$sizes) <- NULL
  out
}

mixed_prefix <- function(family) if (family == "binomial") "logistic" else "linear"

print_mixed <- function(x, what) {
  cat(what, if (!is.null(x$vars$decompose)) " (correlated random effects)" else "", "\n", sep = "")
  cat("  ", nrow(x$data), " rows, ", length(x$coefficients$alpha), " providers, ",
      length(x$coefficients$beta) - 1, " fixed-effect terms\n", sep = "")
  cat("  provider variance ", format(x$variance$alpha, digits = 4),
      if (!is.null(x$sigma)) paste0("; residual standard deviation ", format(x$sigma, digits = 4)),
      "\n", sep = "")
  cat("  log-likelihood ", format(x$loglik, digits = 6), "\n", sep = "")
  invisible(x)
}

#' @export
print.logistic_re <- function(x, ...) {
  print_mixed(x, "Logistic model with random provider effects")
}

#' @export
print.linear_re <- function(x, ...) {
  print_mixed(x, "Linear model with random provider effects")
}

#' @export
coef.logistic_re <- function(object, ...) {
  object$coefficients
}

#' @export
coef.linear_re <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.logistic_re <- function(object, ...) {
  out <- object$variance$beta
  attr(out, "post_var") <- object$post_var
  out
}

#' @export
vcov.linear_re <- function(object, ...) {
  out <- object$variance$beta
  attr(out, "post_var") <- object$post_var
  out
}

check_mixed <- function(fit) {
  if (!inherits(fit, c("logistic_re", "linear_re"))) {
    stop("`fit` must be a random-effect model.", call. = FALSE)
  }
  invisible(TRUE)
}
