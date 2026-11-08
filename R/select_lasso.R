#' Variable selection with unpenalized provider effects
#'
#' Lasso and group-lasso paths for the logistic model with fixed provider
#' effects, \eqn{\mathrm{logit}\, p_{ij} = \gamma_i + Z_{ij}^\top\beta}, where
#' the provider effects (and any covariate named in `penalize`'s complement)
#' are never penalized. The objective at penalty \eqn{\lambda} is
#' \deqn{-\frac{1}{n}\,\ell(\gamma, \beta) + \lambda \sum_j m_j |\beta_j|
#' \qquad \text{or} \qquad -\frac{1}{n}\,\ell(\gamma, \beta) + \lambda
#' \sum_g m_g \|\beta_g\|_2,}
#' solved by coordinate descent on standardized covariates (each penalized
#' group orthonormalized) with one Newton step on the provider effects per
#' cycle, warm starts along a decreasing sequence of \eqn{\lambda}, and an
#' active-set strategy.
#'
#' @param data an object from [pp_data()] with a binary, binomial or survival
#'   outcome.
#' @param penalty `"lasso"` (default) or `"group_lasso"`.
#' @param penalize names of the covariates that carry a penalty; default all.
#'   The others enter the model unpenalized.
#' @param group for the group lasso, the group of each covariate: a vector of
#'   labels in the order of the covariates, or named by covariate. The label
#'   `0`, `NA` or `""` marks an unpenalized covariate.
#' @param multiplier penalty multipliers \eqn{m_j} (named by covariate) or
#'   \eqn{m_g} (named by group); default 1 per covariate and the square root
#'   of the group size per group.
#' @param lambda the penalty sequence; by default `n_lambda` values from the
#'   smallest \eqn{\lambda} with every penalized coefficient zero down to
#'   `lambda_min_ratio` times it, equally spaced on the log scale.
#' @param n_lambda,lambda_min_ratio length and range of the default sequence.
#' @param standardize centre and scale the covariates before penalizing (the
#'   estimates are returned on the original scale). Default `TRUE`.
#' @param bound the provider effects are kept within `median +/- bound`.
#' @param tol convergence tolerance on the largest change of a coefficient
#'   within a cycle. Default 1e-4.
#' @param max_iter maximum number of cycles per value of `lambda`.
#' @param active_set start each `lambda` from the coefficients that are
#'   nonzero and add the covariates that violate the optimality conditions.
#' @param mm for the lasso, use the majorization with the constant working
#'   weight \eqn{N/4} instead of the IRLS weights (the group lasso always
#'   does). Default `FALSE`.
#' @param backtrack line search in the provider-effect step. Default `FALSE`.
#' @param threads OpenMP threads for the provider-effect step.
#' @param model for a survival outcome: `"cox"` (default; the Cox model
#'   stratified by provider, whose baselines absorb the provider effects) or
#'   `"discrete"` (the discrete-time hazard model, lasso only). `NULL` for
#'   the other outcomes.
#'
#' @details Providers with fewer than the cutoff patients are excluded as in
#' every fitter. The path starts from the fit with the unpenalized covariates
#' only (or the provider-only model), whose residuals give the largest
#' \eqn{\lambda}. The deviance reported per \eqn{\lambda} is \eqn{-2\ell}
#' without the constant term.
#'
#' With a survival outcome ([pp_data()] with `time` and `event`) the paths
#' are those of `grplasso`'s `Strat.cox()` and `pp.DiscSurv()`:
#' * `model = "cox"`: the weighted partial likelihood of the Cox model with
#'   one baseline per provider, counting-process rows `(start, stop]`, case
#'   weights and offsets as in [fit_cox()]; the provider effects are not
#'   parameters (`gamma` is `NULL`) and the working weight is 1 (Breheny and
#'   Huang's majorization for the Cox model). The covariates of a group are
#'   orthonormalized as for the logistic model.
#' * `model = "discrete"`: \eqn{\mathrm{logit}\, h_{jk} = \alpha_k + \gamma_i
#'   + Z_j^\top\beta} for subject \eqn{j} of provider \eqn{i} in period
#'   \eqn{k}, the periods being the distinct values of `time`; one row per
#'   subject (no `start`, no weights), the subject at risk in every period up
#'   to its own. The period effects \eqn{\alpha} and the provider effects
#'   \eqn{\gamma} are unpenalized, \eqn{\gamma} identified by a zero median
#'   (`grplasso` fixed the first provider at zero) and kept within
#'   `+/- bound`; the likelihood is computed without expanding the data. The
#'   IRLS weights \eqn{\sum_k h_{jk}(1 - h_{jk})} or, with `mm = TRUE`,
#'   \eqn{t_j / 4}.
#'
#' @return An object of class `pp_lasso` (and `pp_group_lasso` for the group
#'   lasso):
#' \item{beta}{a matrix of covariate coefficients, one column per `lambda`.}
#' \item{gamma}{a matrix of provider effects, one column per `lambda`
#'   (`NULL` for the stratified Cox model).}
#' \item{alpha, time_grid}{the discrete-time model's period effects (one
#'   column per `lambda`) and the periods' times.}
#' \item{model}{`"logistic"`, `"cox"` or `"discrete"`.}
#' \item{lambda, deviance, iter}{per value of `lambda` (the deviance of the
#'   Cox model is minus twice its partial log-likelihood).}
#' \item{n_selected}{number of penalized covariates with a nonzero
#'   coefficient per `lambda`.}
#' \item{penalized}{logical, per covariate.}
#' \item{group}{the group of each covariate (group lasso).}
#' \item{multiplier}{the multipliers used.}
#' \item{data, vars, sizes, rows}{as in [fit_logistic_fe()].}
#'
#' @references
#' Breheny P, Huang J (2015). Group descent algorithms for nonconvex
#' penalized linear and logistic regression models with grouped predictors.
#' *Statistics and Computing*, 25, 173-187.
#'
#' Yang Y, Wu W, He K (2025+). Variable selection for provider profiling
#' with high-dimensional fixed effects (the `grplasso` package).
#' @seealso [cv()], [coef.pp_lasso()], [predict.pp_lasso()]
#' @examples
#' d <- sim_kecc(n_providers = 20, n_covariates = 6, seed = 1)
#' pd <- pp_data(d, y = "y", x = paste0("x", 1:6), provider = "provider",
#'               patient = "patient", message = FALSE)
#' path <- var_select(pd, n_lambda = 20)
#' path
#' coef(path, lambda = path$lambda[10])$beta
#' grp <- var_select(pd, penalty = "group_lasso", group = c(1, 1, 2, 2, 3, 0), n_lambda = 20)
#' @name var_select
NULL

#' @rdname var_select
#' @export
var_select <- function(data, penalty = c("lasso", "group_lasso"), group = NULL, penalize = NULL,
                       multiplier = NULL, lambda = NULL, n_lambda = 100, lambda_min_ratio = 1e-4,
                       standardize = TRUE, bound = 10, tol = 1e-4, max_iter = 10000,
                       active_set = TRUE, mm = FALSE, backtrack = FALSE, threads = 1,
                       model = NULL) {
  penalty <- match.arg(penalty)
  if (penalty == "group_lasso" && is.null(group)) {
    stop("penalty = \"group_lasso\" needs `group`.", call. = FALSE)
  }
  if (penalty == "lasso" && !is.null(group)) {
    stop("`group` is only used with penalty = \"group_lasso\".", call. = FALSE)
  }
  out <- if (penalty == "lasso") {
    fit_lasso_path(data, penalize = penalize, multiplier = multiplier, lambda = lambda,
                   n_lambda = n_lambda, lambda_min_ratio = lambda_min_ratio,
                   standardize = standardize, bound = bound, tol = tol, max_iter = max_iter,
                   active_set = active_set, mm = mm, backtrack = backtrack, threads = threads,
                   model = model)
  } else {
    fit_group_lasso_path(data, group = group, multiplier = multiplier, lambda = lambda,
                         n_lambda = n_lambda, lambda_min_ratio = lambda_min_ratio,
                         standardize = standardize, bound = bound, tol = tol,
                         max_iter = max_iter, active_set = active_set, backtrack = backtrack,
                         threads = threads, model = model)
  }
  out$settings <- list(penalty = penalty, group = group, penalize = penalize,
                       multiplier = multiplier, standardize = standardize, bound = bound, tol = tol,
                       max_iter = max_iter, active_set = active_set, mm = mm,
                       backtrack = backtrack, threads = threads, model = model)
  out$cutoff <- data$cutoff
  out
}

fit_lasso_path <- function(data, penalize = NULL, multiplier = NULL, lambda = NULL, n_lambda = 100,
                         lambda_min_ratio = 1e-4, standardize = TRUE, bound = 10, tol = 1e-4,
                         max_iter = 10000, active_set = TRUE, mm = FALSE, backtrack = FALSE,
                         threads = 1, model = NULL) {
  model <- survival_model(data, model)
  if (identical(model, "cox")) {
    return(cox_lasso(data, penalize = penalize, multiplier = multiplier, lambda = lambda,
                     n_lambda = n_lambda, lambda_min_ratio = lambda_min_ratio,
                     standardize = standardize, tol = tol, max_iter = max_iter,
                     active_set = active_set))
  }
  if (identical(model, "discrete")) {
    return(discrete_lasso(data, penalize = penalize, multiplier = multiplier, lambda = lambda,
                          n_lambda = n_lambda, lambda_min_ratio = lambda_min_ratio,
                          standardize = standardize, bound = bound, tol = tol,
                          max_iter = max_iter, active_set = active_set, mm = mm,
                          backtrack = backtrack))
  }
  m <- selection_matrices(data)
  x <- data$vars$x
  if (is.null(penalize)) penalize <- x
  bad <- setdiff(penalize, x)
  if (length(bad)) stop("Not covariates: ", paste(bad, collapse = ", "), call. = FALSE)
  penalized <- x %in% penalize
  mult <- rep(1, length(x))
  if (!is.null(multiplier)) {
    if (is.null(names(multiplier)) || !all(names(multiplier) %in% penalize)) {
      stop("`multiplier` must be named by penalized covariate.", call. = FALSE)
    }
    mult[match(names(multiplier), x)] <- as.numeric(multiplier)
  }
  std <- standardize_columns(m$z, standardize)
  start <- path_start(data, m, std, which(!penalized))
  if (is.null(lambda)) {
    score <- abs(as.numeric(crossprod(std$z[, penalized, drop = FALSE], start$residual)))
    lambda_max <- max(score / (mult[penalized] * m$n_obs))
    lambda <- lambda_sequence(lambda_max, n_lambda, lambda_min_ratio)
  } else {
    lambda <- sort(as.numeric(lambda), decreasing = TRUE)
  }
  res <- lasso_path(m$n, m$y, std$z, m$n_prov, start$gamma, start$beta, as.integer(penalized),
                    mult, lambda, mm, backtrack, bound, tol, max(1L, as.integer(max_iter)),
                    active_set, max(1L, as.integer(threads)))
  beta <- res$beta / std$scale
  out <- path_object(data, m, std, beta, res, lambda, penalized, mult,
                     control = list(standardize = standardize, bound = bound, tol = tol,
                                    max_iter = max_iter, active_set = active_set, mm = mm,
                                    backtrack = backtrack, threads = threads))
  class(out) <- "pp_lasso"
  out
}

fit_group_lasso_path <- function(data, group, multiplier = NULL, lambda = NULL, n_lambda = 100,
                               lambda_min_ratio = 1e-4, standardize = TRUE, bound = 10,
                               tol = 1e-4, max_iter = 10000, active_set = TRUE,
                               backtrack = FALSE, threads = 1, model = NULL) {
  model <- survival_model(data, model)
  if (identical(model, "discrete")) {
    stop("The discrete-time model has no group lasso; use model = \"cox\".", call. = FALSE)
  }
  if (identical(model, "cox")) {
    return(cox_lasso(data, group = group, multiplier = multiplier, lambda = lambda,
                     n_lambda = n_lambda, lambda_min_ratio = lambda_min_ratio,
                     standardize = standardize, tol = tol, max_iter = max_iter,
                     active_set = active_set))
  }
  m <- selection_matrices(data)
  x <- data$vars$x
  if (!is.null(names(group))) {
    if (!setequal(names(group), x)) stop("`group` must name every covariate.", call. = FALSE)
    group <- group[x]
  } else if (length(group) != length(x)) {
    stop("`group` must have one label per covariate.", call. = FALSE)
  }
  label <- as.character(group)
  unpen <- is.na(label) | label %in% c("0", "")
  levels_g <- unique(label[!unpen])
  if (length(levels_g) == 0) stop("No penalized group.", call. = FALSE)
  g_index <- match(label, levels_g) # NA for the unpenalized covariates
  size <- as.numeric(table(factor(g_index, levels = seq_along(levels_g))))
  mult <- sqrt(size)
  if (!is.null(multiplier)) {
    if (is.null(names(multiplier)) || !all(names(multiplier) %in% levels_g)) {
      stop("`multiplier` must be named by penalized group.", call. = FALSE)
    }
    mult[match(names(multiplier), levels_g)] <- as.numeric(multiplier)
  }
  std <- standardize_columns(m$z, standardize)
  ortho <- orthonormalize_groups(std$z, g_index, levels_g, m$n_obs)
  start <- path_start(data, m, std, which(unpen))
  beta0 <- ortho$to_ortho(start$beta)
  if (is.null(lambda)) {
    score <- vapply(seq_along(levels_g), function(g) {
      cols <- ortho$col_first[g]:ortho$col_last[g]
      sqrt(sum(crossprod(ortho$z[, cols, drop = FALSE], start$residual)^2))
    }, numeric(1))
    lambda_max <- max(score / (mult * m$n_obs))
    lambda <- lambda_sequence(lambda_max, n_lambda, lambda_min_ratio)
  } else {
    lambda <- sort(as.numeric(lambda), decreasing = TRUE)
  }
  res <- group_lasso_path(m$n, m$y, ortho$z, m$n_prov, start$gamma, beta0,
                          seq_along(ortho$unpenalized) - 1L, as.integer(ortho$col_first - 1),
                          as.integer(ortho$col_last - 1), mult, lambda, backtrack, bound, tol,
                          max(1L, as.integer(max_iter)), active_set, max(1L, as.integer(threads)))
  beta <- ortho$from_ortho(res$beta) / std$scale
  out <- path_object(data, m, std, beta, res, lambda, !unpen, stats::setNames(mult, levels_g),
                     control = list(standardize = standardize, bound = bound, tol = tol,
                                    max_iter = max_iter, active_set = active_set, mm = TRUE,
                                    backtrack = backtrack, threads = threads))
  out$group <- stats::setNames(ifelse(unpen, NA_character_, label), x)
  class(out) <- c("pp_group_lasso", "pp_lasso")
  out
}

# The survival model of a selection path ("cox" or "discrete"), or NULL for
# the other outcomes.
survival_model <- function(data, model) {
  if (!inherits(data, "pp_data")) stop("`data` must come from pp_data().", call. = FALSE)
  if (outcome_of(data) != "survival") {
    if (!is.null(model)) stop("`model` is only used with a survival outcome.", call. = FALSE)
    return(NULL)
  }
  if (is.null(model)) return("cox")
  match.arg(model, c("cox", "discrete"))
}

# Included rows of a pp_data object as the matrices the paths read.
selection_matrices <- function(data) {
  check_model_data(data, c("binary", "binomial"), "Variable selection")
  v <- data$vars
  d <- data$data[data$data$included, , drop = FALSE]
  if (nrow(d) == 0) stop("No provider reaches the cutoff.", call. = FALSE)
  rownames(d) <- NULL
  prov <- d[[v$provider]]
  levels <- unique(prov)
  z <- as.matrix(d[, v$x, drop = FALSE])
  storage.mode(z) <- "double"
  list(d = d, y = as.numeric(d[[v$y]]),
       n = if (is.null(v$n)) rep(1, nrow(d)) else as.numeric(d[[v$n]]),
       z = z, prov = prov, levels = levels, code = match(prov, levels),
       n_prov = as.numeric(table(factor(prov, levels = levels))), n_obs = nrow(d),
       sizes = data$sizes[data$sizes$included, , drop = FALSE])
}

# Centre and scale the columns (population standard deviation); constant
# columns are an error.
standardize_columns <- function(z, standardize) {
  p <- ncol(z)
  if (!standardize) {
    return(list(z = z, center = rep(0, p), scale = rep(1, p)))
  }
  center <- colMeans(z)
  zc <- sweep(z, 2, center)
  scale <- sqrt(colMeans(zc^2))
  if (any(scale <= 1e-10)) {
    stop("Constant covariate(s): ", paste(colnames(z)[scale <= 1e-10], collapse = ", "),
         call. = FALSE)
  }
  list(z = sweep(zc, 2, scale, "/"), center = center, scale = scale)
}

# Starting values of the path: the fit with the unpenalized covariates only
# (or the provider-only model), expressed in the standardized coordinates, and
# its residuals Y - N p.
path_start <- function(data, m, std, unpenalized) {
  p <- ncol(m$z)
  beta <- rep(0, p)
  if (length(unpenalized)) {
    sub <- structure(list(data = m$d, vars = utils::modifyList(data$vars,
                                                               list(x = data$vars$x[unpenalized])),
                          sizes = m$sizes, cutoff = data$cutoff), class = "pp_data")
    fit <- fit_logistic_fe(sub, tol = 1e-8)
    beta_orig <- unname(fit$coefficients$beta)
    beta[unpenalized] <- beta_orig * std$scale[unpenalized]
    gamma <- unname(fit$coefficients$gamma) + sum(std$center[unpenalized] * beta_orig)
    residual <- m$y - m$n * fit$fitted
  } else {
    rate <- as.numeric(rowsum(m$y, m$code)) / as.numeric(rowsum(m$n, m$code))
    gamma <- stats::qlogis(rate)
    gamma <- pmin(pmax(gamma, stats::median(gamma[is.finite(gamma)]) - 10),
                  stats::median(gamma[is.finite(gamma)]) + 10)
    residual <- m$y - m$n * rate[m$code]
  }
  list(gamma = gamma, beta = beta, residual = residual)
}

lambda_sequence <- function(lambda_max, n_lambda, lambda_min_ratio) {
  if (n_lambda < 2) stop("`n_lambda` must be at least 2.", call. = FALSE)
  # nudged up so that every penalized coefficient is exactly zero at the first
  # penalty despite the rounding of the score computed in R and in C++
  lambda_max <- lambda_max * (1 + 1e-9)
  exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio), length.out = n_lambda))
}

# Orthonormalize each penalized group of standardized columns (Z_g' Z_g / n =
# I), dropping directions with a zero singular value; returns the reordered
# matrix (unpenalized columns first, then the groups), the column ranges of
# the groups and the maps between the two coordinate systems.
orthonormalize_groups <- function(z, g_index, levels_g, n) {
  unpen <- which(is.na(g_index))
  blocks <- list()
  q_list <- list()
  col_first <- col_last <- integer(length(levels_g))
  pos <- length(unpen)
  for (g in seq_along(levels_g)) {
    cols <- which(g_index == g)
    s <- svd(z[, cols, drop = FALSE], nu = 0)
    keep <- which(s$d > 1e-10 * max(s$d))
    q <- s$v[, keep, drop = FALSE] * rep(sqrt(n) / s$d[keep], each = length(cols))
    blocks[[g]] <- z[, cols, drop = FALSE] %*% q
    q_list[[g]] <- list(cols = cols, q = q)
    col_first[g] <- pos + 1
    col_last[g] <- pos + length(keep)
    pos <- pos + length(keep)
  }
  z_new <- cbind(z[, unpen, drop = FALSE], do.call(cbind, blocks))
  to_ortho <- function(beta) {
    # only used for the start, where every penalized coefficient is zero
    c(beta[unpen], rep(0, pos - length(unpen)))
  }
  from_ortho <- function(b) {
    out <- matrix(0, ncol(z), ncol(b))
    out[unpen, ] <- b[seq_along(unpen), , drop = FALSE]
    for (g in seq_along(levels_g)) {
      out[q_list[[g]]$cols, ] <- q_list[[g]]$q %*% b[col_first[g]:col_last[g], , drop = FALSE]
    }
    out
  }
  list(z = z_new, unpenalized = unpen, col_first = col_first, col_last = col_last,
       to_ortho = to_ortho, from_ortho = from_ortho)
}

# The common path object; `beta` is already on the original scale of the
# covariates but still refers to the centred covariates.
path_object <- function(data, m, std, beta, res, lambda, penalized, multiplier, control) {
  x <- data$vars$x
  gamma <- res$gamma - matrix(as.numeric(crossprod(std$center, beta)), nrow(res$gamma),
                              ncol(res$gamma), byrow = TRUE)
  dimnames(beta) <- list(x, signif(lambda, 6))
  dimnames(gamma) <- list(m$levels, signif(lambda, 6))
  list(beta = beta, gamma = gamma, alpha = NULL, time_grid = NULL, lambda = lambda,
       deviance = as.numeric(res$deviance), iter = as.integer(res$iter),
       n_selected = as.integer(colSums(beta[penalized, , drop = FALSE] != 0)),
       penalized = stats::setNames(penalized, x), multiplier = multiplier,
       data = m$d, vars = data$vars, sizes = m$sizes, rows = m$n_prov, control = control,
       model = "logistic")
}

#' @export
print.pp_lasso <- function(x, ...) {
  cat(if (inherits(x, "pp_group_lasso")) "Group-lasso" else "Lasso",
      switch(path_model(x),
             logistic = " path of the logistic model with unpenalized provider effects\n",
             cox = " path of the Cox model stratified by provider\n",
             discrete = " path of the discrete-time hazard model (period and provider effects)\n"),
      sep = "")
  cat("  ", nrow(x$data), " rows, ", nrow(x$sizes), " providers, ", nrow(x$beta),
      " covariates (", sum(x$penalized), " penalized)",
      if (!is.null(x$alpha)) paste0(", ", nrow(x$alpha), " periods"), "\n", sep = "")
  cat("  ", length(x$lambda), " values of lambda from ", format(max(x$lambda), digits = 4),
      " to ", format(min(x$lambda), digits = 4), "; selected covariates from ",
      min(x$n_selected), " to ", max(x$n_selected), "\n", sep = "")
  invisible(x)
}

#' Coefficients and predictions along a selection path
#'
#' @param object a path from [var_select()] or [var_select()].
#' @param lambda penalty values at which to report; values between two points
#'   of the path are interpolated linearly.
#' @param which alternatively, indices into `object$lambda`; default all.
#' @param drop return vectors when a single value is requested.
#' @param newdata a data frame with the covariate and provider columns of the
#'   fit; default the fitted rows.
#' @param type `"link"` (the linear predictor, including any offset) or
#'   `"response"`: the probability for the logistic model, the relative
#'   hazard \eqn{e^\eta} for the Cox model, and for the discrete-time model
#'   the probability of the event by the row's time, \eqn{1 - \prod_k (1 -
#'   h_{jk})} over the fit's periods up to that time (a time between two
#'   periods counts up to the earlier one, a time before the first period
#'   counts the first).
#' @param unseen for providers of `newdata` that are not in the fit:
#'   `"missing"` (default; `NA`) or `"median"` (the median provider effect).
#'   Not used by the Cox model, which has no provider parameters.
#' @param ... unused.
#' @return `coef()`: a list with `gamma` (`NULL` for the Cox model) and
#'   `beta`, plus `alpha` for the discrete-time model (matrices, or vectors
#'   with `drop = TRUE` and one value of `lambda`). `predict()`: a matrix with one
#'   column per requested `lambda`. `selected_covariates()`: the names of the
#'   penalized covariates with a nonzero coefficient, as a list over
#'   `lambda`.
#' @examples
#' d <- sim_kecc(n_providers = 20, n_covariates = 6, seed = 1)
#' pd <- pp_data(d, y = "y", x = paste0("x", 1:6), provider = "provider",
#'               patient = "patient", message = FALSE)
#' path <- var_select(pd, n_lambda = 20)
#' selected_covariates(path, which = c(5, 15))
#' head(predict(path, which = 15, type = "response"))
#' @name coef.pp_lasso
NULL

#' @rdname coef.pp_lasso
#' @export
coef.pp_lasso <- function(object, lambda = NULL, which = NULL, drop = TRUE, ...) {
  if (!is.null(lambda)) {
    rng <- range(object$lambda)
    if (any(lambda > rng[2] | lambda < rng[1])) {
      stop("`lambda` must lie within the fitted path.", call. = FALSE)
    }
    ind <- stats::approx(object$lambda, seq_along(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    pick <- function(mat) {
      if (is.null(mat)) return(NULL)
      out <- mat[, l, drop = FALSE] * (1 - w) + mat[, r, drop = FALSE] * w
      colnames(out) <- signif(lambda, 6)
      out
    }
  } else {
    if (is.null(which)) which <- seq_along(object$lambda)
    pick <- function(mat) if (is.null(mat)) NULL else mat[, which, drop = FALSE]
  }
  beta <- pick(object$beta)
  gamma <- pick(object$gamma)
  alpha <- pick(object$alpha)
  if (drop && ncol(beta) == 1) {
    beta <- beta[, 1]
    if (!is.null(gamma)) gamma <- gamma[, 1]
    if (!is.null(alpha)) alpha <- alpha[, 1]
  }
  out <- list(gamma = gamma, beta = beta)
  if (!is.null(alpha)) out$alpha <- alpha
  out
}

#' @rdname coef.pp_lasso
#' @export
predict.pp_lasso <- function(object, newdata = NULL, lambda = NULL, which = NULL,
                             type = c("link", "response"), unseen = c("missing", "median"),
                             ...) {
  type <- match.arg(type)
  unseen <- match.arg(unseen)
  cf <- coef(object, lambda = lambda, which = which, drop = FALSE)
  v <- object$vars
  if (is.null(newdata)) newdata <- object$data
  missing_cols <- setdiff(c(v$x, v$provider), names(newdata))
  if (length(missing_cols)) {
    stop("Column(s) not in `newdata`: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  z <- as.matrix(newdata[, v$x, drop = FALSE])
  storage.mode(z) <- "double"
  model <- path_model(object)
  offset <- if (!is.null(v$offset) && v$offset %in% names(newdata)) {
    as.numeric(newdata[[v$offset]])
  } else {
    0
  }
  eta <- z %*% cf$beta + offset
  if (model != "cox") {
    pos <- match(as.character(newdata[[v$provider]]), rownames(cf$gamma))
    g <- cf$gamma[ifelse(is.na(pos), 1L, pos), , drop = FALSE]
    if (anyNA(pos)) {
      g[is.na(pos), ] <- if (unseen == "median") {
        rep(apply(cf$gamma, 2, stats::median), each = sum(is.na(pos)))
      } else {
        NA_real_
      }
    }
    eta <- eta + g
  }
  dimnames(eta) <- list(NULL, colnames(cf$beta))
  if (type == "link") return(eta)
  if (model == "logistic") return(stats::plogis(eta))
  if (model == "cox") return(exp(eta))
  ex <- discrete_expanded(object, newdata, lambda = lambda, which = which)
  prob <- 1 - exp(rowsum(log1p(-ex$hazard), ex$row, reorder = TRUE))
  dimnames(prob) <- dimnames(eta)
  prob
}

#' @rdname coef.pp_lasso
#' @export
selected_covariates <- function(object, lambda = NULL, which = NULL) {
  cf <- coef(object, lambda = lambda, which = which, drop = FALSE)
  pen <- object$penalized
  out <- lapply(seq_len(ncol(cf$beta)), function(k) {
    names(pen)[pen & cf$beta[, k] != 0]
  })
  names(out) <- colnames(cf$beta)
  out
}
