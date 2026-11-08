#' Cross-validation of a selection path
#'
#' Chooses the penalty of a [var_select()] path by
#' K-fold cross-validation of the binomial deviance (for a survival outcome,
#' the partial-likelihood deviance per event of the held-out rows within
#' their own strata and risk sets for the Cox model, the binomial deviance
#' per person-period for the discrete-time model). Folds are formed within
#' every provider over its patients (rows when there is no patient column),
#' so that each provider is present in every training set whenever it has at
#' least two units; a provider held out entirely is predicted with the median
#' provider effect of the training fit.
#'
#' @param object a path from [var_select()]; every fold refits it with the
#'   same settings and penalties.
#' @param ... unused.
#' @param n_folds number of folds. Default 10.
#' @param folds optional fold membership (integers from 1) of the rows of
#'   `data$data`; rows of excluded providers are ignored.
#' @param seed optional seed for the fold assignment; the random-number
#'   state of the session is restored afterwards.
#'
#' @return An object of class `pp_cv`:
#' \item{lambda}{the penalty sequence of the full-data path.}
#' \item{cve, cvse}{the cross-validated deviance per row, event or
#'   person-period (mean over folds of
#'   the fold means) and its standard error across folds.}
#' \item{lambda_min, lambda_1se}{the penalty with the smallest `cve` and the
#'   largest penalty whose `cve` is within one standard error of it.}
#' \item{selected_min, selected_1se}{the covariates selected at those
#'   penalties.}
#' \item{fit}{the path on the full data.}
#' \item{folds}{the fold of each row of the fitted data.}
#' @examples
#' \donttest{
#' d <- sim_kecc(n_providers = 20, n_covariates = 6, seed = 1)
#' pd <- pp_data(d, y = "y", x = paste0("x", 1:6), provider = "provider",
#'               patient = "patient", message = FALSE)
#' cvp <- cv(var_select(pd, n_lambda = 20), n_folds = 5, seed = 1)
#' cvp
#' cvp$selected_min
#' }
#' @name cv
NULL

#' @rdname cv
#' @export
cv <- function(object, ...) {
  UseMethod("cv")
}

#' @rdname cv
#' @export
cv.pp_lasso <- function(object, n_folds = 10, folds = NULL, seed = NULL, ...) {
  s <- object$settings
  if (is.null(s)) stop("`object` must come from var_select().", call. = FALSE)
  data <- structure(list(data = object$data, vars = object$vars, sizes = object$sizes,
                         cutoff = if (is.null(object$cutoff)) 0 else object$cutoff),
                    class = "pp_data")
  cv_selection(data, var_select, s, n_folds, folds, seed, full = object)
}

cv_selection <- function(data, fitter, args, n_folds, folds, seed, full = NULL) {
  if (is.null(full)) full <- do.call(fitter, c(list(data), args))
  m <- if (outcome_of(data) == "survival") survival_matrices(data) else selection_matrices(data)
  included <- data$data$included
  if (is.null(folds)) {
    folds <- make_folds(data, n_folds, seed)
  } else {
    if (length(folds) != nrow(data$data)) {
      stop("`folds` must have one value per row of `data$data`.", call. = FALSE)
    }
    folds <- as.integer(folds)[included]
    n_folds <- max(folds)
  }
  args$lambda <- full$lambda
  loss <- matrix(NA_real_, n_folds, length(full$lambda))
  n_fold <- integer(n_folds)
  for (k in seq_len(n_folds)) {
    train <- which(folds != k)
    test <- which(folds == k)
    if (length(test) == 0) next
    keep <- which(included)[train]
    fit_k <- do.call(fitter, c(list(pp_subset(data, keep)), args))
    res <- cv_fold_loss(fit_k, m, test)
    if (is.null(res)) next
    loss[k, ] <- res$loss
    n_fold[k] <- res$n
  }
  used <- n_fold > 0
  cve <- colSums(loss[used, , drop = FALSE] * n_fold[used]) / sum(n_fold)
  cvse <- apply(loss[used, , drop = FALSE], 2, stats::sd) / sqrt(sum(used))
  i_min <- which.min(cve)
  i_1se <- min(which(cve <= cve[i_min] + cvse[i_min]))
  structure(list(lambda = full$lambda, cve = cve, cvse = cvse,
                 lambda_min = full$lambda[i_min], lambda_1se = full$lambda[i_1se],
                 index_min = i_min, index_1se = i_1se,
                 selected_min = selected_covariates(full, which = i_min)[[1]],
                 selected_1se = selected_covariates(full, which = i_1se)[[1]],
                 fit = full, folds = folds, n_folds = n_folds),
            class = "pp_cv")
}

# The held-out loss of one fold, per value of lambda, and the number of
# units it averages over: the binomial deviance per row (logistic), the
# partial-likelihood deviance per weighted event within the held-out rows'
# own strata and risk sets (Cox), or the binomial deviance per person-period
# (discrete time). NULL when the fold has nothing to score.
cv_fold_loss <- function(fit, m, test) {
  model <- path_model(fit)
  newdata <- m$d[test, , drop = FALSE]
  if (model == "logistic") {
    p <- predict(fit, newdata = newdata, type = "response", unseen = "median")
    p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    dev <- -2 * (m$y[test] * log(p) + (m$n[test] - m$y[test]) * log(1 - p))
    return(list(loss = colMeans(dev), n = length(test)))
  }
  if (model == "cox") {
    n_events <- sum(m$w[test] * m$event[test])
    if (n_events == 0) return(NULL)
    eta <- predict(fit, newdata = newdata, type = "link")
    code <- match(m$prov[test], unique(m$prov[test]))
    st <- cox_structure(m$start[test], m$stop[test], m$event[test], m$w[test], code)
    dev <- cox_deviance(m$event[test], m$w[test], st$k_stop, st$k_start, st$row_first,
                        st$row_last, st$grid_first, st$grid_count, st$d_grid, eta)
    return(list(loss = as.numeric(dev) / n_events, n = n_events))
  }
  ex <- discrete_expanded(fit, newdata)
  h <- pmin(pmax(ex$hazard, 1e-10), 1 - 1e-10)
  dev <- -2 * (ex$y * log(h) + (1 - ex$y) * log(1 - h))
  list(loss = colMeans(dev), n = length(ex$y))
}

# Fold membership of the included rows: within each provider, its units
# (patients, or rows) are permuted and dealt to the folds in turn from a
# random starting fold.
make_folds <- function(data, n_folds, seed) {
  if (!is.numeric(n_folds) || n_folds < 2) stop("`n_folds` must be at least 2.", call. = FALSE)
  v <- data$vars
  d <- data$data[data$data$included, , drop = FALSE]
  unit <- if (is.null(v$patient)) seq_len(nrow(d)) else d[[v$patient]]
  unit_key <- paste(d[[v$provider]], unit, sep = "\r")
  with_seed(seed, {
    folds <- integer(nrow(d))
    for (p in unique(d[[v$provider]])) {
      rows <- which(d[[v$provider]] == p)
      units <- unique(unit_key[rows])
      units <- units[sample.int(length(units))]
      fold_of_unit <- (seq_along(units) - 1 + sample.int(n_folds, 1)) %% n_folds + 1
      folds[rows] <- fold_of_unit[match(unit_key[rows], units)]
    }
    folds
  })
}

# Evaluate `expr` under `seed` (when given) and restore the session's
# random-number state afterwards.
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  env <- globalenv()
  seed_name <- ".Random.seed"
  old <- if (exists(seed_name, envir = env, inherits = FALSE)) {
    get(seed_name, envir = env, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old)) rm(list = seed_name, envir = env) else assign(seed_name, old, envir = env)
  })
  set.seed(seed)
  expr
}

#' @export
print.pp_cv <- function(x, ...) {
  cat("Cross-validated ", if (inherits(x$fit, "pp_group_lasso")) "group-lasso" else "lasso",
      " path (", x$n_folds, " folds)\n", sep = "")
  cat("  lambda_min ", format(x$lambda_min, digits = 4), ": ", length(x$selected_min),
      " covariates selected", if (length(x$selected_min)) paste0(" (",
      paste(x$selected_min, collapse = ", "), ")"), "\n", sep = "")
  cat("  lambda_1se ", format(x$lambda_1se, digits = 4), ": ", length(x$selected_1se),
      " covariates selected\n", sep = "")
  invisible(x)
}
