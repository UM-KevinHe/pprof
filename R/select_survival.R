# Selection paths for survival outcomes: the Cox model stratified by provider
# (lasso or group lasso on the covariates) and the discrete-time hazard model
# with period and provider effects (lasso). Both are reached through
# fit_lasso_path() / fit_group_lasso_path() on a survival pp_data object.

# Rows and columns of a survival pp_data object for the paths.
survival_matrices <- function(data) {
  if (!inherits(data, "pp_data")) stop("`data` must come from pp_data().", call. = FALSE)
  v <- data$vars
  d <- data$data[data$data$included, , drop = FALSE]
  if (nrow(d) == 0) stop("No provider reaches the cutoff.", call. = FALSE)
  rownames(d) <- NULL
  n <- nrow(d)
  prov <- d[[v$provider]]
  levels <- unique(prov)
  z <- as.matrix(d[, v$x, drop = FALSE])
  storage.mode(z) <- "double"
  list(d = d, z = z, prov = prov, levels = levels, code = match(prov, levels),
       n_prov = as.numeric(table(factor(prov, levels = levels))), n_obs = n,
       start = if (is.null(v$start)) rep(0, n) else as.numeric(d[[v$start]]),
       stop = as.numeric(d[[v$time]]), event = as.numeric(d[[v$event]]),
       w = if (is.null(v$weights)) rep(1, n) else as.numeric(d[[v$weights]]),
       offset = if (is.null(v$offset)) rep(0, n) else as.numeric(d[[v$offset]]),
       sizes = data$sizes[data$sizes$included, , drop = FALSE])
}

# The penalty layout of the covariates: which are unpenalized, the penalized
# groups (every penalized covariate its own group for the lasso) and the
# multipliers.
penalty_layout <- function(x, penalize, group, multiplier) {
  if (is.null(group)) {
    if (is.null(penalize)) penalize <- x
    bad <- setdiff(penalize, x)
    if (length(bad)) stop("Not covariates: ", paste(bad, collapse = ", "), call. = FALSE)
    unpen <- !x %in% penalize
    label <- ifelse(unpen, NA_character_, x)
    levels_g <- x[!unpen]
    mult <- rep(1, length(levels_g))
    if (!is.null(multiplier)) {
      if (is.null(names(multiplier)) || !all(names(multiplier) %in% levels_g)) {
        stop("`multiplier` must be named by penalized covariate.", call. = FALSE)
      }
      mult[match(names(multiplier), levels_g)] <- as.numeric(multiplier)
    }
  } else {
    if (!is.null(names(group))) {
      if (!setequal(names(group), x)) stop("`group` must name every covariate.", call. = FALSE)
      group <- group[x]
    } else if (length(group) != length(x)) {
      stop("`group` must have one label per covariate.", call. = FALSE)
    }
    label <- as.character(group)
    unpen <- is.na(label) | label %in% c("0", "")
    label[unpen] <- NA_character_
    levels_g <- unique(label[!unpen])
    size <- as.numeric(table(factor(label[!unpen], levels = levels_g)))
    mult <- sqrt(size)
    if (!is.null(multiplier)) {
      if (is.null(names(multiplier)) || !all(names(multiplier) %in% levels_g)) {
        stop("`multiplier` must be named by penalized group.", call. = FALSE)
      }
      mult[match(names(multiplier), levels_g)] <- as.numeric(multiplier)
    }
  }
  if (length(levels_g) == 0) stop("No penalized covariate.", call. = FALSE)
  list(unpen = unpen, label = stats::setNames(label, x), levels_g = levels_g,
       g_index = match(label, levels_g), mult = stats::setNames(mult, levels_g))
}

# Per-stratum grids of distinct event times and the row indices into them
# (see cox_working in src/survival_penalized.cpp). Rows must be grouped by
# `code`, which runs 1..m.
cox_structure <- function(start, stop, event, w, code) {
  n <- length(stop)
  rows_of <- split(seq_len(n), factor(code, levels = seq_len(max(code))))
  m <- length(rows_of)
  k_stop <- k_start <- integer(n)
  grid_first <- grid_count <- row_first <- row_last <- integer(m)
  d_grid <- vector("list", m)
  pos <- 0L
  for (j in seq_len(m)) {
    rows <- rows_of[[j]]
    row_first[j] <- rows[1]
    row_last[j] <- rows[length(rows)]
    ev <- rows[event[rows] > 0]
    times <- sort(unique(stop[ev]))
    k_stop[rows] <- findInterval(stop[rows], times)
    k_start[rows] <- findInterval(start[rows], times)
    grid_first[j] <- pos
    grid_count[j] <- length(times)
    pos <- pos + length(times)
    d_grid[[j]] <- if (length(times)) {
      as.numeric(rowsum(w[ev] * event[ev], match(stop[ev], times), reorder = TRUE))
    } else {
      numeric(0)
    }
  }
  list(k_stop = k_stop, k_start = k_start, row_first = row_first - 1L,
       row_last = row_last - 1L, grid_first = grid_first, grid_count = grid_count,
       d_grid = unlist(d_grid))
}

cox_lasso <- function(data, penalize = NULL, group = NULL, multiplier = NULL, lambda = NULL,
                      n_lambda = 100, lambda_min_ratio = 1e-4, standardize = TRUE, tol = 1e-4,
                      max_iter = 10000, active_set = TRUE) {
  m <- survival_matrices(data)
  x <- data$vars$x
  gl <- penalty_layout(x, penalize, group, multiplier)
  std <- standardize_columns(m$z, standardize)
  ortho <- orthonormalize_groups(std$z, gl$g_index, gl$levels_g, m$n_obs)
  st <- cox_structure(m$start, m$stop, m$event, m$w, m$code)
  beta0 <- rep(0, ncol(ortho$z))
  n_unpen <- length(ortho$unpenalized)
  if (n_unpen) {
    dd <- data.frame(.start = m$start, .time = m$stop, .event = m$event, .offset = m$offset,
                     .provider = factor(m$code))
    dd$.weights <- m$w
    dd$zu <- ortho$z[, seq_len(n_unpen), drop = FALSE]
    fit <- cox_stage(c("strata(.provider)", "zu", "offset(.offset)"),
                     "Surv(.start, .time, .event)", dd, "breslow")
    beta0[seq_len(n_unpen)] <- as.numeric(stats::coef(fit$fit))
  }
  eta0 <- as.numeric(ortho$z %*% beta0) + m$offset
  if (is.null(lambda)) {
    r <- cox_residuals(m$event, m$w, st$k_stop, st$k_start, st$row_first, st$row_last,
                       st$grid_first, st$grid_count, st$d_grid, eta0)$residual
    score <- vapply(seq_along(gl$levels_g), function(g) {
      cols <- ortho$col_first[g]:ortho$col_last[g]
      sqrt(sum(crossprod(ortho$z[, cols, drop = FALSE], r)^2))
    }, numeric(1))
    lambda_max <- max(score / (gl$mult * m$n_obs))
    lambda <- lambda_sequence(lambda_max, n_lambda, lambda_min_ratio)
  } else {
    lambda <- sort(as.numeric(lambda), decreasing = TRUE)
  }
  res <- cox_lasso_path(m$event, m$w, m$offset, ortho$z, st$k_stop, st$k_start, st$row_first,
                        st$row_last, st$grid_first, st$grid_count, st$d_grid, beta0,
                        seq_along(ortho$unpenalized) - 1L, as.integer(ortho$col_first - 1),
                        as.integer(ortho$col_last - 1), unname(gl$mult), lambda, tol,
                        max(1L, as.integer(max_iter)), active_set)
  beta <- ortho$from_ortho(res$beta) / std$scale
  out <- survival_path_object(data, m, beta, res, lambda, !gl$unpen, gl$mult,
                              control = list(standardize = standardize, tol = tol,
                                             max_iter = max_iter, active_set = active_set),
                              model = "cox")
  if (!is.null(group)) {
    out$group <- gl$label
    class(out) <- c("pp_group_lasso", "pp_lasso")
  } else {
    class(out) <- "pp_lasso"
  }
  out
}

discrete_lasso <- function(data, penalize = NULL, multiplier = NULL, lambda = NULL,
                           n_lambda = 100, lambda_min_ratio = 1e-4, standardize = TRUE,
                           bound = 10, tol = 1e-4, max_iter = 10000, active_set = TRUE,
                           mm = FALSE, backtrack = FALSE) {
  v <- data$vars
  if (!is.null(v$start)) {
    stop("The discrete-time model takes one row per subject, without `start`.", call. = FALSE)
  }
  if (!is.null(v$weights)) {
    stop("The discrete-time model does not use case weights.", call. = FALSE)
  }
  m <- survival_matrices(data)
  x <- v$x
  gl <- penalty_layout(x, penalize, NULL, multiplier)
  grid <- sort(unique(m$stop))
  t_index <- match(m$stop, grid)
  n_period <- length(grid)
  d_period <- as.numeric(tabulate(t_index[m$event > 0], n_period))
  at_risk <- rev(cumsum(rev(tabulate(t_index, n_period))))
  alpha0 <- stats::qlogis(pmin(pmax(d_period / at_risk, 1e-10), 1 - 1e-10))
  std <- standardize_columns(m$z, standardize)
  scale2 <- sqrt(colMeans(std$z^2))
  z <- sweep(std$z, 2, scale2, "/")
  penalized <- as.integer(!gl$unpen)
  mult <- rep(1, length(x))
  mult[!gl$unpen] <- unname(gl$mult)
  n_prov <- as.numeric(m$n_prov)
  null <- discrete_lasso_path(t_index, m$event, z, m$offset, n_prov, d_period, alpha0,
                              rep(0, length(m$levels)), rep(0, length(x)), penalized, mult, Inf,
                              mm, backtrack, bound, tol, max(1L, as.integer(max_iter)), FALSE)
  alpha_n <- null$alpha[, 1]
  gamma_n <- null$gamma[, 1]
  beta_n <- null$beta[, 1]
  if (is.null(lambda)) {
    eta_n <- gamma_n[m$code] + as.numeric(z %*% beta_n) + m$offset
    r <- m$event - discrete_hazards(t_index, alpha_n, eta_n)$expected
    score <- abs(as.numeric(crossprod(z[, !gl$unpen, drop = FALSE], r)))
    lambda_max <- max(score / (mult[!gl$unpen] * m$n_obs))
    lambda <- lambda_sequence(lambda_max, n_lambda, lambda_min_ratio)
  } else {
    lambda <- sort(as.numeric(lambda), decreasing = TRUE)
  }
  res <- discrete_lasso_path(t_index, m$event, z, m$offset, n_prov, d_period, alpha_n, gamma_n,
                             beta_n, penalized, mult, lambda, mm, backtrack, bound, tol,
                             max(1L, as.integer(max_iter)), active_set)
  beta <- res$beta / (std$scale * scale2)
  alpha <- res$alpha - matrix(as.numeric(crossprod(std$center, beta)), n_period, ncol(beta),
                              byrow = TRUE)
  out <- survival_path_object(data, m, beta, res, lambda, !gl$unpen, gl$mult,
                              control = list(standardize = standardize, bound = bound,
                                             tol = tol, max_iter = max_iter,
                                             active_set = active_set, mm = mm,
                                             backtrack = backtrack),
                              model = "discrete", alpha = alpha, gamma = res$gamma,
                              time_grid = grid)
  class(out) <- "pp_lasso"
  out
}

survival_path_object <- function(data, m, beta, res, lambda, penalized, multiplier, control,
                                 model, alpha = NULL, gamma = NULL, time_grid = NULL) {
  x <- data$vars$x
  lam <- signif(lambda, 6)
  dimnames(beta) <- list(x, lam)
  if (!is.null(gamma)) dimnames(gamma) <- list(m$levels, lam)
  if (!is.null(alpha)) dimnames(alpha) <- list(signif(time_grid, 6), lam)
  list(beta = beta, gamma = gamma, alpha = alpha, time_grid = time_grid, lambda = lambda,
       deviance = as.numeric(res$deviance), iter = as.integer(res$iter),
       n_selected = as.integer(colSums(beta[penalized, , drop = FALSE] != 0)),
       penalized = stats::setNames(penalized, x), multiplier = multiplier,
       data = m$d, vars = data$vars, sizes = m$sizes, rows = m$n_prov, control = control,
       model = model)
}

# The model of a path object (older objects are logistic).
path_model <- function(object) {
  if (is.null(object$model)) "logistic" else object$model
}

# Hazards of every row of `newdata` in every period up to its exit period
# (the fit's periods at or before the row's time, at least one), for each
# requested lambda: the expanded response `y`, the row of each expanded
# entry, and a matrix of hazards.
discrete_expanded <- function(object, newdata, lambda = NULL, which = NULL) {
  v <- object$vars
  cf <- coef(object, lambda = lambda, which = which, drop = FALSE)
  eta <- predict(object, newdata = newdata, lambda = lambda, which = which, type = "link",
                 unseen = "median")
  periods <- pmax(findInterval(as.numeric(newdata[[v$time]]), object$time_grid), 1L)
  n <- nrow(newdata)
  y <- numeric(sum(periods))
  y[cumsum(periods)] <- as.numeric(newdata[[v$event]])
  hazard <- matrix(NA_real_, sum(periods), ncol(eta))
  for (l in seq_len(ncol(eta))) {
    hazard[, l] <- discrete_hazards(periods, cf$alpha[, l], eta[, l])$hazard
  }
  colnames(hazard) <- colnames(eta)
  list(y = y, row = rep(seq_len(n), periods), hazard = hazard)
}
