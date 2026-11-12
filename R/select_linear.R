# Selection paths of the linear model with fixed provider effects (a
# continuous outcome), reached through fit_lasso_path() / fit_group_lasso_path().
# The provider effects are profiled out exactly by centring the outcome and
# the covariates within providers; the path then runs on the Gram matrix of
# the centred data (src/linear_penalized.cpp). Standardization and the
# orthonormalization of groups are those of the logistic paths (on the
# covariates as given), so the penalty means the same in both models.
linear_lasso <- function(data, penalize = NULL, group = NULL, multiplier = NULL, lambda = NULL,
                         n_lambda = 100, lambda_min_ratio = 1e-4, standardize = TRUE,
                         tol = 1e-4, max_iter = 10000, active_set = TRUE) {
  m <- selection_matrices(data)
  x <- data$vars$x
  gl <- penalty_layout(x, penalize, group, multiplier)
  std <- standardize_columns(m$z, standardize)
  within <- function(a) {
    a - rowsum(a, m$code, reorder = TRUE)[m$code, , drop = FALSE] / m$n_prov[m$code]
  }
  flat <- colMeans(within(std$z)^2) <= 1e-10 * pmax(colMeans(std$z^2), 1e-300)
  if (any(flat)) {
    stop("Covariate(s) constant within every provider, absorbed by the provider effects: ",
         paste(x[flat], collapse = ", "), call. = FALSE)
  }
  ortho <- orthonormalize_groups(std$z, gl$g_index, gl$levels_g, m$n_obs)
  zw <- within(ortho$z)
  yw <- as.numeric(within(matrix(m$y)))
  n <- m$n_obs
  gram <- crossprod(zw) / n
  cross <- as.numeric(crossprod(zw, yw)) / n
  n_unpen <- length(ortho$unpenalized)
  beta0 <- rep(0, ncol(zw))
  if (n_unpen) {
    u <- seq_len(n_unpen)
    beta0[u] <- solve(gram[u, u, drop = FALSE], cross[u])
  }
  if (is.null(lambda)) {
    grad <- cross - as.numeric(gram %*% beta0)
    score <- vapply(seq_along(gl$levels_g), function(g) {
      sqrt(sum(grad[ortho$col_first[g]:ortho$col_last[g]]^2))
    }, numeric(1))
    lambda <- lambda_sequence(max(score / gl$mult), n_lambda, lambda_min_ratio)
  } else {
    lambda <- sort(as.numeric(lambda), decreasing = TRUE)
  }
  res <- linear_lasso_path(gram, cross, sum(yw^2), n, beta0, seq_len(n_unpen) - 1L,
                           as.integer(ortho$col_first - 1), as.integer(ortho$col_last - 1),
                           unname(gl$mult), lambda, tol, max(1L, as.integer(max_iter)),
                           active_set)
  beta <- ortho$from_ortho(res$beta) / std$scale
  # the provider effect at each lambda: the provider's mean of y - Z beta
  ybar <- as.numeric(rowsum(m$y, m$code, reorder = TRUE)) / m$n_prov
  zbar <- rowsum(m$z, m$code, reorder = TRUE) / m$n_prov
  gamma <- ybar - zbar %*% beta
  lam <- signif(lambda, 6)
  dimnames(beta) <- list(x, lam)
  dimnames(gamma) <- list(m$levels, lam)
  out <- list(beta = beta, gamma = gamma, alpha = NULL, time_grid = NULL, lambda = lambda,
              deviance = as.numeric(res$deviance), iter = as.integer(res$iter),
              n_selected = as.integer(colSums(beta[!gl$unpen, , drop = FALSE] != 0)),
              penalized = stats::setNames(!gl$unpen, x), multiplier = gl$mult,
              data = m$d, vars = data$vars, sizes = m$sizes, rows = m$n_prov,
              control = list(standardize = standardize, tol = tol, max_iter = max_iter,
                             active_set = active_set),
              model = "linear")
  if (!is.null(group)) {
    out$group <- gl$label
    class(out) <- c("pp_group_lasso", "pp_lasso")
  } else {
    class(out) <- "pp_lasso"
  }
  out
}
