#' Data checks for a provider model
#'
#' Missing values, (near-)zero variance of covariates, pairwise correlation
#' among covariates and variance inflation factors. [pp_data()] runs these
#' checks when `check = TRUE`; the function is exported for exploratory use.
#'
#' @inheritParams pp_data
#' @param cor_threshold warn about pairs of covariates whose absolute
#'   correlation exceeds this value. Default 0.9.
#' @param vif_threshold warn about covariates whose variance inflation
#'   factor is at least this value. Default 10.
#'
#' @details Missing values in any named column and covariates with zero
#' variance are errors. Near-zero variance (frequency ratio of the two most
#' common values above 95/5 and at most 10% distinct values), high
#' correlation and large VIF are warnings. The VIF of a covariate is
#' \eqn{1/(1 - R^2)} from regressing it on the other covariates.
#'
#' @return Invisibly, a list: `missing` (missing values per column), `nzv`
#'   (near-zero-variance metrics per covariate), `cor` (correlation matrix),
#'   `high_cor` (pairs above the threshold) and `vif`.
#' @examples
#' d <- sim_kecc(n_providers = 10, seed = 2)
#' chk <- pp_check(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider", message = FALSE)
#' chk$vif
#' @export
pp_check <- function(data, y, x, provider, n = NULL, cor_threshold = 0.9, vif_threshold = 10,
                     message = TRUE) {
  cols <- c(y, x, provider, n)
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols)) {
    stop("Column(s) not in `data`: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  say <- function(...) if (message) message(...)

  say("Checking missing values ...")
  n_missing <- vapply(data[cols], function(v) sum(is.na(v)), integer(1))
  if (any(n_missing > 0)) {
    for (v in names(n_missing)[n_missing > 0]) {
      warning(n_missing[[v]], " of ", nrow(data), " values of '", v, "' missing.",
              immediate. = TRUE, call. = FALSE)
    }
    share <- 100 * mean(!stats::complete.cases(data[cols]))
    stop(format(share, digits = 3),
         "% of the rows have a missing value; remove or impute them first.", call. = FALSE)
  }

  say("Checking variation in covariates ...")
  nzv <- nzv_metrics(data[x])
  if (any(nzv$zero_var)) {
    stop("Covariate(s) with zero variance: ", paste(rownames(nzv)[nzv$zero_var], collapse = ", "),
         call. = FALSE)
  }
  if (any(nzv$nzv)) {
    warning("Covariate(s) with near-zero variance: ",
            paste(rownames(nzv)[nzv$nzv], collapse = ", "), immediate. = TRUE, call. = FALSE)
  }

  say("Checking pairwise correlation among covariates ...")
  cor_mat <- stats::cor(as.matrix(data[x]))
  upper <- cor_mat
  upper[lower.tri(upper, diag = TRUE)] <- 0
  high <- which(abs(upper) > cor_threshold, arr.ind = TRUE)
  high_cor <- data.frame(a = rownames(cor_mat)[high[, 1]], b = colnames(cor_mat)[high[, 2]],
                         cor = cor_mat[high], stringsAsFactors = FALSE)
  if (nrow(high_cor)) {
    warning(nrow(high_cor), " pair(s) of covariates with |correlation| > ", cor_threshold, ": ",
            paste0(high_cor$a, "/", high_cor$b, collapse = ", "), immediate. = TRUE, call. = FALSE)
  }

  say("Checking variance inflation factors ...")
  vif <- vif_values(data[x])
  if (any(vif >= vif_threshold)) {
    warning("Covariate(s) with VIF >= ", vif_threshold, ": ",
            paste(names(vif)[vif >= vif_threshold], collapse = ", "),
            immediate. = TRUE, call. = FALSE)
  }
  say("Checks completed.")
  invisible(list(missing = n_missing, nzv = nzv, cor = cor_mat, high_cor = high_cor, vif = vif))
}

# Frequency ratio and percent-unique metrics for near-zero variance (the
# usual caret rule: ratio > 95/5 and at most 10% unique values).
nzv_metrics <- function(x, freq_cut = 95 / 5, unique_cut = 10) {
  x <- as.data.frame(x)
  rows <- lapply(x, function(v) {
    n_all <- length(v)
    v <- v[!is.na(v)]
    tab <- table(v)
    freq_ratio <- if (length(tab) <= 1) {
      0
    } else {
      w <- which.max(tab)
      max(tab) / max(tab[-w])
    }
    n_unique <- length(tab)
    percent_unique <- 100 * n_unique / n_all
    zero_var <- n_unique <= 1
    c(freq_ratio = freq_ratio, percent_unique = percent_unique, zero_var = zero_var,
      nzv = (freq_ratio > freq_cut && percent_unique <= unique_cut) || zero_var)
  })
  out <- as.data.frame(do.call(rbind, rows))
  out$zero_var <- as.logical(out$zero_var)
  out$nzv <- as.logical(out$nzv)
  rownames(out) <- names(x)
  out
}

# Variance inflation factors: 1 / (1 - R^2_j), regressing covariate j on the
# others.
vif_values <- function(x) {
  x <- as.matrix(x)
  p <- ncol(x)
  vif <- stats::setNames(rep(1, p), colnames(x))
  if (p < 2) return(vif)
  for (j in seq_len(p)) {
    fit <- stats::lm.fit(cbind(1, x[, -j, drop = FALSE]), x[, j])
    tss <- sum((x[, j] - mean(x[, j]))^2)
    vif[j] <- 1 / (1 - (1 - sum(fit$residuals^2) / tss))
  }
  vif
}
