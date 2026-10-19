#' Data quality checks for the fixed provider effects model
#'
#' Checks missingness, (near-)zero variance of covariates, pairwise
#' correlation among covariates and variance inflation factors. This is the
#' check performed by [fe_data_prep()] when `check = TRUE`, available on its
#' own for exploratory use.
#'
#' @param Y the outcome: a numeric vector or one-column data frame.
#' @param Z a numeric matrix or data frame of covariates.
#' @param ProvID the provider identifier: a vector or one-column data frame.
#' @param N optional number of trials (binomial outcome).
#' @param message a logical; print progress messages.
#'
#' @details Near-zero variance uses the usual rule (frequency ratio of the two
#' most common values > 95/5 and at most 10% unique values). A covariate with
#' zero variance or data with missing values stop with an error; near-zero
#' variance, correlation above 0.9 and VIF of 10 or more give warnings.
#'
#' @return Invisibly, a list with the near-zero-variance metrics, the
#' correlation matrix and the VIFs.
#'
#' @export
data_check <- function(Y, Z, ProvID, N = NULL, message = TRUE) {
  Z <- as.data.frame(Z)
  data <- data.frame(Y = as.numeric(unlist(Y)), ProvID = unlist(ProvID), Z, check.names = FALSE)
  Z.char <- colnames(Z)
  if (!is.null(N)) data$N <- as.numeric(unlist(N))
  check_vars(data, N.char = if (is.null(N)) NULL else "N", Y.char = "Y", Z.char = Z.char,
             prov.char = "ProvID", message = message)
}

# Frequency ratio / percent unique metrics for near-zero variance.
nzv_metrics <- function(x, freqCut = 95 / 5, uniqueCut = 10) {
  x <- as.data.frame(x)
  metrics <- lapply(x, function(v) {
    n_all <- length(v)
    v <- v[!is.na(v)]
    tab <- table(v)
    freqRatio <- if (length(tab) <= 1) 0 else {
      w <- which.max(tab)
      max(tab) / max(tab[-w])
    }
    lunique <- length(unique(v))
    percentUnique <- 100 * lunique / n_all
    zeroVar <- (lunique == 1) || length(v) == 0
    c(freqRatio = freqRatio, percentUnique = percentUnique, zeroVar = zeroVar,
      nzv = (freqRatio > freqCut && percentUnique <= uniqueCut) || zeroVar)
  })
  out <- as.data.frame(do.call(rbind, metrics))
  out$zeroVar <- as.logical(out$zeroVar)
  out$nzv <- as.logical(out$nzv)
  rownames(out) <- colnames(x)
  out
}

# Variance inflation factors: 1 / (1 - R^2_j) regressing covariate j on the
# other covariates.
vif_values <- function(Z) {
  Z <- as.matrix(Z)
  p <- ncol(Z)
  vif <- stats::setNames(rep(1, p), colnames(Z))
  if (p < 2) return(vif)
  for (j in seq_len(p)) {
    fit <- stats::lm.fit(cbind(1, Z[, -j, drop = FALSE]), Z[, j])
    tss <- sum((Z[, j] - mean(Z[, j]))^2)
    r2 <- 1 - sum(fit$residuals^2) / tss
    vif[j] <- 1 / (1 - r2)
  }
  vif
}

# The checks run by fe_data_prep(check = TRUE).
check_vars <- function(data, N.char, Y.char, Z.char, prov.char, message = TRUE) {
  vars <- c(N.char, Y.char, Z.char, prov.char)
  ## check missingness of variables
  if (message == TRUE) message("Checking missingness of variables ... ")
  if (sum(stats::complete.cases(data[, vars])) == NROW(data)) {
    if (message == TRUE) message("Missing values NOT found. Checking missingness of variables completed!")
  } else {
    for (name in vars) {
      if (sum(is.na(data[, name])) > 0) {
        warning(sum(is.na(data[, name])), " out of ", NROW(data[, name]), " in '", name, "' missing!",
                immediate. = TRUE, call. = FALSE)
      }
    }
    missingness <- (1 - sum(stats::complete.cases(data[, vars])) / NROW(data)) * 100
    stop(paste(round(missingness, 2), "% of all observations are missing!", sep = ""), call. = FALSE)
  }
  ## check variation in covariates
  if (message == TRUE) message("Checking variation in covariates ... ")
  nzv <- nzv_metrics(data[, Z.char, drop = FALSE])
  if (sum(nzv$zeroVar) > 0) {
    stop("Covariate(s) '", paste(rownames(nzv)[nzv$zeroVar], collapse = "', '"),
         "' with zero variance(s)!", call. = FALSE)
  } else if (sum(nzv$nzv) > 0) {
    warning("Covariate(s) '", paste(rownames(nzv)[nzv$nzv], collapse = "', '"),
            "' with near zero variance(s)!", immediate. = TRUE, call. = FALSE)
  }
  if (message == TRUE) message("Checking variation in covariates completed!")
  ## check correlation
  if (message == TRUE) message("Checking pairwise correlation among covariates ... ")
  cor_mat <- stats::cor(as.matrix(data[, Z.char, drop = FALSE]))
  threshold.cor <- 0.9
  upper <- cor_mat
  upper[lower.tri(upper, diag = TRUE)] <- 0
  high <- which(abs(upper) > threshold.cor, arr.ind = TRUE)
  if (NROW(high) > 0) {
    warning("The following ", NROW(high), " pair(s) of covariates are highly correlated (correlation > ",
            threshold.cor, "): ", immediate. = TRUE, call. = FALSE)
    for (k in seq_len(NROW(high))) {
      message('("', rownames(cor_mat)[high[k, 1]], '", "', colnames(cor_mat)[high[k, 2]], '")')
    }
  }
  if (message == TRUE) message("Checking pairwise correlation among covariates completed!")
  ## check VIF
  if (message == TRUE) message("Checking VIF of covariates ... ")
  vif <- vif_values(data[, Z.char, drop = FALSE])
  if (sum(vif >= 10) > 0) {
    warning("Covariate(s) '", paste(names(vif)[vif >= 10], collapse = "', '"),
            "' with serious multicollinearity!", immediate. = TRUE, call. = FALSE)
  }
  if (message == TRUE) message("Checking VIF of covariates completed!")
  invisible(list(nzv = nzv, cor = cor_mat, vif = vif))
}
