#' Tests and confidence intervals for covariate coefficients
#'
#' @param fit a model fitted by [logis_fe()].
#' @param parm optional subset of covariates, by name or position.
#' @param level confidence level. Default 0.95.
#' @param test the test:
#'   \itemize{
#'   \item `"wald"` (default): model-based Wald test and confidence interval;
#'   \item `"robust_wald"`: Wald test with the patient-clustered sandwich
#'     variance (requires `PatientID` in [fe_data_prep()]);
#'   \item `"lr"`: likelihood ratio test (refits the model without the covariate);
#'   \item `"score"`: score test (refits the model without the covariate).
#'   }
#' @param null the null value of the coefficient. `"lr"` and `"score"` only
#'   allow 0.
#'
#' @return A data frame with one row per covariate. For `"wald"` and
#'   `"robust_wald"`: `beta`, `se.beta`, `Z` (test statistic), `p`,
#'   `CI.lower`, `CI.upper`. For `"lr"` and `"score"`: `beta`, `stat`
#'   (chi-square statistic, 1 df) and `p`.
#'
#' @export
summary_fe_covar <- function(fit, parm, level = 0.95, test = "wald", null = 0) {
  if (missing(fit)) stop("Argument 'fit' is required!", call. = FALSE)
  check_fit(fit)
  if (!test %in% c("wald", "robust_wald", "lr", "score")) stop("Argument 'test' NOT as required!", call. = FALSE)
  N.char <- fit$char_list$N.char
  Y.char <- fit$char_list$Y.char
  Z.char <- fit$char_list$Z.char
  prov.char <- fit$char_list$prov.char
  alpha <- 1 - level

  if (missing(parm)) {
    ind <- 1:length(Z.char)
  } else if (is.character(parm)) {
    ind <- which(Z.char %in% parm)
  } else if (is.numeric(parm) & max(abs(as.integer(parm) - parm)) == 0 & !(0 %in% parm)) {
    ind <- parm
  } else {
    stop("Argument 'parm' includes invalid elements!", call. = FALSE)
  }

  data <- fit$data_include
  N <- trials_of(data, N.char)
  n.prov <- sapply(split(data[, Y.char], data[, prov.char]), length)
  gamma <- fit$df.prov$gamma_est
  beta <- as.numeric(fit$beta)
  Z <- as.matrix(data[, Z.char])

  if (test %in% c("wald", "robust_wald")) {
    if (test == "wald") {
      ls <- wald_covar(N, data[, Y.char], Z, n.prov, gamma, beta, ind, null, alpha)
    } else {
      if (is.null(fit$char_list$PatientID.char)) {
        stop("Argument 'test' = 'robust_wald' requires 'PatientID' in the data!", call. = FALSE)
      }
      residual <- data[, Y.char] - N * fit$pred
      meat <- robust_beta_meat(residual, Z, data[, prov.char], data[, fit$char_list$PatientID.char],
                               names(n.prov))
      ls <- robust_wald_covar(N, data[, Y.char], Z, n.prov, gamma, beta, ind, null, alpha,
                              meat$A0, meat$A1_t, meat$A2)
    }
    df <- data.frame(beta = fit$beta[ind],
                     se.beta = c(ls$se.beta),
                     Z = c(ls$stat),
                     p = c(ls$p),
                     CI.lower = c(ls$beta.lower),
                     CI.upper = c(ls$beta.upper))
    row.names(df) <- Z.char[ind]
    return(df)
  }

  # likelihood ratio and score tests refit the model without one covariate
  if (null != 0) stop("Argument 'null' is invalid!", call. = FALSE)
  if (length(Z.char) < 2) stop("The 'lr' and 'score' tests need at least two covariates.", call. = FALSE)
  refit_without <- function(index) {
    char_list.null <- fit$char_list
    char_list.null$Z.char <- Z.char[-index]
    keep <- setdiff(colnames(data), Z.char[index])
    data.prep.null <- structure(list(data = data[, keep], char_list = char_list.null), class = "data_prep")
    ctrl <- fit$control
    logis_fe(data.prep.null, algorithm = ctrl$algorithm, max.iter = ctrl$max.iter, tol = ctrl$tol,
             bound = ctrl$bound, backtrack = ctrl$backtrack, threads = ctrl$threads, stop = ctrl$stop)
  }
  neg2lkd <- function(g, b, Zm) {
    g.obs <- rep(pmax(pmin(g, stats::median(g) + 10), stats::median(g) - 10), n.prov)
    eta <- g.obs + Zm %*% b
    -2 * sum(eta * data[, Y.char] - N * log(1 + exp(eta)))
  }

  if (test == "lr") {
    full <- neg2lkd(gamma, beta, Z)
    stat <- sapply(ind, function(index) {
      fe.null <- refit_without(index)
      neg2lkd(fe.null$df.prov$gamma_est, as.numeric(fe.null$beta), Z[, -index, drop = FALSE]) - full
    })
  } else {
    stat <- sapply(ind, function(index) {
      fe.null <- refit_without(index)
      Z.null <- Z[, -index, drop = FALSE]
      z.add <- Z[, index]
      probs.null <- as.numeric(stats::plogis(rep(fe.null$df.prov$gamma_est, n.prov) + Z.null %*% as.numeric(fe.null$beta)))
      probs.null <- pmin(pmax(probs.null, 1e-10), 1 - 1e-10)
      w <- N * probs.null * (1 - probs.null)
      prov <- provider_codes(data[, prov.char], names(n.prov))
      info.gamma.inv.null <- 1 / as.numeric(rowsum(w, prov, reorder = TRUE))
      info.betagamma.null <- rowsum(w * Z.null, prov, reorder = TRUE) # m x (p-1)
      info.beta.null <- t(Z.null) %*% (w * Z.null)
      mat.tmp <- info.gamma.inv.null * info.betagamma.null            # m x (p-1)
      schur.inv.null <- solve(info.beta.null - t(info.betagamma.null) %*% mat.tmp)
      info.inv.11 <- mat.tmp %*% schur.inv.null %*% t(mat.tmp)
      diag(info.inv.11) <- diag(info.inv.11) + info.gamma.inv.null
      info.inv.21 <- -schur.inv.null %*% t(mat.tmp)
      info.beta.add.gamma <- as.numeric(rowsum(w * z.add, prov, reorder = TRUE))
      info.beta.add.beta <- t(z.add) %*% (w * Z.null)
      info.beta.add <- sum(z.add * w * z.add) -
        t(info.beta.add.gamma) %*% info.inv.11 %*% info.beta.add.gamma -
        2 * info.beta.add.beta %*% info.inv.21 %*% info.beta.add.gamma -
        info.beta.add.beta %*% schur.inv.null %*% t(info.beta.add.beta)
      score.beta.add <- sum(z.add * (data[, Y.char] - N * probs.null))
      as.numeric(score.beta.add^2 / info.beta.add)
    })
  }
  df <- data.frame(beta = fit$beta[ind], stat = stat,
                   p = stats::pchisq(stat, 1, lower.tail = FALSE))
  rownames(df) <- Z.char[ind]
  df
}
