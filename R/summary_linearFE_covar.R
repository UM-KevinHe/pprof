#=== SUMMARY STATISTICS FOR COVARIATE ESTIMATES =================================
#' Provide the summary statistics for covariate estimates
#'
#' @param fit an object as output of \code{linear_fe} function.
#'
#' @param parm a character vector specifies a subset of covariates. All covariates are included by default.
#'
#' @param level confidence level used for constructing confidence intervals. Default is 0.95.
#'
#' @param null the null value of the covariate estimate that requires testing. (e.g. test \eqn{H_0: \beta = 0})
#'
#' @param alternative a character string specifies the alternative hypothesis.
#'  Must be one of \code{"two.sided"}, \code{"less"}, or \code{"greater"}.
#'  Default is \code{"two.sided"}.
#'
#' @return a dataframe containing summary statistics for covariate estimates
#'
#' @importFrom stats pnorm qnorm pt qt
#'
#' @export
summary_linearFE_covar <- function(fit, parm, level = 0.95, null = 0, alternative = "two.sided") {
  alpha <- 1 - level

  if (missing(fit)) stop ("Argument 'fit' is required!",call.=F)
  if (!class(fit) %in% c("linear_fe")) stop("Object fit is not of the classes 'linear_fe'!",call.=F)

  Z.char <- fit$char_list$Z.char
  beta <- fit$coefficient$beta
  se.beta <- sqrt(diag(fit$variance$beta))
  m <- length(fit$coefficient$gamma)
  p <- length(fit$coefficient$beta)
  n <- nrow(fit$data_include)
  model.method <- fit$method

  # Test Statistics
  stat <- (beta - null) / se.beta

  if (alternative == "two.sided") {
    p_value <- switch(model.method,
                      "Profile Likelihood" = 2 * (1 - pnorm(abs(stat))),
                      "Dummy" = 2 * (1 - pt(abs(stat), df = n - p - m)))
    crit_value <- ifelse(model.method == "Profile Likelihood", qnorm(1 - alpha / 2),
                         qt(1 - alpha / 2, df = n - p - m))
    lower_bound <- beta - crit_value * se.beta
    upper_bound <- beta + crit_value * se.beta
  }
  else if (alternative == "greater") {
    p_value <- switch(model.method,
                      "Profile Likelihood" = 1 - pnorm(stat),
                      "Dummy" = 1 - pt(stat, df = n - p - m))
    crit_value <- ifelse(fit$method == "Profile Likelihood", qnorm(level),
                           qt(level, df = n - p - m))

    lower_bound <- beta - crit_value * se.beta
    upper_bound <- Inf
  }
  else if (alternative == "less") {
    p_value <- switch(model.method,
                      "Profile Likelihood" = pnorm(stat),
                      "Dummy" = pt(stat, df = n - p - m))
    crit_value <- ifelse(fit$method == "Profile Likelihood", qnorm(level),
                         qt(level, df = n - p - m))

    lower_bound <- -Inf
    upper_bound <- beta + crit_value * se.beta
  }
  else {
    stop("Argument 'alternative' should be one of 'two.sided', 'less', 'greater'.")
  }

  p_value <- format.pval(p_value, digits = 7, eps = .Machine$double.eps)

  result <- data.frame(beta = beta, se.beta = se.beta, stat = stat, p_value = p_value,
                       lower_bound = lower_bound, upper_bound = upper_bound)
  colnames(result) <- c("Estimate", "Std.Error", "Stat", "p value", "CI.Lower", "CI.Upper")

  if (missing(parm)) {
    ind <- 1:length(Z.char)
  } else if (is.character(parm)) {
    ind <- which(Z.char %in% parm)
  } else if (is.numeric(parm) & max(abs(as.integer(parm) - parm)) == 0 & !(0 %in% parm)) {
    ind <- parm
  } else {
    stop("Argument 'parm' includes invalid elements!")
  }

  result <- result[ind, ]
  return(result)
}
