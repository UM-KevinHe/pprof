#=== SUMMARY STATISTICS FOR COVARIATE ESTIMATES =================================
#' Provide the summary statistics for covariate estimates
#'
#' @param fit an object as output of \code{linear_fe} function.
#'
#' @param parm a character vector specifies a subset of covariates. All covariates are included by default.
#'
#' @param level confidence level used for constructing confidence intervals. Defaulting to 0.95.
#'
#' @param null the null value of the covariate estimate that requires testing. (e.g. test \eqn{H_0: \beta = 0})
#'
#' @param ...
#'
#'
#' @return a dataframe containing summary statistics for covariate estimates
#'
#' @importFrom stats pnorm qnorm pt qt
#'
#'
#' @export
summary_linearFE_covar <- function(fit, parm, level = 0.95, method = "pl", null = 0) {
  alpha <- 1 - level

  if (missing(fit)) stop ("Argument 'fit' is required!",call.=F)
  if (!class(fit) %in% c("linear_fe")) stop("Object fit is not of the classes 'linear_fe'!",call.=F)

  beta <- fit$coefficient$beta
  se.beta <- sqrt(diag(fit$variance$beta))
  m <- length(fit$coefficient$gamma)
  p <- length(fit$coefficient$beta)
  n <- nrow(fit$data_include)

  # Test Statistics
  stat <- (beta - null) / se.beta

  if (method == "pl") {
    p_value <- 2 * (1 - pnorm(abs(stat)))
    crit_value <- qnorm(1 - alpha / 2)
  } else if (method == "lm") {
    df <- n - p - m
    p_value <- 2 * (1 - pt(abs(stat), df))
    crit_value <- qt(1 - alpha / 2, df)
  } else {
    stop("Reference distribution must be either 'normal' or 't'")
  }

  p_value <- format.pval(p_value, digits = 7, eps = .Machine$double.eps)

  # Confidence Interval
  lower_bound <- beta - crit_value * se.beta
  upper_bound <- beta + crit_value * se.beta

  result <- data.frame(beta = beta, se.beta = se.beta, stat = stat, p_value = p_value,
                       lower_bound = lower_bound, upper_bound = upper_bound)
  colnames(result) <- c("Estimate", "Std.Error", "Stat", "p", "CI.Lower", "CI.Upper")

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
