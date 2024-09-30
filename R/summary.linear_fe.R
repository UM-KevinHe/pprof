#' Result Summaries of Covariate Estimates
#'
#' Provide the summary statistics for the covariate estimates for a fixed/random effect linear model.
#'
#' @param fit a model fitted from \code{linear_fe} or \code{linear_re}.
#' @param parm Specifies a subset of covariates for which the result summaries should be output.
#' By default, all covariates are included.
#' @param level the confidence level during the hypothesis test, meaning a significance level of \eqn{1 - \text{level}}.
#' The default value is 0.95.
#' @param null a number defining the null hypothesis for the covariate estimates. The default value is \code{0}.
#' @param alternative a character string specifying the alternative hypothesis, must be one of
#' \code{"two.sided"} (default), \code{"greater"}, or \code{"less"}.
#'
#' @return A data frame containing summary statistics for covariate estimates, with the following columns:
#' \item{Estimate}{the estimates of covariate coefficients.}
#' \item{Std.Error}{the standard error of the estimate.}
#' \item{Stat}{the test statistic.}
#' \item{p value}{the p-value for the hypothesis test.}
#' \item{CI.upper}{the lower bound of the confidence interval.}
#' \item{CI.lower}{the upper bound of the confidence interval.}
#'
#' @examples
#' data(ExampleDataLinear)
#' Y <- ExampleDataLinear$Y
#' Z <- ExampleDataLinear$Z
#' ID <- ExampleDataLinear$ID
#'
#' fit_fe <- linear_fe(Y = Y, Z = Z, ID = ID)
#' summary(fit_fe)
#'
#' @importFrom stats pnorm qnorm pt qt
#'
#' @exportS3Method summary linear_fe

summary.linear_fe <- function(fit, parm, level = 0.95, null = 0, alternative = "two.sided") {
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
