#' Conduct hypothesis testing for provider effects from a fitted `linear_re` object
#'
#' Conduct hypothesis tests on provider effects and identify outlying providers for a random effect linear model.
#'
#' @param fit a model fitted from \code{linear_re}.
#' @param parm specifies a subset of providers for which confidence intervals are to be given.
#' By default, all providers are included. The class of `parm` should match the class of the provider IDs.
#' @param level the confidence level during the hypothesis test, meaning a significance level of \eqn{1 - \text{level}}.
#' The default value is 0.95.
#' @param null a number defining the null hypothesis for the provider effects.
#' The default value is 0.
#' @param alternative a character string specifying the alternative hypothesis, must be one of
#' \code{"two.sided"} (default), \code{"greater"}, or \code{"less"}.
#' @param \dots additional arguments that can be passed to the function.
#'
#' @return A data frame containing the results of the hypothesis test, with the following columns:
#' \item{flag}{a flagging indicator where \code{1} means statistically higher than expected
#' and \code{-1} means statistically lower than expected.}
#' \item{p-value}{the p-value of the hypothesis test.}
#' \item{stat}{the test statistic.}
#' \item{Std.Error}{the standard error of the provider effect estimate.}
#'
#' @details
#' The function identifies outlying providers based on hypothesis test results.
#' For two-sided tests, \code{1} indicates performance significantly higher than expected, \code{-1} indicates lower,
#' For one-sided tests, \code{1} (right-tailed) or \code{-1} (left-tailed) flags are used.
#' Providers whose performance falls within the central range are flagged as \code{0}.
#' Outlying providers are determined by the test statistic falling beyond the threshold based on the significance level \eqn{1 - \text{level}}.
#'
#' @examples
#' data(ExampleDataLinear)
#' outcome <- ExampleDataLinear$Y
#' ProvID <- ExampleDataLinear$ProvID
#' covar <- ExampleDataLinear$Z
#' fit_re <- linear_re(Y = outcome, Z = covar, ProvID = ProvID)
#' test(fit_re)
#'
#' @importFrom stats pnorm qnorm pt qt
#'
#' @exportS3Method test linear_re

test.linear_re <- function(fit, parm, level = 0.95, null = 0, alternative = "two.sided", ...) {
  alpha <- 1 - level

  data <- fit$data_include
  Y.char <- fit$char_list$Y.char
  ProvID.char <- fit$char_list$ProvID.char
  Z.char <- fit$char_list$Z.char
  FEcoef <- fit$coefficient$FE

  var_alpha <- fit$variance$alpha
  sigma_sq <- fit$sigma^2
  n.prov <- sapply(split(data[, Y.char], data[, ProvID.char]), length)
  R_i <- as.vector(var_alpha) / (as.vector(var_alpha) + as.vector(sigma_sq) / n.prov)

  # Y_bar_i <- sapply(split(data[, Y.char], data[, ProvID.char]), mean)
  # Z_bar_i <- t(matrix(sapply(split(data[,Z.char, drop = FALSE], data[,ProvID.char]), colMeans),
  #                             ncol=length(Y_bar_i), nrow = length(beta)))
  #
  # estimate_alpha <- Y_bar_i - Z_bar_i%*%beta - mu

  denom <- sqrt(R_i * sigma_sq / n.prov)
  Z_score <- (fit$coefficient$RE - null)/denom
  p <- pnorm(Z_score, lower.tail=F)

  if (alternative == "two.sided") {
    p_value <- 2 * pmin(p, 1-p)
    flag <- ifelse(p < alpha/2, 1, ifelse(p > 1 - alpha/2, -1, 0))
  }
  else if (alternative == "greater") {
    p_value <- p
    flag <- ifelse(p < alpha, 1, 0)
  }
  else if (alternative == "less") {
    p_value <- 1 - p
    flag <- ifelse(1 - p < alpha, -1, 0)
  }
  else {
    stop("Argument 'alternative' should be one of 'two.sided', 'less', 'greater'")
  }

  result <- data.frame(flag = factor(flag), p = p_value, stat = Z_score, Std.Error = denom)
  colnames(result) <- c("flag", "p value", "stat", "Std.Error")

  if (missing(parm)) {
    attr(result, "provider size") <- n.prov
    return(result)
  }
  else {
    if (is.numeric(parm)) {  #avoid "integer" class
      parm <- as.numeric(parm)
    }
    if (class(parm) == class(data[, ProvID.char])) {
      attr(result, "provider size") <- n.prov[names(n.prov) %in% parm]
      result <- result[row.names(result) %in% parm, ]
      return(result)
    } else {
      stop("Argument 'parm' includes invalid elements!")
    }
  }
}
