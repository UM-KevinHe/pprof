#' Get Confidence Intervals for Provider Effects or Standardized Measures
#'
#' Provide confidence interval for provider effects \eqn{\hat{\gamma}} or standardized measures from a fixed effect linear model.
#'
#' @param fit a model fitted from \code{linear_fe}.
#' @param parm specifies a subset of providers for which confidence intervals are to be given.
#' By default, all providers are included. The class of `parm` should match the class of the provider IDs.
#' @param level the confidence level. The default value is 0.95.
#' @param option 	a character string specifying whether the confidence intervals
#' should be provided for provider effects, standardized measures:
#'   \itemize{
#'   \item {\code{"gamma"}} provider effect (only supports \code{"two.sided"} confidence interval).
#'   \item {\code{"SM"}} standardized measures.
#'   }
#' @param stdz a character string or a vector specifying the standardization method
#' if `option` includes \code{"SM"}. See `stdz` argument in \code{\link{SM_output.linear_fe}}.
#' @param null a character string or a number specifying the population norm for calculating standardized measures
#' if `option` includes \code{"SM"}. See `null` argument in \code{\link{SM_output.linear_fe}}.
#' @param alternative a character string specifying the alternative hypothesis, must be one of
#' \code{"two.sided"} (default), \code{"greater"}, or \code{"less"}.
#' Note that \code{"gamma"} for argument `option` only supports \code{"two.sided"}.
#'
#' @return A list of data frames containing the confidence intervals based on the values of `option` and `stdz`.
#' \item{CI.gamma}{Confidence intervals for provider effects if `option` includes \code{"gamma"}.}
#' \item{CI.indirect}{Confidence intervals for indirect standardized difference if `option` includes \code{"SM"} and `stdz` includes \code{"indirect"}.}
#' \item{CI.direct}{Confidence intervals for direct standardized difference if `option` includes \code{"SM"} and `stdz` includes \code{"direct"}.}
#'
#' @examples
#' data(ExampleDataLinear)
#' Y <- ExampleDataLinear$Y
#' ID <- ExampleDataLinear$ID
#' Z <- ExampleDataLinear$Z
#'
#' fit_fe <- linear_fe(Y = Y, Z = Z, ID = ID)
#' confint(fit_fe)
#'
#' @importFrom stats pnorm qnorm pt qt
#'
#' @exportS3Method confint linear_fe

confint.linear_fe <- function(fit, parm, level = 0.95, option = "SM", stdz = "indirect",
                              null = "median", alternative = "two.sided") {
  return_ls <- list()

  alpha <- 1 - level

  if (missing(fit)) stop ("Argument 'fit' is required!",call.=F)
  if (!class(fit) %in% c("linear_fe")) stop("Object fit is not of the classes 'linear_fe'!",call.=F)
  if (! "gamma" %in% option & !"SM" %in% option) stop("Argument 'option' NOT as required!", call.=F)
  if (!"indirect" %in% stdz & !"direct" %in% stdz) stop("Argument 'stdz' NOT as required!", call.=F)
  if ("gamma" %in% option && alternative != "two.sided")
    stop("Provider effect (option = 'gamma') only supports two-sided confidence intervals.", call. = FALSE)

  data <- fit$data_include
  prov <- data[ ,fit$char_list$ID.char]
  prov.name <- rownames(fit$coefficient$gamma)
  m <- length(fit$coefficient$gamma)
  p <- length(fit$coefficient$beta)
  n <- nrow(fit$data_include)
  n.prov <- sapply(split(data[, fit$char_list$Y.char], data[, fit$char_list$ID.char]), length)

  # CI of Gamma
  gamma <- fit$coefficient$gamma
  se.gamma <- sqrt(fit$variance$gamma)

  if (alternative == "two.sided") {
    crit_value <- ifelse(fit$method == "Profile Likelihood", qnorm(1 - alpha / 2),
                         crit_value <- qt(1 - alpha / 2, df = n - m - p))
    U_gamma <- gamma + crit_value * se.gamma
    L_gamma <- gamma - crit_value * se.gamma
  }
  else if (alternative == "greater") {
    crit_value <- ifelse(fit$method == "Profile Likelihood", qnorm(1 - alpha),
                         crit_value <- qt(1 - alpha, df = n - m - p))
    U_gamma <- Inf
    L_gamma <- gamma - crit_value * se.gamma
  }
  else if (alternative == "less") {
    crit_value <- ifelse(fit$method == "Profile Likelihood", qnorm(1 - alpha),
                         crit_value <- qt(1 - alpha, df = n - m - p))
    U_gamma <- gamma + crit_value * se.gamma
    L_gamma <- -Inf
  }
  else {
    stop("Argument 'alternative' should be one of 'two.sided', 'less', 'greater'.")
  }

  CI_gamma <- data.frame(gamma = gamma, gamma.Lower = L_gamma, gamma.Upper = U_gamma)
  colnames(CI_gamma) <- c("gamma", "gamma.Lower", "gamma.Upper")

  if (missing(parm)) {
    ind <- 1:length(prov.name)
  } else {
    if (is.numeric(parm)) {  #avoid "integer" class
      parm <- as.numeric(parm)
    }
    if (class(parm) == class(data[, fit$char_list$ID.char])) {
      ind <- which(prov.name %in% parm)
    } else {
      stop("Argument 'parm' includes invalid elements.")
    }
  }

  if (option == "gamma") return (CI_gamma[ind, ])

  # CI of SM
  if (option == "SM") {
    if ("indirect" %in% stdz) {
      SM <- SM_output(fit, stdz = "indirect", null = null)

      if (alternative == "two.sided") {
        L.obs <- rep(L_gamma, n.prov) + fit$linear_pred
        L.prov <- sapply(split(L.obs, prov), sum)
        L_indirect <- (L.prov - SM$OE$OE_indirect$Exp)/n.prov

        U.obs <- rep(U_gamma, n.prov) + fit$linear_pred
        U.prov <- sapply(split(U.obs, prov), sum)
        U_indirect <- (U.prov - SM$OE$OE_indirect$Exp)/n.prov
      }
      else if (alternative == "greater") {
        L.obs <- rep(L_gamma, n.prov) + fit$linear_pred
        L.prov <- sapply(split(L.obs, prov), sum)
        L_indirect <- (L.prov - SM$OE$OE_indirect$Exp)/n.prov

        U_indirect <- Inf
      }
      else if (alternative == "less") {
        U.obs <- rep(U_gamma, n.prov) + fit$linear_pred
        U.prov <- sapply(split(U.obs, prov), sum)
        U_indirect <- (U.prov - SM$OE$OE_indirect$Exp)/n.prov

        L_indirect <- -Inf
      }
      else {
        stop("Argument 'alternative' should be one of 'two.sided', 'less', 'greater'.")
      }

      CI_indirect <- data.frame(SM = SM$indirect.difference, indirect.Lower = L_indirect, indirect.Upper = U_indirect)
      colnames(CI_indirect) <- c("Indirect.Difference", "indirect.Lower", "indirect.Upper")
      return_ls$CI.indirect <- CI_indirect[ind, ]
    }

    if ("direct" %in% stdz) {
      SM <- SM_output(fit, stdz = "direct", null = null)

      Exp.direct <- function(gamma){
        sum(gamma + fit$linear_pred)
      }

      if (alternative == "two.sided") {
        L.prov <- sapply(L_gamma, Exp.direct)
        L_direct <- (L.prov - SM$OE$OE_direct$Obs)/n
        U.prov <- sapply(U_gamma, Exp.direct)
        U_direct <- (U.prov - SM$OE$OE_direct$Obs)/n
      }
      else if (alternative == "greater") {
        L.prov <- sapply(L_gamma, Exp.direct)
        L_direct <- (L.prov - SM$OE$OE_direct$Obs)/n
        U_direct <- Inf
      }
      else if (alternative == "less") {
        U.prov <- sapply(U_gamma, Exp.direct)
        U_direct <- (U.prov - SM$OE$OE_direct$Obs)/n
        L_direct <- -Inf
      }
      else {
        stop("Argument 'alternative' should be one of 'two.sided', 'less', 'greater'.")
      }

      CI_direct <- data.frame(SM = SM$direct.difference, direct.Lower = L_direct, direct.Upper = U_direct)
      colnames(CI_direct) <- c("Direct.Difference", "direct.Lower", "direct.Upper")
      return_ls$CI.direct <- CI_direct[ind, ]
    }
  }

  return(return_ls)
}
