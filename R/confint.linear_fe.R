#=== CONFIDENCE INTERVAL FUNCTION of Linear FE Model=================================
#' Provide confidence interval for provider effects or standardization ratios/rates
#'
#' @param fit an object as output of \code{linear_fe} function.
#'
#' @param parm specify a subset of which providers are to be given confidence intervals. All providers are included by default.
#'
#' @param level confidence level used for constructing confidence intervals. Defaulting to 0.95.
#'
#' @param option the confidence interval for the function's output, whether it is for gamma or standardization ratios/rates.
#'   \itemize{
#'   \item "gamma": provider effect
#'   \item "SR": standardization ratios/rates
#'   }
#'
#' @param stdz if option = 'SR', a character string specifying the standardization method. Defaulting to "indirect".
#'   \itemize{
#'   \item "indirect": using indirect standardized method
#'   \item "direct": using direct standardized method
#'   }
#'
#' @param ...
#'
#' @return A dataframe containing the point estimate, and lower and upper bounds of the estimate.
#'
#' @importFrom stats pnorm qnorm pt qt
#'
#' @exportS3Method confint linear_fe
#'

confint.linear_fe <- function(fit, parm, level = 0.95, ref.dis = "normal",
                              option = c("gamma", "SR"), null = "median", stdz = "indirect") {
  return_ls <- list()

  alpha <- 1 - level

  if (missing(fit)) stop ("Argument 'fit' is required!",call.=F)
  if (!class(fit) %in% c("linear_fe")) stop("Object fit is not of the classes 'linear_fe'!",call.=F)
  if (! "gamma" %in% option & !"SR" %in% option) stop("Argument 'option' NOT as required!", call.=F)

  data <- fit$data_include
  prov.name <- rownames(fit$coefficient$gamma)
  m <- length(fit$coefficient$gamma)
  p <- length(fit$coefficient$beta)
  n <- nrow(fit$data_include)
  n.prov <- sapply(split(data[, fit$char_list$Y.char], data[, fit$char_list$prov.char]), length)

  # CI of Gamma
  gamma <- fit$coefficient$gamma
  sd.gamma <- sqrt(fit$variance$gamma)

  if (ref.dis == "normal") {
    crit_value <- qnorm(1 - alpha / 2)
  } else if (ref.dis == "t") {
    df <- n - m - p
    crit_value <- qt(1 - alpha / 2, df)
  } else {
    stop("Reference distribution must be either 'normal' or 't'")
  }

  L_gamma <- gamma - crit_value * sd.gamma
  U_gamma <- gamma + crit_value * sd.gamma

  CI_gamma <- data.frame(gamma = gamma, gamma.Lower = L_gamma, gamma.Upper = U_gamma)
  colnames(CI_gamma) <- c("gamma", "gamma.Lower", "gamma.Upper")

  if (missing(parm)) {
    ind <- 1:length(prov.name)
  } else {
    if (is.numeric(parm)) {  #avoid "integer" class
      parm <- as.numeric(parm)
    }

    if (class(parm) == class(data[, fit$char_list$prov.char])) {
      ind <- which(prov.name %in% parm)
    } else {
      stop("Argument 'parm' includes invalid elements!")
    }
  }

  if ("gamma" %in% option) {return_ls$CI.gamma <- CI_gamma[ind, ]}

  # CI of SR
  if ("SR" %in% option) {
    if ("indirect" %in% stdz) {
      SR <- SR_linear(fit, stdz = "indirect", null)

      L.obs <- rep(L_gamma, n.prov) + fit$linear_pred
      L.prov <- sapply(split(L.obs, fit$prov), sum)
      L_indirect <- (L.prov - SR$OE$OE_indirect$Exp)/n.prov

      U.obs <- rep(U_gamma, n.prov) + fit$linear_pred
      U.prov <- sapply(split(U.obs, fit$prov), sum)
      U_indirect <- (U.prov - SR$OE$OE_indirect$Exp)/n.prov

      CI_indirect <- data.frame(SR = SR$indirect.difference, indirect.Lower = L_indirect, indirect.Upper = U_indirect)
      colnames(CI_indirect) <- c("Indirect.Difference", "indirect.Lower", "indirect.Upper")

      return_ls$CI.indirect <- CI_indirect[ind, ]
    }

    if ("direct" %in% stdz) {
      SR <- SR_linear(fit, stdz = "direct", null)

      Exp.direct <- function(gamma){
        sum(gamma + fit$linear_pred)
      }

      L.prov <- sapply(L_gamma, Exp.direct)
      L_direct <- (L.prov - SR$OE$OE_direct$Obs)/n

      U.prov <- sapply(U_gamma, Exp.direct)
      U_direct <- (U.prov - SR$OE$OE_direct$Obs)/n

      CI_direct <- data.frame(SR = SR$direct.difference, direct.Lower = L_direct, direct.Upper = U_direct)
      colnames(CI_direct) <- c("Direct.Difference", "direct.Lower", "direct.Upper")

      return_ls$CI.direct <- CI_direct[ind, ]
    }
  }

  return(return_ls)
}
