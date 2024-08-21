#=== CONFIDENCE INTERVAL FUNCTION of Linear FE Model=================================
#' Provide confidence interval for provider effects or standardization ratios/rates
#'
#' @param fit an object as output of \code{linear_re} function.
#'
#' @param parm specify a subset of which providers are to be given confidence intervals. All providers are included by default.
#'
#' @param level confidence level used for constructing confidence intervals. Defaulting to 0.95.
#'
#' @param option the confidence interval for the function's output, whether it is for gamma or standardization ratios/rates.
#'   \itemize{
#'   \item "alpha": provider effect
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
#' @importFrom stats pnorm pt
#'
#' @exportS3Method confint linear_re
#'

confint.linear_re <- function(fit, parm, level = 0.95, option = c("alpha", "SR"),
                              null = "median", stdz = "indirect") {
  return_ls <- list()

  alpha <- 1 - level

  if (missing(fit)) stop ("Argument 'fit' is required!",call.=F)
  if (!class(fit) %in% c("linear_re")) stop("Object fit is not of the classes 'linear_re'!",call.=F)
  if (! "alpha" %in% option & !"SR" %in% option) stop("Argument 'option' NOT as required!", call.=F)

  data <- fit$data_include
  n <- nrow(data)
  Y.char <- fit$char_list$Y.char
  ID.char <- fit$char_list$ID.char
  Z.char <- fit$char_list$Z.char
  prov.name <- rownames(fit$coefficient$alpha)
  mu <- fit$coefficient$mu

  var_alpha <- fit$variance$alpha
  sigma_sq <- fit$sigma^2
  n.prov <- sapply(split(data[, Y.char], data[, ID.char]), length)
  R_i <- as.vector(var_alpha) / (as.vector(var_alpha) + as.vector(sigma_sq) / n.prov)

  se.alpha <- sqrt(R_i * sigma_sq / n.prov)
  crit_value <- qnorm(1 - alpha / 2)

  L_alpha <- fit$coefficient$alpha - crit_value * se.alpha
  U_alpha <- fit$coefficient$alpha + crit_value * se.alpha

  CI_alpha <- data.frame(alpha = fit$coefficient$alpha, alpha.Lower = L_alpha, alpha.Upper = U_alpha)
  colnames(CI_alpha) <- c("Estimate", "alpha.Lower", "alpha.Upper")

  if (missing(parm)) {
    ind <- 1:length(prov.name)
  } else {
    if (is.numeric(parm)) {  #avoid "integer" class
      parm <- as.numeric(parm)
    }

    if (class(parm) == class(data[, ID.char])) {
      ind <- which(prov.name %in% parm)
    } else {
      stop("Argument 'parm' includes invalid elements!")
    }
  }

  if ("alpha" %in% option) {return_ls$CI.alpha <- CI_alpha[ind, ]}

  # CI of SR
  if ("SR" %in% option) {
    if ("indirect" %in% stdz) {
      SR <- SR_linear(fit, stdz = "indirect")

      L.obs <- rep(mu,n) + rep(L_alpha, n.prov) + fit$linear_pred
      L.prov <- sapply(split(L.obs, fit$prov), sum)
      L_indirect <- (L.prov - SR$OE$OE_indirect$Exp)/n.prov

      U.obs <- rep(mu,n) + rep(U_alpha, n.prov) + fit$linear_pred
      U.prov <- sapply(split(U.obs, fit$prov), sum)
      U_indirect <- (U.prov - SR$OE$OE_indirect$Exp)/n.prov

      CI_indirect <- data.frame(SR = SR$indirect.difference, indirect.Lower = L_indirect, indirect.Upper = U_indirect)
      colnames(CI_indirect) <- c("Indirect.Difference", "indirect.Lower", "indirect.Upper")

      return_ls$CI.indirect <- CI_indirect[ind, ]
    }

    if ("direct" %in% stdz) {
      SR <- SR_linear(fit, stdz = "direct")

      mu = as.vector(mu)
      Exp.direct <- function(alpha){
        sum(mu + alpha + fit$linear_pred)
      }

      L.prov <- sapply(L_alpha, Exp.direct)
      L_direct <- (L.prov - SR$OE$OE_direct$Obs)/n

      U.prov <- sapply(U_alpha, Exp.direct)
      U_direct <- (U.prov - SR$OE$OE_direct$Obs)/n

      CI_direct <- data.frame(SR = SR$direct.difference, direct.Lower = L_direct, direct.Upper = U_direct)
      colnames(CI_direct) <- c("Direct.Difference", "direct.Lower", "direct.Upper")

      return_ls$CI.direct <- CI_direct[ind, ]
    }
  }

  return(return_ls)
}
