#' Get Confidence Intervals for Provider Effects or Standardized Measures
#'
#' Provide confidence interval for provider effects \eqn{\hat{\gamma}} or standardized measures from a fixed effect logistic model.
#'
#' @param fit a model fitted from \code{logis_fe}.
#' @param parm specifies a subset of providers for which confidence intervals are to be given.
#' By default, all providers are included. The class of `parm` should match the class of the provider IDs.
#' @param level the confidence level. The default value is 0.95.
#' @param test a character string specifying the type of testing method. The default is "exact".
#'   \itemize{
#'   \item {\code{"exact"}} exact test.
#'   \item {\code{"wald"}} wald test.
#'   \item {\code{"score"}} score test.
#'   }
#' @param option 	a character string specifying whether the confidence intervals
#' should be provided for provider effects, standardized measures:
#'   \itemize{
#'   \item {\code{"gamma"}} provider effect.
#'   \item {\code{"SM"}} standardized measures.
#'   }
#' @param stdz a character string or a vector specifying the standardization method
#' if \code{option = "SM"}. See `stdz` argument in \code{\link{SM_output.logis_fe}}.
#' @param null a character string or a number defining the population norm if \code{option = "SM"}.
#' @param measure a character string or a vector indicating whether the output measure is "ratio" or "rate" if \code{option = "SM"}.
#' Both "rate" and "ratio" will be provided by default.
#'   \itemize{
#'   \item {\code{"rate"}} output the standardized rate. The "rate" has been restricted to 0% - 100%.
#'   \item {\code{"ratio"}}  output the standardized ratio.
#'   \item {\code{c("ratio", "rate")}} output both the standardized rate and ratio.
#'   }
#' @param alternative a character string specifying the alternative hypothesis, must be one of
#' \code{"two.sided"} (default), \code{"greater"}, or \code{"less"}.
#'
#' @details
#' The wald test is invalid for extreme providers (i.e. when provider effect goes to infinity).
#' We suggest using score or exact test to generate the confidence interval.
#'
#' @return A dataframe (\code{option = "gamma"}) or a list of data frames (\code{option = "SM"}) containing the point estimate, and lower and upper bounds of the estimate.
#'
#' @examples
#' data(data_FE)
#' fit_fe <- logis_fe(Y = data_FE$Y, Z = data_FE$Z, ID = data_FE$ID, message = FALSE)
#' confint(fit_fe, option = "gamma")
#' confint(fit_fe, option = "SR")
#'
#' @importFrom stats plogis
#'
#' @exportS3Method confint logis_fe

confint.logis_fe <- function(fit, parm, level = 0.95, test = "exact",
                             option = "SM", stdz = "indirect", null = "median",
                             measure = c("rate", "ratio"), alternative = "two.sided") {
  if (missing(fit)) stop ("Argument 'fit' is required!",call.=F)
  if (!class(fit) %in% c("logis_fe")) stop("Object fit is not of the classes 'logis_fe'!",call.=F)
  if (! "gamma" %in% option & !"SM" %in% option) stop("Argument 'option' NOT as required!", call.=F)
  if (!(test %in% c("exact", "score", "wald"))) stop("Argument 'test' NOT as required!", call.=F)
  if (!"indirect" %in% stdz & !"direct" %in% stdz) stop("Argument 'stdz' NOT as required!", call.=F)

  alpha <- 1 - level

  Y.char <- fit$char_list$Y.char
  Z.char <- fit$char_list$Z.char
  ID.char <- fit$char_list$ID.char
  gamma <- fit$coefficient$gamma
  beta <- fit$coefficient$beta
  # df.prov <- fit$df.prov
  # names(gamma) <- rownames(df.prov)
  prov.order <- rownames(gamma)

  if (!missing(parm)) {
    if (is.numeric(parm)) {  #avoid "integer" class
      parm <- as.numeric(parm)
    }
  }

  #confidence of gamma
  confint_fe_gamma <- function(fit, test, parm, alpha, alternative) {
    data <- fit$data_include
    Obs_provider <- sapply(split(data[,Y.char],data[,ID.char]),sum)
    if (missing(parm)) {
      # pass
    } else if (class(parm)==class(data[,ID.char]) & test!="wald") {
      data <- data[data[,ID.char] %in% parm,]
    } else if (class(parm)==class(data[,ID.char]) & test=="wald") {
      indices <- which(unique(data[,ID.char]) %in% parm)
    } else {
      stop("Argument 'parm' includes invalid elements!")
    }

    if (test %in% c("score", "exact")) {
      if (test=="score") {
        qnorm.halfalpha <- qnorm(alpha/2, lower=F)
        qnorm.alpha <- qnorm(alpha, lower=F)
        CL.finite <- function(df, alternative) {
          prov <- ifelse(length(unique(df[,ID.char]))==1, unique(df[,ID.char]),
                         stop("Number of providers involved NOT equal to one!"))
          UL.gamma <- function(Gamma) { #increasing function w.r.t gamma_null
            p <- plogis(Gamma+Z.beta)
            if (alternative == "two.sided") {
              return((Obs - sum(p)) / sqrt(sum(p*(1-p))) + qnorm.halfalpha)
            }
            else if (alternative == "less") {
              return((Obs - sum(p)) / sqrt(sum(p*(1-p))) + qnorm.alpha)
            }
          }
          LL.gamma <- function(Gamma) {
            p <- plogis(Gamma+Z.beta)
            if (alternative == "two.sided") {
              return((Obs-sum(p)) / sqrt(sum(p*(1-p))) - qnorm.halfalpha)
            }
            else if (alternative == "greater") {
              return((Obs-sum(p)) / sqrt(sum(p*(1-p))) - qnorm.alpha)
            }
          }
          #Obs <- df.prov[as.character(prov), "Obs_provider"]  #Number of events for "prov"
          Obs <- Obs_provider[[as.character(prov)]]
          Z.beta <- as.matrix(df[,Z.char])%*%beta
          # gamma.lower <- uniroot(LL.gamma, gamma[as.character(prov),]+c(-5,0))$root
          # gamma.upper <- uniroot(UL.gamma, gamma[as.character(prov),]+c(0,5))$root
          max_attempts <- 3
          gamma.lower <- -Inf
          gamma.upper <- Inf
          if (!(alternative %in% c('two.sided', 'less', 'greater')))
            {stop("Argument 'alternative' should be one of 'two.sided', 'less', 'greater'.", call.=F)}
          if (alternative == "two.sided" | alternative == "less") {
            for (i in 0:(max_attempts-1)) {
              result_upper <- try(uniroot(UL.gamma, gamma[as.character(prov),]+c(5*i,5*(i+1))), silent = TRUE)
              if (class(result_upper)[1] == "try-error") {

              } else {
                gamma.upper <- result_upper$root
                break # Exit loop upon successful root finding
              }
            }
          }
          if (alternative == "two.sided" | alternative == "greater") {
            for (i in 0:(max_attempts-1)) {
              result_lower <- try(uniroot(LL.gamma, gamma[as.character(prov),]+c((-5)*(i+1),(-5)*i)), silent = TRUE)
              if (class(result_lower)[1] == "try-error") {

              } else {
                gamma.lower <- result_lower$root
                break # Exit loop upon successful root finding
              }
            }
          }
          return_mat <- c(gamma[as.character(prov),], gamma.lower, gamma.upper)
          return(return_mat)
        }
        CL.no.events <- function(df) { #only upper bound
          prov <- ifelse(length(unique(df[,ID.char]))==1, unique(df[,ID.char]),
                         stop("Number of providers involved NOT equal to one!"))
          Z.beta <- as.matrix(df[,Z.char])%*%beta
          max.Z.beta <- norm(Z.beta, "I")
          UL.gamma <- function(Gamma) {
            p <- plogis(Gamma+Z.beta)
            return(qnorm.alpha-sum(p)/sqrt(sum(p*(1-p))))
          }
          #gamma.upper <- uniroot(UL.gamma,(10+max.Z.beta)*c(-1,1))$root
          max_attempts <- 3
          gamma.upper <- Inf
          for (i in 1:max_attempts){
            result_upper <- try(uniroot(UL.gamma,(10+max.Z.beta)*c(-i,i)), silent = T)
            if (class(result_upper)[1] == "try-error") {

            } else {
              gamma.upper <- result_upper$root
              break # Exit loop upon successful root finding
            }
          }
          return_mat <- c(gamma[as.character(prov),], -Inf, gamma.upper)
          return(return_mat)
        }
        CL.all.events <- function(df) { #only lower bound
          prov <- ifelse(length(unique(df[,ID.char]))==1, unique(df[,ID.char]),
                         stop("Number of providers involved NOT equal to one!"))
          Z.beta <- as.matrix(df[,Z.char])%*%beta
          max.Z.beta <- norm(Z.beta, "I")
          LL.gamma <- function(Gamma) {
            p <- plogis(Gamma+Z.beta)
            pq <- p*(1-p)
            pq[pq == 0] <- 1e-20
            return(sum(1-p)/sqrt(sum(pq))-qnorm.alpha)
          }
          #gamma.lower <- uniroot(LL.gamma,(10+max.Z.beta)*c(-1,1))$root
          max_attempts <- 3
          gamma.lower <- -Inf
          for (i in 1:max_attempts) {
            result_lower <- try(uniroot(LL.gamma,(10+max.Z.beta)*c(-i,i)), silent = TRUE)
            if (class(result_lower)[1] == "try-error") {

            } else {
              gamma.lower <- result_lower$root
              break # Exit loop upon successful root finding
            }
          }
          return_mat <- c(gamma[as.character(prov),], gamma.lower, Inf)
          return(return_mat)
        }
      } else {
        CL.finite <- function(df, alternative) {
          UL.gamma <- function(Gamma) {
            if (alternative == "two.sided") {
              poibin::ppoibin(Obs-1,plogis(Gamma+Z.beta))+0.5*poibin::dpoibin(Obs,plogis(Gamma+Z.beta))-alpha/2
            }
            else if (alternative == "less") {
              poibin::ppoibin(Obs,plogis(Gamma+Z.beta))-alpha
            }
          }
          LL.gamma <- function(Gamma) {
            if (alternative == "two.sided") {
              1-poibin::ppoibin(Obs,plogis(Gamma+Z.beta))+0.5*poibin::dpoibin(Obs,plogis(Gamma+Z.beta))-alpha/2
            }
            else if (alternative == "greater") {
              1-poibin::ppoibin(Obs-1,plogis(Gamma+Z.beta))-alpha
            }
          }
          prov <- ifelse(length(unique(df[,ID.char]))==1, unique(df[,ID.char]),
                         stop("Number of providers involved NOT equal to one!"))
          # Obs <- df.prov[as.character(prov), "Obs_provider"]
          Obs <- Obs_provider[[as.character(prov)]]
          Z.beta <- as.matrix(df[,Z.char])%*%beta
          # gamma.lower <- uniroot(LL.gamma, gamma[as.character(prov),]+c(-5,0))$root
          # gamma.upper <- uniroot(UL.gamma, gamma[as.character(prov),]+c(0,5))$root
          max_attempts <- 3
          gamma.lower <- -Inf
          gamma.upper <- Inf
          if (!(alternative %in% c('two.sided', 'less', 'greater')))
            {stop("Argument 'alternative' should be one of 'two.sided', 'less', 'greater'.", call.=F)}
          if (alternative == "two.sided" | alternative == "less") {
            for (i in 0:(max_attempts-1)) {
              result_upper <- try(uniroot(UL.gamma, gamma[as.character(prov),]+c(5*i,5*(i+1))), silent = TRUE)
              if (class(result_upper)[1] == "try-error") {

              } else {
                gamma.upper <- result_upper$root
                break # Exit loop upon successful root finding
              }
            }
          }
          if (alternative == "two.sided" | alternative == "greater") {
            for (i in 0:(max_attempts-1)) {
              result_lower <- try(uniroot(LL.gamma, gamma[as.character(prov),]+c((-5)*(i+1),(-5)*i)), silent = TRUE)
              if (class(result_lower)[1] == "try-error") {

              } else {
                gamma.lower <- result_lower$root
                break # Exit loop upon successful root finding
              }
            }
          }
          return_mat <- c(gamma[as.character(prov),], gamma.lower, gamma.upper)
          return(return_mat)
        }
        CL.no.events <- function(df) {
          prov <- ifelse(length(unique(df[,ID.char]))==1, unique(df[,ID.char]),
                         stop("Number of providers involved NOT equal to one!"))
          Z.beta <- as.matrix(df[,Z.char])%*%beta
          max.Z.beta <- norm(Z.beta, "I")
          # gamma.upper <- uniroot(function(x) prod(plogis(-x-Z.beta))/2-alpha,
          #                        (10+max.Z.beta)*c(-1,1))$root
          max_attempts <- 3
          gamma.upper <- Inf
          for (i in 1:max_attempts){
            result_upper <- try(uniroot(function(x) prod(plogis(-x-Z.beta))/2-alpha,
                                        (10+max.Z.beta)*c(-i,i)), silent = T)
            if (class(result_upper)[1] == "try-error") {

            } else {
              gamma.upper <- result_upper$root
              break # Exit loop upon successful root finding
            }
          }
          return_mat <- c(gamma[as.character(prov),], -Inf, gamma.upper)
          return(return_mat)
        }
        CL.all.events <- function(df) {
          prov <- ifelse(length(unique(df[,ID.char]))==1, unique(df[,ID.char]),
                         stop("Number of providers involved NOT equal to one!"))
          Z.beta <- as.matrix(df[,Z.char])%*%beta
          max.Z.beta <- norm(Z.beta, "I")
          # gamma.lower <- uniroot(function(x) prod(plogis(x+Z.beta))/2-alpha,
          #                        (10+max.Z.beta)*c(-1,1))$root
          max_attempts <- 3
          gamma.lower <- -Inf
          for (i in 1:max_attempts) {
            result_lower <- try(uniroot(function(x) prod(plogis(x+Z.beta))/2-alpha,
                                        (10+max.Z.beta)*c(-i,i)), silent = TRUE)
            if (class(result_lower)[1] == "try-error") {

            } else {
              gamma.lower <- result_lower$root
              break # Exit loop upon successful root finding
            }
          }
          return_mat <- c(gamma[as.character(prov),], gamma.lower, Inf)
          return(return_mat)
        }
      }
      confint.finite <- sapply(by(data[(data$no.events==0) & (data$all.events==0),],
                                  data[(data$no.events==0) & (data$all.events==0),ID.char],identity),
                               FUN=function(df) CL.finite(df, alternative))
      prov_finite <- unique(data[(data$no.events==0) & (data$all.events==0),ID.char])
      confint.no.events <- sapply(by(data[data$no.events==1,], data[data$no.events==1,ID.char],identity),
                                  FUN=function(df) CL.no.events(df))
      prov_no.events <- unique(data[data$no.events==1,ID.char])
      confint.all.events <- sapply(by(data[data$all.events==1,], data[data$all.events==1,ID.char],identity),
                                   FUN=function(df) CL.all.events(df))
      prov_all.events <- unique(data[data$all.events==1,ID.char])
      confint_df <- as.numeric(cbind(confint.finite, confint.no.events, confint.all.events))
      confint_df <- as.data.frame(matrix(confint_df, ncol = 3, byrow = T))
      colnames(confint_df) <- c("gamma", "gamma.lower", "gamma.upper")
      rownames(confint_df) <- c(prov_finite, prov_no.events, prov_all.events)
      # return(confint_df[order(match(rownames(confint_df), prov.order)),])
      return(confint_df[order(as.numeric(rownames(confint_df))),])
    } else if (test=="wald") {
      warning("The Wald test fails for datasets with providers having all or no events. Score test or exact test are recommended.")
      if (missing(parm)) {
        indices <- 1:length(prov.order)
      }

      se.gamma <- sqrt(fit$variance$gamma)
      if (alternative == "two.sided") {
        crit_value <- qnorm(1 - alpha / 2)
        U_gamma <- gamma + crit_value * se.gamma
        L_gamma <- gamma - crit_value * se.gamma
      }
      else if (alternative == "greater") {
        crit_value <- qnorm(1 - alpha)
        U_gamma <- Inf
        L_gamma <- gamma - crit_value * se.gamma
      }
      else if (alternative == "less") {
        crit_value <- qnorm(1 - alpha)
        U_gamma <- gamma + crit_value * se.gamma
        L_gamma <- -Inf
      }
      else {
        stop("Argument 'alternative' should be one of 'two.sided', 'less', 'greater'.")
      }

      return_mat <- data.frame(gamma, L_gamma, U_gamma,
                               row.names=rownames(gamma))
      colnames(return_mat) <- c("gamma", "gamma.lower", "gamma.upper")
      return(return_mat[indices,])
    }
  }
  if (option == "gamma"){
    if (alternative != "two.sided")
      stop("Provider effect (option = 'gamma') only supports two-sided confidence intervals.", call. = FALSE)
    return_mat <- confint_fe_gamma(fit, test, parm, alpha, alternative)
    return(return_mat)
  }
  else if (option == "SM"){
    data.ori <- fit$data_include
    population_rate <- sum(data.ori[,Y.char])/nrow(data.ori) * 100  #sum(O_i)/N *100%
    return_ls <- list()
    if ("indirect" %in% stdz) {
      SR.indirect <- SM_output(fit, stdz = c("indirect"), measure = c("ratio", "rate"), null = null)
      if (missing(parm)) {
        OE_df.indirect <- SR.indirect$OE$OE_indirect
        indirect.ratio_df <- SR.indirect$indirect.ratio
        indirect.rate_df <- SR.indirect$indirect.rate
        data <- data.ori
      } else if (class(parm) == class(data.ori[,ID.char])) {
        OE_df.indirect <- SR.indirect$OE$OE_indirect[rownames(SR.indirect$OE$OE_indirect) %in% parm, , drop = FALSE]
        data <- data.ori[data.ori[,ID.char] %in% parm,]
        indirect.ratio_df <- SR.indirect$indirect.ratio[rownames(SR.indirect$indirect.ratio) %in% parm, , drop = FALSE]
        indirect.rate_df <- SR.indirect$indirect.rate[rownames(SR.indirect$indirect.rate) %in% parm, , drop = FALSE]
      } else {
        stop("Argument 'parm' includes invalid elements!")
      }
      #functions for calculate CI of SRs
      qnorm.halfalpha <- qnorm(alpha/2, lower=F)
      qnorm.alpha <- qnorm(alpha, lower=F)
      SR_indirect.finite <- function(df) {
        prov <- ifelse(length(unique(df[,ID.char]))==1, unique(df[,ID.char]),
                       stop("Number of providers involved NOT equal to one!"))
        Z.beta <- as.matrix(df[,Z.char])%*%beta
        confint_gamma <- confint_fe_gamma(fit, test = test, parm = unique(df$ID), alpha = alpha, alternative = alternative)
        gamma.lower <- confint_gamma$gamma.lower
        gamma.upper <- confint_gamma$gamma.upper
        EXP.i <- OE_df.indirect[rownames(OE_df.indirect) == unique(df[,ID.char]), "Exp.indirect_provider"]
        SR.lower <- sum(plogis(gamma.lower+Z.beta)) / EXP.i
        SR.upper <- sum(plogis(gamma.upper+Z.beta)) / EXP.i
        return(c(SR.lower, SR.upper))
      }
      SR_indirect.no.events <- function(df) {
        prov <- ifelse(length(unique(df[,ID.char]))==1, unique(df[,ID.char]),
                       stop("Number of providers involved NOT equal to one!"))
        Z.beta <- as.matrix(df[,Z.char])%*%beta
        confint_gamma <- confint_fe_gamma(fit, test = test, parm = unique(df$ID), alpha = alpha, alternative = alternative)
        gamma.upper <- confint_gamma$gamma.upper
        EXP.i <- OE_df.indirect[rownames(OE_df.indirect) == unique(df[,ID.char]), "Exp.indirect_provider"]
        if (alternative == "greater") {
          SR.upper <- nrow(df) / EXP.i
        }
        else {
          SR.upper <- sum(plogis(gamma.upper+Z.beta)) / EXP.i
        }
        return(c(0, SR.upper))
      }
      SR_indirect.all.events <- function(df) {
        prov <- ifelse(length(unique(df[,ID.char]))==1, unique(df[,ID.char]),
                       stop("Number of providers involved NOT equal to one!"))
        Z.beta <- as.matrix(df[,Z.char])%*%beta
        confint_gamma <- confint_fe_gamma(fit, test = test, parm = unique(df$ID), alpha = alpha, alternative = alternative)
        gamma.lower <- confint_gamma$gamma.lower
        EXP.i <- OE_df.indirect[rownames(OE_df.indirect) == unique(df[,ID.char]), "Exp.indirect_provider"]
        if (alternative == "less") {
          SR.lower <- 0
        }
        else {
          SR.lower <- sum(plogis(gamma.lower+Z.beta)) / EXP.i
        }
        SR.upper <- nrow(df) / EXP.i
        return(c(SR.lower, SR.upper))
      }

      confint.finite <- sapply(by(data[(data$no.events==0) & (data$all.events==0),],
                                  data[(data$no.events==0) & (data$all.events==0),ID.char],identity),
                               FUN=function(df) SR_indirect.finite(df))
      prov_finite <- unique(data[(data$no.events==0) & (data$all.events==0),ID.char])
      confint.no.events <- sapply(by(data[data$no.events==1,], data[data$no.events==1,ID.char],identity),
                                  FUN=function(df) SR_indirect.no.events(df))
      prov_no.events <- unique(data[data$no.events==1,ID.char])
      confint.all.events <- sapply(by(data[data$all.events==1,], data[data$all.events==1,ID.char],identity),
                                   FUN=function(df) SR_indirect.all.events(df))
      prov_all.events <- unique(data[data$all.events==1,ID.char])
      CI.combined <- cbind(confint.finite, confint.no.events, confint.all.events)
      if (ncol(CI.combined) > 1) {
        CI.combined <- CI.combined[, order(as.numeric(colnames(CI.combined)))]
      }

      CI.indirect_ratio <- as.numeric(rbind(t(indirect.ratio_df), CI.combined))
      CI.indirect_ratio <- as.data.frame(matrix(CI.indirect_ratio, ncol = 3, byrow = T))
      colnames(CI.indirect_ratio) <- c("indirect_ratio", "CI_ratio.lower", "CI_ratio.upper")
      rownames(CI.indirect_ratio) <- rownames(indirect.ratio_df)

      if ("ratio" %in% measure){
        attr(CI.indirect_ratio, "confidence_level") <- paste(level * 100, "%")
        attr(CI.indirect_ratio, "type") <- ifelse(alternative == "greater", "upper one-sided",
                                            ifelse(alternative == "less", "lower one-sided",
                                                   "two-sided"))
        attr(CI.indirect_ratio, "description") <- "Indirect Standardized Ratio"
        attr(CI.indirect_ratio, "model") <- "FE logis"
        return_ls$CI.indirect_ratio <- CI.indirect_ratio
      }

      if ("rate" %in% measure){
        rate.lower <- pmax(pmin(CI.indirect_ratio$CI_ratio.lower * population_rate, 100), 0)
        rate.upper <- pmax(pmin(CI.indirect_ratio$CI_ratio.upper * population_rate, 100), 0)
        CI.indirect_rate <- as.data.frame(cbind(indirect.rate_df, rate.lower, rate.upper))
        colnames(CI.indirect_rate) <- c("indirect_rate", "CI_rate.lower", "CI_rate.upper")

        attr(CI.indirect_rate, "confidence_level") <- paste(level * 100, "%")
        attr(CI.indirect_rate, "type") <- ifelse(alternative == "greater", "upper one-sided",
                                                  ifelse(alternative == "less", "lower one-sided",
                                                         "two-sided"))
        attr(CI.indirect_rate, "description") <- "Indirect Standardized Rate"
        attr(CI.indirect_rate, "model") <- "FE logis"
        attr(CI.indirect_rate, "population_rate") <- population_rate
        return_ls$CI.indirect_rate <- CI.indirect_rate
      }
    }


    if ("direct" %in% stdz) {
      SR.direct <- SM_output(fit, stdz = c("direct"), measure = c("ratio", "rate"), null = null)
      if (missing(parm)) {
        OE_df.direct <- SR.direct$OE$OE_direct[1,1]
        direct.ratio_df <- SR.direct$direct.ratio
        direct.rate_df <- SR.direct$direct.rate
        data <- data.ori
      } else if (class(parm) == class(data[,ID.char])) {
        OE_df.direct <- SR.direct$OE$OE_direct[1,1]
        direct.ratio_df <- SR.direct$direct.ratio[rownames(SR.direct$direct.ratio) %in% parm, , drop = FALSE]
        direct.rate_df <- SR.direct$direct.rate[rownames(SR.direct$direct.rate) %in% parm, , drop = FALSE]
        data <- data.ori[data.ori[,ID.char] %in% parm,]
      } else {
        stop("Argument 'parm' includes invalid elements!")
      }
      #funcitons for calculate CI of SRs
      qnorm.halfalpha <- qnorm(alpha/2, lower=F)
      qnorm.alpha <- qnorm(alpha, lower=F)

      SR_direct.finite <- function(ID) {
        Z.beta.all <- as.matrix(data.ori[,Z.char])%*%beta
        confint_gamma <- confint_fe_gamma(fit, test = test, parm = ID, alpha = alpha, alternative = alternative)
        gamma.lower <- confint_gamma$gamma.lower
        gamma.upper <- confint_gamma$gamma.upper
        SR.lower <- sum(plogis(gamma.lower+Z.beta.all)) / OE_df.direct
        SR.upper <- sum(plogis(gamma.upper+Z.beta.all)) / OE_df.direct
        return(c(SR.lower, SR.upper))
      }
      SR_direct.no.events <- function(ID) {
        Z.beta.all <- as.matrix(data.ori[,Z.char])%*%beta
        confint_gamma <- confint_fe_gamma(fit, test = test, parm = ID, alpha = alpha, alternative = alternative)
        gamma.upper <- confint_gamma$gamma.upper
        if (alternative == "greater") {
          SR.upper <- nrow(data.ori)/sum(data.ori[,Y.char])
        }
        else {
          SR.upper <- sum(plogis(gamma.upper+Z.beta.all)) / OE_df.direct
        }
        return(c(0, SR.upper))
      }
      SR_direct.all.events <- function(ID) {
        Z.beta.all <- as.matrix(data.ori[,Z.char])%*%beta
        confint_gamma <- confint_fe_gamma(fit, test = test, parm = ID, alpha = alpha, alternative = alternative)
        gamma.lower <- confint_gamma$gamma.lower
        if (alternative == "less") {
          SR.lower <- 0
        }
        else {
          SR.lower <- sum(plogis(gamma.lower+Z.beta.all)) / OE_df.direct
        }
        SR.upper <- nrow(data.ori) / OE_df.direct
        return(c(SR.lower, SR.upper))
      }

      confint.finite <- sapply(unique(data[(data$no.events==0) & (data$all.events==0),]$ID),
                               FUN = function(ID) SR_direct.finite(ID))
      prov_finite <- unique(data[(data$no.events==0) & (data$all.events==0),ID.char])
      confint.no.events <- sapply(unique(data[(data$no.events==1) & (data$all.events==0),]$ID),
                                  FUN = function(ID) SR_direct.no.events(ID))
      prov_no.events <- unique(data[data$no.events==1,ID.char])
      confint.all.events <- sapply(unique(data[(data$no.events==0) & (data$all.events==1),]$ID),
                                   FUN = function(ID) SR_direct.all.events(ID))
      prov_all.events <- unique(data[data$all.events==1,ID.char])
      CI.combined <- cbind(confint.finite, confint.no.events, confint.all.events)
      colnames(CI.combined) <- c(prov_finite, prov_no.events, prov_all.events)
      if (ncol(CI.combined) > 1) {
        CI.combined <- CI.combined[, order(as.numeric(colnames(CI.combined)))]
      }

      CI.direct_ratio <- as.numeric(rbind(t(direct.ratio_df), CI.combined))
      CI.direct_ratio <- as.data.frame(matrix(CI.direct_ratio, ncol = 3, byrow = T))
      colnames(CI.direct_ratio) <- c("direct_ratio", "CI_ratio.lower", "CI_ratio.upper")
      rownames(CI.direct_ratio) <- rownames(direct.ratio_df)

      if ("ratio" %in% measure){
        attr(CI.direct_ratio, "confidence_level") <- paste(level * 100, "%")
        attr(CI.direct_ratio, "type") <- ifelse(alternative == "greater", "upper one-sided",
                                                  ifelse(alternative == "less", "lower one-sided",
                                                         "two-sided"))
        attr(CI.direct_ratio, "description") <- "Direct Standardized Ratio"
        attr(CI.direct_ratio, "model") <- "FE logis"
        attr(CI.direct_ratio, "population_rate") <- population_rate
        return_ls$CI.direct_ratio <- CI.direct_ratio
      }

      if ("rate" %in% measure){
        rate.lower <- pmax(pmin(CI.direct_ratio$CI_ratio.lower * population_rate, 100), 0)
        rate.upper <- pmax(pmin(CI.direct_ratio$CI_ratio.upper * population_rate, 100), 0)
        CI.direct_rate <- as.data.frame(cbind(direct.rate_df, rate.lower, rate.upper))
        colnames(CI.direct_rate) <- c("direct_rate", "CI_rate.lower", "CI_rate.upper")

        attr(CI.direct_rate, "confidence_level") <- paste(level * 100, "%")
        attr(CI.direct_rate, "type") <- ifelse(alternative == "greater", "upper one-sided",
                                                 ifelse(alternative == "less", "lower one-sided",
                                                        "two-sided"))
        attr(CI.direct_rate, "description") <- "Direct Standardized Rate"
        attr(CI.direct_rate, "model") <- "FE logis"
        attr(CI.direct_rate, "population_rate") <- population_rate
        return_ls$CI.direct_rate <- CI.direct_rate
      }
    }
    return(return_ls)
  }


}
