#=== FUNCTIONS FOR FITTING THE Fixed Effects MODEL =================================
#' Main function for fitting fixed effects model
#'
#' @param Y a numerical vector, with values of 0 or 1, indicating the outcome variable.
#'
#' @param Z a matrix or data frame containing covariates.
#'
#' @param ID a vector representing the provider id. Its elements can be either numeric values or characters.
#'
#' @param algorithm a string specifying the algorithm to be used. Defaulting to "SerBIN".
#'   \itemize{
#'   \item "SerBIN": using the Serial blockwise inversion Newton algorithm to fit the model (See [Wu et al. (2022)](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.9387)).
#'   \item "BAN": using the block ascent Newton algorithm to fit the model (See [He et al. (2013)](https://link.springer.com/article/10.1007/s10985-013-9264-6)).
#'   }
#' @param max.iter maximum number of iterations. Defaulting to 10,000.
#'
#' @param tol a small positive number specifying stopping criterion of Newton-Raphson algorithm. Defaulting to 1e-5.
#'
#' @param bound a positive number to avoid inflation of provider effect. Defaulting to 10.
#'
#' @param backtrack a Boolean indicating whether backtracking line search is implemented. Defaulting to FALSE.
#'
#' @param Rcpp a Boolean indicating whether the Rcpp function is used. Defaulting to TRUE.
#'
#' @param AUC a Boolean indicating whether report AUC. Defaulting to FALSE.
#'
#' @param message a Boolean indicating whether print out the information of the data preparation process and track the fitting process. Defaulting to TRUE.
#'
#' @param cutoff an integer as cutoff of provider size with 10 as default. Providers with observations fewer than the "cutoff" value will be labeled as "include = 0".
#'
#' @param stop a character string specifying the stopping rule to determine convergence.
#' `"incre"` means we stop the algorithm when the infinity norm of the difference between current and previous beta coefficients is less than the `tol`.
#' `"relch"` means we stop the algorithm when the \eqn{(loglik(m)-loglik(m-1))/(loglik(m))} is less than the `tol`,
#'  where \eqn{loglik(m)} denotes the log-partial likelihood at iteration step m.
#' `"ratch"` means we stop the algorithm when \eqn{(loglik(m)-loglik(m-1))/(loglik(m)-loglik(0))} is less than the `tol`.
#' `"all"` means we stop the algorithm when all the stopping rules (`"beta"`, `"relch"`, `"ratch"`) are met.
#' `"or` means we stop the algorithm if any one of the rules (`"beta"`, `"relch"`, `"ratch"`) is met.
#' The default value is `or`.
#' If `iter.max` is achieved, it overrides any stop rule for algorithm termination.
#'
#' @param check a Boolean indicating whether to the `data_check` function and to check the data quality. Defaulting to FALSE
#'
#' @param ...
#'
#'
#' @details
#'
#' The default algorithm is based on Serial blockwise inversion Newton (SerBIN) proposed by
#' [Wu et al. (2022)](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.9387),
#' but users can also choose to use the block ascent Newton (BAN) algorithm proposed by
#' [He et al. (2013)](https://link.springer.com/article/10.1007/s10985-013-9264-6) to fit the model.
#' Both methodologies build upon the Newton-Raphson method, yet SerBIN simultaneously updates both the provider effect and covariate coefficient.
#' This concurrent update necessitates the inversion of the complete information matrix at each iteration.
#' In contrast, BAN adopts a two-layer updating approach, where the covariate coefficient is sequentially fixed to update the provider effect,
#' followed by fixing the provider effect to update the covariate coefficient.
#'
#' We suggest using the default `"SerBIN"` option as it typically converges much faster for most datasets.
#' However, in rare cases where the SerBIN algorithm encounters second-order derivative irreversibility leading to an error,
#' users can consider using the `"BAN"` option as an alternative.
#'
#' For a deeper understanding, please consult the original article for detailed insights.
#'
#'
#'
#' @return An object with S3 class \code{logis_fe}.
#'
#' \item{beta}{a vector of fixed effects estimates of covariates}
#'
#' \item{gamma}{a vector of estimates of provider effects}
#'
#' \item{linear_pred}{a vector of linear predictors}
#'
#' \item{pred}{a vector of predicted probabilities}
#'
#' \item{neg2Loglkd}{minus two times log likelihood}
#'
#' \item{AIC}{Akaike info criterion}
#'
#' \item{BIC}{Bayesian info criterion}
#'
#' \item{AUC}{area under the ROC curve}
#'
#'
#' @examples
#' data(data_FE)
#' fit_fe <- logis_fe(data_FE$Y, data_FE$Z, data_FE$ID)
#'
#' @importFrom Rcpp evalCpp
#' @importFrom pROC auc
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @references
#' \itemize{
#' \item Wu, W, Yang, Y, Kang, J, He, K. (2022) Improving large-scale estimation and inference for profiling health care providers.
#' \emph{Statistics in Medicine}, \strong{41(15)}: 2840-2853.
#'
#' \item He K, Kalbfleisch, J, Li, Y, and et al. (2013) Evaluating hospital readmission rates in dialysis providers; adjusting for hospital effects.
#' \emph{Lifetime Data Analysis}, \strong{19}: 490-512.
#' }
#'
#' @keywords Cost-Efficient Newton-Raphson (CENR) Algorithm, Fixed Provider Effects
#'
#' @export
#'
#' @useDynLib pprof, .registration = TRUE

logis_fe <- function(Y, Z, ID, algorithm = "SerBIN", max.iter = 10000, tol = 1e-5, bound = 10,
                     backtrack = TRUE, Rcpp = TRUE, AUC = FALSE, message = TRUE, cutoff = 10,
                     stop = "or", check = FALSE){

  # Check input
  if (missing(Y) || missing(Z) || missing(ID))
    stop("Arguments 'Y', 'Z', and 'ID' are all required!", call.=FALSE)

  if (!is.logical(backtrack)) stop("Argument 'backtrack' NOT as required!", call.=F)

  #check dimensions of the input data
  if (length(Y) != length(ID) | length(ID) != nrow(Z)){
    stop("Dimensions of the input data do not match!!", call.=F)
  }

  if (check == TRUE)
    data_check(Y, Z, ID)

  # Data Preparation
  data <- as.data.frame(cbind(Y, ID, Z))
  Y.char <- colnames(data)[1]
  prov.char <- colnames(data)[2]
  Z.char <- colnames(Z)

  data <- data[order(factor(data[,prov.char])),] # sort data by provider ID
  prov.size <- as.integer(table(data[,prov.char])) # provider sizes
  prov.size.long <- rep(prov.size,prov.size) # provider sizes assigned to patients
  data$included <- 1 * (prov.size.long > cutoff) # create variable 'included' as an indicator
  if (message == TRUE) warning(sum(prov.size<=cutoff)," out of ",length(prov.size),
                               " providers considered small and filtered out!",immediate.=T,call.=F)

  prov.list <- unique(data[data$included==1,prov.char])   # a reduced list of provider IDs
  prov.no.events <-      # providers with no events
    prov.list[sapply(split(data[data$included==1,Y.char], factor(data[data$included==1,prov.char])),sum)==0]
  data$no.events <- 0
  data$no.events[data[,prov.char]%in%c(prov.no.events)] <- 1
  if (message == TRUE) message(paste(length(prov.no.events),"out of",length(prov.list),
                                     "remaining providers with no events."))
  prov.all.events <-     # providers with all events
    prov.list[sapply(split(1-data[data$included==1,Y.char],factor(data[data$included==1,prov.char])),sum)==0]
  data$all.events <- 0
  data$all.events[data[,prov.char]%in%c(prov.all.events)] <- 1
  if (message == TRUE) message(paste(length(prov.all.events),"out of",length(prov.list),
                                     "remaining providers with all events."))
  if (message == TRUE) message(paste0("After screening, ", round(sum(data[data$included==1,Y.char])/length(data[data$included==1,Y.char])*100,2),
                                      "% of all records exhibit occurrences of events (Y = 1)"))


  # for the remaining parts, only use the data with "included==1" ("cutoff" of provider size)
  data <- data[data$included==1,]
  n.prov <- sapply(split(data[, Y.char], data[, prov.char]), length) # provider-specific number of discharges
  n.events.prov <- sapply(split(data[, Y.char], data[, prov.char]), sum) # provider-specific number of events
  Z <- as.matrix(data[,Z.char])
  gamma.prov <- rep(log(mean(data[,Y.char])/(1-mean(data[,Y.char]))), length(n.prov))
  beta <- rep(0, NCOL(Z))

  Loglkd <- function(gamma.obs, beta) {
    sum((gamma.obs+Z%*%beta)*data[,Y.char]-log(1+exp(gamma.obs+Z%*%beta)))
  }
  loglkd_init = Loglkd(rep(gamma.prov, n.prov), beta)

  # Model
  if (algorithm == "SerBIN") {
    if (Rcpp) { #Rcpp always use "backtrack"
      ls <- logis_BIN_fe_prov(as.matrix(data[,Y.char]),Z,n.prov,gamma.prov,beta,
                              0,1,tol,max.iter, bound, message, backtrack, stop)
      gamma.prov <- as.numeric(ls$gamma)
      beta <- as.numeric(ls$beta)
    } else {
      iter <- 0
      crit <- 100 # initialize stop criterion
      if (message){
        message("Implementing SerBIN algorithm for fixed provider effects model ...")
      }

      if (backtrack){ # initialize parameters for backtracking line search
        s <- 0.01
        t <- 0.6
      }

      while (iter<=max.iter & crit>=tol) {
        iter <- iter + 1
        gamma.obs <- rep(gamma.prov, n.prov)
        loglkd = Loglkd(gamma.obs, beta)
        p <- c(plogis(gamma.obs+Z%*%beta))
        pq <- p*(1-p)
        pq[pq == 0] <- 1e-20
        score.gamma <- sapply(split(data[,Y.char]-p, data[,prov.char]), sum)
        score.beta <- t(Z)%*%(data[,Y.char]-p)
        info.gamma.inv <- 1/sapply(split(pq, data[,prov.char]),sum) #I_11^(-1)
        info.betagamma <- sapply(by(pq*Z,data[,prov.char],identity),colSums) #I_21
        info.beta <- t(Z)%*%(pq*Z) #I_22
        mat.tmp1 <- info.gamma.inv*t(info.betagamma) #J_1^T
        schur.inv <- solve(info.beta-info.betagamma%*%mat.tmp1) #S^-1
        mat.tmp2 <- mat.tmp1%*%schur.inv #J_2^T

        d.gamma.prov <- info.gamma.inv*score.gamma +
          mat.tmp2%*%(t(mat.tmp1)%*%score.gamma-score.beta)
        d.beta <- -t(mat.tmp2)%*%score.gamma+schur.inv%*%score.beta
        v <- 1 # initialize step size
        if (backtrack) {
          d.loglkd <- Loglkd(rep(gamma.prov+v*d.gamma.prov, n.prov), beta+v*d.beta) - loglkd
          lambda <- c(score.gamma,score.beta)%*%c(d.gamma.prov,d.beta)
          while (d.loglkd < s*v*lambda) {  #update step size
            v <- t * v
            d.loglkd <- Loglkd(rep(gamma.prov+v*d.gamma.prov, n.prov), beta+v*d.beta) - loglkd
          }
        }
        gamma.prov <- gamma.prov + v * d.gamma.prov
        gamma.prov <- pmin(pmax(gamma.prov, median(gamma.prov)-bound), median(gamma.prov)+bound)
        beta.new <- beta + v * d.beta

        d.loglkd = Loglkd(rep(gamma.prov, n.prov), beta.new) - loglkd

        # stopping criterion
        if (stop=="beta"){
          crit <- norm(matrix(beta-beta.new),"I")
          if (message){
            cat(paste0("Iter ",iter,": Inf norm of running diff in est reg parm is ",
                       formatC(crit,digits=3,format="e"),";\n"))
          }
        }
        else if (stop=="relch"){
          crit <- abs(d.loglkd/(d.loglkd+loglkd))
          if (message) {
            cat(paste0("Iter ",iter,": Relative change in est log likelihood is ",
                       formatC(crit,digits=3,format="e"),";\n"))
          }
        }
        else if (stop=="ratch"){
          crit <- abs(d.loglkd/(d.loglkd+loglkd-loglkd_init))
          if (message) {
            cat(paste0("Iter ",iter,": Adjusted relative change in est log likelihood is ",
                       formatC(crit,digits=3,format="e"),";\n"))
          }
        }
        else if (stop=="all"){
          crit_beta <- norm(matrix(beta-beta.new),"I")
          crit_relch <- abs(d.loglkd/(d.loglkd+loglkd))
          crit_ratch <- abs(d.loglkd/(d.loglkd+loglkd-loglkd_init))
          crit <- max(crit_beta, crit_relch, crit_ratch)
          if (message) {
            cat(sprintf("Iter %d: Maximum criterion across all checks is %.3e;\n", iter, crit))
          }
        }
        else if (stop=="or"){
          crit_beta <- norm(matrix(beta-beta.new),"I")
          crit_relch <- abs(d.loglkd/(d.loglkd+loglkd))
          crit_ratch <- abs(d.loglkd/(d.loglkd+loglkd-loglkd_init))
          crit <- min(crit_beta, crit_relch, crit_ratch)
          if (message) {
            cat(sprintf("Iter %d: Minimum criterion across all checks is %.3e;\n", iter, crit))
          }
        }

        beta <- beta.new

      }
      if (message){
        message("\n SerBIN algorithm converged after ",iter," iterations!")
      }
    }
  } else if (algorithm == "BAN"){
    if (Rcpp) {
      ls <- logis_fe_prov(as.matrix(data[,Y.char]),Z,n.prov,gamma.prov,beta,backtrack,max.iter,bound,tol,message,stop)
      gamma.prov <- as.numeric(ls$gamma); beta <- as.numeric(ls$beta)
    } else {
      iter <- 0
      crit <- 100 # initialize stop criterion
      if (message){
        message("Implementing BAN algorithm for fixed provider effects model ...")
      }
      if (backtrack){ # initialize parameters for backtracking line search
        s <- 0.01
        t <- 0.8
      }

      while (iter<=max.iter & crit>=tol) {
        iter <- iter + 1
        # provider effect update
        gamma.obs <- rep(gamma.prov, n.prov)
        loglkd.old = Loglkd(gamma.obs, beta)
        Z.beta <- Z%*%beta
        p <- c(plogis(gamma.obs+Z.beta)); pq <- p*(1-p)
        pq[pq == 0] <- 1e-20
        score.gamma.prov <- sapply(split(data[,Y.char]-p, data[,prov.char]), sum)
        d.gamma.prov <- score.gamma.prov / sapply(split(pq, data[,prov.char]), sum)
        v <- 1 # initialize step size
        if (backtrack) {
          loglkd <- Loglkd(rep(gamma.prov, n.prov), beta)
          d.loglkd <- Loglkd(rep(gamma.prov+v*d.gamma.prov, n.prov), beta) - loglkd
          lambda <- score.gamma.prov%*%d.gamma.prov
          while (d.loglkd < s*v*lambda) {
            v <- t * v
            d.loglkd <- Loglkd(rep(gamma.prov+v*d.gamma.prov, n.prov), beta) - loglkd
          }
        }
        gamma.prov <- gamma.prov + v * d.gamma.prov
        gamma.prov <- pmin(pmax(gamma.prov, median(gamma.prov)-bound), median(gamma.prov)+bound)
        gamma.obs <- rep(gamma.prov, n.prov)

        # regression parameter update
        p <- c(plogis(gamma.obs+Z.beta)); pq <- p*(1-p)
        score.beta <- t(Z)%*%(data[,Y.char]-p)
        info.beta <- t(Z)%*%(c(pq)*Z)
        d.beta <- as.numeric(solve(info.beta)%*%score.beta)
        v <- 1 # initialize step size
        if (backtrack) {
          loglkd <- Loglkd(gamma.obs, beta)
          d.loglkd <- Loglkd(gamma.obs, beta+v*d.beta) - loglkd
          lambda <- c(score.beta)%*%d.beta
          while (d.loglkd < s*v*lambda) {
            v <- t * v
            d.loglkd <- Loglkd(gamma.obs, beta+v*d.beta) - loglkd
          }
        }
        beta.new <- beta + v * d.beta
        d.loglkd = Loglkd(rep(gamma.prov, n.prov), beta.new) - loglkd.old

        # stopping criterion
        if (stop=="beta"){
          crit <- norm(matrix(beta-beta.new),"I")
          if (message){
            cat(paste0("Iter ",iter,": Inf norm of running diff in est reg parm is ",
                       formatC(crit,digits=3,format="e"),";\n"))
          }
        }
        else if (stop=="relch"){
          crit <- abs(d.loglkd/(d.loglkd+loglkd.old))
          if (message) {
            cat(paste0("Iter ",iter,": Relative change in est log likelihood is ",
                       formatC(crit,digits=3,format="e"),";\n"))
          }
        }
        else if (stop=="ratch"){
          crit <- abs(d.loglkd/(d.loglkd+loglkd.old-loglkd_init))
          if (message) {
            cat(paste0("Iter ",iter,": Adjusted relative change in est log likelihood is ",
                       formatC(crit,digits=3,format="e"),";\n"))
          }
        }
        else if (stop=="all"){
          crit_beta <- norm(matrix(beta-beta.new),"I")
          crit_relch <- abs(d.loglkd/(d.loglkd+loglkd.old))
          crit_ratch <- abs(d.loglkd/(d.loglkd+loglkd.old-loglkd_init))
          crit <- max(crit_beta, crit_relch, crit_ratch)
          if (message) {
            cat(sprintf("Iter %d: Maximum criterion across all checks is %.3e;\n", iter, crit))
          }
        }
        else if (stop=="or"){
          crit_beta <- norm(matrix(beta-beta.new),"I")
          crit_relch <- abs(d.loglkd/(d.loglkd+loglkd.old))
          crit_ratch <- abs(d.loglkd/(d.loglkd+loglkd.old-loglkd_init))
          crit <- min(crit_beta, crit_relch, crit_ratch)
          if (message) {
            cat(sprintf("Iter %d: Minimum criterion across all checks is %.3e;\n", iter, crit))
          }
        }

        beta <- beta.new
      }
      if (message){
        message("\n BAN algorithm converged after ",iter," iterations!")
      }
    }
  } else {
    stop("Argument 'algorithm' NOT as required!")
  }

  gamma.obs <- rep(gamma.prov, n.prov)
  neg2Loglkd <- -2*sum((gamma.obs+Z%*%beta)*data[,Y.char]-log(1+exp(gamma.obs+Z%*%beta)))
  AIC <- neg2Loglkd + 2 * (length(gamma.prov)+length(beta))
  BIC <- neg2Loglkd + log(nrow(data)) * (length(gamma.prov)+length(beta))

  df.prov <- data.frame(Obs_provider = sapply(split(data[,Y.char],data[,prov.char]),sum),
                        gamma_est = gamma.prov) #original gamma-hat, for internal using
  linear_pred <- Z %*% beta
  pred <- as.numeric(plogis(gamma.obs + linear_pred))

  #change output format
  beta <- matrix(beta)
  gamma.prov <- matrix(gamma.prov)
  dimnames(beta) <- list(Z.char, "beta")
  dimnames(gamma.prov) <- list(names(n.prov), "gamma")

  char_list <- list(Y.char = Y.char,
                    prov.char = prov.char,
                    Z.char = Z.char)

  return_ls <- structure(list(beta = beta,
                              gamma = gamma.prov, #provider effect
                              linear_pred = linear_pred, #linear predictor
                              pred = pred, #predicted probability
                              neg2Loglkd = neg2Loglkd,
                              AIC = AIC,
                              BIC = BIC,
                              obs = data[, Y.char], #patient-level obs
                              prov = data[, prov.char]),
                         class = "logis_fe")
  if (AUC) {
    AUC <- pROC::auc(data[,Y.char], pred)
    return_ls$AUC <- AUC[1]
  }
  return_ls$df.prov <- df.prov
  return_ls$char_list <- char_list
  return_ls$data_include <- data
  return(return_ls)
}
