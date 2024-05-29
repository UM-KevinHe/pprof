#=== FUNCTIONS FOR Data Quality Check =================================
#' Data Quality Check including
#'
#' @param Y a numerical vector, with values of 0 or 1, indicating the outcome variable.
#'
#' @param Z a matrix or data frame containing covariates.
#'
#' @param ID a vector representing the provider id. Its elements can be either numeric values or characters.
#'
#' @param message a Boolean indicating whether printing out the information of the data preparation process. Defaulting to "TRUE".
#'
#' @param ...
#'
#'
#' @details The function performs the following checks:
#'   \itemize{
#'     \item \strong{Missingness}: Checks for any missing values in the dataset and provides a summary of missing data.
#'     \item \strong{Variation}: Identifies covariates with zero or near-zero variance which might affect model stability.
#'     \item \strong{Correlation}: Analyzes pairwise correlation among covariates and highlights highly correlated pairs.
#'     \item \strong{VIF}: Computes the Variable Inflation Factors to identify covariates with potential multicollinearity issues.
#'   }
#'
#' The `fe_data_prep()` function returns data sorted by the provider identifiers,
#' accompanied by additional provider-related information indicating whether the provider's size exceeds a specified "cutoff"
#' and whether the respective provider has experienced either zero or all events.
#' The reason behind introducing a "cutoff" lies in findings from both simulated and real data studies,
#' revealing the instability of coefficient estimates for providers with small sizes.
#' Consequently, we recommend excluding small providers during the model fitting process.
#' It is important to note that the resultant data frame retains all providers, including small ones,
#' but utilizes the "included = 0" label to signify the small providers. Subsequently, during the model fitting stage,
#' the `logis_fe()` function disregards records marked with "included = 0".
#'
#'
#' @return
#'
#' \item{data}{a sorted data frame including response, provider identifiers, covariates, and additional provider information.}
#'
#' \item{char_list}{a list including variable names.}
#'
#'
#' @examples
#' data(data_FE)
#' results <- data_check(data_FE$Y, data_FE$Z, data_FE$ID)
#' head(data.prep$data)
#' data.prep$char_list
#'
#' @importFrom caret nearZeroVar
#' @importFrom olsrr ols_vif_tol
#'
#' @keywords data quality, data check
#'
#' @export

data_check <- function(Y, Z, ID, message = TRUE) {
  data <- as.data.frame(cbind(Y, ID, Z))
  check_results <- list()

  ## check missingness of variables
  if (message == TRUE) message("Checking missingness of variables ... ")
  if (sum(complete.cases(data[,c(Y.char,Z.char,prov.char)]))==NROW(data)) {
    if (message == TRUE) message("Missing values NOT found. Checking missingness of variables completed!")
  } else {
    check.na <- function(name) {
      if (sum(is.na(data[,name])) > 0) {
        warning(sum(is.na(data[,name]))," out of ",NROW(data[,name])," in '",name,"' missing!",immediate.=T,call.=F)
      }
    }
    invisible(sapply(c(Y.char,Z.char,prov.char), check.na))
    missingness <- (1 - sum(complete.cases(data[,c(Y.char,Z.char,prov.char)])) / NROW(data)) * 100
    stop(paste(round(missingness,2), "% of all observations are missing!",sep=""),call.=F)
  }

  ## check variation in covariates
  if (message == TRUE) message("Checking variation in covariates ... ")
  nzv <- caret::nearZeroVar(data[,Z.char], saveMetrics=T)
  if (sum(nzv$zeroVar==T) > 0) {
    stop("Covariate(s) '", paste(row.names(nzv[nzv$zeroVar==T,]), collapse="', '"),
         "' with zero variance(s)!", call.=F)
  } else if (sum(nzv$nzv==T) > 0) {
    warning("Covariate(s) '",paste(row.names(nzv[nzv$nzv==T,]), collapse="', '"),
            "' with near zero variance(s)!",immediate.=T,call.=F)
  }
  if (message == TRUE) message("Checking variation in covariates completed!")

  ## check correlation
  if (message == TRUE) message("Checking pairwise correlation among covariates ... ")
  cor <- cor(data[,Z.char])
  threshold.cor <- 0.9
  if (sum(abs(cor[upper.tri(cor)])>threshold.cor) > 0) {
    cor[lower.tri(cor,diag=T)] <- 0
    ind <- which(abs(cor)>threshold.cor)
    pairs <- sapply(ind, function(ind) c(rownames(cor)[ind%%NROW(cor)],
                                         colnames(cor)[ind%/%NROW(cor)+1]))
    warning("The following ", NCOL(pairs),
            " pair(s) of covariates are highly correlated (correlation > ",
            threshold.cor,"): ", immediate.=T, call.=F)
    invisible(apply(pairs,2,function(col) message('("',paste(col, collapse='", "'),'")')))
  }
  if (message == TRUE) message("Checking pairwise correlation among covariates completed!")

  ## check VIF
  if (message == TRUE) message("Checking VIF of covariates ... ")
  m.lm <- lm(as.formula(paste(Y.char,"~",paste(Z.char, collapse="+"))), data=data)
  vif <- olsrr::ols_vif_tol(m.lm)
  if(sum(vif$VIF >= 10) > 0){
    warning("Covariate(s) '",
            paste(as.data.frame(vif)[vif$VIF>=10,"Variables"], collapse="', '"),
            "' with serious multicollinearity!",immediate.=T,call.=F)
  }
  if (message == TRUE) message("Checking VIF of covariates completed!")

}









