#' Main Function for Fitting the Fixed Effect Linear Model
#'
#' Fit a fixed effect linear model via profile likelihood or dummy encoding.
#'
#' @param formula a two-sided formula object describing the model to be fitted,
#' with the response variable on the left of a ~ operator and covariates on the right,
#' separated by + operators. The fixed effect of the grouping identifier is specified using \code{id()}.
#' @param data a data frame containing the variables named in the `formula`,
#' or the columns specified by `Y.char`, `Z.char`, and `ID.char`.
#' @param Y.char a character string specifying the column name of the response variable in the `data`.
#' @param Z.char a character vector specifying the column names of the covariates in the `data`.
#' @param ID.char a character string specifying the column name of the grouping identifier in the `data`.
#' @param Y a numeric vector representing the response variable.
#' @param Z a matrix or data frame representing the covariates, which can include both numeric and categorical variables.
#' @param ID a numeric vector representing the grouping identifier.
#' @param method a character string specifying the method to fit the model.
#' `"pl"` (default) uses profile likelihood to fit the model,
#' while `"dummy"` calls \code{\link{lm}} to fit the model using dummy variables for the grouping identifier.
#'
#' @return A list of objects with S3 class \code{"linear_fe"}:
#' \item{coefficient}{a list containing the estimated coefficients:
#'   \code{beta}, the fixed effects for each predictor, and \code{gamma}, the effect for each group.}
#' \item{variance}{a list containing the variance estimates:
#'   \code{beta}, the variance-covariance matrix of the predictor coefficients, and \code{gamma}, the variance of the group effects.}
#' \item{sigma}{the residual standard error.}
#' \item{fitted}{the fitted values of each individual.}
#' \item{observation}{the original response of each individual.}
#' \item{residuals}{the residuals of each individual, that is response minus fitted values}
#' \item{linear_pred}{the linear predictor of each individual.}
#' \item{data_include}{the data used to fit the model, sorted by the group identifier.
#' For categorical covariates, this includes the dummy variables created for
#' all categories except the reference level.}
#' \item{char_list}{a list of the character vectors representing the column names for
#' the response variable, covariates, and group identifier.
#' For categorical variables, the names reflect the dummy variables created for each category.}
#' \item{method}{the method used for model fitting, either \code{"Profile Likelihood"} or \code{"Dummy"}.}
#'
#' @details
#' This function is used to fit a fixed effect linear model of the form:
#' \deqn{Y_{ij} = \gamma_i + \mathbf{Z}_{ij}^\top\boldsymbol\beta + \epsilon_{ij}}
#' where \eqn{Y_{ij}} is the outcome for individual \eqn{j} in group \eqn{i}, \eqn{\gamma_i} is the group-specific effect, \eqn{\mathbf{Z}_{ij}} are the covariates, and \eqn{\boldsymbol\beta} is the vector of coefficients for the covariates.
#' The default method for fitting the model is profile likelihood, but dummy encoding can also be used by specifying the appropriate method.
#' When the number of groups is very large, we recommend using the profile likelihood method, as it is significantly faster than dummy encoding.
#'
#' The function accepts three different input formats:
#' a formula and dataset, where the formula is of the form \code{response ~ covariates + id(group)}, with \code{group} representing the group identifier;
#' a dataset along with the column names of the response, covariates, and group identifier;
#' or the outcome vector \eqn{\boldsymbol{Y}}, the covariate matrix or data frame \eqn{\mathbf{Z}}, and the group identifier vector \eqn{\boldsymbol{\gamma}}.
#'
#'
#' @importFrom Matrix bdiag
#'
#' @export
#'
#' @examples
#' data(ExampleDataLinear)
#' Y <- ExampleDataLinear$Y
#' Z <- ExampleDataLinear$Z
#' ID <- ExampleDataLinear$ID
#' data <- data.frame(Y, ID, Z)
#' Z.char <- colnames(Z)
#' Y.char <- "Y"
#' ID.char <- "ID"
#' formula <- as.formula(paste("Y ~", paste(Z.char, collapse = " + "), "+ id(ID)"))
#'
#' # Fit fixed linear effect model using three input formats
#' fit_fe1 <- linear_fe(Y = Y, Z = Z, ID = ID)
#' fit_fe2 <- linear_fe(data = data, Y.char = Y.char, Z.char = Z.char, ID.char = ID.char)
#' fit_fe3 <- linear_fe(formula, data)
#'
#' @references
#' R Core Team (2023). \emph{The R Stats Package: lm}.
#' Available at: \url{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/lm.html}
#' \cr

linear_fe <- function(formula = NULL, data = NULL,
                      Y = NULL, Z = NULL, ID = NULL,
                      Y.char = NULL, Z.char = NULL, ID.char = NULL,
                      method = "pl"){
  if (method == "pl") {
    if (!is.null(formula) && !is.null(data)) {
      message("Input format: formula and data.")

      formula_terms <- terms(formula)
      response <- as.character(attr(formula_terms, "variables"))[2]
      predictors <- attr(formula_terms, "term.labels")
      # id_var <- NULL
      # for (term in predictors) {
      #   if (grepl("as.factor", term)) {
      #     id_var <- gsub("as.factor\\((.*)\\)", "\\1", term)
      #     predictors <- predictors[predictors != term]
      #     break
      #   }
      # }
      # if (is.null(id_var)) {
      #   for (term in predictors) {
      #     if (is.factor(data[[term]])) {
      #       id_var <- term
      #       predictors <- predictors[predictors != id_var]
      #       break
      #     }
      #   }
      # }

      ID.char <- gsub(".*id\\(([^)]+)\\).*", "\\1", predictors[grepl("id\\(", predictors)])
      Z.char <- predictors[!grepl("id\\(", predictors)]

      if (!all(c(Y.char, Z.char, ID.char) %in% colnames(data)))
        stop("The formula contains variables that are not in the data.", call.=F)

      #mf <- model.frame(formula, data)
      #Y <- model.response(mf)
      Y <- data[,response, drop = F]
      Z <- model.matrix(reformulate(Z.char), data)[, -1, drop = F]
      # Z <- model.matrix(~ data[[predictors]] - 1)
      ID <- data[,ID.char, drop = F]
    }
    else if (!is.null(data) && !is.null(Y.char) && !is.null(Z.char) && !is.null(ID.char)) {
      message("Input format: data, Y.char, Z.char, and ID.char.")

      if (!all(c(Y.char, Z.char, ID.char) %in% colnames(data)))
        stop("Some of the specified columns are not in the data!", call.=FALSE)

      Y <- data[, Y.char]
      Z <- model.matrix(reformulate(Z.char), data)[, -1, drop = FALSE]
      ID <- data[, ID.char, drop = F]
    }
    else if (!is.null(Y) && !is.null(Z) && !is.null(ID)) {
      message("Input format: Y, Z, and ID.")

      if (length(Y) != length(ID) | (length(ID) != nrow(Z))) {
        stop("Dimensions of the input data do not match!!", call.=F)
      }

      Z.char <- colnames(Z)
      Z <- model.matrix(reformulate(Z.char), Z)[, -1, drop = FALSE]
    }
    else {
      stop("Insufficient or incompatible arguments provided. Please provide either (1) formula and data, (2) data, Y.char, Z.char, and ID.char, or (3) Y, Z, and ID.", call.=FALSE)
    }

    data <- data.frame(Y, ID, Z)
    Y.char <- colnames(data)[1]
    ID.char <- colnames(data)[2]
    Z.char <- colnames(Z)
    data <- data[order(factor(data[,ID.char])),]

    n.prov <- sapply(split(data[, Y.char], data[, ID.char]), length)
    m <- length(n.prov) # number of providers
    n <- sum(n.prov) # number of observations
    p <- length(Z.char) # number of covariates

    Z <- as.matrix(data[,Z.char], drop = F)
    Y <- as.matrix(data[, Y.char, drop = F])
    ID <- as.matrix(data[, ID.char, drop = F])

    Q <- lapply(n.prov, function(n) diag(n)-matrix(1, nrow = n, ncol = n)/n)

    # Coefficients
    beta <- matrix(solve(t(Z)%*%bdiag(Q)%*%Z)%*%t(Z)%*%bdiag(Q)%*%Y, ncol = 1)
    colnames(beta) <- "beta"
    rownames(beta) <- Z.char

    y_bar <- sapply(split(data[,Y.char], data[,ID.char]),mean)
    Z_bar <- t(matrix(sapply(split(data[,Z.char, drop = FALSE], data[,ID.char]), colMeans),
                      ncol=length(y_bar), nrow = length(beta)))
    gamma.prov <- as.matrix(y_bar - Z_bar %*% beta)
    colnames(gamma.prov) <- "gamma"
    rownames(gamma.prov)<- names(n.prov)

    coefficient <- list()
    coefficient$beta <- beta
    coefficient$gamma <- gamma.prov

    # Prediction
    linear_pred <- Z %*% beta
    colnames(linear_pred) <- "Linear Predictor"
    rownames(linear_pred) <- seq_len(nrow(linear_pred))
    gamma.obs <- rep(gamma.prov, n.prov)
    pred <- gamma.obs + linear_pred
    colnames(pred) <- "Prediction"
    rownames(pred) <- seq_len(nrow(pred))
    residuals <- matrix(Y - pred, ncol = 1)
    colnames(residuals) <- "Residuals"
    rownames(residuals) <- seq_len(nrow(residuals))
    sigma_hat_sq <- sum(residuals^2)/(n - m - p)

    # Variance
    varcov_beta <- matrix(sigma_hat_sq * solve(t(Z)%*%bdiag(Q)%*%Z), ncol = p, nrow = p)
    rownames(varcov_beta) <- Z.char
    colnames(varcov_beta) <- Z.char

    var_gamma <- matrix(sigma_hat_sq/n.prov, ncol = 1)
    rownames(var_gamma) <- names(n.prov)
    colnames(var_gamma) <- "Variance.Gamma"

    variance <- list()
    variance$beta <- varcov_beta
    variance$gamma <- var_gamma

    char_list <- list(Y.char = Y.char,
                      ID.char = ID.char,
                      Z.char = Z.char)

    model_method <- "Profile Likelihood"
  }
  else if (method == "dummy") {
    if (!is.null(formula) && !is.null(data)){
      message("Input format: formula and data.")

      formula_terms <- terms(formula)
      Y.char <- as.character(attr(formula_terms, "variables"))[2]
      predictors <- attr(formula_terms, "term.labels")
      ID.char <- gsub(".*id\\(([^)]+)\\).*", "\\1", predictors[grepl("id\\(", predictors)])
      Z.char <- predictors[!grepl("id\\(", predictors)]
      original_ID <- data[, ID.char, drop = F]
      data[,ID.char] <- as.factor(data[,ID.char])

      if (!all(c(Y.char, Z.char, ID.char) %in% colnames(data)))
        stop("The formula contains variables that are not in the data.", call.=F)

      # ID.char is always in the first position
      new_formula <- as.formula(paste(Y.char, "~", ID.char, "+",
                                      paste(Z.char, collapse = " + ")))

      data <- data[order(factor(data[,ID.char])),]
      formula <- update(new_formula, . ~ . - 1)
      fit_lm <- lm(formula, data = data)
    }
    else if (!is.null(data) && !is.null(Y.char) && !is.null(Z.char) && !is.null(ID.char)) {
      message("Input format: data, Y.char, Z.char, and ID.char.")

      if (!all(c(Y.char, Z.char, ID.char) %in% colnames(data)))
        stop("Some of the specified columns are not in the data!", call.=FALSE)

      original_ID <- data[, ID.char, drop = F]
      data[,ID.char] <- as.factor(data[,ID.char])
      data <- data[order(factor(data[,ID.char])),]
      formula <- as.formula(paste(Y.char, "~", ID.char, "+", paste(Z.char, collapse = " + "), "-1"))
      fit_lm <- lm(formula, data = data)
    }
    else if (!is.null(Y) && !is.null(Z) && !is.null(ID)) {
      message("Input format: Y, Z, and ID.")

      if (length(Y) != length(ID) | (length(ID) != nrow(Z))) {
        stop("Dimensions of the input data do not match!!", call.=F)
      }
      data <- as.data.frame(cbind(Y, ID, Z))
      Y.char <- colnames(data)[1]
      ID.char <- colnames(data)[2]
      Z.char <- colnames(Z)
      original_ID <- data[, ID.char, drop = F]
      data[,ID.char] <- as.factor(data[,ID.char])
      data <- data[order(factor(data[,ID.char])),]
      formula <- as.formula(paste(Y.char, "~", ID.char, "+", paste(Z.char, collapse = " + "), "-1"))
      fit_lm <- lm(formula, data = data)
    }
    else {
      stop("Insufficient or incompatible arguments provided. Please provide either (1) formula and data, (2) data, Y.char, Z.char, and ID.char, or (3) Y, Z, and ID.", call.=FALSE)
    }

    sum <- summary(fit_lm)
    X.model <- model.matrix(fit_lm)
    Y <- as.matrix(data[, Y.char, drop = F])
    # ID <- as.matrix(data[, ID.char, drop = F])

    n.prov <- sapply(split(data[, Y.char], data[, ID.char]), length)
    m <- length(n.prov) # number of providers
    n <- sum(n.prov) # number of observations
    Z.char <- colnames(X.model)[(m+1):length(colnames(X.model))]
    Z <- X.model[,Z.char, drop = F]
    p <- length(Z.char) # number of covariates
    data <- data.frame(Y, original_ID, Z)

    # Coefficients
    beta <- matrix(sum$coefficients[(m+1):(m+p), 1], ncol = 1)
    colnames(beta) <- "beta"
    rownames(beta) <- Z.char
    gamma.prov <- matrix(sum$coefficients[1:m, 1], ncol = 1)
    colnames(gamma.prov) <- "gamma"
    rownames(gamma.prov)<- names(n.prov)

    coefficient <- list()
    coefficient$beta <- beta
    coefficient$gamma <- gamma.prov

    # Variance
    sigma_hat_sq <- sum$sigma^2
    varcov <- sigma_hat_sq * solve(t(X.model)%*%X.model)
    varcov_beta <- matrix(varcov[(m+1):(m+p), (m+1):(m+p)], ncol = p, nrow = p)
    rownames(varcov_beta) <- Z.char
    colnames(varcov_beta) <- Z.char

    var_gamma <- matrix(diag(varcov)[1:m], ncol = 1)
    rownames(var_gamma) <- names(n.prov)
    colnames(var_gamma) <- "Variance.Gamma"

    variance <- list()
    variance$beta <- varcov_beta
    variance$gamma <- var_gamma

    # Prediction
    residuals <- matrix(sum$residuals, ncol = 1)
    colnames(residuals) <- "Residuals"
    rownames(residuals) <- seq_len(nrow(residuals))

    pred <- Y - sum$residuals
    pred <- matrix(pred, ncol = 1)
    colnames(pred) <- "Prediction"
    rownames(pred) <- seq_len(nrow(pred))
    linear_pred <- Z %*% beta
    colnames(linear_pred) <- "Linear Predictor"
    rownames(linear_pred) <- seq_len(nrow(linear_pred))

    char_list <- list(Y.char = Y.char,
                      ID.char = ID.char,
                      Z.char = Z.char)

    model_method <- "Dummy"
  }
  else stop("Method should be either 'pl' or 'dummy'.")


  result <- structure(list(coefficient = coefficient,
                           variance = variance,
                           sigma = sqrt(sigma_hat_sq),
                           fitted = pred,
                           observation = Y,
                           residuals = residuals,
                           linear_pred = linear_pred
                           ),
                      class = "linear_fe")  #define a list for prediction

  result$data_include <- data
  result$char_list <- char_list
  result$method <- model_method

  return(result)
}














