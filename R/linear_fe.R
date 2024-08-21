#' Main Function for Fitting the Fixed Effect Linear Model using profile likelihood
#'
#' @importFrom Matrix bdiag
#' @importFrom Rcpp evalCpp
#' @importFrom pROC auc
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @export
#'
#' @useDynLib pprof, .registration = TRUE

linear_fe.complete <- function(formula = NULL, data = NULL,
                               Y = NULL, Z = NULL, ID = NULL,
                               Y.char = NULL, Z.char = NULL, ID.char = NULL,
                               method = "pl"){
  if (method == "pl") {
    if (!is.null(formula) && !is.null(data) && !is.null(ID.char)) {
      formula_terms <- terms(formula)
      response <- as.character(attr(formula_terms, "variables"))[2]
      predictors <- attr(formula_terms, "term.labels")
      id_var <- NULL
      for (term in predictors) {
        if (grepl("as.factor", term)) {
          id_var <- gsub("as.factor\\((.*)\\)", "\\1", term)
          predictors <- predictors[predictors != term]
          break
        }
      }
      if (is.null(id_var)) {
        for (term in predictors) {
          if (is.factor(data[[term]])) {
            id_var <- term
            predictors <- predictors[predictors != id_var]
            break
          }
        }
      }

      mf <- model.frame(formula, data)
      Y <- model.response(mf)
      Z <- model.matrix(reformulate(predictors), data)[, -1, drop = FALSE]
      ID <- data[,id_var]
    }
    else if (!is.null(data) && !is.null(Y.char) && !is.null(Z.char) && !is.null(ID.char)) {
      if (!all(c(Y.char, Z.char, ID.char) %in% colnames(data)))
        stop("Some of the specified columns are not in the data!", call.=FALSE)

      Y <- data[, Y.char]
      Z <- as.matrix(data[, Z.char, drop = FALSE])
      ID <- data[, ID.char]
    }
    else if (!is.null(Y) && !is.null(Z) && !is.null(ID)) {
      if (length(Y) != length(ID) | (length(ID) != nrow(Z))) {
        stop("Dimensions of the input data do not match!!", call.=F)
      }
    }
    else {
      stop("Insufficient or incompatible arguments provided. Please provide either (1) formula and data, (2) data, Y.char, Z.char, and ID.char, or (3) Y, Z, and ID.", call.=FALSE)
    }

    data <- as.data.frame(cbind(Y, ID, Z))
    Y.char <- colnames(data)[1]
    prov.char <- colnames(data)[2]
    Z.char <- colnames(Z)
    data <- data[order(factor(data[,prov.char])),]

    n.prov <- sapply(split(data[, Y.char], data[, prov.char]), length)
    m <- length(n.prov) # number of providers
    n <- sum(n.prov) # number of observations
    p <- length(Z.char) # number of covariates

    Z <- as.matrix(data[,Z.char], drop = F)
    Y <- as.matrix(data[, Y.char, drop = F])
    ID <- as.matrix(data[, prov.char, drop = F])

    Q <- lapply(n.prov, function(n) diag(n)-matrix(1, nrow = n, ncol = n)/n)


    # Coefficients
    beta <- matrix(solve(t(Z)%*%bdiag(Q)%*%Z)%*%t(Z)%*%bdiag(Q)%*%Y, ncol = 1)
    colnames(beta) <- "beta"
    rownames(beta) <- Z.char

    y_bar <- sapply(split(data[,Y.char], data[,prov.char]),mean)
    Z_bar <- t(matrix(sapply(split(data[,Z.char, drop = FALSE], data[,prov.char]), colMeans),
                      ncol=length(y_bar), nrow = length(beta)))
    gamma.prov <- as.matrix(y_bar - Z_bar %*% beta)
    colnames(gamma.prov) <- "gamma"
    rownames(gamma.prov)<- names(n.prov)

    coefficient <- list()
    coefficient$beta <- beta
    coefficient$gamma <- gamma.prov

    # Prediction
    linear_pred <- Z %*% beta
    gamma.obs <- rep(gamma.prov, n.prov)
    pred = gamma.obs + linear_pred
    residuals <- Y - pred
    sigma_hat_sq <- sum(residuals^2)/(n - m - p)

    # Variance
    varcor_beta <- matrix(sigma_hat_sq * solve(t(Z)%*%bdiag(Q)%*%Z), ncol = p, nrow = p)
    rownames(varcor_beta) <- Z.char
    colnames(varcor_beta) <- Z.char

    var_gamma <- matrix(sigma_hat_sq/n.prov, ncol = 1)
    rownames(var_gamma) <- names(n.prov)
    colnames(var_gamma) <- "Variance.Gamma"

    variance <- list()
    variance$beta <- varcor_beta
    variance$gamma <- var_gamma

    char_list <- list(Y.char = Y.char,
                      prov.char = prov.char,
                      Z.char = Z.char)

  }
  else if (method == "lm") {
    if (!is.null(formula) && !is.null(data)){
      formula <- update(formula, . ~ . - 1)
      fit_lm <- lm(formula, data = data)
    }
    else if (!is.null(data) && !is.null(Y.char) && !is.null(Z.char) && !is.null(ID.char)) {
      formula <- as.formula(paste(Y.char, "~", paste(Z.char, collapse = " + "), "+", ID.char, "-1"))
      fit_lm <- lm(formula, data = data)
    }
    else if (!is.null(Y) && !is.null(Z) && !is.null(ID)) {
      if (length(Y) != length(ID) | (length(ID) != nrow(Z))) {
        stop("Dimensions of the input data do not match!!", call.=F)
      }
      Y.char <- deparse(substitute(Y))
      ID.char <- deparse(substitute(prov.ID))
      formula <- as.formula(paste(Y.char, "~", paste(Z.char, collapse = " + "), "+", ID.char, "-1"))
      fit_lm <- lm(formula, data = data)

    }
    else {
      stop("Insufficient or incompatible arguments provided. Please provide either (1) formula and data, (2) data, Y.char, Z.char, and ID.char, or (3) Y, Z, and ID.", call.=FALSE)
    }

  }
  else stop("Method should be either 'pl' or 'lm'")


  result <- structure(list(coefficient = coefficient,
                           variance = variance,
                           sigma = sqrt(sigma_hat_sq),
                           linear_pred = linear_pred,
                           prediction = pred,
                           observation = Y,
                           prov = data[, prov.char]),
                      class = "linear_fe")  #define a list for prediction

  result$data_include <- data
  result$char_list <- char_list

  return(result)
}














