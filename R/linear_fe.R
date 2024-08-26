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


linear_fe <- function(Y, Z, ID, Rcpp = TRUE){
  if (missing(Y) || missing(Z) || missing(ID))
    stop("Arguments 'Y', 'Z', and 'ID' are all required!", call.=FALSE)

  if (length(Y) != length(ID) | (length(ID) != nrow(Z))){
    stop("Dimensions of the input data do not match!!", call.=F)
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

  if (Rcpp == F) {
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

    char_list <- list(Y.char = Y.char,
                      prov.char = prov.char,
                      Z.char = Z.char)

    result <- structure(list(coefficient = coefficient,
                             linear_pred = linear_pred,
                             prediction = pred,
                             observation = Y,
                             prov = data[, prov.char]),
                        class = "linear_fe")  #define a list for prediction

    result$data_include <- data
    result$char_list <- char_list
  }

  if (Rcpp == T) {
    Y <- as.vector(Y)
    ID <- as.vector(ID)
    ls <- compute_profilkd_linear(Y, Z, ID, n.prov)
    gamma.prov <- as.numeric(ls$gamma); beta <- as.numeric(ls$beta)
    return_ls$coefficient$gamma <- gamma.prov
    return_ls$coefficient$beta <- beta
  }

  return(result)
}



linear_fe.fit <- function(formula, data, ID.char, message = TRUE){
  if (missing(formula) || missing(data))
    stop("Arguments 'formula' and 'data' are required!", call.=FALSE)

  if (!inherits(formula, "formula"))
    stop("Argument 'formula' must be a formula!", call.=FALSE)

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

  result <- linear_fe(Y, Z, ID, message)

  return(result)
}


linear_fe.data <- function(data, Y.char, Z.char, ID.char, message) {
  if (missing(data) || missing(Y.char) || missing(Z.char) || missing(ID.char))
    stop("Arguments 'data', 'Y.char', 'Z.char', and 'ID.char' are all required!", call.=FALSE)

  if (!all(c(Y.char, Z.char, ID.char) %in% colnames(data)))
    stop("Some of the specified columns are not in the data!", call.=FALSE)

  Y = data[, Y.char]
  Z = data[, Z.char, drop = FALSE]
  ID = data[, ID.char]

  result <- linear_fe(Y, Z, ID, message)

  return(result)
}


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
      # Z <- model.matrix(~ data[[predictors]] - 1)
      ID <- data[,id_var, drop = F]
    }
    else if (!is.null(data) && !is.null(Y.char) && !is.null(Z.char) && !is.null(ID.char)) {
      if (!all(c(Y.char, Z.char, ID.char) %in% colnames(data)))
        stop("Some of the specified columns are not in the data!", call.=FALSE)

      Y <- data[, Y.char]
      Z <- model.matrix(reformulate(Z.char), data)[, -1, drop = FALSE]
      ID <- data[, ID.char, drop = F]
    }
    else if (!is.null(Y) && !is.null(Z) && !is.null(ID)) {
      if (length(Y) != length(ID) | (length(ID) != nrow(Z))) {
        stop("Dimensions of the input data do not match!!", call.=F)
      }

      Z.char <- colnames(Z)
      Z <- model.matrix(reformulate(Z.char), Z)[, -1, drop = FALSE]
    }
    else {
      stop("Insufficient or incompatible arguments provided. Please provide either (1) formula and data, (2) data, Y.char, Z.char, and ID.char, or (3) Y, Z, and ID.", call.=FALSE)
    }

    data <- as.data.frame(cbind(Y, ID, Z))
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
    gamma.obs <- rep(gamma.prov, n.prov)
    pred <- gamma.obs + linear_pred
    colnames(pred) <- "Prediction"
    rownames(pred) <- seq_len(nrow(pred))
    residuals <- matrix(Y - pred, ncol = 1)
    colnames(residuals) <- "Residuals"
    rownames(residuals) <- seq_len(nrow(residuals))
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
                      ID.char = ID.char,
                      Z.char = Z.char)

  }
  else if (method == "lm") {
    if (!is.null(formula) && !is.null(data) && !is.null(ID.char)){
      original_ID <- data[, ID.char, drop = F]
      data[,ID.char] <- as.factor(data[,ID.char])
      formula_terms <- all.vars(formula)
      formula_terms <- formula_terms[formula_terms != ID.char]
      Z.char <- formula_terms[2:length(formula_terms)]
      Y.char <- formula_terms[1]
      # new_formula_terms <- c(formula_terms, ID.char)
      # ID.char is always in the last position
      new_formula <- as.formula(paste(Y.char, "~", ID.char, "+",
                                      paste(Z.char, collapse = " + ")))

      data <- data[order(factor(data[,ID.char])),]
      formula <- update(new_formula, . ~ . - 1)
      fit_lm <- lm(formula, data = data)
    }
    else if (!is.null(data) && !is.null(Y.char) && !is.null(Z.char) && !is.null(ID.char)) {
      original_ID <- data[, ID.char, drop = F]
      data[,ID.char] <- as.factor(data[,ID.char])
      data <- data[order(factor(data[,ID.char])),]
      formula <- as.formula(paste(Y.char, "~", ID.char, "+", paste(Z.char, collapse = " + "), "-1"))
      fit_lm <- lm(formula, data = data)
    }
    else if (!is.null(Y) && !is.null(Z) && !is.null(ID)) {
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
    ID <- as.matrix(data[, ID.char, drop = F])

    n.prov <- sapply(split(data[, Y.char], data[, ID.char]), length)
    m <- length(n.prov) # number of providers
    n <- sum(n.prov) # number of observations
    Z.char <- colnames(X.model)[(m+1):length(colnames(X.model))]
    Z <- X.model[,Z.char, drop = F]
    p <- length(Z.char) # number of covariates
    data <- as.data.frame(cbind(Y, ID, Z))

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
    X <- model.matrix(fit_lm)
    varcor <- sigma_hat_sq * solve(t(X)%*%X)
    varcor_beta <- matrix(varcor[(m+1):(m+p), (m+1):(m+p)], ncol = p, nrow = p)
    rownames(varcor_beta) <- Z.char
    colnames(varcor_beta) <- Z.char

    var_gamma <- matrix(diag(varcor)[1:m], ncol = 1)
    rownames(var_gamma) <- names(n.prov)
    colnames(var_gamma) <- "Variance.Gamma"

    variance <- list()
    variance$beta <- varcor_beta
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

    char_list <- list(Y.char = Y.char,
                      ID.char = ID.char,
                      Z.char = Z.char)

    # Restore the ID column to keep the original class
    data[, ID.char] <- original_ID
  }
  else stop("Method should be either 'pl' or 'lm'")


  result <- structure(list(coefficient = coefficient,
                           variance = variance,
                           sigma = sqrt(sigma_hat_sq),
                           linear_pred = linear_pred,
                           prediction = pred,
                           observation = Y,
                           residuals = residuals,
                           prov = data[, ID.char]),
                      class = "linear_fe")  #define a list for prediction

  result$data_include <- data
  result$char_list <- char_list

  return(result)
}














