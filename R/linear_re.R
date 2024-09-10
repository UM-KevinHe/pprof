#' @importFrom lme4 lmer
#' @importFrom nlme lme
linear_re <- function(Y, Z, ID) {
  if (missing(Y) || missing(Z) || missing(ID))
    stop("Arguments 'Y', 'Z', and 'ID' are all required!", call.=FALSE)

  if (length(Y) != length(ID) | length(ID) != nrow(Z)){
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

  Z <- as.matrix(data[,Z.char])
  Y <- as.matrix(data[, Y.char, drop = F])
  ID <- as.matrix(data[, prov.char, drop = F])

  formula <- as.formula(paste(Y.char, "~ (1|", prov.char, ")+", paste0(Z.char, collapse = "+")))
  fit_re <- lmer(formula, data)

  # Coefficients
  beta <- matrix(fixef(fit_re)[2:length(fixef(fit_re))], ncol=1)
  colnames(beta) <- "beta"
  rownames(beta) <- Z.char

  mu <- matrix(fixef(fit_re)[1])
  colnames(mu) <- "mu"
  rownames(mu) <- "intercept"

  alpha <- as.matrix(ranef(fit_re)$ID, ncol = 1)
  colnames(alpha) <- "alpha"
  rownames(alpha) <- names(n.prov)

  coefficient <- list()
  coefficient$beta <- beta
  coefficient$mu <- mu
  coefficient$alpha <- alpha

  # Variance
  sum = summary(fit_re)
  var_alpha = matrix(as.data.frame(sum$varcor)[1,"sdcor"]^2)
  colnames(var_alpha) <- "Variance.Alpha"
  rownames(var_alpha) <- "ID"

  var_beta = matrix(as.data.frame(sum$coefficients)[2:(length(Z.char)+1),2]^2)
  colnames(var_beta) <- "Variance.Beta"
  rownames(var_beta) <- Z.char

  var_mu = matrix(as.data.frame(sum$coefficients)[1,2]^2)
  colnames(mu) = "Variance.Mu"
  rownames(mu) = "Intercept"

  # prediction
  linear_pred <- Z %*% beta
  alpha.obs = rep(alpha, n.prov)
  pred = rep(mu,n) + alpha.obs + linear_pred

  char_list <- list(Y.char = Y.char,
                    prov.char = prov.char,
                    Z.char = Z.char)

  result <- structure(list(coefficient = coefficient,
                           prediction = pred,
                           linear_pred = linear_pred,
                           observation = Y,
                           prov = data[, prov.char]),
                      class = "linear_re")  #define a list for prediction

  result$data_include <- data
  result$char_list <- char_list

  return(result)
}


linear_re.complete <- function(formula = NULL, data = NULL,
                      Y = NULL, Z = NULL, ID = NULL,
                      Y.char = NULL, Z.char = NULL, ID.char = NULL, ...) {
  if (!is.null(formula) && !is.null(data)) {
    terms <- terms(formula)
    response <- as.character(attr(terms, "variables"))[2]
    predictors <- attr(terms, "term.labels")

    RE_term <- predictors[grepl("\\|", predictors)]
    id_var <- trimws(gsub(".*\\|", "", RE_term))

    Y.char <- response
    ID.char <- id_var
    Z.char <- predictors[!grepl("\\|", predictors)]
    data <- data[order(factor(data[, id_var])),]
    # Y <- data[, Y.char, drop = F]
    # Z <- as.matrix(data[, Z.char, drop = F])
    # ID <- data[, ID.char, drop = F]
    fit_re <- lmer(formula, data, ...)
  }
  else if (!is.null(data) && !is.null(Y.char) && !is.null(Z.char) && !is.null(ID.char)) {
    if (!all(c(Y.char, Z.char, ID.char) %in% colnames(data)))
      stop("Some of the specified columns are not in the data!", call.=FALSE)
    data <- data[order(factor(data[, ID.char])),]
    # Y <- data[, Y.char, drop = F]
    # Z <- as.matrix(data[, Z.char, drop = F])
    # ID <- data[, ID.char, drop = F]

    formula <- as.formula(paste(Y.char, "~ (1|", ID.char, ") +", paste(Z.char, collapse = " + ")))
    fit_re <- lmer(formula, data, ...)
  }
  else if (!is.null(Y) && !is.null(Z) && !is.null(ID)) {
    if (length(Y) != length(ID) | length(ID) != nrow(Z)){
      stop("Dimensions of the input data do not match!!", call.=F)}

    data <- as.data.frame(cbind(Y, ID, Z))
    Y.char <- colnames(data)[1]
    ID.char <- colnames(data)[2]
    Z.char <- colnames(Z)
    data <- data[order(factor(data[,ID.char])),]

    # Z <- as.matrix(data[,Z.char], drop = F)
    # Y <- as.matrix(data[, Y.char, drop = F])
    # ID <- as.matrix(data[, ID.char, drop = F])

    formula <- as.formula(paste(Y.char, "~ (1|", ID.char, ")+", paste0(Z.char, collapse = "+")))

    fit_re <- lmer(formula, data, ...)
  }

  X.model <- model.matrix(fit_re)
  Y <- as.matrix(data[, Y.char, drop = F])
  ID <- as.matrix(data[, ID.char, drop = F])
  data <- as.data.frame(cbind(Y, ID, X.model))

  n.prov <- sapply(split(data[, Y.char], data[, ID.char]), length)
  m <- length(n.prov) # number of providers
  n <- sum(n.prov) # number of observations
  # p <- length(Z.char) # number of covariates

  # Coefficients
  FE_coefficient <- matrix(fixef(fit_re))
  colnames(FE_coefficient) <- "Coefficient"
  rownames(FE_coefficient) <- names(fixef(fit_re))

  RE_coefficient <- as.matrix(ranef(fit_re)[[ID.char]], ncol = 1)
  colnames(RE_coefficient) <- "alpha"
  rownames(RE_coefficient) <- names(n.prov)

  coefficient <- list()
  coefficient$FE <- FE_coefficient
  coefficient$RE <- RE_coefficient

  # Variance
  sum <- summary(fit_re)
  var_alpha <- matrix(as.data.frame(sum$varcor)[1,"sdcor"]^2)
  colnames(var_alpha) <- "Variance.Alpha"
  rownames(var_alpha) <- "ID"

  varcov_FE <- matrix(sum$vcov, ncol = length(FE_coefficient))
  colnames(varcov_FE) <- colnames(sum$vcov)
  rownames(varcov_FE) <- rownames(sum$vcov)

  variance <- list()
  variance$alpha <- var_alpha
  variance$FE <- varcov_FE

  sigma <- sum$sigma

  # prediction
  linear_pred <- X.model %*% FE_coefficient
  colnames(linear_pred) <- "Fixed Fitted"
  rownames(linear_pred) <- seq_len(nrow(linear_pred))

  pred <- matrix(fitted(fit_re), ncol = 1)
  colnames(pred) <- "Prediction"
  rownames(pred) <- seq_len(nrow(pred))

  res <- matrix(residuals(fit_re), ncol = 1)
  colnames(res) <- "Residuals"
  rownames(res) <- seq_len(nrow(res))

  char_list <- list(Y.char = Y.char,
                    ID.char = ID.char,
                    Z.char = Z.char)

  result <- structure(list(coefficient = coefficient,
                           variance = variance,
                           sigma = sigma,
                           prediction = pred,
                           observation = Y,
                           residuals = res,
                           linear_pred = linear_pred,
                           prov = data[, ID.char],
                           model = fit_re),
                      class = "linear_re")  #define a list for prediction

  result$data_include <- data
  result$char_list <- char_list

  return(result)
}


linear_re.nlme <- function(formula = NULL, data = NULL,
                           Y = NULL, Z = NULL, ID = NULL,
                           Y.char = NULL, Z.char = NULL, ID.char = NULL, ...) {
  if (!is.null(formula) && !is.null(data)) {
    terms <- terms(formula)
    response <- as.character(attr(terms, "variables"))[2]
    predictors <- attr(terms, "term.labels")

    RE_term <- predictors[grepl("\\|", predictors)]
    id_var <- trimws(gsub(".*\\|", "", RE_term))

    Y.char <- response
    ID.char <- id_var
    Z.char <- predictors[!grepl("\\|", predictors)]
    data <- data[order(factor(data[, id_var])),]
    # Y <- data[, Y.char, drop = F]
    # Z <- as.matrix(data[, Z.char, drop = F])
    # ID <- data[, ID.char, drop = F]
    fe_formula <- reformulate(Z.char, response = Y.char)
    re_formula <- as.formula(paste("~ 1 |", id_var))
    fit_re <- lme(fixed = fe_formula, random = re_formula, data = data, ...)
  }
  else if (!is.null(data) && !is.null(Y.char) && !is.null(Z.char) && !is.null(ID.char)) {
    if (!all(c(Y.char, Z.char, ID.char) %in% colnames(data)))
      stop("Some of the specified columns are not in the data!", call.=FALSE)
    data <- data[order(factor(data[, ID.char])),]
    # Y <- data[, Y.char, drop = F]
    # Z <- as.matrix(data[, Z.char, drop = F])
    # ID <- data[, ID.char, drop = F]

    fe_formula <- reformulate(Z.char, response = Y.char)
    re_formula <- as.formula(paste("~ 1 |", ID.char))
    fit_re <- lme(fixed = fe_formula, random = re_formula, data = data, ...)
  }
  else if (!is.null(Y) && !is.null(Z) && !is.null(ID)) {
    if (length(Y) != length(ID) | length(ID) != nrow(Z)){
      stop("Dimensions of the input data do not match!!", call.=F)}

    data <- as.data.frame(cbind(Y, ID, Z))
    Y.char <- colnames(data)[1]
    ID.char <- colnames(data)[2]
    Z.char <- colnames(Z)
    data <- data[order(factor(data[,ID.char])),]

    # Z <- as.matrix(data[,Z.char], drop = F)
    # Y <- as.matrix(data[, Y.char, drop = F])
    # ID <- as.matrix(data[, ID.char, drop = F])

    fe_formula <- reformulate(Z.char, response = Y.char)
    re_formula <- as.formula(paste("~ 1 |", ID.char))
    fit_re <- lme(fixed = fe_formula, random = re_formula, data = data, ...)
  }

  X.model <- model.matrix(formula(fit_re), data = data)
  Y <- as.matrix(data[, Y.char, drop = F])
  ID <- as.matrix(data[, ID.char, drop = F])
  data <- as.data.frame(cbind(Y, ID, X.model))

  n.prov <- sapply(split(data[, Y.char], data[, ID.char]), length)
  m <- length(n.prov) # number of providers
  n <- sum(n.prov) # number of observations
  # p <- length(Z.char) # number of covariates

  # Coefficients
  FE_coefficient <- matrix(fixef(fit_re))
  colnames(FE_coefficient) <- "Coefficient"
  rownames(FE_coefficient) <- names(fixef(fit_re))

  RE_coefficient <- as.matrix(ranef(fit_re), ncol = 1)
  colnames(RE_coefficient) <- "alpha"
  rownames(RE_coefficient) <- names(n.prov)

  coefficient <- list()
  coefficient$FE <- FE_coefficient
  coefficient$RE <- RE_coefficient

  # Variance
  sum <- summary(fit_re)
  var_alpha <- matrix(as.numeric(VarCorr(fit_re)[1]))
  colnames(var_alpha) <- "Variance.Alpha"
  rownames(var_alpha) <- "ID"

  varcov_FE <- matrix(sum$varFix, ncol = length(FE_coefficient))
  colnames(varcov_FE) <- colnames(sum$varFix)
  rownames(varcov_FE) <- rownames(sum$varFix)

  variance <- list()
  variance$alpha <- var_alpha
  variance$FE <- varcov_FE

  sigma <- sum$sigma

  # prediction
  linear_pred <- sum$fitted[, 1, drop = FALSE]
  colnames(linear_pred) <- "Fixed"
  pred <- sum$fitted[, 2, drop = FALSE]
  colnames(pred) <- "prediction"
  res <- matrix(residuals(fit_re), ncol = 1)
  colnames(res) <- "Residuals"
  rownames(res) <- seq_len(nrow(res))


  char_list <- list(Y.char = Y.char,
                    ID.char = ID.char,
                    Z.char = Z.char)

  result <- structure(list(coefficient = coefficient,
                           variance = variance,
                           sigma = sigma,
                           prediction = pred,
                           observation = Y,
                           residuals = res,
                           linear_pred = linear_pred,
                           prov = data[, ID.char],
                           model = fit_re),
                      class = "linear_re")  #define a list for prediction

  result$data_include <- data
  result$char_list <- char_list

  return(result)
}
