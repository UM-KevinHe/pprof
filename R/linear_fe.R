#' Fit the linear model using profile likelihood
#'
#' Fitting a linear model for data with high-dimensional center effect through profile likelihood estimation

linear_fe <- function(Y, Z, ID){
  if (missing(Y) || missing(Z) || missing(ID))
    stop("Arguments 'Y', 'Z', and 'ID' are all required!", call.=FALSE)

  if (length(Y) != length(ID) | length(ID) != nrow(Z)){
    stop("Dimensions of the input data do not match!!", call.=F)
  }

  data <- as.data.frame(cbind(Y, ID, Z))
  Y.char <- colnames(data)[1]
  prov.char <- colnames(data)[2]
  Z.char <- colnames(Z)
  Z <- as.matrix(Z)

  data <- data[order(factor(data[,prov.char])),]
  n.prov <- sapply(split(data[, Y.char], data[, prov.char]), length)
  m <- length(n.prov)
  sum.first.term <- matrix(0, nrow = length(Z.char), ncol = length(Z.char))
  sum.second.term <- matrix(0, nrow = length(Z.char), ncol = 1)
  for (j in names(n.prov)){
    temp.X <- as.matrix(data[which(data$ID == j), Z.char])
    temp.Y <- as.matrix(data[which(data$ID == j), Y.char, drop = F])
    n <- length(temp.Y)
    Qn <- diag(1, nrow = n, ncol = n) - matrix(1, nrow = n, ncol = n)/n
    sum.first.term <- sum.first.term + t(temp.X) %*% Qn %*% temp.X
    sum.second.term <- sum.second.term + t(temp.X) %*% Qn %*% temp.Y
  }
  beta <- solve(sum.first.term) %*% sum.second.term
  colnames(beta) <- "beta"

  gamma.prov <- c()
  for (j in names(n.prov)){
    temp.y_bar <- mean(as.matrix(data[which(data$ID == j), Y.char, drop = F]))
    temp.x_bar <- as.matrix(colMeans(as.matrix(data[which(data$ID == j), Z.char])))
    gamma.prov <- c(gamma.prov, temp.y_bar - t(temp.x_bar) %*% beta)
  }
  gamma.prov <- as.matrix(gamma.prov)
  colnames(gamma.prov) <- "gamma"
  rownames(gamma.prov)<- names(n.prov)


  linear_pred <- Z %*% beta
  gamma.obs <- rep(gamma.prov, n.prov)
  pred = gamma.obs + linear_pred
  char_list <- list(Y.char = Y.char,
                    prov.char = prov.char,
                    Z.char = Z.char)

  result <- structure(list(beta = beta,
                           gamma = gamma.prov,
                           linear_pred = linear_pred,
                           pred = pred,
                           obs = data[, Y.char],
                           prov = data[, prov.char]),
                      class = "linear_fe")  #define a list for prediction

  result$data_include <- data
  result$char_list <- char_list

  return(result)
}
