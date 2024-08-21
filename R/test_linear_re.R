test_linear_re <- function(fit, parm, level = 0.95, null = 0) {
  alpha <- 1 - level

  data <- fit$data_include
  Y.char <- fit$char_list$Y.char
  ID.char <- fit$char_list$ID.char
  Z.char <- fit$char_list$Z.char
  beta <- fit$coefficient$beta
  mu <- fit$coefficient$mu

  var_alpha <- fit$variance$alpha
  sigma_sq <- fit$sigma^2
  n.prov <- sapply(split(data[, Y.char], data[, ID.char]), length)
  R_i <- as.vector(var_alpha) / (as.vector(var_alpha) + as.vector(sigma_sq) / n.prov)

  # Y_bar_i <- sapply(split(data[, Y.char], data[, ID.char]), mean)
  # Z_bar_i <- t(matrix(sapply(split(data[,Z.char, drop = FALSE], data[,ID.char]), colMeans),
  #                             ncol=length(Y_bar_i), nrow = length(beta)))
  #
  # estimate_alpha <- Y_bar_i - Z_bar_i%*%beta - mu

  denom <- sqrt(R_i * sigma_sq / n.prov)
  Z_score <- (fit$coefficient$alpha - null)/denom
  p <- pnorm(Z_score, lower=F)
  p_value <- 2 * pmin(p, 1-p)
  flag <- ifelse(p < alpha/2, 1, ifelse(p > 1 - alpha/2, -1, 0))

  result <- data.frame(flag = factor(flag), p = p_value, stat = Z_score, Std.Error = denom)
  colnames(result) <- c("flag", "p", "stat", "Std.Error")

  if (missing(parm)) {return(result)}
  else {
    if (is.numeric(parm)) {  #avoid "integer" class
      parm <- as.numeric(parm)
    }
    if (class(parm) == class(data[, ID.char])) {
      result = result[row.names(result) %in% parm, ]
      return(result)
    } else {
      stop("Argument 'parm' includes invalid elements!")
    }
  }
}
