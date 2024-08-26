test_linear_fe <- function(fit, parm, level = 0.95, null = "median", method = "pl") {
  alpha <- 1 - level

  data <- fit$data_include
  ID.char <- fit$char_list$ID.char
  gamma <- fit$coefficient$gamma
  se.gamma <- sqrt(fit$variance$gamma)
  gamma.null <- ifelse(null=="median", median(gamma),
                       ifelse(is.numeric(null), null[1],
                              stop("Argument 'null' NOT as required!", call.=F)))

  m <- length(fit$coefficient$gamma)
  p <- length(fit$coefficient$beta)
  n <- nrow(fit$data_include)

  # test statistics
  stat <- (gamma - gamma.null)/se.gamma

  if (method == "pl") {
    p <- pnorm(stat, lower=F)
  }
  else if (method == "lm") {
    df <- n - m - p
    p <- pt(stat, df, lower = F)
  }
  else {
    stop("Method should be either 'pl' or 'lm'")
  }

  flag <- ifelse(p < alpha/2, 1, ifelse(p <= 1-alpha/2, 0, -1))
  p_value <- 2 * pmin(p, 1-p)

  result <- data.frame(flag = factor(flag), p = p_value, stat = stat, Std.Error = se.gamma)
  colnames(result) <- c("flag", "p", "stat", "Std.Error")

  if (missing(parm)) {return(result)}
  else {
    if (is.integer(parm)) {  #avoid "integer" class
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
