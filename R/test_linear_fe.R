test_linear_fe <- function(fit, parm, level = 0.95, null = "median", alternative = "two.sided") {
  alpha <- 1 - level

  data <- fit$data_include
  ID.char <- fit$char_list$ID.char
  gamma <- fit$coefficient$gamma
  se.gamma <- sqrt(fit$variance$gamma)
  gamma.null <- ifelse(null=="median", median(gamma),
                       ifelse(null=="mean", sum(n.prov*gamma)/n,
                              ifelse(class(null)=="numeric", null[1],
                                     stop("Argument 'null' NOT as required!",call.=F))))

  n.prov <- sapply(split(data[, fit$char_list$Y.char], data[, ID.char]), length)
  m <- length(fit$coefficient$gamma)
  p <- length(fit$coefficient$beta)
  n <- nrow(fit$data_include)

  # test statistics
  stat <- (gamma - gamma.null)/se.gamma

  prob <- switch(fit$method,
              "Profile Likelihood" = pnorm(stat, lower=F),
              "Dummy" = pt(stat, df = n - m - p, lower = F))

  if (alternative == "two.sided") {
    flag <- ifelse(prob < alpha/2, 1, ifelse(prob <= 1-alpha/2, 0, -1))
    p_value <- 2 * pmin(prob, 1-prob)
  }
  else if (alternative == "greater") {
    flag <- ifelse(prob < alpha, 1, 0)
    p_value <- prob
  }
  else if (alternative == "less") {
    flag <- ifelse(1 - prob < alpha, -1, 0)
    p_value <- 1 - prob
  }
  else {
    stop("Argument 'alternative' should be one of 'two.sided', 'less', 'greater'")
  }

  result <- data.frame(flag = factor(flag), p = p_value, stat = stat, Std.Error = se.gamma)
  colnames(result) <- c("flag", "p value", "stat", "Std.Error")

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
