test_linear_fe <- function(fit, parm, level = 0.95, null = "median",
                           tail = "two", direction = "smaller") {
  alpha <- 1 - level

  data <- fit$data_include
  ID.char <- fit$char_list$ID.char
  gamma <- fit$coefficient$gamma
  se.gamma <- sqrt(fit$variance$gamma)
  n.prov <- sapply(split(data[, fit$char_list$Y.char], data[, ID.char]), length)
  n <- nrow(fit$data_include)
  gamma.null <- ifelse(null=="median", median(gamma),
                       ifelse(null=="mean", sum(n.prov*gamma)/n,
                              ifelse(class(null)=="numeric", null[1],
                                     stop("Argument 'null' NOT as required!",call.=F))))

  m <- length(fit$coefficient$gamma)
  p <- length(fit$coefficient$beta)
  n <- nrow(fit$data_include)

  # test statistics
  stat <- (gamma - gamma.null)/se.gamma

  if (fit$method == "Profile Likelihood") {
    p <- pnorm(stat, lower=F)
  }
  else if (fit$method == "Dummy") {
    df <- n - m - p
    p <- pt(stat, df, lower = F)
  }


  if (tail == "two") {
    flag <- ifelse(p < alpha/2, 1, ifelse(p <= 1-alpha/2, 0, -1))
    p_value <- 2 * pmin(p, 1-p)
  }
  else if (tail == "one") {
    if (direction == "larger") {
      p_value <- p
      flag <- ifelse(p < alpha, 1, 0)
    }
    else if (direction == "smaller") {
      p_value <- 1 - p
      flag <- ifelse(1-p < alpha, -1, 0)
    }
    else {
      stop("Argument 'direction' must be either 'smaller' or 'larger' for one-tail tests.")
    }
  }
  else {
    stop("Argument 'tail' must be either 'one' or 'two'.")
  }

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
