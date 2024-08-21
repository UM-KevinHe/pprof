test_linear_fe <- function(fit, parm, level = 0.95, null = "median", ref.dis = "normal") {
  alpha <- 1 - level

  data <- fit$data_include
  prov.char <- fit$char_list$prov.char
  gamma <- fit$coefficient$gamma
  sd.gamma <- sqrt(fit$variance$gamma)
  gamma.null <- ifelse(null=="median", median(gamma),
                       ifelse(is.numeric(null), null[1],
                              stop("Argument 'null' NOT as required!", call.=F)))

  m <- length(fit$coefficient$gamma)
  p <- length(fit$coefficient$beta)
  n <- nrow(fit$data_include)

  # test statistics
  stat <- (gamma - gamma.null)/sd.gamma

  if (ref.dis == "normal") {
    p <- pnorm(stat, lower=F)
  }
  else if (ref.dis == "t") {
    df <- n - m - p
    p <- pt(stat, df, lower = F)
  }
  else {
    stop("Reference distribution must be either normal distribution or t distribution")
  }

  flag <- ifelse(p < alpha/2, 1, ifelse(p <= 1-alpha/2, 0, -1))
  p_value <- 2 * pmin(p, 1-p)

  result <- data.frame(flag = factor(flag), p = p_value, stat = stat, Std.Error = sd.gamma)
  colnames(result) <- c("flag", "p", "stat", "Std.Error")

  if (missing(parm)) {return(result)}
  else {
    if (is.integer(parm)) {  #avoid "integer" class
      parm <- as.numeric(parm)
    }
    if (class(parm) == class(data[, prov.char])) {
      result = result[row.names(result) %in% parm, ]
      return(result)
    } else {
      stop("Argument 'parm' includes invalid elements!")
    }
  }
}
