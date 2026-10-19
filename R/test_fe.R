#' Hypothesis tests for provider effects
#'
#' Tests \eqn{H_0: \gamma_i = \gamma_0} for every provider (or a subset) and
#' flags providers that are significantly higher or lower than expected.
#'
#' @param fit a model fitted by [logis_fe()].
#' @param parm optional subset of providers (same class as the provider IDs).
#' @param level confidence level; the significance level is `1 - level`.
#'   Default 0.95.
#' @param test the test:
#'   \itemize{
#'   \item `"exact.poisbinom"` (default): exact test using the Poisson-binomial
#'     distribution of the provider's total events under the null (each row
#'     contributes `N` Bernoulli trials);
#'   \item `"exact.bootstrap"`: the same exact test by simulation (`n` draws);
#'   \item `"wald"`: model-based Wald test;
#'   \item `"robust_wald"`: Wald test with a patient-clustered sandwich
#'     standard error (requires `PatientID` in [fe_data_prep()]);
#'   \item `"stabilized_wald"`: Wald test with the stabilized sandwich
#'     standard error (requires `PatientID`; see Details);
#'   \item `"score"`: score test with \eqn{\hat\beta} held fixed;
#'   \item `"robust_score"`: score test with a patient-clustered variance
#'     (requires `PatientID`);
#'   \item `"stabilized_score"`: score test with the stabilized variance
#'     (requires `PatientID`).
#'   }
#' @param null the null value \eqn{\gamma_0}: `"median"` (default) or a number.
#' @param n number of simulated draws for `"exact.bootstrap"`. Default 10,000.
#' @param threads number of OpenMP threads (`"wald"`).
#' @param alternative `"two.sided"` (default), `"greater"` or `"less"`.
#' @param ... unused. The standard score test of version 0.1.0
#'   (`score_modified = FALSE`) was removed.
#'
#' @details The Wald tests are unreliable for providers whose effect is at
#' the bound (no or all events).
#'
#' The patient-clustered tests work with the row residuals
#' \eqn{r_j = Y_j - N_j p_j} and variances \eqn{v_j = N_j p_j (1 - p_j)},
#' at the fitted \eqn{\hat p_j} for the Wald tests and at
#' \eqn{p_j = expit(\gamma_0 + Z_j^T\hat\beta)} for the score tests. With
#' \eqn{U_i = \sum_{j \in i} r_j} over provider \eqn{i}'s rows:
#' \itemize{
#' \item `"robust_wald"`: \eqn{SE = I_{\gamma\gamma}^{-1}\sqrt{V_i}} with
#'   \eqn{V_i = \sum_{k \in i} (\sum_{j \in k} r_j)^2}, summing residuals
#'   within patient \eqn{k};
#' \item `"robust_score"`: \eqn{z = U_i / \sqrt{V_i}};
#' \item `"stabilized_wald"` and `"stabilized_score"` replace \eqn{V_i} by
#'   \eqn{\sum_{k \in i} a_k^T \hat R a_k}, where \eqn{a_k} holds
#'   \eqn{\sqrt{v_j}} for patient \eqn{k}'s rows and \eqn{\hat R} is one
#'   matrix pooled over all patients of all providers: entry \eqn{(s, t)}
#'   averages \eqn{e_s e_t}, \eqn{e_j = r_j / \sqrt{v_j}}, over the patients
#'   with at least \eqn{\max(s, t)} rows.
#' }
#' \eqn{I_{\gamma\gamma}} is the provider's own information
#' \eqn{\sum_j v_j}. The stabilized tests take a patient's rows in their order
#' in the data: sort by time within patient before [fe_data_prep()] (which
#' keeps that order). All providers enter \eqn{\hat R}; `parm` only selects
#' the output. With \eqn{N = 1} every formula is the binary version.
#'
#' @return A data frame, one row per tested provider (row names = provider
#' IDs), with columns `flag` (factor: 1 = higher than expected, -1 = lower,
#' 0 = as expected), `p`, `stat`, and `se.gamma` for the Wald tests. The
#' attribute `"provider size"` holds the number of rows per provider.
#'
#' @references
#' Wu W, Yang Y, Kang J, He K. (2022) Improving large-scale estimation and
#' inference for profiling health care providers. \emph{Statistics in
#' Medicine}, 41(15): 2840-2853.
#'
#' @export
test_fe <- function(fit, parm, level = 0.95, test = "exact.poisbinom", null = "median", n = 10000,
                    threads = 1, alternative = "two.sided", ...) {
  if (missing(fit)) stop("Argument 'fit' is required!", call. = FALSE)
  check_fit(fit)
  if (isFALSE(list(...)$score_modified)) # would otherwise vanish into `...`
    stop("The standard score test (score_modified = FALSE) was removed in pprofR 0.2.0; ",
         "test = 'score' keeps beta fixed.", call. = FALSE)
  clustered <- c("robust_wald", "stabilized_wald", "robust_score", "stabilized_score")
  if (!(test %in% c("exact.poisbinom", "exact.bootstrap", "score", "wald", clustered)))
    stop("Argument 'test' NOT as required!", call. = FALSE)
  if (test %in% clustered && is.null(fit$char_list$PatientID.char))
    stop(sprintf("Argument 'test' = '%s' requires 'PatientID' in the data!", test), call. = FALSE)
  alternative <- match_alternative(alternative)
  alpha <- 1 - level

  N.char <- fit$char_list$N.char
  Y.char <- fit$char_list$Y.char
  Z.char <- fit$char_list$Z.char
  prov.char <- fit$char_list$prov.char
  data <- fit$data_include
  N <- trials_of(data, N.char)
  prov <- data[, prov.char]
  n.prov <- sapply(split(data[, Y.char], prov), length)
  prov_levels <- names(n.prov)

  gamma <- fit$df.prov$gamma_est
  beta <- fit$beta
  gamma.null <- gamma_null_value(null, gamma)
  indices <- provider_indices(if (missing(parm)) NULL else parm, prov, prov_levels)
  selected <- prov_levels[indices]

  if (test %in% c("exact.poisbinom", "exact.bootstrap", "score")) {
    # tests evaluated provider by provider under the null
    keep <- as.character(prov) %in% selected
    Z_null <- unname(as.matrix(data[keep, Z.char]))
    probs <- as.numeric(stats::plogis(gamma.null + Z_null %*% beta))
    obs_row <- data[keep, Y.char]
    N_row <- N[keep]
    grp <- factor(as.character(prov[keep]), levels = selected)

    if (test == "score") {
      z.score <- sapply(split(obs_row - N_row * probs, grp), sum) /
        sqrt(sapply(split(N_row * probs * (1 - probs), grp), sum))
      out <- z_test_output(z.score, alpha, alternative)
      res <- data.frame(flag = factor(out$flag), p = out$p, stat = z.score, row.names = selected)
    } else {
      if (test == "exact.bootstrap" && (n <= 0 | as.integer(n) != n))
        stop("Argument 'n' NOT a positive integer ! ", call. = FALSE)
      one_provider <- function(idx) {
        pr <- probs[idx]
        nt <- N_row[idx]
        obs <- sum(obs_row[idx])
        if (test == "exact.poisbinom") {
          pr_trials <- if (all(nt == 1)) pr else rep(pr, nt) # each row = N Bernoulli trials
          exact_poisbinom(obs, pr_trials, alpha, alternative)
        } else {
          exact_bootstrap(obs, pr, nt, n, alpha, alternative)
        }
      }
      results <- sapply(split(seq_along(probs), grp), one_provider)
      res <- data.frame(flag = factor(results[1, ]), p = results[2, ], stat = results[3, ],
                        row.names = selected)
    }
    attr(res, "provider size") <- n.prov[indices]
    return(res)
  }

  if (test %in% c("stabilized_wald", "robust_score", "stabilized_score")) {
    # every provider enters the stabilized pooling; parm only selects the output
    m <- length(prov_levels)
    prov_row <- provider_codes(prov, prov_levels)
    cl <- cluster_codes(prov, data[, fit$char_list$PatientID.char])
    prov_cl <- prov_row[cl$first]
    p <- if (test == "stabilized_wald") fit$pred else
      as.numeric(stats::plogis(gamma.null + fit$linear_pred)) # under the null
    residual <- data[, Y.char] - N * p
    if (test == "robust_score") {
      s <- cluster_sums(residual, cl$code, prov_row, prov_cl, m)
      stat <- s$U[indices] / sqrt(s$meat[indices])
    } else {
      rows <- max(tabulate(cl$code))
      if (rows > 5000)
        stop("The stabilized tests pool a matrix over the rows of a patient; the largest patient has ",
             rows, " rows.", call. = FALSE)
      meat <- stabilized_meat(residual, N * p * (1 - p), cl$code, prov_cl, m)[indices]
      if (test == "stabilized_score") {
        stat <- cluster_sums(residual, cl$code, prov_row, prov_cl, m)$U[indices] / sqrt(meat)
      } else {
        se.gamma <- robust_wald_gamma(N, as.matrix(data[, Z.char]), fit$pred, n.prov, indices)
        se.gamma <- sqrt(se.gamma^4 * meat)
        stat <- (gamma[indices] - gamma.null) / se.gamma
      }
    }
    stat[!is.finite(stat)] <- NA
    out <- z_test_output(stat, alpha, alternative)
    res <- if (test == "stabilized_wald") {
      data.frame(flag = factor(out$flag), se.gamma = se.gamma, p = out$p, stat = stat, row.names = selected)
    } else {
      data.frame(flag = factor(out$flag), p = out$p, stat = stat, row.names = selected)
    }
    attr(res, "provider size") <- n.prov[indices]
    return(res)
  }

  if (test == "wald") {
    se.gamma <- wald_gamma(N, as.matrix(data[, Z.char]), fit$pred, n.prov, indices, threads)
  } else { # robust_wald
    residual <- data[, Y.char] - N * fit$pred # O - E at the row level
    meat <- robust_gamma_meat(residual, prov, data[, fit$char_list$PatientID.char], prov_levels)
    se.gamma <- robust_wald_gamma(N, as.matrix(data[, Z.char]), fit$pred, n.prov, indices)
    se.gamma <- sqrt(se.gamma^4 * meat[indices])
  }
  stat <- (gamma[indices] - gamma.null) / se.gamma
  out <- z_test_output(stat, alpha, alternative)
  res <- data.frame(flag = factor(out$flag), se.gamma = se.gamma, p = out$p, stat = stat,
                    row.names = selected)
  attr(res, "provider size") <- n.prov[indices]
  res
}

# Exact test for one provider: `probs` holds one success probability per
# Bernoulli trial. Returns c(flag, p-value, z-score).
exact_poisbinom <- function(obs, probs, alpha, alternative) {
  if (alternative == "two.sided") {
    p <- 1 - poibin::ppoibin(obs, probs) + 0.5 * poibin::dpoibin(obs, probs) # ppoibin: P(#events <= obs)
    z.score <- stats::qnorm(p, lower.tail = FALSE)
    flag <- flag_two_sided(p, alpha)
    p.val <- 2 * min(p, 1 - p)
  } else if (alternative == "greater") {
    p <- 1 - poibin::ppoibin(obs - 1, probs)
    z.score <- stats::qnorm(p, lower.tail = FALSE)
    p.val <- p
    flag <- ifelse(p < alpha, 1, 0)
  } else {
    p <- poibin::ppoibin(obs, probs)
    z.score <- stats::qnorm(p, lower.tail = TRUE)
    p.val <- p
    flag <- ifelse(p < alpha, -1, 0)
  }
  c(flag, p.val, z.score)
}

# Simulation version of the exact test for one provider.
exact_bootstrap <- function(obs, probs, trials, n, alpha, alternative) {
  size <- if (all(trials == 1)) 1 else rep(trials, times = n)
  sums <- colSums(matrix(stats::rbinom(n = length(probs) * n, size = size, prob = rep(probs, times = n)), ncol = n))
  if (alternative == "two.sided") {
    p <- (sum(sums > obs) + 0.5 * sum(sums == obs)) / n
    z.score <- stats::qnorm(p, lower.tail = FALSE)
    p.val <- 2 * min(p, 1 - p)
    flag <- flag_two_sided(p, alpha)
  } else if (alternative == "greater") {
    p <- sum(sums >= obs) / n
    z.score <- stats::qnorm(p, lower.tail = FALSE)
    p.val <- p
    flag <- ifelse(p < alpha, 1, 0)
  } else {
    p <- sum(sums <= obs) / n
    z.score <- stats::qnorm(p, lower.tail = TRUE)
    p.val <- p
    flag <- ifelse(p < alpha, -1, 0)
  }
  c(flag, p.val, z.score)
}
