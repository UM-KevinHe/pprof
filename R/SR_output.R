#' Direct and indirect standardized ratios and rates
#'
#' @param fit a model fitted by [logis_fe()]. `fit$gamma` may be replaced by
#'   a longer vector of provider effects (e.g. with the +/-17 values assigned
#'   to providers with 0% / 100% crude rates) before calling, as the PPPW
#'   measure does for direct standardization.
#' @param stdz `"indirect"` (default), `"direct"`, or both.
#' @param measure `"rate"`, `"ratio"`, or both (default).
#' @param null for indirect standardization, the population norm of the
#'   provider effect: `"median"` (default) or a number.
#' @param threads number of OpenMP threads for direct standardization.
#' @param se a logical; for direct standardization also return the
#'   delta-method numerator of the standard error (`numer_se`) and the
#'   expected count scaled by `N_pop` (`numer_measure`). Default `TRUE`.
#'   Ignored when `stdz` does not include `"direct"`.
#' @param N_pop the denominator of the rates: the total number of trials in
#'   the standard population. Defaults to the number of trials among the rows
#'   used to fit the model.
#' @param parm optional subset of providers.
#' @param ... unused.
#'
#' @details With \eqn{p_{ij}(\gamma) = expit(\gamma + Z_{ij}^T\hat\beta)} and
#' \eqn{N_{ij}} trials per row:
#' \itemize{
#' \item indirect: \eqn{E_i = \sum_{j} N_{ij} p_{ij}(\gamma_0)} over provider
#'   i's rows, ratio \eqn{O_i / E_i}, rate = ratio x \eqn{\sum O / N_{pop}};
#' \item direct: \eqn{E_i = \sum_{k} N_k p_k(\hat\gamma_i)} over the whole
#'   standard population, ratio \eqn{E_i / \sum O}, rate = ratio x
#'   \eqn{\sum O / N_{pop}}.
#' }
#' Rates are truncated to \[0, 1\]. The direct-rate standard error used by
#' the measures is `numer_se / N_pop * se(gamma_i)`.
#'
#' @return A list with the requested `indirect.ratio`, `indirect.rate`,
#'   `direct.ratio`, `direct.rate` (one-column matrices), `OE` (observed and
#'   expected counts) and, for direct standardization with `se = TRUE`,
#'   `numer_se`, `N_pop` and `numer_measure`.
#'
#' @references He K. (2019) Indirect and direct standardization for evaluating
#' transplant centers. \emph{Journal of Hospital Administration}, 8(1): 9-14.
#'
#' @export
SR_output <- function(fit, stdz = "indirect", measure = c("rate", "ratio"), null = "median",
                      threads = 1, se = TRUE, N_pop, parm, ...) {
  if (missing(fit)) stop("Argument 'fit' is required!", call. = FALSE)
  check_fit(fit)
  if (!"indirect" %in% stdz & !"direct" %in% stdz) {
    stop("Argument 'stdz' NOT as required!", call. = FALSE)
  }
  if (!"rate" %in% measure & !"ratio" %in% measure) {
    stop("Argument 'measure' NOT as required!", call. = FALSE)
  }

  gamma.prov <- fit$gamma
  Z_beta <- fit$linear_pred
  N <- trials_of(fit$data_include, fit$char_list$N.char)
  if (missing(N_pop)) {
    N_pop <- sum(N)
  }
  ind <- if (missing(parm)) seq_len(NROW(gamma.prov)) else {
    if (is.numeric(parm)) parm <- as.numeric(parm)
    which(rownames(gamma.prov) %in% as.character(parm))
  }

  return_ls <- list()
  OE_list <- list()

  if ("indirect" %in% stdz) {
    gamma.null <- gamma_null_value(null, gamma.prov)
    p0 <- as.numeric(stats::plogis(gamma.null + Z_beta)) # expected prob of events under null
    Exp <- N * p0
    df.prov <- data.frame(Obs_provider = fit$df.prov$Obs_provider,
                          Exp.indirect_provider = sapply(split(Exp, fit$prov), sum),
                          Var.indirect_provider = sapply(split(N * p0 * (1 - p0), fit$prov), sum))
    OE_list$OE_indirect <- df.prov[ind, ]
    df.prov$IS_Ratio <- df.prov$Obs_provider / df.prov$Exp.indirect_provider # indirect standardized ratio: O_i/E_i
    if ("ratio" %in% measure) {
      indirect_stdz.ratio <- matrix(df.prov$IS_Ratio)
      dimnames(indirect_stdz.ratio) <- list(rownames(fit$gamma), "Indirect_standardized.ratio")
      return_ls$indirect.ratio <- indirect_stdz.ratio[ind, , drop = FALSE]
    }
    if ("rate" %in% measure) {
      population_rate <- sum(fit$obs) / N_pop # sum(O_i)/N
      df.prov$IS_Rate <- pmax(pmin(df.prov$IS_Ratio * population_rate, 1), 0) # restricted to 0%-100%
      indirect_stdz.rate <- matrix(df.prov$IS_Rate)
      dimnames(indirect_stdz.rate) <- list(rownames(fit$gamma), "Indirect_standardized.rate")
      return_ls$indirect.rate <- indirect_stdz.rate[ind, , drop = FALSE]
    }
  }

  if ("direct" %in% stdz) {
    res <- computeDirectExp(as.numeric(gamma.prov), as.numeric(Z_beta), N, max(1L, as.integer(threads)))
    Exp <- res[, 1, drop = FALSE]
    numer_se <- res[, 2, drop = FALSE]

    df.prov <- data.frame(Obs_all = rep(sum(fit$obs), length(gamma.prov)), # denominator
                          Exp.direct_all = Exp)                           # numerator
    OE_list$OE_direct <- df.prov[ind, ]
    df.prov$DS_ratio <- df.prov$Exp.direct_all / df.prov$Obs_all # direct standardized ratio: E/sum(O)
    if ("ratio" %in% measure) {
      direct_stdz.ratio <- matrix(df.prov$DS_ratio)
      dimnames(direct_stdz.ratio) <- list(rownames(fit$gamma), "Direct_standardized.ratio")
      return_ls$direct.ratio <- direct_stdz.ratio[ind, , drop = FALSE]
    }
    if ("rate" %in% measure) {
      population_rate <- sum(fit$obs) / N_pop # sum(O_i)/N
      df.prov$DS_Rate <- pmax(pmin(df.prov$DS_ratio * population_rate, 1), 0) # restricted to 0%-100%
      direct_stdz.rate <- matrix(df.prov$DS_Rate)
      dimnames(direct_stdz.rate) <- list(rownames(fit$gamma), "Direct_standardized.rate")
      return_ls$direct.rate <- direct_stdz.rate[ind, , drop = FALSE]
    }
  }
  return_ls$OE <- OE_list

  if (se == TRUE && "direct" %in% stdz) {
    return_ls$numer_se <- numer_se[ind, , drop = FALSE]
    return_ls$N_pop <- N_pop
    if ("rate" %in% measure) return_ls$numer_measure <- return_ls$direct.rate * N_pop
  }
  return_ls
}

#' Direct standardization over an arbitrary standard population
#'
#' Low-level building block of direct standardization. For each provider
#' effect \eqn{\gamma_i} it sums, over the rows of a standard population with
#' linear predictors `Z_beta` and trials `N`,
#' \eqn{E_i = \sum_k N_k\, expit(\gamma_i + Z\beta_k)} and
#' \eqn{V_i = \sum_k N_k\, p_k (1 - p_k)}.
#'
#' The direct standardized rate is \eqn{E_i / \sum_k N_k}; its delta-method
#' standard error (treating \eqn{\hat\beta} as fixed) is
#' \eqn{V_i / \sum_k N_k \times se(\hat\gamma_i)}.
#'
#' @param gamma provider effects.
#' @param Z_beta linear predictors of the standard population.
#' @param N number of trials per row of the standard population
#'   (default 1: binary outcome).
#' @param threads number of OpenMP threads.
#'
#' @return A data frame with columns `Exp` and `Var`, one row per element of
#'   `gamma`.
#'
#' @export
direct_exp <- function(gamma, Z_beta, N = 1, threads = 1) {
  Z_beta <- as.numeric(Z_beta)
  if (length(N) == 1) N <- rep(N, length(Z_beta))
  if (length(N) != length(Z_beta)) stop("'N' and 'Z_beta' must have the same length.", call. = FALSE)
  res <- computeDirectExp(as.numeric(gamma), Z_beta, as.numeric(N), max(1L, as.integer(threads)))
  data.frame(Exp = res[, 1], Var = res[, 2])
}
