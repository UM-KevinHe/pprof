SR_linear <- function(fit, stdz = "indirect", null = "median") {
  data = fit$data_include
  Z_beta <- fit$linear_pred
  n <- nrow(fit$data_include)

  return_ls <- list()
  OE_list <- list()

  # linear FE
  if (class(fit) == "linear_fe") {
    gamma.prov <- fit$coefficient$gamma
    n.prov <- sapply(split(data[, fit$char_list$Y.char], data[, fit$char_list$prov.char]), length)
    gamma.null <- ifelse(null=="median", median(gamma.prov),
                         ifelse(null=="mean", sum(n.prov*gamma.prov)/n,
                                ifelse(class(null)=="numeric", null[1],
                                       stop("Argument 'null' NOT as required!",call.=F))))
    if ("indirect" %in% stdz) {
      Exp <- gamma.null + Z_beta
      Exp.indirect_provider <- sapply(split(Exp, fit$prov), sum)
      Obs.indirect_provider <- sapply(split(fit$observation, fit$prov), sum)

      indirect_stdz.diff <- matrix((Obs.indirect_provider - Exp.indirect_provider)/n.prov)
      dimnames(indirect_stdz.diff) <- list(rownames(gamma.prov), "Indirect_standardized.difference")
      return_ls$indirect.difference <- indirect_stdz.diff

      OE.df <- data.frame(Obs = Obs.indirect_provider, Exp = Exp.indirect_provider)
      OE_list$OE_indirect <- OE.df
    }

    if ("direct" %in% stdz) {
      Obs.direct_provider <- sum(gamma.null + fit$linear_pred)
      Exp.direct <- function(gamma){
        sum(gamma + Z_beta)
      }
      Exp.direct_provider <- sapply(gamma.prov, Exp.direct)
      direct_stdz.diff <- matrix((Exp.direct_provider - Obs.direct_provider)/n)
      dimnames(direct_stdz.diff) <- list(rownames(gamma.prov), "Direct_standardized.difference")
      return_ls$direct.difference <- direct_stdz.diff

      OE.df <- data.frame(Obs = Obs.direct_provider, Exp = Exp.direct_provider)
      OE_list$OE_direct <- OE.df
    }
  }

  # linear RE
  if (class(fit) == "linear_re") {
    alpha.prov = fit$coefficient$alpha
    if ("indirect" %in% stdz) {
      n.prov <- sapply(split(data[, fit$char_list$Y.char], data[, fit$char_list$ID.char]), length)

      Exp <- rep(fit$coefficient$mu, n) + Z_beta
      Exp.indirect_provider <- sapply(split(Exp, fit$prov), sum)
      Obs.indirect_provider <- sapply(split(fit$prediction, fit$prov), sum)

      indirect_stdz.diff <- matrix((Obs.indirect_provider - Exp.indirect_provider)/n.prov)
      dimnames(indirect_stdz.diff) <- list(rownames(alpha.prov), "Indirect_standardized.difference")
      return_ls$indirect.difference <- indirect_stdz.diff

      OE.df <- data.frame(Obs = Obs.indirect_provider, Exp = Exp.indirect_provider)
      OE_list$OE_indirect <- OE.df
    }

    if ("direct" %in% stdz) {
      mu = as.vector(fit$coefficient$mu)
      Obs.direct_provider <- sum(fit$observation)
      Exp.direct <- function(alpha){
        sum(mu + alpha + Z_beta)
      }
      Exp.direct_provider <- sapply(alpha.prov, Exp.direct)
      direct_stdz.diff <- matrix((Exp.direct_provider - Obs.direct_provider)/n)
      dimnames(direct_stdz.diff) <- list(rownames(alpha.prov), "Direct_standardized.difference")
      return_ls$direct.difference <- direct_stdz.diff

      OE.df <- data.frame(Obs = Obs.direct_provider, Exp = Exp.direct_provider)
      OE_list$OE_direct <- OE.df
    }
  }

  return_ls$OE <- OE_list
  return(return_ls)
}
