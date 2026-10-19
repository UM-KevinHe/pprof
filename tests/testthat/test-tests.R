test_that("Wald tests of beta and gamma use the glm standard errors", {
  d <- sim_fe(seed = 11)
  fit <- logis_fe(prep_fe(d), tol = 1e-12)
  ref <- glm_fe(d)
  se_ref <- sqrt(diag(vcov(ref)))
  sb <- summary_fe_covar(fit, test = "wald")
  expect_equal(unname(sb$se.beta), unname(se_ref[c("x1", "x2", "x3")]), tolerance = 1e-5)
  expect_equal(sb$Z, sb$beta / sb$se.beta, tolerance = 1e-10)
  tg <- test_fe(fit, test = "wald", null = 0)
  expect_equal(unname(tg$se.gamma), unname(se_ref[grep("^prov", names(se_ref))]), tolerance = 1e-5)
})

test_that("likelihood ratio and score tests of beta agree with glm", {
  d <- sim_fe(seed = 12)
  fit <- logis_fe(prep_fe(d), tol = 1e-12)
  ref <- glm_fe(d)
  lr <- summary_fe_covar(fit, test = "lr")
  drop <- stats::drop1(ref, test = "LRT")
  expect_equal(unname(lr$stat), unname(drop[c("x1", "x2", "x3"), "LRT"]), tolerance = 1e-5)
  sc <- summary_fe_covar(fit, test = "score")
  rao <- stats::drop1(ref, test = "Rao")
  expect_equal(unname(sc$stat), unname(rao[c("x1", "x2", "x3"), "Rao score"]), tolerance = 1e-4)
})

test_that("robust Wald tests need a patient ID and reduce to sandwich formulas", {
  d <- sim_fe(seed = 13)
  fit <- logis_fe(prep_fe(d))
  rg <- test_fe(fit, test = "robust_wald")
  expect_true(all(rg$se.gamma > 0))
  expect_named(rg, c("flag", "se.gamma", "p", "stat"))
  rb <- summary_fe_covar(fit, test = "robust_wald")
  expect_true(all(rb$se.beta > 0))
  dp <- suppressWarnings(suppressMessages(fe_data_prep(Y = d$Y, Z = d$Z, ProvID = d$ProvID, check = FALSE, message = FALSE)))
  for (tst in c("robust_wald", "robust_score", "stabilized_wald", "stabilized_score"))
    expect_error(test_fe(logis_fe(dp), test = tst), "PatientID")
})

# Straightforward R versions of the patient-clustered variances, rows of a
# patient taken in data order.
clustered_ref <- function(fit, p) {
  din <- fit$data_include
  r <- din$Y - din$N * p
  v <- din$N * p * (1 - p)
  pat <- paste(din$prov, din$pat)
  prov_of <- tapply(as.character(din$prov), pat, `[`, 1)
  cs <- tapply(r, pat, sum)
  robust <- tapply(cs^2, prov_of[names(cs)], sum)
  e <- ifelse(v > 0, r / sqrt(v), 0)
  a <- sqrt(v)
  by_pat <- split(seq_along(r), factor(pat, levels = unique(pat)))
  dim <- max(lengths(by_pat))
  S <- matrix(0, dim, dim)
  n_ge <- numeric(dim)
  for (rows in by_pat) {
    s <- length(rows)
    S[1:s, 1:s] <- S[1:s, 1:s] + tcrossprod(e[rows])
    n_ge[1:s] <- n_ge[1:s] + 1
  }
  R <- S / outer(seq_len(dim), seq_len(dim), function(k, l) n_ge[pmax(k, l)])
  W <- vapply(by_pat, function(rows) {
    s <- length(rows)
    drop(t(a[rows]) %*% R[1:s, 1:s, drop = FALSE] %*% a[rows])
  }, numeric(1))
  stable <- tapply(W, prov_of[names(by_pat)], sum)
  list(U = tapply(r, din$prov, sum), robust = robust, stable = stable,
       info = tapply(v, din$prov, sum))
}

test_that("robust and stabilized tests follow their definitions (binary and binomial)", {
  for (binomial in c(FALSE, TRUE)) {
    fit <- logis_fe(prep_fe(sim_fe(seed = 31, binomial = binomial), binomial = binomial), tol = 1e-12)
    g0 <- stats::median(fit$gamma)
    at_null <- clustered_ref(fit, stats::plogis(g0 + as.numeric(fit$linear_pred)))
    at_fit <- clustered_ref(fit, fit$pred)
    rs <- test_fe(fit, test = "robust_score")
    expect_equal(rs$stat, as.numeric(at_null$U / sqrt(at_null$robust)), tolerance = 1e-10)
    ss <- test_fe(fit, test = "stabilized_score")
    expect_equal(ss$stat, as.numeric(at_null$U / sqrt(at_null$stable)), tolerance = 1e-10)
    sw <- test_fe(fit, test = "stabilized_wald")
    se <- as.numeric(sqrt(at_fit$stable) / at_fit$info)
    expect_equal(sw$se.gamma, se, tolerance = 1e-10)
    expect_equal(sw$stat, (as.numeric(fit$gamma) - g0) / se, tolerance = 1e-10)
    expect_named(sw, c("flag", "se.gamma", "p", "stat"))
    expect_named(rs, c("flag", "p", "stat"))
  }
})

test_that("provider tests return one row per provider, subsets and alternatives work", {
  d <- sim_fe(seed = 14)
  fit <- logis_fe(prep_fe(d))
  provs <- rownames(fit$gamma)
  for (tst in c("exact.poisbinom", "wald", "robust_wald", "score",
                "robust_score", "stabilized_wald", "stabilized_score")) {
    res <- test_fe(fit, test = tst)
    expect_equal(rownames(res), provs)
    sub <- test_fe(fit, test = tst, parm = provs[c(2, 5)])
    expect_equal(rownames(sub), provs[c(2, 5)])
    expect_equal(sub$p, res$p[c(2, 5)], tolerance = 1e-12)
  }
  expect_error(test_fe(fit, test = "score", score_modified = FALSE), "removed")
  set.seed(1)
  bt <- test_fe(fit, test = "exact.bootstrap", n = 4000)
  ex <- test_fe(fit, test = "exact.poisbinom")
  expect_lt(max(abs(bt$p - ex$p)), 0.05)
  g <- test_fe(fit, test = "wald", alternative = "greater")
  expect_true(all(g$flag %in% c(0, 1)))
  l <- test_fe(fit, test = "exact.poisbinom", alternative = "less")
  expect_true(all(l$flag %in% c(0, -1)))
})

test_that("binomial exact test equals the Bernoulli-expanded data", {
  d <- sim_fe(m = 6, n_range = c(10, 15), binomial = TRUE, seed = 15)
  fit <- logis_fe(prep_fe(d, binomial = TRUE), tol = 1e-12)
  res <- test_fe(fit, test = "exact.poisbinom")
  # expand every binomial row into N Bernoulli rows and use the same fitted beta
  din <- fit$data_include
  idx <- rep(seq_len(nrow(din)), din$N)
  yb <- unlist(lapply(seq_len(nrow(din)), function(i) c(rep(1, din$Y[i]), rep(0, din$N[i] - din$Y[i]))))
  p0 <- stats::plogis(stats::median(fit$df.prov$gamma_est) + as.matrix(din[idx, c("x1", "x2", "x3")]) %*% fit$beta)
  for (j in seq_len(nrow(res))) {
    rows <- din$prov[idx] == rownames(res)[j]
    obs <- sum(yb[rows])
    p <- 1 - poibin::ppoibin(obs, p0[rows]) + 0.5 * poibin::dpoibin(obs, p0[rows])
    expect_equal(res$p[j], 2 * min(p, 1 - p), tolerance = 1e-10)
  }
})

test_that("confidence intervals cover the estimates and are ordered", {
  fit <- logis_fe(prep_fe(sim_fe(seed = 16)))
  for (tst in c("exact", "score", "wald", "robust_wald")) {
    ci <- confint(fit, option = "gamma", test = tst)
    expect_true(all(ci$gamma.lower <= ci$gamma & ci$gamma <= ci$gamma.upper))
  }
  sm <- confint(fit, stdz = c("indirect", "direct"), test = "score")
  expect_named(sm, c("CI.indirect_ratio", "CI.indirect_rate", "CI.direct_ratio", "CI.direct_rate"))
  r <- sm$CI.direct_rate
  expect_true(all(r$CI_rate.lower <= r$direct_rate & r$direct_rate <= r$CI_rate.upper))
})
