test_that("indirect and direct standardization follow their definitions", {
  fit <- logis_fe(prep_fe(sim_fe(seed = 21)))
  sr <- SR_output(fit, stdz = c("indirect", "direct"), measure = c("ratio", "rate"))
  zb <- as.numeric(fit$linear_pred)
  g0 <- stats::median(fit$gamma)
  E_ind <- tapply(stats::plogis(g0 + zb), fit$prov, sum)
  O <- tapply(fit$obs, fit$prov, sum)
  expect_equal(as.numeric(sr$indirect.ratio), as.numeric(O / E_ind), tolerance = 1e-12)
  E_dir <- sapply(as.numeric(fit$gamma), function(g) sum(stats::plogis(g + zb)))
  expect_equal(as.numeric(sr$direct.ratio), E_dir / sum(fit$obs), tolerance = 1e-12)
  expect_equal(as.numeric(sr$direct.rate), pmin(E_dir / length(zb), 1), tolerance = 1e-12)
  expect_equal(as.numeric(sr$numer_se),
               sapply(as.numeric(fit$gamma), function(g) { p <- stats::plogis(g + zb); sum(p * (1 - p)) }),
               tolerance = 1e-12)
})

test_that("SR_output default call works (indirect only, se ignored)", {
  fit <- logis_fe(prep_fe(sim_fe(seed = 22)))
  sr <- SR_output(fit)
  expect_named(sr, c("indirect.ratio", "indirect.rate", "OE"))
})

test_that("binomial direct standardization weights rows by their trials", {
  fit <- logis_fe(prep_fe(sim_fe(binomial = TRUE, seed = 23), binomial = TRUE))
  zb <- as.numeric(fit$linear_pred)
  N <- fit$data_include$N
  sr <- SR_output(fit, stdz = "direct", measure = "rate")
  E <- sapply(as.numeric(fit$gamma), function(g) sum(N * stats::plogis(g + zb)))
  expect_equal(as.numeric(sr$direct.rate), E / sum(N), tolerance = 1e-12)
  de <- direct_exp(fit$gamma, zb, N = N)
  expect_equal(de$Exp, E, tolerance = 1e-12)
})

test_that("direct standardization is exact with repeated linear predictors and extreme effects", {
  set.seed(24)
  n <- 3000
  zb <- as.numeric(cbind(sample(0:4, n, TRUE), rbinom(n, 1, 0.4)) %*% c(0.3, -0.5)) # 10 distinct values
  N <- sample(1:30, n, TRUE)
  g <- c(-2, 0, 1.5, -17, 17, -800, 800)
  de <- direct_exp(g, zb, N = N, threads = 2)
  expect_equal(de$Exp, sapply(g, function(x) sum(N * stats::plogis(x + zb))), tolerance = 1e-12)
  expect_equal(de$Var, sapply(g, function(x) { p <- stats::plogis(x + zb); sum(N * p * (1 - p)) }),
               tolerance = 1e-12)
  expect_equal(direct_exp(g, zb, N = N, threads = 1), de)
})

test_that("test_SM computes calibrated z-tests and intervals", {
  m <- c(0.20, 0.35, 0.50)
  se <- c(0.05, 0.02, 0.20)
  r <- test_SM(m, se, ref = 0.30, crit = 1.96)
  expect_equal(r$z, (m - 0.30) / se)
  expect_equal(r$flag, c(-1, 1, 0))
  expect_equal(r$lower, m - 1.96 * se)
  expect_equal(r$p, 2 * pnorm(-abs(r$z)), tolerance = 1e-12)
  lg <- test_SM(m, se, ref = 0.30, transform = "logit", location = 0.1, scale = 1.81, crit = 1.96)
  z <- (qlogis(m) - qlogis(0.30)) / (se / (m * (1 - m)))
  expect_equal(lg$z, z, tolerance = 1e-12)
  expect_equal(lg$prob, pnorm(z, 0.1, 1.81), tolerance = 1e-12)
  expect_equal(lg$upper, plogis(qlogis(m) + (1.96 * 1.81 - 0.1) * se / (m * (1 - m))), tolerance = 1e-12)
  expect_true(is.na(test_SM(0.2, NA, ref = 0.3)$flag))
})
