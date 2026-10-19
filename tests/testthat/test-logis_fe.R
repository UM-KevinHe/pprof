test_that("fe_data_prep sorts by provider and flags small providers", {
  d <- sim_fe()
  d$ProvID$prov[1:5] <- "SMALL"
  dp <- prep_fe(d)
  expect_s3_class(dp, "data_prep")
  expect_false(is.unsorted(as.character(dp$data$prov)))
  expect_true(all(dp$data$included[dp$data$prov == "SMALL"] == 0))
  expect_equal(dp$char_list$N.char, "N")
  expect_equal(unique(dp$data$N), 1)
})

test_that("logis_fe reproduces glm with provider dummies (binary)", {
  d <- sim_fe()
  fit <- logis_fe(prep_fe(d), tol = 1e-12)
  ref <- glm_fe(d)
  expect_equal(as.numeric(fit$beta), unname(coef(ref)[c("x1", "x2", "x3")]), tolerance = 1e-7)
  expect_equal(as.numeric(fit$gamma), unname(coef(ref)[grep("^prov", names(coef(ref)))]), tolerance = 1e-7)
  expect_equal(fit$neg2Loglkd, ref$deviance, tolerance = 1e-6) # binary: deviance = -2 loglik
})

test_that("logis_fe reproduces glm for a binomial outcome (N > 1)", {
  d <- sim_fe(binomial = TRUE, seed = 2)
  fit <- logis_fe(prep_fe(d, binomial = TRUE), tol = 1e-12)
  ref <- glm_fe(d)
  expect_equal(as.numeric(fit$beta), unname(coef(ref)[c("x1", "x2", "x3")]), tolerance = 1e-7)
  expect_equal(as.numeric(fit$gamma), unname(coef(ref)[grep("^prov", names(coef(ref)))]), tolerance = 1e-7)
})

test_that("SerBIN and BAN agree, and every stopping rule converges", {
  dp <- prep_fe(sim_fe(seed = 3))
  f1 <- logis_fe(dp, tol = 1e-12)
  f2 <- logis_fe(dp, algorithm = "BAN", tol = 1e-12)
  expect_equal(as.numeric(f1$beta), as.numeric(f2$beta), tolerance = 1e-6)
  for (st in c("relch", "ratch", "all", "or")) {
    f <- logis_fe(dp, stop = st, tol = 1e-12)
    expect_equal(as.numeric(f$beta), as.numeric(f1$beta), tolerance = 1e-5)
  }
  expect_error(logis_fe(dp, stop = "bad"), "stop")
})

test_that("results do not depend on the number of threads", {
  dp <- prep_fe(sim_fe(seed = 4))
  expect_identical(logis_fe(dp, threads = 1)$beta, logis_fe(dp, threads = 1)$beta)
  expect_equal(logis_fe(dp, threads = 1)$beta, logis_fe(dp, threads = 3)$beta, tolerance = 1e-12)
})

test_that("coef and vcov methods", {
  d <- sim_fe(seed = 5)
  fit <- logis_fe(prep_fe(d), tol = 1e-12)
  expect_named(coef(fit), c("gamma", "beta"))
  ref <- glm_fe(d)
  v <- vcov(fit)
  expect_equal(unname(sqrt(diag(v))), unname(sqrt(diag(vcov(ref)))[c("x1", "x2", "x3")]), tolerance = 1e-5)
  vg <- attr(v, "var.gamma")
  expect_equal(unname(vg), unname(diag(vcov(ref))[grep("^prov", names(coef(ref)))]), tolerance = 1e-5)
})

test_that("AUC is reported for binary outcomes only", {
  fit <- logis_fe(prep_fe(sim_fe(seed = 6)), AUC = TRUE)
  expect_true(fit$AUC > 0.5 && fit$AUC < 1)
  expect_warning(logis_fe(prep_fe(sim_fe(binomial = TRUE, seed = 6), binomial = TRUE), AUC = TRUE), "binary")
})
