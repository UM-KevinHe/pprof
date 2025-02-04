test_that("test.confint_fe function behaves correctly", {
  data(ExampleDataLinear)
  Y <- ExampleDataLinear$Y
  Z <- ExampleDataLinear$Z
  ProvID <- ExampleDataLinear$ProvID
  data <- data.frame(Y, ProvID, Z)
  Z.char <- colnames(Z)
  Y.char <- "Y"
  ProvID.char <- "ProvID"

  fit_fe <- linear_fe(data = data, Y.char = Y.char, Z.char = Z.char, ProvID.char = ProvID.char)

  CI <- confint(fit_fe, stdz = c("indirect", "direct"))
  CI_parm <- confint(fit_fe, parm = c(2, 4, 5, 7), stdz = c("indirect", "direct"))

  expect_true(all(all(CI_parm$CI.gamma == CI$CI.gamma[c(2, 4, 5, 7),]),
                  all(CI_parm$CI.indirect == CI$CI.indirect[c(2, 4, 5, 7),]),
                  all(CI_parm$CI.direct == CI$CI.direct[c(2, 4, 5, 7),])),
              info = "Argument 'parm' behaves correctly.")

  CI_one.side <- confint(fit_fe, stdz = c("indirect", "direct"), alternative = "greater")
  expect_true(all(all(CI_one.side$CI.gamma$gamma.Upper == Inf),
                  all(CI_one.side$CI.indirect$indirect.Upper == Inf)),
              info = "Argument 'alternative' behaves correctly.")
})
