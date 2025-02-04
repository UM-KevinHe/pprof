test_that("SM_output.linear_fe function behaves correctly", {
  data(ExampleDataLinear)
  Y <- ExampleDataLinear$Y
  Z <- ExampleDataLinear$Z
  ProvID <- ExampleDataLinear$ProvID

  fit_fe1 <- linear_fe(Y = Y, Z = Z, ProvID = ProvID)
  SM_median <- SM_output(fit_fe1, stdz = c("indirect", "direct"))
  dimnames(SM_median$indirect.difference) <- NULL
  dimnames(SM_median$direct.difference) <- NULL

  SM_mean <- SM_output(fit_fe1, parm = c(1,3,5), stdz = c("indirect", "direct"), null = "mean")

  expect_true(all.equal(SM_median$indirect.difference, SM_median$direct.difference),
              info = "ISD and DSD are equal.")

  expect_true(all(nrow(SM_mean$indirect.difference) == 3,
                  nrow(SM_mean$direct.difference) == 3),
              info = "Argument 'parm' performs correctly.")
})
