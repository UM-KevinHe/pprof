test_that("test.linear_re function behaves correctly", {
  data(ExampleDataLinear)
  Y <- ExampleDataLinear$Y
  Z <- ExampleDataLinear$Z
  ProvID <- ExampleDataLinear$ProvID

  fit_re <- linear_re(Y = Y, Z = Z, ProvID = ProvID)

  test_all <- test(fit_re)
  test_parm <- test(fit_re, parm = c(1, 3:5, 9))

  expect_true(all(test_parm == test_all[c(1,3:5,9),]),
              info = "Argument 'parm' behaves correctly.")
})
