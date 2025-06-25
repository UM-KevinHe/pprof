test_that("test.linear_fe function behaves correctly", {
  data(ExampleDataLinear)
  Y <- ExampleDataLinear$Y
  Z <- ExampleDataLinear$Z
  ProvID <- ExampleDataLinear$ProvID
  data <- data.frame(Y, ProvID, Z)
  Z.char <- colnames(Z)
  Y.char <- "Y"
  ProvID.char <- "ProvID"

  fit_fe <- linear_fe(data = data, Y.char = Y.char, Z.char = Z.char, ProvID.char = ProvID.char)

  test_all <- test(fit_fe)
  test_parm <- test(fit_fe, parm = c(1, 3, 9))

  expect_true(all(test_parm == test_all[c(1,3,9),]),
              info = "Argument 'parm' behaves correctly.")
})
