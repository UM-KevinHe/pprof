test_that("test.linear_fe function behaves correctly", {
  data(ExampleDataLinear)
  Y <- ExampleDataLinear$Y
  Z <- ExampleDataLinear$Z
  ID <- ExampleDataLinear$ID
  data <- data.frame(Y, ID, Z)
  Z.char <- colnames(Z)
  Y.char <- "Y"
  ID.char <- "ID"

  fit_fe <- linear_fe(data = data, Y.char = Y.char, Z.char = Z.char, ID.char = ID.char)

  test_all <- test(fit_fe)
  test_parm <- test(fit_fe, parm = c(1, 3, 9))

  expect_true(all(test_parm == test_all[c(1,3,9),]),
              info = "Argument 'parm' behaves correctly.")
})
