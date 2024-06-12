test_that("logis_fe function behaves correctly", {
  data(data_FE)
  Y = data_FE$Y
  Z = data_FE$Z
  ID = data_FE$ID

  fit <- logis_fe(Y,Z,ID)

  # Check if "fit" has the correct class "logis_fe"
  expect_true("logis_fe" == class(fit), "'fit' should be an object of class 'logis_fe'.")

  # Check if "fit" includes expected components
  expect_true(all(c("beta", "gamma", "linear_pred", "data_include") %in% names(fit)),
              "'fit' from logis_fe should include 'beta', 'gamma', 'linear_pred', and 'data_include' components.")
})
