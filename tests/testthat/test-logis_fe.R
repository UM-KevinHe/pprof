test_that("logis_fe function behaves correctly", {
  data(data_FE)
  Y = data_FE$Y
  Z = data_FE$Z
  ID = data_FE$ID

  fit <- logis_fe(Y,Z,ID)

  # Check if "fit" has the correct class "logis_fe"
  expect_true("logis_fe" == class(fit), "'fit' should be an object of class 'logis_fe'.")

  # Check if "fit$data_included" includes expected components
  expect_true(all(c("included", "no.events", "all.events") %in% names(fit$data_include)),
              "'fit$data_included' should contain 'included', 'no.events' and 'all.events' columns")

  # Check if "fit" contains the correct number of components
  expect_length(fit, 12)
})
