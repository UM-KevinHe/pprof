test_that("SR_output funtion returns correct structure", {
  data(data_FE)
  Y = data_FE$Y
  Z = data_FE$Z
  ID = data_FE$ID

  fit <- logis_fe(Y,Z,ID)

  sr <- SR_output(fit)

  # Check if "sr" includes expected components
  expect_true(all(c("indirect.ratio", "indirect.rate", "OE") %in% names(sr)),
              "The output from SR_output should include indirect.ratio, indirect.rate, and OE components.")
})
