

# Define a test case
test_that("statsCpG calculates statistics and correlations", {
  # Create a sample cleaned object
  cleaned_obj <- list(
      Coverage_matrix = data.frame(CpG1 = c(10, NA, 30), CpG2 = c(NA, 6, 40)),
      Met_matrix = data.frame(CpG1 = c(8, 2, 10), CpG2 = c(NA, 4, 21)),
      Unmet_matrix = data.frame(CpG1 = c(2, NA, 20), CpG2 = c(NA, 2, 19)),
      Beta_matrix = data.frame(CpG1 = c(0.8, NA, 0.3), CpG2 = c(NA, 0.3, 0.48)),
      M_matrix = data.frame(CpG1 = c(2, NA, 3), CpG2 = c(NA, 2, 8))
    )


  # Test the function
  result <- statsCpG(cleaned_obj, varselect = 5, plot = FALSE)

  # Perform assertions using expect_* functions
  expect_true(is.list(result))
  expect_true("Matrix" %in% names(result))
  expect_true("Individuals" %in% names(result))
  expect_true("CpGs" %in% names(result))
  expect_true("Table_missing" %in% names(result))
  expect_true("Fraction_missingnessp" %in% names(result))
  expect_true("Linear_correlation" %in% names(result))
  expect_true("long_correlation_matrix" %in% names(result))
  expect_true("means" %in% names(result))
  expect_true("sds" %in% names(result))

})
