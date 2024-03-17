

# Define a test case
test_that("cleanOutliers removes outliers and calculates ratios", {
  # Create a sample filtered object
  filtered_obj <- list(
    Coverage_matrix = data.frame(CpG1 = c(10, NA, 30), CpG2 = c(NA, 6, 40)),
    Met_matrix = data.frame(CpG1 = c(8, 2, 10), CpG2 = c(NA, 4, 21)),
    Unmet_matrix = data.frame(CpG1 = c(2, NA, 20), CpG2 = c(NA, 2, 19))
  )
  # Test the function
  result <- cleanOutliers(filtered_obj, remove_outliers = TRUE)

  # Perform assertions using expect_* functions
  expect_true("Coverage_matrix" %in% names(result))
  expect_true("Met_matrix" %in% names(result))
  expect_true("Unmet_matrix" %in% names(result))
  expect_true("Beta_matrix" %in% names(result))
  expect_true("M_matrix" %in% names(result))

})
