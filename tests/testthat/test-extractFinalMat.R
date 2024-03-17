

# Define a test case
test_that("extractFinalMat removes CpGs with high missing values", {
  # Create a sample input object
  imputed_list <- list(
    list(matrix(1:6, nrow = 2)),
    list(matrix(7:12, nrow = 2)))


  # Test the function
  result <- extractFinalMat(imputed_list)

  # Perform assertions using expect_* functions
  expect_type(result, "list")

})
