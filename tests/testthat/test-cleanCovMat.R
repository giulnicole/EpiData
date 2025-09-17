

# Define a test case
test_that("cleanCovMat removes CpGs with high missing values", {
  # Create a sample input object
  # Creating the matrix with chosen values
  data("meth_data")
  input.list<- assays(meth_data)


  # Test the function
  library(dplyr)
  library(SummarizedExperiment)
  result <- cleanCovMat(input.list)

  # Perform assertions using expect_* functions
  expect_true("Coverage_matrix" %in% names(result))
  expect_true("Met_matrix" %in% names(result))
  expect_true("Unmet_matrix" %in% names(result))

})

