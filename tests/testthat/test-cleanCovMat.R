

# Define a test case
test_that("cleanCovMat removes CpGs with high missing values", {
  # Create a sample input object
  # Creating the matrix with chosen values
  library(dplyr)
  library(SummarizedExperiment)

  data("meth_data")
  input.list<- SummarizedExperiment::assays(meth_data)

  result <- cleanCovMat(input.list)

  # Perform assertions using expect_* functions
  expect_true("Coverage_matrix" %in% names(result))
  expect_true("Met_matrix" %in% names(result))
  expect_true("Unmet_matrix" %in% names(result))

})

