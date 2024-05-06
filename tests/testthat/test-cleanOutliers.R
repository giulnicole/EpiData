

# Define a test case


test_that("cleanCovMat removes CpGs with high missing values", {
  # Create a sample input object
  # Creating the matrix with chosen values
  data("rangedObject")
  input.list<- assays(rangedObject)

  # Test the function
  library(tidyverse)
  result0 <- cleanCovMat(input.list)
  result <- cleanOutliers(result0, remove_outliers = FALSE)

  # Perform assertions using expect_* functions
  expect_true("coverage" %in% names(result))
  expect_true("methylated" %in% names(result))
  expect_true("unmethylated" %in% names(result))
  expect_true("beta" %in% names(result))
  expect_true("m" %in% names(result))

})
