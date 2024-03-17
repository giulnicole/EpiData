

# Define a test case
test_that("adjustCounts handles NAs in coverage, methylated, and unmethylated matrices", {

  # Create a sample input
    coverage = data.frame(CpG1 = c(10, NA, 30), CpG2 = c(NA, 6, 40))
    methylated = data.frame(CpG1 = c(8, 2, 10), CpG2 = c(NA, 4, 21))
    unmethylated = data.frame(CpG1 = c(2, NA, 20), CpG2 = c(NA, 2, 19))



  # Test the function
  result <- adjustCounts(coverage, methylated, unmethylated)


  expect_true(is.list(result))


})

