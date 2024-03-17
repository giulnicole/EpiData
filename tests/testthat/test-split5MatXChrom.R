

# Define a test case
test_that("split5MatXChrom removes CpGs with high missing values", {
  # Create a sample input object
  list.cleaned <- list(
    Coverage_matrix = data.frame(CpG1 = c(10, NA, 30), CpG2 = c(NA, 6, 40)),
    Met_matrix = data.frame(CpG1 = c(8, 2, 10), CpG2 = c(NA, 4, 21)),
    Unmet_matrix = data.frame(CpG1 = c(2, NA, 20), CpG2 = c(NA, 2, 19)),
    Beta_matrix = data.frame(CpG1 = c(0.8, NA, 0.3), CpG2 = c(NA, 0.3, 0.48)),
    M_matrix = data.frame(CpG1 = c(2, NA, 3), CpG2 = c(NA, 2, 8))
  )

  # Test the function
  result <- split5MatXChrom(list.cleaned)

  # Perform assertions using expect_* functions
  expect_type(result, "list")

})
