

# Define a test case
test_that("statsCpG calculates statistics and correlations", {

 library(SummarizedExperiment)
 data("rangedObject")
 input.list<- assays(rangedObject)
 clean.coverage <- cleanCovMat(input.obj = input.list,
                          max_na_cpg = 0.5, max_na_ind = 0.2,
                            cpg_removal_threshold = 10)

# Cleaning outliers in methylated and unmethylated counts
 clean.out <- cleanOutliers(filtered.obj=clean.coverage,
                              outlier_threshold=5,
                               remove_outliers = TRUE)
  # Test the function
  result <- statsCpG(clean.out, varselect = 5, plot = FALSE)

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
