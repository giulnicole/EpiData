

# Define a test case
test_that("statsCpG calculates statistics and correlations", {

 library(SummarizedExperiment)
 data("meth_data")
 input.list<- assays(meth_data)
 clean.coverage <- cleanCovMat(input.obj = input.list,
                          max_na_cpg = 0.5, max_na_ind = 0.2,
                            cpg_removal_threshold = 10)

# Cleaning outliers in methylated and unmethylated counts
 clean.out <- cleanOutliers(filtered.obj=clean.coverage,
                              outlier_threshold=5,
                               remove_outliers = TRUE)


# Calculating beta and M values
 mat.values <- methValues(clean.coverage)

# Test the function
 stat.result <- missingStats(list(mat.values))

  # Perform assertions using expect_* functions
  expect_true(is.list(stat.result[[1]]))
  expect_true("Matrix" %in% names(stat.result[[1]]))
  expect_true("Individuals" %in% names(stat.result[[1]]))
  expect_true("CpGs" %in% names(stat.result[[1]]))
  expect_true("Table_missing" %in% names(stat.result[[1]]))
  expect_true("Fraction_missingnessp" %in% names(stat.result[[1]]))
  expect_true("Linear_correlation" %in% names(stat.result[[1]]))
  expect_true("long_correlation_matrix" %in% names(stat.result[[1]]))
  expect_true("means" %in% names(stat.result[[1]]))
  expect_true("sds" %in% names(stat.result[[1]]))

})
