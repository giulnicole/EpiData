
# Define a test case
test_that("missingCpG calculates statistics and correlations with parallelization", {
  # Create a sample cleaned object
  library(SummarizedExperiment)
  library(BiocParallel)
  library(matrixcalc)

  data("rangedObject")
  input.list<- assays(rangedObject)
  clean.coverage <- cleanCovMat(input.obj = input.list,
                                max_na_cpg = 0.5, max_na_ind = 0.2,
                                cpg_removal_threshold = 10)

  # Cleaning outliers in methylated and unmethylated counts
  clean.out <- cleanOutliers(filtered.obj=clean.coverage,
                             outlier_threshold=5,
                             remove_outliers = TRUE)

  list.cleaned<- list(clean.out, clean.out)
  result <- missingStats(list.cleaned)


  # Perform assertions using expect_* functions
  expect_true(is.list(result))

})

