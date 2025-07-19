#' @title eigenImpute
#'
#' @description A function for computing the non-negative definite matrix of covariance matrix, by eigenvectors' decompisition and
#' and simulate new values drawn from the same distribution to impute remaining missing values in the original cleaned dataset (Parallelized version)
#' \code{} computes the imputation based on correlation pattern between CpGs
#'
#' @import matrixcalc
#' @importFrom BiocParallel bplapply
#'
#' @param list.stats list of SummarizedExperimet object cleaned and after statistics and NAs' pattern computation (ParallStats output)
#'
#'
#' @name eigenImpute
#'
#' @return
#'  List with SummarizedExperiment objects and imputed datasets in metadata:
#' \item{imputed}{Matrix of CpG imputed}
#'
#' @examples
#' \dontrun{
#'
#'  data("meth_data")
#'  print(meth_data)
#'
#'  # Cleaning low counts for coverage
#'  input.list <- assays(meth_data)
#'  clean.coverage <- cleanCovMat(input.obj = input.list,
#'                               max_na_cpg = 0.5, max_na_ind = 0.2,
#'                               cpg_removal_threshold = 10)
#'
#' # Cleaning outliers in methylated and unmethylated counts
#'  clean.out <- cleanOutliers(filtered.obj=clean.coverage,
#'                             outlier_threshold=5, remove_outliers = TRUE)ù
#'
#'  list.cleaned<- list(clean.out, clean.out)
#'
#'  stats <- missingStats(list.cleaned)
#'  simu.mix <- eigenImpute(list.stats = stats)
#' }
#'
#' @export
eigenImpute <- function(list.stats){


  meta.simu <- bplapply(list.stats, imputeCorr)
  simulated.mat <- extractFinalMat(meta.simu)

  return(simulated.mat)

}# eigenImpute



