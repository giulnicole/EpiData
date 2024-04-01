#' @title ParallSimu
#'
#' @description A function for computing the non-negative definite matrix of covariance matrix, by eigenvectors' decompisition and
#' and simulate new values drawn from the same distribution to impute remaining missing values in the original cleaned dataset (Parallelized version)
#' \code{} computes the imputation based on correlation pattern between CpGs
#'
#' @import matrixcalc
#'
#' @param stats.obj list of SummarizedExperimet object cleaned and after statistics and NAs' pattern computation (ParallStats output)
#' @param correlation.type Type of correlation to be imputed (default is with meanCpG imputed per rows, otherwise pairwise correlation).
#'
#'
#' @name ParallSimu
#'
#' @return
#'  List with SummarizedExperiment objects and imputed datasets in metadata:
#' \item{imputed}{Matrix of CpG imputed}
#'
#' @examples
#' \dontrun{
#'
#'  data("rangedObject")
#'
#'  # Cleaning low counts for coverage
#'  input.list <- assays(data)
#'  clean.coverage <- cleanCovMat(input.obj = input.list,
#'                               max_na_cpg = 0.5, max_na_ind = 0.2,
#'                              cpg_removal_threshold = 10)
#'
#' # Cleaning outliers in methylated and unmethylated counts
#'  clean.out <- cleanOutliers(filtered.obj=clean.coverage,
#'                             outlier_threshold=5, remove_outliers = TRUE)ù
#'
#'  list.cleaned1<- clean.out
#'  list.cleaned2<- clean.out
#'  list.cleaned<- list(list.cleaned1$Output_outliers@assays@data@listData,
#'                      list.cleaned2$Output_outliers@assays@data@listData)
#'  stats <- ParallStats(list.cleaned)
#'  simu.mix <- ParallSimu(list.stats = stats)
#' }
#'
#' @export
ParallSimu <- function(list.stats, correlation.type="meanCpG"){


  meta.simu <- bplapply(list.stats, imputeCorr, correlation.type)
  simulated.mat <- extractFinalMat(meta.simu)

  return(simulated.mat)

}# parallSimu



