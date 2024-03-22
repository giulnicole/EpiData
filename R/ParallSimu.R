#' @title ParallSimu
#'
#' @description A function for computing the non-negative definite matrix of covariance matrix, by eigenvectors' decompisition and
#' and simulate new values drawn from the same distribution to impute remaining missing values in the original cleaned dataset (Parallelized version)
#' \code{} computes the imputation based on correlation pattern between CpGs
#'
#' @import matrixcalc
#'
#' @param stats.obj list of SummarizedExperimet object cleaned and after statistics and NAs' pattern computation (ParallStats output)
#' @param matrix matrix on which computing the correlation matrix; "M" as default
#' @param varselect index of the dataset to be used (numeric value 1-5)
#'
#' @name ParallSimu
#'
#' @return
#'  List with SummarizedExperiment objects and imputed datasets in metadata:
#' \item{imputed}{Matrix of CpG imputed}
#'
#' @examples
#' \dontrun{
#'  data('clean.out')
#'  list.cleaned1<- clean.out
#'  list.cleaned2<- clean.out
#'  list.cleaned<- list(list.cleaned1$Output_outliers@assays@data@listData, list.cleaned2$Output_outliers@assays@data@listData)
#'  stats2 <- ParallStats(list.cleaned)
#'  simu.mix2 <- ParallSimu(list.stats = stats2)
#' }
#'
#' @export
ParallSimu <- function(list.stats, correlation.type="meanCpG"){


  meta.simu <- bplapply(list.stats, imputeCorr, correlation.type)
  simulated.mat <- extractFinalMat(meta.simu)

  return(simulated.mat)

}# parallSimu



