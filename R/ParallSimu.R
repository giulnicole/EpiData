#' @title ParallSimu
#'
#' @description A function for computing the non-negative definite matrix of covariance matrix, by eigenvectors' decompisition and
#' and simulate new values drawn from the same distribution to impute remaining missing values in the original cleaned dataset (Parallelized version)
#' \code{\link{ParallSimu}} computes the imputation based on correlation pattern between CpGs
#'
#'
#' @param stats.obj list of SummarizedExperimet object cleaned and after statistics and NAs' pattern computation (ParallStats output)
#' @param matrix matrix on which computing the correlation matrix; "M" as default
#' @param varselect index of the dataset to be used (numeric value 1-5)
#'
#'
#'
#' @name ParallSimu
#'
#' @return
#'  list with SummarizedExperiment objects and imputed datasets in metadata:
#' \item{imputed}{Matrix of CpG imputed}
#'
#' @examples
#' \dontrun{
#'  data('clean.out')
#'  list.cleaned1<- clean.out
#'  list.cleaned2<- clean.out
#'  list.cleaned<- list(list.cleaned1$Output_outliers@assays@data@listData, list.cleaned2$Output_outliers@assays@data@listData)
#'  input.stats<- list()
#'  for (i in 1:length(list.cleaned)){
#'
#'   input.stats[[i]] <- GRconversion2(list.cleaned[[i]])
#'
#'   names(input.stats)[[i]] <- names(list.cleaned)[[i]]
#'   }
#'
#'  stats2 <- ParallStats(input.stats)
#'  simu.mix2 <- ParallSimu(list.stats = stats2)
#' }
#'
#' @export
ParallSimu <- function(list.stats){

  library(stats)
  meta.simu <- bplapply(list.stats, patternCpG)

  return(meta.simu)

}# parallSimu



