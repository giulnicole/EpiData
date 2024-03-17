#' @title ParallStats
#'
#' @description A function for calculating statistics and missing values patterns
#'
#' \code{} computes the statistics on missing values per dataset, highlighting the pattern of missing values per each CpG (Parallelized version)
#'
#' @import kableExtra
#'
#' @param cleaned.obj SummarizedExperiment list object from filtering coverage (cleanCovMat) + cleaning outliers (cleanOutliers)
#' @param varselect is the index of the dataset to be used (numeric value 1-5) to extract the information related to the CpGs' pattern. 1 = coverage counts, 2 = methylated counts, 3 = unmethylated counts, 4 = beta values, 5 = M values.
#'
#' @name ParallStats
#'
#'
#' @return
#'  List with SummarizedExperiment objects and metadata statistics:
#'  \item{Rows}{Subject passing filters}
#'  \item{Columns}{Filtered CpGs}
#'  \item{Table_missing}{Table with NAs statistics per CpG with 5 columns: CpG's id, nObs, percObs, nNA, percNA}
#'  \item{Total_NA}{Number of NAs in the dataset}
#'  \item{NA_per_variable}{Number of NAs per variable}
#'  \item{Fraction_missingness}{Fraction of NAs in total}
#'  \item{K_table}{Table K table extra for table_missing}
#'  \item{Linear_correlation}{Matrix with liean correlation}
#'  \item{long_correlation_matrix}{Matrix  with linear correlation in longitudinal format}
#'
#' @examples
#' \dontrun{
#'  data('clean.out')
#'  list.cleaned1<- clean.out
#'  list.cleaned2<- clean.out
#'  list.cleaned<- list(list.cleaned1$Output_outliers@assays@data@listData, list.cleaned2$Output_outliers@assays@data@listData)
#'  stats2 <- ParallStats(list.cleaned)
#' }
#'
#'
#'
#'
#' @export
ParallStats <- function(list.cleaned){

  meta.stat <- bplapply(list.cleaned, statsCpG)

  return(meta.stat)

}# parallStats



