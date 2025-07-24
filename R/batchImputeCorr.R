#' @title batchImputeCorr
#'
#' @description #' This function applies the \code{\link{imputeCorr}} function to a list of chromosome-specific objects previously processed with \code{\link{statsCpG}}.
#' Each element in the input list should represent a chromosome and must contain the output from \code{statsCpG()}, which includes statistics
#' such as the number of individuals, CpGs, means, standard deviations, and the CpG matrix.
#'
#'
#' @importFrom matrixcalc is.positive.semi.definite
#' @import stats
#'
#' @name batchImputeCorr
#'
#' @param stats_list A named list where each element is the output of statsCpG() per chromosome
#' @param matrix Type of matrix for imputeCorr (default: "M")
#' @param varselect Index for statsCpG selection, passed to imputeCorr (default: 5, the M values matrixs)
#'
#' @return A named list where each element is the output of imputeCorr()
#'
#'
#' @examples
#' \dontrun{
#' # Data
#'  data("meth_data")
#'  print(meth_data)
#'
#'  input.list<- assays(meth_data)
#'
#'  # Cleaning low counts for coverage
#'
#'  clean.coverage <- cleanCovMat(input.obj = input.list,
#'                               max_na_cpg = 0.5, max_na_ind = 0.2,
#'                              cpg_removal_threshold = 10)
#'
#' # Cleaning outliers in methylated and unmethylated counts
#'  clean.out <- cleanOutliers(filtered.obj=clean.coverage,
#'                             outlier_threshold=5, remove_outliers = TRUE)
#'
#'  stat.result <- missingStats(list(mat.values))
#'  imputed <- batchImputeCorr(stat.result)
#'
#'  }
#'
#' @export
batchImputeCorr <- function(stats_list,
                            matrix = "M",
                            varselect = 5) {

  # Check input is a list
  if (!is.list(stats_list)) {
    stop("Input must be a list of statsCpG outputs (one per chromosome).")
  }

  # Apply imputeCorr to each element of the list
  imputed_results <- lapply(stats_list, function(chr_stats) {
    imputeCorr(cleaned.obj = chr_stats, matrix = matrix, varselect = varselect)
  })

  return(imputed_results)
}

