#' @title missingStats
#'
#' @description A function for calculating statistics and missing values patterns
#'
#' \code{} computes the statistics on missing values per dataset, highlighting the pattern of missing values per each CpG (Parallelized version)
#'
#'
#' @param list.cleaned SummarizedExperiment list object from filtering coverage (cleanCovMat) + cleaning outliers (cleanOutliers)
#'
#' @name missingStats
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
#'
#'  data("meth_data")
#'
#'  # Cleaning low counts for coverage
#'  input.list <- assays(meth_data)
#'
#'  clean.coverage <- cleanCovMat(input.obj = input.list,
#'                               max_na_cpg = 0.5, max_na_ind = 0.2,
#'                              cpg_removal_threshold = 10)
#'
#' # Cleaning outliers in methylated and unmethylated counts
#'   clean.out <- cleanOutliers(filtered.obj=clean.coverage,
#'                             outlier_threshold=5, remove_outliers = TRUE)
#'
#'   mat.values <- methValues(clean.coverage)
#'   stat.result <- missingStats(list(mat.values))
#'
#'
#'
#' }
#'
#'
#'
#'
#' @export
missingStats <- function(input.list) {
  if (!is.list(input.list)) {
    stop("Input must be a list of cleaned objects.")
  }

  results <- lapply(seq_along(input.list), function(i) {
    message("Processing element ", i, "...")
    obj <- input.list[[i]]

    if (!is.list(obj) || length(obj) < 5) {
      warning(sprintf("Element %d skipped: expected a list of 5 elements.", i))
      return(NULL)
    }

    tryCatch({
      statsCpG(obj, varselect = 5, plot = FALSE)
    }, error = function(e) {
      warning(sprintf("Element %d failed: %s", i, e$message))
      NULL
    })
  })

  # Preserve names if present
  if (!is.null(names(input.list))) {
    names(results) <- names(input.list)
  }

  return(results)
}
