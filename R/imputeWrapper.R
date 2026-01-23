#' @title imputeWrapper
#'
#' @description Apply imputation method(s) to a BSseq object and return it with imputed M-values.
#'
#' @param bs `BSseq` object containing methylation data with NAs.
#' @param methods Vector of method numbers to apply:
#'   - 1: Mean imputation (BSImpute with mean)
#'   - 2: KNN imputation (BSImpute with knn)
#'   - 3: KNN imputation (knnBSImpute)
#'   Can be a single value (e.g., 2) or multiple (e.g., c(1, 2)).
#'   Default is 2.
#' @param dist_threshold Maximum genomic distance (bp) between CpGs in same cluster. Default 1000.
#'
#' @return BSseq object with imputed M-values stored in assay slots:
#'   - "M_imputed_mean" (if method 1 is used)
#'   - "M_imputed_knn_BSImpute" (if method 2 is used)
#'   - "M_imputed_knn_knnBSImpute" (if method 3 is used)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(bsseq)
#' data(BS.chr22)
#' bs <- BS.chr22
#' bs_filt <- filterCov(bs = bs,
#'                      min_cov = 10,
#'                      row_min = 0.5,
#'                      col_min = 0.5)
#' # Apply single method (KNN via BSImpute)
#' bs_result <- imputeWrapper(bs, methods = 2)
#'
#' # Apply multiple methods
#' bs_result <- imputeWrapper(bs, methods = c(1, 2, 3))
#'
#' # Extract imputed data
#' m_knn <- SummarizedExperiment::assay(bs_result, "M_imputed_knn_BSImpute")
#' m_mean <- SummarizedExperiment::assay(bs_result, "M_imputed_mean")
#' }
imputeWrapper <- function(bs, methods = 2, dist_threshold = 1000) {
  #TODO: you don't need the SummarizedExperiment any more, please check the new BSImpute function
  #TODO: this function needs more work, let's discuss over a short meeting
  #TODO: remove "\dontrun" once this function is fixed!
  # Validate methods
  if (!all(methods %in% c(1, 2, 3))) {
    stop("methods must contain only: 1 (mean), 2 (KNN via BSImpute), or 3 (knnBSImpute)")
  }

  # Method names for assay slots
  assay_names <- c(
    "M_imputed_mean",
    "M_imputed_knn_BSImpute",
    "M_imputed_knn_knnBSImpute"
  )

  method_labels <- c(
    "Mean (BSImpute)",
    "KNN (BSImpute)",
    "KNN (knnBSImpute)"
  )

  # Apply each requested method
  for (method in methods) {

    cat("Applying method", method, ":", method_labels[method], "...\n")

    if (method == 1) {
      # Method 1: Mean imputation via BSImpute
      bs_temp <- BSImpute(
        bs = bs,
        dist_threshold = dist_threshold,
        impute_method = "mean"
      )
      # Extract imputed M-values and store with descriptive name
      m_imputed <- SummarizedExperiment::assay(bs_temp, "M_values_imputed")
      SummarizedExperiment::assays(bs, withDimnames = FALSE)[[assay_names[1]]] <- m_imputed

    } else if (method == 2) {
      # Method 2: KNN imputation via BSImpute
      bs_temp <- BSImpute(
        bs = bs,
        dist_threshold = dist_threshold,
        impute_method = "knn"
      )
      # Extract imputed M-values and store with descriptive name
      m_imputed <- SummarizedExperiment::assay(bs_temp, "M_values_imputed")
      SummarizedExperiment::assays(bs, withDimnames = FALSE)[[assay_names[2]]] <- m_imputed

    } else if (method == 3) {
      # Method 3: KNN imputation via knnBSImpute
      bs_temp <- knnBSImpute(
        bs = bs,
        dist_threshold = dist_threshold
      )
      # Extract imputed M-values and store with descriptive name
      m_imputed <- SummarizedExperiment::assay(bs_temp, "M_values_imputed")
      SummarizedExperiment::assays(bs, withDimnames = FALSE)[[assay_names[3]]] <- m_imputed
    }
  }

  cat("\nImputation complete!\n")
  cat("Available assays in BSseq object:\n")
  print(names(SummarizedExperiment::assays(bs)))

  return(bs)
}

