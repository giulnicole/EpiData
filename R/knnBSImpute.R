#' @title knnBSImpute
#'
#' @description Cluster CpGs by genomic proximity and impute missing M-values using KNN.
#'
#' @import impute
#' @import SummarizedExperiment
#'
#' @param bs `BSseq` object.
#' @param dist_threshold Maximum genomic distance (bp) between CpGs in same cluster. Default 1000.
#'
#' @return BSseq object with imputed M-values in assay slot "M_values_imputed".
#'
#' @export
knnBSImpute <- function(bs, dist_threshold = 1000) {
  #TODO: please check the function BSImpute to see how to generate the imputed BSseq object and revise this function (if necessary) accordingly
  # Check if impute package is available
  if (!requireNamespace("impute", quietly = TRUE)) {
    stop("Package 'impute' is required. Install it with: BiocManager::install('impute')")
  }

  # Extract methylation data
  bs_mat <- extractMethData(bs)
  m_mat <- bs_mat[["m"]]
  gr <- bs_mat$gr

  # Create positions dataframe
  positions <- data.frame(
    CpG = rownames(m_mat),
    chr = GenomicRanges::seqnames(gr),
    pos = GenomicRanges::start(gr),
    stringsAsFactors = FALSE
  )

  # Get M-values matrix
  m_values <- as.matrix(m_mat[positions$CpG, , drop = FALSE])

  # Cluster CpGs by chromosome and distance
  cluster_ids <- integer(nrow(positions))
  cluster_num <- 1
  uniq_chr <- unique(positions$chr)

  for (chr in uniq_chr) {
    idx <- which(positions$chr == chr)
    if (length(idx) == 0) next
    cluster_ids[idx[1]] <- cluster_num
    for (i in 2:length(idx)) {
      dist_bp <- positions$pos[idx[i]] - positions$pos[idx[i - 1]]
      if (dist_bp > dist_threshold) cluster_num <- cluster_num + 1
      cluster_ids[idx[i]] <- cluster_num
    }
    cluster_num <- cluster_num + 1
  }

  positions$cluster <- cluster_ids
  m_imp <- m_values

  # KNN imputation per cluster
  for (cl in unique(positions$cluster)) {
    idx <- which(positions$cluster == cl)
    block <- m_values[idx, , drop = FALSE]
    if (any(is.na(block))) {
      m_imp[idx, ] <- impute::impute.knn(block, colmax = 1)$data
    }
  }

  # Store imputed M-values in BSseq object
  rownames(m_imp) <- positions$CpG
  SummarizedExperiment::assays(bs, withDimnames = FALSE)$M_values_imputed <- m_imp

  return(bs)
}
