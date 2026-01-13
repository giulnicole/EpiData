#' @title positionClusterImpute
#'
#' @description Cluster CpGs by genomic proximity and impute missing M-values.
#'
#' @import impute
#'
#' @param bs `BSseq` object.
#' @param dist_threshold Maximum genomic distance (bp) between CpGs in same cluster. Default 1000.
#' @param impute_method "mean" or "knn" (requires `impute` package).
#'
#' @return Imputed M-value matrix (rows = CpGs, columns = samples).
#'
#' @export
positionClusterImpute <- function(bs,
                                  dist_threshold = 1000,
                                  impute_method = "mean") {

  # Check if bs is processed
  if(is.null(rownames(bs))){
    bs <- processBSseq(bs)
  }

  # M-values
  m_values <- bsseq::getCoverage(bs, type = "M")

  # Calculate positions
  positions <- as.data.frame(bsseq::granges(bs))
  positions <- positions[,c("seqnames","start")]
  rownames(positions) <- rownames(m_values)
  colnames(positions) <- c("chr", "pos")

  # Clusters
  cluster_ids <- integer(nrow(positions))
  cluster_num <- 1
  uniq_chr <- unique(positions$chr)

  # Cluster CpGs per chromosome
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

  # Run imputations
  m_imp <- m_values

  for (cl in unique(positions$cluster)) {
    idx <- which(positions$cluster == cl)
    block <- m_values[idx, , drop = FALSE]
    if (!any(is.na(block))) next

    if (impute_method == "knn" && requireNamespace("impute", quietly = TRUE)) {
      block <- impute::impute.knn(block, colmax = 1)$data
    } else {
      # mean-based imputation
      for (r in seq_len(nrow(block))) {
        rowv <- block[r, ]
        if (all(is.na(rowv))) {
          block[r, ] <- colMeans(block, na.rm = TRUE)
        } else {
          nas <- is.na(rowv)
          if (any(nas)) block[r, nas] <- mean(rowv, na.rm = TRUE)
        }
      }
    }
    m_imp[idx, ] <- block
  }

  rownames(m_imp) <- positions$CpG
  return(as.data.frame(m_imp))

  #TODO: Why don't we have a `BSseq` object to return? Then, the user can extract the imputed values by `extractMethData`
}
