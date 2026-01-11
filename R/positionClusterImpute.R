#' @title positionClusterImpute
#'
#' @description Cluster CpGs by genomic proximity and impute missing M-values.
#'
#' @import impute
#'
#' @param m_values Data frame or matrix of M-values (rows = CpGs, columns = samples).
#' @param positions Data frame with columns: CpG, chr, pos.
#' @param dist_threshold Maximum genomic distance (bp) between CpGs in same cluster. Default 1000.
#' @param impute_method "mean" or "knn" (requires `impute` package).
#'
#' @return Imputed M-value matrix (rows = CpGs, columns = samples).
#'
#' @export
positionClusterImpute <- function(m_values, positions,
                                  dist_threshold = 1000,
                                  impute_method = "mean") {

  if (!all(c("CpG", "chr", "pos") %in% colnames(positions))) {
    stop("positions must have columns: CpG, chr, pos")
  }

  if (is.null(rownames(m_values))) stop("m_values must have CpG IDs as rownames")

  # keep only CpGs in both
  positions <- positions[positions$CpG %in% rownames(m_values), , drop = FALSE]
  positions <- positions[order(positions$chr, positions$pos), , drop = FALSE]
  m_values <- as.matrix(m_values[positions$CpG, , drop = FALSE])

  cluster_ids <- integer(nrow(positions))
  cluster_num <- 1
  uniq_chr <- unique(positions$chr)

  # cluster CpGs per chromosome
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
}
