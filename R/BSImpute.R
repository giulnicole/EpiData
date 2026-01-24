#' @title BSImpute
#'
#' @description Cluster CpGs by genomic proximity and impute missing M-values.
#'
#' @import impute
#' @import SummarizedExperiment
#'
#' @param bs `BSseq` object.
#' @param dist_threshold Maximum genomic distance (bp) between CpGs in same cluster. Default 1000.
#' @param impute_method "mean" or "knn" (requires `impute` package).
#' @param min_cpgs Minimum number of CpGs per cluster. Default 10.
#'
#' @return BSseq object with imputed M-values.
#'
#' @export
BSImpute <- function(bs,
                      dist_threshold = 1000,
                      impute_method = "mean",
                      min_cpgs = 10) {

  bs_mat <- extractMethData(bs)
  m_mat <- bs_mat[["m"]]
  gr <- bs_mat$gr

  positions <- data.frame(
    CpG = rownames(m_mat),
    chr = GenomicRanges::seqnames(gr),
    pos = GenomicRanges::start(gr),
    stringsAsFactors = FALSE
  )

  m_values <- as.matrix(m_mat[positions$CpG, , drop = FALSE])

  ## Initial clustering by distance
  cluster_ids <- integer(nrow(positions))
  cluster_num <- 1

  for (chr in unique(positions$chr)) {
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

  ## Introducing minimum CpGs per cluster
  for (chr in unique(positions$chr)) {
    chr_idx <- which(positions$chr == chr)
    chr_clusters <- unique(positions$cluster[chr_idx])

    i <- 1
    while (i < length(chr_clusters)) {
      cl <- chr_clusters[i]
      next_cl <- chr_clusters[i + 1]

      cl_idx <- chr_idx[positions$cluster[chr_idx] == cl]

      if (length(cl_idx) < min_cpgs) {
        positions$cluster[cl_idx] <- next_cl
        chr_clusters <- unique(positions$cluster[chr_idx])
      } else {
        i <- i + 1
      }
    }
  }

  ## Imputation
  m_imp <- m_values

  for (cl in unique(positions$cluster)) {
    idx <- which(positions$cluster == cl)
    block <- m_values[idx, , drop = FALSE]
    if (!any(is.na(block))) next

    if (impute_method == "knn" && requireNamespace("impute", quietly = TRUE)) {
      block <- impute::impute.knn(block, colmax = 1)$data
    } else {
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
  SummarizedExperiment::assays(bs, withDimnames = FALSE)$M_values_imputed <- m_imp

  return(bs)
}
