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
#'
#' @return Imputed M-value matrix (rows = CpGs, columns = samples).
#'
#' @export
BSImpute <- function(bs, dist_threshold = 1000, impute_method = "mean") {

  bs_mat <- extractMethData(bs)

  # Methylation values
  #TODO: The user should have the option to select between m-value or beta-value below
  m_mat <- bs_mat[["m"]]

  # gr object
  gr <- bs_mat$gr

  # keep only CpGs in both
  positions <- data.frame(
    CpG = rownames(m_mat),
    chr = gr@seqnames,
    pos = gr@ranges@start,
    stringsAsFactors = FALSE)

  # m-values
  m_values <- as.matrix(m_mat[positions$CpG, , drop = FALSE])

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

  # Run imputations (for m-value)
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

  # Run imputations (for coverage)
  #TODO: We need the imputed coverage matrix as well. Since this is a repeatation with above, we need to build a function later
  cov_values <- bsseq::getCoverage(bs, type = "Cov")
  cov_imp <- cov_values

  for (cl in unique(positions$cluster)) {
    idx <- which(positions$cluster == cl)
    block <- cov_values[idx, , drop = FALSE]
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
    cov_imp[idx, ] <- block
  }

  rownames(cov_imp) <- positions$CpG

  # reconstruct BS
  bs_imp <- reconstructBSseq(val_matrix = m_imp,
                             cov_matrix = cov_imp,
                             gr = bsseq::granges(bs),
                             type = c("m", "beta")[1],
                             log2_offset = 2,
                             sampleNames = NULL)

  return(bs_imp)
}
