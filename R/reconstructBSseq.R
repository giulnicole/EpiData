#' @title reconstructBSseq
#' @description
#' Back-calculates counts and wraps them into a BSseq object
#'
#' @param val_matrix A numeric matrix of M-values or Beta values.
#' @param cov_matrix A numeric matrix of total coverage.
#' @param gr A \code{GRanges} object representing the CpG sites.
#' @param type Character, either "m" or "beta".
#' @param log2_offset The numeric offset used for M-values (default 2).
#' @param sampleNames Optional vector of sample names.
#'
#' @importFrom bsseq BSseq
#'
#' @examples
#' library(bsseq)
#' data(BS.chr22)
#' bs <- BS.chr22
#'
#' dat <- extractMethData(bs)
#' m_value <- dat$m
#' coverage <- bsseq::getCoverage(bs, type = "Cov")
#' gr <- bsseq::granges(bs)
#'
#' bs_recons <- reconstructBSseq(val_matrix = m_value,
#'                               cov_matrix = coverage,
#'                               gr = gr,
#'                               type = c("m"),
#'                               log2_offset = 2,
#'                               sampleNames = NULL)
#'
#' @export
reconstructBSseq <- function(val_matrix,
                             cov_matrix,
                             gr,
                             type = c("m", "beta"),
                             log2_offset = 2,
                             sampleNames = NULL) {

  type <- match.arg(type)

  # 1. Math Logic
  if (type == "m") {
    R <- 2^val_matrix
    meth <- (R * (cov_matrix + log2_offset) - log2_offset) / (R + 1)
  } else {
    meth <- val_matrix * cov_matrix
  }

  # 2. Cleanup: Handle NAs and Boundaries
  meth[is.na(meth)] <- 0
  meth[meth < 0] <- 0
  meth[meth > cov_matrix] <- cov_matrix[meth > cov_matrix]

  # 3. Integerize properly
  meth <- round(meth)
  storage.mode(meth) <- "integer"

  cov_matrix <- round(cov_matrix)
  storage.mode(cov_matrix) <- "integer"

  # 4. CRITICAL: Synchronize Names to avoid the SummarizedExperiment Error
  # If sampleNames wasn't provided, try to grab from the matrix
  if (is.null(sampleNames)) {
    sampleNames <- colnames(val_matrix)
    if (is.null(sampleNames)) sampleNames <- paste0("Sample_", seq_len(ncol(meth)))
  }

  # Match names exactly
  colnames(meth) <- colnames(cov_matrix) <- sampleNames

  # Usually, BSseq prefers NULL rownames or names that match the GRanges exactly.
  # We'll set them to NULL to let the GRanges handle the indexing.
  rownames(meth) <- rownames(cov_matrix) <- NULL

  # 5. Construct BSseq object
  bs <- bsseq::BSseq(
    M = meth,
    Cov = cov_matrix,
    gr = gr,
    sampleNames = sampleNames
  )

  return(bs)
}
