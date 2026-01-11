#' @title filterCov
#' @description
#' \code{} This function performs basic quality filtering on bisulfite sequencing data stored in a **BSseq** object.
#' CpG loci with insufficient read coverage across samples are removed, followed by optional filtering of CpGs and/or samples based on the proportion of missing or low-coverage values.
#' The intent is to retain only the most informative sites and samples for downstream analysis (e.g., methylation modelling or differential methylation testing),
#' while dropping loci with sparse data that may bias or reduce statistical power.
#'
#' @importFrom bsseq getCoverage
#'
#' @param bs A BSseq object.
#' @param min_cov Minimum coverage required for a CpG in a sample.
#' @param row_min The minimum percent missing data allowed in any row (CpG).
#' @param col_min The minimum percent missing data allowed in any column (Samples).
#'
#' @name filterCov
#'
#' @return A filtered BSseq object with CpGs and samples removed.
#'
#'
#' @examples
#' library(bsseq)
#' data(BS.chr22)
#' bs <- BS.chr22
#' bs_filt <- filterCov(bs = bs,
#'                      min_cov = 10,
#'                      row_min = 0.5,
#'                      col_min = 0.5)
#' bs
#' bs_filt
#' @export
filterCov <- function(bs,
                      min_cov = 10,
                      row_min = 0.5,
                      col_min = 0.5){

  if (!inherits(bs, "BSseq"))
    stop("Input must be a BSseq object")

  cov <- getCoverage(bs, type = "Cov")

  #### 1) Filter CpGs (rows) ----
  pass_per_cpg <- rowMeans(cov >= min_cov)
  keep_cpg <- pass_per_cpg >= row_min

  message(sum(keep_cpg), " of ", length(keep_cpg),
          " CpGs passed CpG-level coverage filter.")

  bs_filtered <- bs[keep_cpg, ]

  #### 2) Filter samples (columns) ----
  cov_filt <- getCoverage(bs_filtered, type = "Cov")

  pass_per_sample <- colMeans(cov_filt >= min_cov)
  keep_samples <- pass_per_sample >= col_min

  message(sum(keep_samples), " of ", length(keep_samples),
          " samples passed sample-level coverage filter.")

  bs_filtered <- bs_filtered[, keep_samples]

  return(bs_filtered)
}
