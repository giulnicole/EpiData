#' Process `BSseq` object
#'
#' @description
#' \itemize{
#'   \item Zero-coverage sites are converted to \code{NA} since they are considered as missing data
#'   \item Adds rownames to the \code{BSseq} object
#' }
#'
#' @param bs A \code{BSseq} object.
#' @return Returns the same \code{BSseq} object with rownames
#'
#' @importFrom bsseq granges getCoverage
#' @importFrom GenomicRange seqnames ranges
#'
#' @examples
#'
#' library(bsseq)
#' data(BS.chr22)
#'
#' bs1 <- BS.chr22
#' cov1 <- bsseq::getCoverage(BSseq = bs1, type = "Cov")
#' head(cov1, n = 10)
#'
#' bs2 <- processBSseq(bs1)
#' cov2 <- bsseq::getCoverage(BSseq = bs2, type = "Cov")
#' head(cov2, n = 10)
#'
#' @export
processBSseq <- function(bs, verbose = TRUE){
  if(verbose == TRUE){
    message("Adding rownames to the BSseq object")
  }
  gr <- bsseq::granges(bs)
  row_id <- paste(GenomicRanges::seqnames(gr), GenomicRanges::ranges(gr), sep = "_")
  rownames(bs) <- row_id

  if(verbose == TRUE){
    message("Zero-coverage sites are converted to 'NA'")
  }
  cov <- bsseq::getCoverage(BSseq = bs, type = "Cov")
  na_id <- cov == 0

  bs@assays@data$Cov[na_id] <- NA
  bs@assays@data$M[na_id] <- NA

  return(bs)
}
