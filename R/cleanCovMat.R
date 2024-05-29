#' @title clenCovMat
#' @description
#' \code{}First part of CpGs' cleaning. `cleanCovMat` helps in the conversion of zeroes into missing values (0->NA) in the coverage matrix,
#'  removing rows (CpGs) and columns (individuals) above pre-specified missingness threshold.
#'
#' @import tidyverse
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import stats
#' @import htmltools
#' @import ggplot2
#' @importFrom dplyr mutate_all
#'
#'
#' @param input.obj SummarizedExperiment assay object (list) as input with list with coverage, methylated counts and unmethylated counts in assays/data/listData (rows are CpGs and individuals are columns).
#' @param max_na_cpg Threshold of missing values per each CpG in the coverage counts matrix.
#' @param max_na_ind Threshold of missing values per each individual.
#' @param cpg_removal_threshold Minimum threshold of coverage counts (across all individuals) for a CpG to be kept.
#'
#' @name cleanCovMat
#'
#' @return
#'  \item{List with 3 elements}{coverage matrix, methylated counts matrix, unmethylated counts matrix, cleaned at the specified missingness thresholds removed.}
#'
#'
#' @examples
#' library(dplyr)
#' library(tidyverse)
#' library(SummarizedExperiment)
#' input.list<- assays(rangedObject)
#' clean.coverage <- cleanCovMat(input.obj = input.list, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)
#'
#'
#' @export
#'
#'
cleanCovMat <- function(input.obj, max_na_cpg=0.5, max_na_ind=0.2, cpg_removal_threshold=10) {

  X<- input.obj$Coverage_matrix
  Y<- input.obj$Met_matrix
  Z<- input.obj$Unmet_matrix

  sites <- rownames(X)
  stopifnot("Input must be numeric dataframe" =is.data.frame(X), all(sapply(X, is.numeric)))
  cat("Converting 0s in NAs for coverage counts ...", "\n")
  X <- X %>% dplyr::mutate_all(~na_if(., 0))


  # Convert the object to a matrix (if needed)
  is_matrix <- is.matrix(X)
  if (!is_matrix) {
    X <- as.matrix(X)
  }


  # FIRST PART: discard CpGs with high percentage of missing (coverage)

  # Calculate the maximum allowed missing values for rows and columns
  max_missing_rows <- ncol(X) * max_na_cpg
  max_missing_cols <- nrow(X) * max_na_ind

  # Discard CpGs with > max_na_cpg % missing values
  row_indices_to_keep <- rowSums(is.na(X)) <= max_missing_rows

  X2 <- X[row_indices_to_keep, ]
  dim(X2)

  eliminated <- nrow(X) - nrow(X2)
  index1 <- which(row_indices_to_keep==TRUE)
  sites_filtered1 <- sites[index1]
  rownames(X2) <- sites_filtered1

  if (eliminated != 0) {
    message(paste(eliminated, " CpG(s) removed due to exceeding the pre-defined removal threshold (>",
                  max_na_cpg * 100, "%) for missingness.", sep = "")) }


  # SECOND PART: discard individuals with high percentage of missing (coverage)

  # Discard individuals with > max_na_ind % missing values
  col_indices_to_keep <- colSums(is.na(X2)) <= max_missing_cols
  X3 <- X2[, col_indices_to_keep]
  eliminated2 <- ncol(X2) - ncol(X3)
  index2 <- which(col_indices_to_keep==TRUE)


  if (eliminated2 != 0) {
    message(paste(eliminated2, " individual(s) removed due to exceeding the pre-defined removal threshold (>",
                  max_na_ind * 100, "%) for missingness.", sep = "")) }

  # THIRD PART: discard CpGs with less than 10 counts per individuals
  # Remove counts below 10 each
  index_above_threshold <- which(apply(X3, 1, median, na.rm=T)  > cpg_removal_threshold)
  X4 <- X3[index_above_threshold, ]
  eliminated3 <- nrow(X3) - nrow(X4)

  sites_filtered2 <- sites_filtered1[index_above_threshold]
  rownames(X4) <- sites_filtered2


  if (eliminated3 != 0) {
    message(paste(eliminated3, " Cpg(s) removed due to not exceeding the pre-defined removal threshold (<",
                  cpg_removal_threshold, ") for counts.", sep = "")) }



  # Filter methylated counts and unmethylated count matrices
  cat("Computing selection of CpGs in methylated counts matrix ...", "\n")

  # Return meth and unmeth counts
  Y <- Y[rownames(Y) %in% rownames(X4),]
  Y <- Y[, colnames(Y) %in% colnames(X4)]
  #Y <- Y %>% mutate_all(~na_if(., 0))

  cat("Computing selection of CpGs in unmethylated counts matrix ...", "\n")
  Z <- Z[rownames(Z) %in% rownames(X4),]
  Z <- Z[, colnames(Z) %in% colnames(X4)]
  #Z <- Z %>% mutate_all(~na_if(., 0))

  cat("Adjusting the NAs in methylated and unmethylated counts' matrices ...", "\n")
  # Adjusting the zeros in methylated and unmethylated counts
  res <- adjustCounts(coverage = X4, methylated = Y, unmethylated = Z)

  Y <- res[[1]]
  Z <- res[[2]]
  X4 <- as.data.frame(X4)

  # Results
  res<- list(Coverage_matrix = X4, Met_matrix = Y, Unmet_matrix = Z)

  cat("Coverting results to matrices ...", "\n")
  res2<- lapply(res, as.data.frame)


  objGR<- res2

  #objGR<- GRconversion2(cleaned.list = results$Cleaned1)

  return(objGR)

}  # (main function)





