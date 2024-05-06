#' @title adjustCounts
#' @description
#' \code{}Filtering step of methylated and unmethylated matrices. `adjustCounts` is a function which helps in adjusting methylated matrix and unmethylated matrix when filtering the CpGs according to the already filtered coverage
#'
#'
#' @param coverage filtered matrix of coverage counts from the list of assays of the SummarizedExperiment object.
#' @param methylated matrix of methylated counts not already adjusted from the list of assays of the SummarizedExperiment object.
#' @param unmethylated matrix of unmethylated counts not already adjusted from the list of assays of the SummarizedExperiment object.
#'
#' @name adjustCounts
#'
#' @return
#'  \item{List with 2 elements}{filtered methylated counts matrix and filtered unmethylated counts matrix.}
#'
#'
#' @examples
#' \dontrun{
#'  data("rangedObject")
#'  C<- assay(data, 1)
#'  M<- assay(data, 2)
#'  U<- assay(data, 3)
#'  res <- adjustCounts(coverage = C, methylated = M, unmethylated = U)
#' }
#'
#' @noRd
#'
adjustCounts <- function(coverage, methylated, unmethylated){


  # NA
  # Create a new matrix based on conditions

  # Initialize
  met2 <- methylated
  co<- coverage
  # Identify the indices where the condition is met
  indices <- which(is.na(co) & unmethylated == 0, arr.ind = TRUE)

  # Replace values at identified indices with NA
  met2[indices] <- NA


  # Initialize
  unmet2 <- unmethylated

  indices2 <- which(is.na(co) & is.na(met2), arr.ind = TRUE)
  unmet2[indices2] <- NA


  methylated <- met2
  unmethylated <- unmet2

  res <- list(methylated, unmethylated)

  return(res)



}  # function
