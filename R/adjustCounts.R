#' @title adjustCounts
#' @description
#' \code{\link{adjustCounts}} internal function which helps in adjusting methylated matrix and unmethylated matrix when filtering the CpGs according to the already filtered coverage
#'
#'
#' @param coverage filtered matrix from the SummarizedExperiment object.
#' @param methylated counts' matrix to be adjusted from the SummarizedExperiment object.
#' @param unmethylated counts'matrix to be adjusted from the SummarizedExperiment object.
#'
#' @name adjustCounts
#'
#' @return
#' a list with 2 elements:
#'  \item{methylated}{filtered methylated counts matrix}
#'  \item{unmethylated}{filtered unmethylated counts matrix}
#'
#'
#' @examples
#' \dontrun{
#'  data("matrices")
#'  C<- assay(dati, 1)
#'  M<- assay(dati, 2)
#'  U<- assay(dati, 3)
#'  res <- adjustCounts(coverage = C, methylated = M, unmethylated = U)
#' }
#'
#' @export
#'
adjustCounts <- function(coverage, methylated, unmethylated){


  # NA
  # Create a new matrix based on conditions

  # Initialize
  met2 <- methylated
  co <- coverage
  un <- unmethylated

  for (i in 1:nrow(met2)) {
    for (j in 1:ncol(met2)) {
      if (is.na(co[i, j]) && un[i, j] == 0) {
        met2[i, j] <- NA
      } else {
        met2[i, j] <-  methylated[i, j]
      }
    }
  }


  # Initialize
  unmet2 <- unmethylated

  for (i in 1:nrow(unmet2)) {
    for (j in 1:ncol(unmet2)) {
      if (is.na(co[i, j]) && is.na(met2[i, j])) {
        unmet2[i, j] <- NA
      } else {
        unmet2[i, j] <-  unmethylated[i, j]
      }
    }
  }


  # Initialize
  met3 <- met2

  for (i in 1:nrow(met3)) {

    for (j in 1:ncol(met3)) {

      if (!is.na(co[i, j]) && co[i, j] == unmet2[i, j]) {
        #if (met3[i, j] == 0) {

        cpg <- as.numeric(met3[i,])
        wil <- wilcox.test(cpg, mu = 0, alternative = "greater")
        pval <- wil$p.value


        if (pval< 0.01) {
          # cat("NA \n")
          met3[i, j] <- NA


        }   else {
          # cat("not NA \n")
          met3[i, j] <- 0

        }


      } else {
        met3[i, j] <-  met2[i, j]
      }
    }
  }


  # Initialize
  unmet3 <- unmet2

  for (i in 1:nrow(unmet3)) {

    for (j in 1:ncol(unmet3)) {

      if (!is.na(co[i, j]) && !is.na(met3[i, j]) && co[i, j] == met3[i, j]) {
        #if (met3[i, j] == 0) {

        cpg <- as.numeric(unmet3[i,])
        wil <- wilcox.test(cpg, mu = 0, alternative = "greater")
        pval <- wil$p.value


        if (pval< 0.01) {
          # cat("NA \n")
          unmet3[i, j] <- NA


        }   else {
          # cat("not NA \n")
          unmet3[i, j] <- 0

        }


      } else {
        unmet3[i, j] <-  unmet2[i, j]
      }
    }
  }

  methylated <- met3
  unmethylated <- unmet3

  res <- list(methylated, unmethylated)

  return(res)



}  # function
