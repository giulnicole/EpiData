#' @title Split3MatXChrom
#' @description
#' \code{\link{Split3MatXChrom}} is a function which computes the Beta and M matrices from the filtered filtered coverage, methylated and unmethylated counts.
#'
#'
#' @import GenomicFeatures
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom IRanges IRanges
#'
#' @param filtered.obj with the three experimental matrices from the SummarizedExperiment object: coverage counts filtered matrix, methylated counts matrix, unmethylated counts matrices.
#'
#' @name Split3MatXChrom
#'
#' @return
#' a list of 5 matrices:
#'  \item{coverage counts}{methylated counts matrix}
#'  \item{methylated counts}{methylated counts matrix}
#'  \item{unmethylated counts}{unmethylated counts matrix}
#'  \item{Beta values}{Beta values counts matrix}
#'  \item{M values}{M values counts matrix}
#'
#' @export
#'
#'
MethilationValues <- function(filtered.obj){

  clean.met2 <- filtered.obj$Met_matrix
  clean.unmet2 <- filtered.obj$Unmet_matrix


  # Calculating beta values matrix and M values on cleaned data
  cat("Calculating beta values and M values ...", "\n")
  p <- 1  # offset parameter


  # B values' matrix
  B <- matrix(0, nrow(clean.met2), ncol(clean.met2))

  # M values' matrix
  M <- matrix(0, nrow(clean.met2), ncol(clean.met2))

  for (j in 1:nrow(clean.met2)){

    for (i in 1:ncol(clean.met2)) {

      meth <- as.numeric(clean.met2[[j,i]])
      unmeth <- as.numeric(clean.unmet2[[j,i]])

      beta <- max(meth, 0)/(max(meth,0) + max(unmeth,0) + p)
      B[j,i] <- beta

      m <- log2((max(meth,0) + p)/ (max(unmeth,0) +p))
      M[j,i] <- m

    } # for i

  } # for j


  colnames(B)<- colnames(clean.met2)
  rownames(B)<- rownames(clean.met2)
  cat("Beta matrix computed ...", "\n")

  colnames(M)<- colnames(clean.met2)
  rownames(M)<- rownames(clean.met2)
  cat("M matrix computed ...", "\n")

  #list.cleaned[[4]] <- B
  #list.cleaned[[5]] <- M


  clean.cov <- filtered.obj$Coverage_matrix
  clean.cov <- as.data.frame(clean.cov[rownames(clean.cov) %in% rownames(clean.met2),])


  list.cleaned2 <- list(Coverage_matrix = as.data.frame(clean.cov),
                        Met_matrix = as.data.frame(clean.met2),
                        Unmet_matrix = as.data.frame(clean.unmet2),
                        Beta_matrix = as.data.frame(B),
                        M_matrix = as.data.frame(M))

  list.cleaned3 <- split5MatXChrom(list.cleaned2)

  # obj<- GRconversion(list.cleaned2)

   return(list.cleaned3)

}  # (main function)

