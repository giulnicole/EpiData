#' @title extractFinalMat
#'
#' @description A function for binding all the CpG of all the chromosomes together.
#'
#' @name extractFinalMat
#'
#' @param list.imputed Imputed list.
#'
#' @return
#' List with all the imputed matrices per each chromosome bound together.
#'
#'
#' @examples
#' \dontrun{
#'  data("matrices")
#'  clean.coverage2 <- cleanCovMat(input.obj=dati, max_na_cpg = 0.5, max_na_ind = 0.2,
#'  cpg_removal_threshold = 10)
#'  clean.out <- cleanOutliers(filtered.obj=clean.coverage2,
#'  outlier_threshold=5, remove_outliers = TRUE)
#'  splitted.list <- split5MatXChrom(clean.out)
#'  stats<- statsCpG(cleaned.obj=clean.out, varselect = 5)
#'  simu3 <- imputeCorr(cleaned.obj=stats, matrix="M", varselect)
#'  imputed.mat <- extractFinalMat(list.imputed=simu3)
#'  }
#'
#' @export
extractFinalMat <- function(list.imputed) {

  # Extract the imputed dataset from each sublist
  elements.to.bind<- lapply(list.imputed, function(sublist) sublist[[1]])

  # Combine the padded elements into a matrix
  big_matrix <- do.call(cbind, elements.to.bind)
  big_matrix <- as.data.frame(big_matrix)

  return(big_matrix)
}
