#' @title patternCpG
#'
#' @description A function for computing the non-negative definite matrix of covariance matrix, by eigenvectors' decompisition and
#' and simulate new values drawn from the same distribution to impute remaining missing values in the original cleaned dataset.
#' \code{\link{patternCpG}} computes the imputation based on correlation pattern between CpGs derived from eigen decomposition of the covariance matrix.
#'
#' @name patternCpG
#'
#' @param cleaned.obj SummarizedExperimet object cleaned and after statistics and NAs' pattern computation (statsCpG output).
#' @param matrix Matrix on which computing the correlation matrix; "M" as default.
#' @param varselect Index of the dataset to be used (numeric value 1-5). 1 = coverage counts, 2 = methylated counts, 3 = unmethylated counts, 4 = beta values, 5 = M values.
#'
#' @return
#' dataset imputed:
#' \item{Simulated_matrix}{Matrix of imputed CpGs via eigenvalue decomposition of covariance matrix and generation of values from a distribution similar per each CpG.}
#'
#'
#' @examples
#' # data("matrices")
#' # clean.coverage2 <- cleanCovMat(input.obj=dati, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)
#' # clean.out <- cleanOutliers(filtered.obj=clean.coverage2, outlier_threshold=5, remove_outliers = TRUE)
#' # stats<- statsCpG(cleaned.obj=clean.out, varselect = 5)
#' # simu3 <- patternCpG(cleaned.obj=stats, matrix="M", varselect)
#'
#' @export
#'
#'
patternCpG <- function(cleaned.obj,
                       matrix="M",
                       varselect) {


  rownum <- cleaned.obj[[1]]@metadata[["statistics"]][["Individuals"]]
  colnum <- cleaned.obj[[1]]@metadata[["statistics"]][["CpGs"]]
  meanval <- cleaned.obj[[1]]@metadata[["statistics"]][["means"]]
  sdval <- cleaned.obj[[1]]@metadata[["statistics"]][["sds"]]

  dataset <-  cleaned.obj[["Output_outliers"]]@assays@data[[varselect]]

  pd_corr_matrix <- bigcor(t(dataset), fun="cor")
  pd_corr_matrix[which(is.na(pd_corr_matrix))] <- 0
  pd_corr_matrix <- pd_corr_matrix[1:nrow(pd_corr_matrix), 1:ncol(pd_corr_matrix)]

  #for (i in 1:colnum){
  mu <- meanval
  stddev <- sdval

  cat("Computing covariance matrix ...\n")
  covMat <- as.data.frame(as.vector(stddev) %*% t(as.vector(stddev))) * pd_corr_matrix
  covMat[is.na(covMat)] <- 0
  epsilon <- 1e-6
  covMat <- covMat + epsilon * diag(length(mu))

  # Transforming in a semi-positive definite matrix
  eigen_decomp <- eigen(as.matrix(covMat), symmetric = TRUE)
  eigenvalues <- eigen_decomp$values
  eigenvalues[eigenvalues < 0] <- 0

  cat("Generating values ...\n")
  covMat2 <- eigen_decomp$vectors %*% diag(eigenvalues) %*% t(eigen_decomp$vectors)
  covMat2 <- as.matrix(covMat2)

  X_hat <- MASS::mvrnorm(n = rownum, mu = mu, Sigma = covMat2)  # Simulated values

  #original_sample <- pd_corr_matrix[1:colnum, 1:colnum]

  #nearPD_sample <- stats::cor(X_hat)[1:colnum, 1:colnum]

  rownames(X_hat) <- 1:nrow(X_hat)

  if (matrix != "M"){
    X_hat[which(X_hat <0)] <- 0
    #X_hat <- round(X_hat)

  }

  cat("Converting results ...\n")
  simulated <- list(Simulated_matrix = X_hat)


  #names(simulated) <- char
  mat.m <- as.data.frame(t(dataset))
  na_positions <- is.na(mat.m)
  simu.m <- simulated$Simulated_matrix
  mat.m[na_positions] <- t(simu.m[na_positions])

  simu.mix2 <- as.data.frame(mat.m)


  return(simu.mix2)

}


