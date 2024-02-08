#' @title ParallSimu
#'
#' @description A function for computing the non-negative definite matrix of covariance matrix, by eigenvectors' decompisition and
#' and simulate new values drawn from the same distribution to impute remaining missing values in the original cleaned dataset (Parallelized version)
#' \code{\link{ParallSimu}} computes the imputation based on correlation pattern between CpGs
#'
#'
#' @param stats.obj list of SummarizedExperimet object cleaned and after statistics and NAs' pattern computation (ParallStats output)
#' @param matrix matrix on which computing the correlation matrix; "M" as default
#' @param varselect index of the dataset to be used (numeric value 1-5)
#'
#'
#'
#' @name ParallSimu
#'
#' @return
#'  list with SummarizedExperiment objects and imputed datasets in metadata:
#' \item{imputed}{Matrix of CpG imputed}
#'
#' @examples
#' # data('clean.out')
#' # list.cleaned1<- clean.out
#' # list.cleaned2<- clean.out
#' # list.cleaned<- list(list.cleaned1$Output_outliers@assays@data@listData, list.cleaned2$Output_outliers@assays@data@listData)
#' # input.stats<- list()
#' # for (i in 1:length(list.cleaned)){
#' #
#' #  input.stats[[i]] <- GRconversion2(list.cleaned[[i]])
#' #
#' #  names(input.stats)[[i]] <- names(list.cleaned)[[i]]
#' #  }
#'
#' # stats2 <- ParallStats(input.stats)
#' # simu.mix2 <- ParallSimu(list.stats = stats2)
#'
#'
#' @export
ParallSimu <- function(list.stats){

  library(stats)
  meta.simu <- bplapply(list.stats, patternCpG)

  return(meta.simu)

}# parallSimu



#' @noRd
patternCpG <- function(cleaned.obj,
                       matrix="M",
                       varselect=5) {


  rownum <- cleaned.obj@metadata[["statistics"]][["Individuals"]]
  colnum <- cleaned.obj@metadata[["statistics"]][["CpGs"]]
  meanval <- cleaned.obj@metadata[["statistics"]][["means"]]
  sdval <- cleaned.obj@metadata[["statistics"]][["sds"]]


  dataset <- as.data.frame(cleaned.obj@assays@data@listData[[varselect]])   # add

  #Row.means <- rowMeans(dataset, na.rm = T)

  #with rowmeans
  #for (i in 1:nrow(dataset)) {
  #  dataset[i, is.na(dataset[i, ])] <- Row.means[i]
  # }


  cat("Computing correlation matrix ...\n")
  pd_corr_matrix <- propagate::bigcor(t(dataset), fun="cor")
  pd_corr_matrix <- round(pd_corr_matrix[1:nrow(pd_corr_matrix), 1:ncol(pd_corr_matrix)], digits = 3)


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
  mat.m[na_positions] <- round(t(simu.m[na_positions]), digits = 3)

  simu.mix2 <- as.data.frame(mat.m)


  cleaned.obj@metadata$imputed <- simu.mix2


  return(cleaned.obj)

}



