#' @title
#'
#' @description
#' \code{\link{patternCpG}} computes the correlation patter on the statistics computed by statsCpG
#'
#'
#' @param rownum
#' @param colnum
#' @param meanval = rep(mean(dataset, na.rm=TRUE), colnum)
#' @param sdval = rep(sd(dataset, na.rm=TRUE), colnum)
#' @param dataset
#' @param matrix default="M"
#'
#'
#' @name patternCpG
#'
#' @return
#' #' list with 3 datasets:
#'
#'  \item{Simulated_matrix}{}
#'  \item{Original_correlation_sample }{}
#'
#'
#'
#' @export
#'
#'
patternCpG <- function(rownum, colnum,
                        meanval = rep(mean(dataset, na.rm=TRUE), colnum),
                        sdval = rep(sd(dataset, na.rm=TRUE), colnum),
                        dataset, matrix="M") {


  pd_corr_matrix <- bigcor(t(dataset), fun="cor")
  pd_corr_matrix[which(is.na(pd_corr_matrix))] <- 0
  pd_corr_matrix <- pd_corr_matrix[1:nrow(pd_corr_matrix), 1:ncol(pd_corr_matrix)]

  #for (i in 1:colnum){

  mu <- meanval
  stddev <- sdval

  cat("Computing covariance matrix ...\n")
  covMat <- as.data.frame(as.vector(stddev) %*% t(as.vector(stddev))) * pd_corr_matrix
  covMat[is.na(covMat)] <- 0


  # Transforming in a semi-positive definite matrix
  eigen_decomp <- eigen(as.matrix(covMat), symmetric = TRUE)
  eigenvalues <- eigen_decomp$values
  eigenvalues[eigenvalues < 0] <- 0

  cat("Generating values ...\n")
  covMat2 <- eigen_decomp$vectors %*% diag(eigenvalues) %*% t(eigen_decomp$vectors)
  covMat2 <- as.matrix(covMat2)

  X_hat <- MASS::mvrnorm(n = rownum, mu = mu, Sigma = covMat2)  # Simulated values

  original_sample <- pd_corr_matrix[1:colnum, 1:colnum]

  #nearPD_sample <- stats::cor(X_hat)[1:colnum, 1:colnum]

  rownames(X_hat) <- 1:nrow(X_hat)

  if (matrix != "M"){
    X_hat[which(X_hat <0)] <- 0
    #X_hat <- round(X_hat)

  }

  cat("Converting results ...\n")
  simulated <- list(Simulated_matrix = X_hat, Original_correlation_sample = original_sample)

  return(simulated)

}


