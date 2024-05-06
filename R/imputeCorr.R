#' @title imputeCorr
#'
#' @description A function for computing the non-negative definite matrix of covariance matrix, by eigenvectors' decompisition and
#' and simulate new values drawn from the same distribution to impute remaining missing values in the original cleaned dataset.
#' \code{} computes the imputation based on correlation pattern between CpGs derived from eigen decomposition of the covariance matrix.
#'
#' @import latentcor
#' @import propagate
#' @importFrom matrixcalc is.positive.semi.definite
#'
#' @name imputeCorr
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
#' \dontrun{
#' # Data
#'  data("rangedObject")
#'  input.list<- assays(rangedObject)
#'
#'  # Cleaning low counts for coverage
#'
#'  clean.coverage <- cleanCovMat(input.obj = input.list,
#'                               max_na_cpg = 0.5, max_na_ind = 0.2,
#'                              cpg_removal_threshold = 10)
#'
#' # Cleaning outliers in methylated and unmethylated counts
#'  clean.out <- cleanOutliers(filtered.obj=clean.coverage,
#'                             outlier_threshold=5, remove_outliers = TRUE)
#'
#'  stats <- statsCpG(cleaned.obj=clean.out, varselect = 5)
#'  simu <- imputeCorr(cleaned.obj=stats, matrix="M", varselect)
#'  }
#'
#'
#'
#'
#'
#'
#' @export
#'
#'
imputeCorr <- function(cleaned.obj,
                       matrix="M",
                       varselect=5) {


  rownum <- cleaned.obj[["Individuals"]]
  colnum <- cleaned.obj[["CpGs"]]
  meanval <- cleaned.obj[["means"]]
  sdval <- cleaned.obj[["sds"]]


  dataset <- as.data.frame(t(cleaned.obj[["Matrix"]]))  # add

      Row.means <- rowMeans(dataset, na.rm = T)

      #with rowmeans
      for (i in 1:nrow(dataset)) {

       dataset[i, is.na(dataset[i, ])] <- Row.means[i]

        }


      cat("Computing correlation matrix based on imputed mean per CpG\n")

      pd_corr_matrix <- cor(dataset)

      # pd_corr_matrix <- pd_corr_matrix[1:nrow(pd_corr_matrix), 1:ncol(pd_corr_matrix)]


    cat("Computing correlation matrix with pairwise observations\n")



  #for (i in 1:colnum){
  mu <- meanval
  stddev <- sdval

  cat("Computing covariance matrix ...\n")
  covMat <- as.data.frame(as.vector(stddev) %*% t(as.vector(stddev))) * pd_corr_matrix
  covMat[is.na(covMat)] <- 0
  #epsilon <- 1e-6
  #covMat <- covMat + epsilon * diag(length(mu))


  if (matrixcalc::is.positive.semi.definite(as.matrix(covMat)) !=TRUE | isSymmetric.matrix(as.matrix(covMat)) != TRUE){

    cat("Eigen decomposition...\n")
    # Transforming in a semi-positive definite matrix
    eigen_decomp <- eigen(as.matrix(covMat), symmetric = TRUE)
    eigenvalues <- eigen_decomp$values
    eigenvalues[eigenvalues < 0] <- 0

    cat("Generating values ...\n")
    covMat2 <- eigen_decomp$vectors %*% diag(eigenvalues) %*% t(eigen_decomp$vectors)
    covMat2 <- as.matrix(covMat2)
    X_hat <- MASS::mvrnorm(n = rownum, mu = mu, Sigma = covMat2)

    } else{

  cat("Generating values ...\n")
  X_hat <- MASS::mvrnorm(n = rownum, mu = mu, Sigma = covMat)

  } # Simulated values

  #original_sample <- pd_corr_matrix[1:colnum, 1:colnum]

  #nearPD_sample <- stats::cor(X_hat)[1:colnum, 1:colnum]

  rownames(X_hat) <- 1:nrow(X_hat)

  if (matrix != "M"){
    X_hat[which(X_hat <0)] <- 0
    #X_hat <- round(X_hat)

  }

  cat("Converting results ...\n")

  simulated <- X_hat



  #names(simulated) <- char
  mat.m <- as.data.frame(dataset)
  na_positions <- is.na(mat.m)
  mat.m[na_positions] <- simulated[na_positions]



  simu.mix2 <- list(Imputed = as.data.frame(mat.m), Simulated = as.data.frame(simulated))

  cat("done...\n")
  return(simu.mix2)


}
