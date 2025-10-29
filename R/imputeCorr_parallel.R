#' @title imputeCorr_parallel
#' @description
#' \code{} Parallelized version of imputeCorr that processes multiple chromosomes. It takes a list where each element contains the chromosome data with "Matrix", "Individuals", "CpGs", "means", "sds".
#'
#' @importFrom matrixcalc is.positive.semi.definite
#' @importFrom parallel makeCluster stopCluster parLapply detectCores clusterExport clusterEvalQ
#' @import stats
#'
#' @param cleaned.obj List of lists. Each element contains: Individuals, CpGs, means, sds, and Matrix.
#' @param matrix Matrix on which computing the correlation matrix; "M" as default.
#' @param varselect Index of the dataset to be used (numeric value 1-5).
#' @param n.cores Number of cores to use. Default is NULL (uses detectCores() - 1).
#' @name imputeCorr_parallel
#'
#'
#' @return List of imputed datasets, one per chromosome
#'
#' @examples
#' \dontrun{
#'  data("bs_list")
#'  clean.coverage <- cleanCovMat(input.obj = bs_list, max_na_cpg = 0.5,
#'  max_na_ind = 0.2, cpg_removal_threshold = 10)
#'  clean.out <- cleanOutliers(filtered.obj=clean.coverage,
#'  outlier_threshold=0, remove_outliers = TRUE)
#'  mat.values <- methValues(clean.coverage)
#'  stat.result <- statsCpG_parallel(mat.values)
#'  imputed <- imputeCorr_parallel(stat.result)
#' }
#'
#'
#' @export
imputeCorr_parallel <- function(cleaned.obj,
                                matrix = "M",
                                varselect = 5,
                                n.cores = NULL) {

  # Detect vignette or non-interactive mode
  if (!interactive() || !is.null(knitr::opts_knit$get("rmarkdown.pandoc.to"))) {
    message("Non-interactive or vignette build detected, running in serial mode.")
    n.cores <- 1
  } else {
    if (is.null(n.cores)) {
      n.cores <- parallel::detectCores() - 1
    }
    n.cores <- max(1, min(n.cores, length(cleaned.obj)))
  }

  cat(paste("Processing", length(cleaned.obj), "chromosomes using", n.cores, "cores (parLapply)\n"))

  matrix_type <- matrix
  var_select <- varselect

  # --- SERIAL EXECUTION (for vignette / non-interactive) ---
  if (n.cores == 1) {
    results <- lapply(cleaned.obj, function(chr_data) {

      rownum <- chr_data[["Individuals"]]
      colnum <- chr_data[["CpGs"]]
      meanval <- chr_data[["means"]]
      sdval <- chr_data[["sds"]]

      dataset <- as.data.frame(t(chr_data[["Matrix"]]))

      Row.means <- rowMeans(dataset, na.rm = TRUE)

      for (i in 1:nrow(dataset)) {
        dataset[i, is.na(dataset[i, ])] <- Row.means[i]
      }

      pd_corr_matrix <- cor(dataset)

      mu <- meanval
      stddev <- sdval

      covMat <- as.data.frame(as.vector(stddev) %*% t(as.vector(stddev))) * pd_corr_matrix
      covMat[is.na(covMat)] <- 0

      if (matrixcalc::is.positive.semi.definite(as.matrix(covMat)) != TRUE |
          isSymmetric.matrix(as.matrix(covMat)) != TRUE) {

        eigen_decomp <- eigen(as.matrix(covMat), symmetric = TRUE)
        eigenvalues <- eigen_decomp$values
        eigenvalues[eigenvalues < 0] <- 0

        covMat2 <- eigen_decomp$vectors %*% diag(eigenvalues) %*% t(eigen_decomp$vectors)
        covMat2 <- as.matrix(covMat2)
        X_hat <- MASS::mvrnorm(n = rownum, mu = mu, Sigma = covMat2)

      } else {
        X_hat <- MASS::mvrnorm(n = rownum, mu = mu, Sigma = covMat)
      }

      rownames(X_hat) <- 1:nrow(X_hat)

      if (matrix_type != "M") {
        X_hat[which(X_hat < 0)] <- 0
      }

      simulated <- X_hat

      mat.m <- as.data.frame(dataset)
      na_positions <- is.na(mat.m)
      mat.m[na_positions] <- simulated[na_positions]

      return(list(Imputed = as.data.frame(mat.m), Simulated = as.data.frame(simulated)))
    })

    return(results)
  }

  # --- PARALLEL EXECUTION (interactive mode only) ---
  cl <- parallel::makeCluster(n.cores)

  tryCatch({
    parallel::clusterEvalQ(cl, {
      library(matrixcalc)
      library(MASS)
      library(stats)
    })

    parallel::clusterExport(cl, varlist = c("matrix_type", "var_select"), envir = environment())

    cat("Starting parallel processing...\n")

    results <- parallel::parLapply(cl, cleaned.obj, function(chr_data) {

      rownum <- chr_data[["Individuals"]]
      colnum <- chr_data[["CpGs"]]
      meanval <- chr_data[["means"]]
      sdval <- chr_data[["sds"]]

      dataset <- as.data.frame(t(chr_data[["Matrix"]]))

      Row.means <- rowMeans(dataset, na.rm = TRUE)

      for (i in 1:nrow(dataset)) {
        dataset[i, is.na(dataset[i, ])] <- Row.means[i]
      }

      pd_corr_matrix <- cor(dataset)

      mu <- meanval
      stddev <- sdval

      covMat <- as.data.frame(as.vector(stddev) %*% t(as.vector(stddev))) * pd_corr_matrix
      covMat[is.na(covMat)] <- 0

      if (matrixcalc::is.positive.semi.definite(as.matrix(covMat)) != TRUE |
          isSymmetric.matrix(as.matrix(covMat)) != TRUE) {

        eigen_decomp <- eigen(as.matrix(covMat), symmetric = TRUE)
        eigenvalues <- eigen_decomp$values
        eigenvalues[eigenvalues < 0] <- 0

        covMat2 <- eigen_decomp$vectors %*% diag(eigenvalues) %*% t(eigen_decomp$vectors)
        covMat2 <- as.matrix(covMat2)
        X_hat <- MASS::mvrnorm(n = rownum, mu = mu, Sigma = covMat2)

      } else {
        X_hat <- MASS::mvrnorm(n = rownum, mu = mu, Sigma = covMat)
      }

      rownames(X_hat) <- 1:nrow(X_hat)

      if (matrix_type != "M") {
        X_hat[which(X_hat < 0)] <- 0
      }

      simulated <- X_hat

      mat.m <- as.data.frame(dataset)
      na_positions <- is.na(mat.m)
      mat.m[na_positions] <- simulated[na_positions]

      return(list(Imputed = as.data.frame(mat.m), Simulated = as.data.frame(simulated)))
    })

    cat("Parallel processing completed!\n")

  }, finally = {
    parallel::stopCluster(cl)
  })

  return(results)
}
