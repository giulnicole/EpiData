#' @title cleanOutliers
#' @description
#' \code{} Second part of CpGs' cleaning. `cleanOutliers` helps in discarding those CpGs with high rate of missing values of methylated and unmethylated counts after
#' coverage cleaning and outlier values in the conversion.
#'
#' @import stats
#'
#' @param filtered.obj List object deriving from cleanMatCov.
#' @param outlier_threshold Threshold for considering a value of the counts' matrices as an outlier.
#' @param remove_outliers TRUE/FALSE whether removing of outliers should be performed or not.
#'
#' @name cleanOutliers
#'
#' @return
#' SummarizedExperiment object with 2 elements:
#'  \item{List with 3 elements}{cleaned matrices (coverage, methylated and unmethylated counts)}
#'
#'
#' @examples
#' library(dplyr)
#' library(tidyverse)
#' data("meth_data")
#' print(meth_data)
#' input.list<- assays(meth_data)
#' subset_list <- lapply(input.list, function(df) df[1:100, 1:10])
#' clean.coverage <- cleanCovMat(input.obj = subset_list, max_na_cpg = 0.5,
#' max_na_ind = 0.2,  cpg_removal_threshold = 10)
#' clean.out <- cleanOutliers(filtered.obj=clean.coverage,
#' outlier_threshold=0, remove_outliers = TRUE)
#'
#'
#'
#' @export
#'
#'
cleanOutliers<- function(filtered.obj, outlier_threshold=5,  remove_outliers=F){

  mat1<- filtered.obj$Met_matrix
  mat2<- filtered.obj$Unmet_matrix


  list.mat <- list(mat1, mat2)
  names(list.mat) <- c("Methylated", "Unmethylated")

  # Adjust this threshold as needed

  ##########

  # FIRST PART: cleaning CpGs that have still high percentage of missing (met, unmet)
  m2 <- list()

  cat("Excluding CpGs that still have too many NAs from in methylated and unmethylated counts' matrices\n")
  for (k in 1:length(list.mat)){

    m <- list.mat[[k]]
    data_with_na <- m

    n <- nrow(m)
    X <- m


    # Convert the object to a matrix (if needed)
    is_matrix <- is.matrix(X)
    if (!is_matrix) {
      X <- as.matrix(X)
    }


    throw<- which(is.na(rowMeans(X, na.rm = T)))
    throw <- as.matrix(throw)

    if(dim(throw)[1]==0) {
      X2 <- X

    } else {
      X2 <- X[-throw,]

    }

    cat(paste(dim(throw)[1], " rows have been removed because all entries were NAs in", names(list.mat)[[k]],"\n"))


    m2[[k]] <- X2

  }


  #######

  # SECOND PART: cleaning outliers (met, unmet)

  if (remove_outliers ==TRUE) {

    cat("Converting outliers into NAs ...\n")

    list.mat2 <- m2

    m3 <- list()


    for (k in 1:length(list.mat2)){

      m <- list.mat2[[k]]

      remove_outlier <- function(row) {
        mean_row <- mean(row)
        sd_row <- sd(row)
        is_outlier <- abs(row - mean_row) > 5 * sd_row
        row[is_outlier] <- NA  # Imposta gli outlier a NA
        return(row)

      }

      # Apply the function
      data_without_outliers <- as.data.frame(t(apply(m, 1, remove_outlier)))

      m3[[k]] <- data_without_outliers


    }



    clean.met<- m3[[1]]
    clean.unmet <- m3[[2]]

    clean.met2 <- clean.met[rownames(clean.met) %in% rownames(clean.unmet) ,]
    clean.unmet2 <-  clean.unmet[rownames(clean.unmet) %in% rownames(clean.met2) ,]


  } else{

    clean.met2 <- m2[[1]]
    clean.unmet2 <- m2[[2]]

  }



  clean.cov <- filtered.obj$Coverage_matrix
  clean.cov <- as.data.frame(clean.cov[rownames(clean.cov) %in% rownames(clean.met2),])


  list.cleaned2 <- list(Coverage_matrix = as.data.frame(clean.cov),
                        Met_matrix = as.data.frame(clean.met2),
                        Unmet_matrix = as.data.frame(clean.unmet2) )


  list.cleaned3 <- split3MatXChrom(list.cleaned2)

  lapply(list.cleaned3, as.matrix)

  return(list.cleaned3)


}  # (main function)



