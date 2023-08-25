#' @title Data matrix cleaning for missing values excesses after first cleaning and for removing outliers in CpGs' values
#'
#' @description
#' \code{\link{clean_counts_outliers}} helps in the conversion of missing values (0->NA),
#' variable types and removes rows and columns above pre-specified missingness threshold
#'
#'
#' @param list.cleaned cleaned list of dataframes from clean_matrix
#' @param outlier_threshold how many standard deviations to be considered in the standardized distribution for removing outliers
#' @param unreliable_outlier_percentage percentage of outliers per CpGs to be considered the threshold to remove the CpG
#' @name clean_counts_outliers
#'
#' @return
#' #' list with 5 datasets:
#'  \item{Coverage_matrix}{cleaned coverage dataset with NAs as missing values and rows/columns above the pre-specified missingness thresholds removed}
#'  \item{Met_matrix }{methylated counts dataset with the filtered CpGs}
#'  \item{Unmet_matrix}{unmethylated counts dataset with the filtered CpGs}
#'  \item{Beta_matrix}{beta values matrix}
#'  \item{M_matrix}{M values matrix}
#'
#'
#'@examples
#'
#'
#'
#'
#'
#' @export
#'
#'
clean_counts_outliers<- function(list.cleaned, outlier_threshold=5,  unreliable_outlier_percentage=2) {


  mat1 <- list.cleaned[[2]]
  mat2 <- list.cleaned[[3]]

  list.mat <- list(mat1, mat2)
  # Adjust this threshold as needed

  ##########

  # FIRST PART: cleaning CpGs that have still high percentage of missing
  m2 <- list()

  cat("Excluding CpGs that still have too many NAs from in methylated and unmethylated copunts' matrices\n")
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

    cat(paste(dim(throw)[1], " rows have been removed because all entries were NAs.\n"))


    clean.met <- X2

    m2[[k]] <- X2

  }


  #######

  # SECOND PART: cleaning outliers

  cat("Converting outliers into NAs ...\n")

  list.mat2 <- m2
  unreliable_outlier_percentage <- 2  # Adjust this threshold as needed

  m3 <- list()

  for (k in 1:length(list.mat2)){

    m <- list.mat2[[k]]
    data_with_na <- m

    for (i in 1:nrow(m)) {

      row_values <- m[i,]


      # Calculate the mean and standard deviation of the column
      mean_value <- mean(row_values, na.rm=T)
      sd_value <- sd(row_values, na.rm=T)

      # The threshold for outliers is 5 by default (discarded values beyond 5 standard deviations)

      # Calculate the percentage of rows with values beyond the threshold
      outlier_percentage <- mean(abs(row_values - mean_value) > outlier_threshold * sd_value, na.rm = T) * 100


      for (j in 1:ncol(m)) {

        if(!is.na(m[i,j])){

          if(abs(m[i,j] - mean_value) > outlier_threshold * sd_value) {


            if(outlier_percentage <= unreliable_outlier_percentage) {

              m[i,j] <- NA
            }    #  3 if

          }   #  2 if
        }    #  1 if


      }

    }

    m3[[k]] <- data_with_na
  }


  clean.met<- m3[[1]]
  clean.unmet <- m3[[2]]

  clean.met2 <- clean.met[rownames(clean.met) %in% rownames(clean.unmet) ,]
  clean.unmet2 <-  clean.unmet[rownames(clean.unmet) %in% rownames(clean.met2) ,]


  #dim(clean.unmet)

  clean.beta <- list.cleaned[[4]]
  clean.beta <- clean.beta[rownames(clean.beta) %in% rownames(clean.met2),]
  #dim(clean.beta)

  clean.m <- list.cleaned[[5]]
  clean.m <- clean.m[rownames(clean.m) %in% rownames(clean.met2),]
  #dim(clean.m)

  clean.cov <- list.cleaned[[1]]
  clean.cov <- clean.cov[rownames(clean.cov) %in% rownames(clean.met2),]

  list.cleaned2 <- list(Coverage_matrix2 = clean.cov, Met_matrix2 = clean.met2,
                        Unmet_matrix2 = clean.unmet2,
                        Beta_matrix2 = clean.beta, M_matrix2 = clean.m)






  return(list.cleaned2)

}




