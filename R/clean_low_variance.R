#' @title Adjusting methylated/unmethylated counts matrix or beta value matrix after cleaning coverage matrix
#'
#' @description
#' \code{\link{clean_low_variance}} helps in the identification and throw out CpGs with very low variance due to missing values in the methylated counts
#'
#'
#' @param X Cleaned matrix list from coverage cleaning (clean_matrix)
#' @param sd_quantile is the threshold for throwing out CpGs not passing this threshold in sd quantile distribution
#' @param max_na in the maximum proportion we want to include in our dataset after
#' @param varselect is the index of the dataset to be used in the output from clean_matrix (numeric value 1-5)
#'
#' @name clean_low_variance
#'
#' @return
#' #' list with 7 datasets:
#'  \item{X3}{cleaned matrix according to max_na and sd_quantile parameters}
#
#'
#'
#'@examples
#'
#'
#'
#' @export
#'
#'
clean_low_variance <- function(X, sd_quantile=0.05, max_na = 0.25, varselect=2) {

  library(matrixStats)
  library(impute)
  n <- nrow(X[[varselect]])

  X <- X[[varselect]]
  # Convert the object to a matrix (if needed)
  is_matrix <- is.matrix(X)
  if (!is_matrix) {
    X <- as.matrix(X)
  }

  row_sd = rowSds(X, na.rm = T)
  l = which(abs(row_sd) <= 1e-10)
  X2<- X[-l,]
  filtered1<- rownames(X2)

  #X2 = X[!l, , drop = T]
  #row_sd = row_sd[!l]

  l <- as.matrix(l)

  #dim(l)[1]

  if(dim(l)[1] !=0) {
    cat(paste(dim(l)[1], " rows have been removed with zero variance.\n"))
  }

  qa = quantile(unique(row_sd), sd_quantile, na.rm = TRUE)
  l2 = which(row_sd < qa)
  l2 <- as.matrix(l2)

  if(dim(l2)[1] !=0) {
    cat(paste(dim(l2)[1], " rows have been removed with too low variance (sd <= ", sd_quantile, " quantile).\n"))
  }


  X3 <- X2[-l2,]
  filtered2<- rownames(X3)

  throw<- which(is.na(rowMeans(X3, na.rm = T)))
  X4 <- X3[-throw,]

  throw <- as.matrix(throw)
  if(dim(throw)[1] !=0) {
    cat(paste(dim(throw)[1], " rows have been removed because all entries were NAs.\n"))
   }

  #rownames(X3) <- sites_filtered2
  return(X4)

}
