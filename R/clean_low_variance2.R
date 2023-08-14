#' @title Adjusting methylated counts matrix for low variance in the dataset
#' subsequently subsetting unmethylated counts, beta and M values
#'
#' @description
#' \code{\link{clean_low_variance2}} helps in the identification and throw out CpGs with very low variance due to missing values in the methylated counts
#'
#'
#' @param X Cleaned matrix list from coverage cleaning (clean_matrix)
#' @param sd_quantile is the threshold for throwing out CpGs not passing this threshold in sd quantile distribution
#' @param max_na in the maximum proportion we want to include in our dataset after
#' @param varselect is the index of the dataset to be used in the output from clean_matrix (numeric value default = 2, because we want to clean generally the meth counts)
#'
#' @name clean_low_variance2
#'
#' @return
#' #' list with 4 datasets:
#'  \item{list.cleaned2}{cleaned matrices list according to max_na and sd_quantile parameters, containing cleaned methylated counts, unmethylated counts, beta and M values matrices}
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

clean_low_variance2 <- function(list.cleaned, sd_quantile=0.05, max_na = 0.25, varselect=2) {

  n <- nrow(list.cleaned[[varselect]])

  # FIRST PART: cleaning low variance for meth counts

  X <- list.cleaned[[varselect]]

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
  #return(X4)

  clean.met <- X4


  # SECOND PART: Subsetting


  # Subset the same unmet and beta and m
  clean.unmet <- list.cleaned[[3]]
  clean.unmet <- clean.unmet[rownames(clean.unmet) %in% rownames(clean.met),]
  #dim(clean.unmet)

  clean.beta <- list.cleaned[[4]]
  clean.beta <- clean.beta[rownames(clean.beta) %in% rownames(clean.met),]
  #dim(clean.beta)

  clean.m <- list.cleaned[[5]]
  clean.m <- clean.m[rownames(clean.m) %in% rownames(clean.met),]
  #dim(clean.m)

  clean.cov <- list.cleaned[[1]]
  clean.cov <- clean.cov[rownames(clean.cov) %in% rownames(clean.met),]

  list.cleaned2 <- list(Coverage_matrix2 = clean.cov, Met_matrix2 = clean.met,
                        Unmet_matrix2 = clean.unmet,
                        Beta_matrix2 = clean.beta, M_matrix2 = clean.m)




  return(list.cleaned2)


}
