#' @title errors
#' @description
#' \code{} internal function which perform the accuracy measure with RMSE, MAE
#'
#'
#' @param origmat1 filtered matrix from the SummarizedExperiment object but without missing
#' @param masked_data.mat1 matrix which has been imputed
#' @param na_positiond position in the matrix with replicate NAs
#'
#' @name errors
#'
#' @return
#' a vector with errors per each iteration:
#'  \item{S3}{Root mean square error (RMSE)}
#'  \item{M3}{Mean absolute error (MAE).}
#'
#'
#' @examples
#'
#' \dontrun{
#'  err <- errors(origmat1, masked_data.mat1, na_positions)
#'
#' }
#'
#'
#' @noRd
#'
#'
errors <- function(origmat1, masked_data.mat1, na_positions) {


  # Calculate RMSE, MAE, KL distance, KL divergence for mat1
  error.rmse.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
  error.mae.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))


  for (j in 1:nrow(origmat1)){

    for (l in 1:ncol(origmat1)) {

      if (na_positions[j,l] == TRUE){
        # computing RMSE
        diff3 <- as.numeric(origmat1[j,l] - masked_data.mat1[j,l])
        sqdiff3 <- (diff3)^2
        error.rmse.mat1[j,l] <- sqdiff3

        #computing MAE
        abs3 <- abs(as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))
        error.mae.mat1[j,l] <- abs3


        } # if


      } # for l


    } # for j

  x3 = mean(error.rmse.mat1)
  S3 = sqrt(x3)/(max(origmat1)-min(origmat1))

  M3 = mean(error.mae.mat1)/(max(origmat1)-min(origmat1))

  cat("RMSE at iteration:\n")
  print(S3)
  cat("\n")

  cat("MAE at iteration:\n")
  print(M3)
  cat("\n")


  errors = c(S3, M3)

  return(errors)




}





