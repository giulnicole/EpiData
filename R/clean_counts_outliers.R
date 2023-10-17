#' @title clean_counts_outliers
#' @description
#' \code{\link{clean_counts_outliers}} helps in discarding those CpGs with high rate of missing values of methylated and unmethylated counts after
#' coverage cleaning and outlier values in the conversion
#'
#'
#' @param list.cleaned list deriving from clean_matrix_cov
#' @param outlier_threshold threshold for considering a value of the counts an outlier
#' @param unreliable_outlier_percentage that we want to tolerate (also considering the sample size we have)
#'
#' @name clean_counts_outliers
#'
#' @return
#' #' list with 2 elements:
#'  \item{Cleaned2}{list with cleaned matrices: coverage matrix, methylated counts matrix, unmethylated counts matrix, cleaned at the specified missingness thresholds removed}
#'  \item{Plots}{Barplot with missing values rates after second cleaning}
#'
#'
#'
#'
#'
#'
#' @export
#'
#'
which_matrix2 <- function(list.cleaned) {

  mat1  <- list.cleaned[[1]]
  mat2 <- list.cleaned[[2]]
  mat3 <-  list.cleaned[[3]]
  mat4<-  list.cleaned[[4]]
  mat5 <- list.cleaned[[5]]

  # Calculate missing value counts
  na_count_mat1 <- colSums(is.na(mat1))
  na_count_mat2 <- colSums(is.na(mat2))
  na_count_mat3 <- colSums(is.na(mat3))
  na_count_mat4 <- colSums(is.na(mat4))
  na_count_mat5 <- colSums(is.na(mat5))


  tot <- dim(mat1)[1]*dim(mat1)[2]


  # Compare missing value counts
  na_summary <- data.frame(
    Dataset = c("Coverage_counts", "Methylated_counts", "Unmethylated_counts",
                "Beta_values", "M_values"),
    TotalMissingValues = c(sum(na_count_mat1), sum(na_count_mat2), sum(na_count_mat3),
                           sum(na_count_mat4), sum(na_count_mat5)),
    PercentageNA = c(sum(na_count_mat1), sum(na_count_mat2), sum(na_count_mat3),
                     sum(na_count_mat4), sum(na_count_mat5))/tot)



  # Or create a bar plot to visualize missing value counts
  library(ggplot2)

  na_plot <- ggplot(na_summary, aes(x = Dataset, y = PercentageNA, fill = Dataset)) +
    geom_bar(stat = "identity") +
    labs(title = "Missing Value Comparison", y = "Pecrentage Missing Values") +
    theme_minimal() + scale_fill_manual(values = c("Coverage_counts" = "#fdbf6f",
                                                   "Methylated_counts" = "#c5b0d5",
                                                   "Unmethylated_counts" = "#a1d99b",
                                                   "Beta_values" = "#FFFF66", "M_values"= "#B0E0E6"))

  res <- list(NA_summary = na_summary, NA_plot = na_plot)


  return(res)


}


clean_counts_outliers<- function(list.cleaned, outlier_threshold=5,  unreliable_outlier_percentage=2){


  mat1 <- list.cleaned[[2]]
  mat2 <- list.cleaned[[3]]

  list.mat <- list(mat1, mat2)
  # Adjust this threshold as needed

  ##########

  # FIRST PART: cleaning CpGs that have still high percentage of missing
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

    cat(paste(dim(throw)[1], " rows have been removed because all entries were NAs.\n"))


    m2[[k]] <- X2

  }


  #######

  # SECOND PART: cleaning outliers

  cat("Converting outliers into NAs ...\n")

  list.mat2 <- m2

  m3 <- list()


  for (k in 1:length(list.mat2)){

    m <- list.mat2[[k]]

    remove_outlier <- function(row) {
      mean_row <- mean(row)
      sd_row <- sd(row)
      is_outlier <- abs(row - mean_row) > 5 * sd_row
      row[is_outlier] <- NA  # outliers as NA
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


  # THIRD PART: calculating ratios

  # Calculating beta values matrix and M values on cleaned data
  cat("Calculating beta values and M values ...", "\n")
  p <- 100  # offset parameter


  # B values' matrix
  B <- matrix(0, nrow(clean.met2), ncol(clean.met2))

  for (j in 1:nrow(clean.met2)){

    for (i in 1:ncol(clean.met2)) {

      meth <- as.numeric(clean.met2[[j,i]])
      unmeth <- as.numeric(clean.unmet2[[j,i]])

      beta <- max(meth, 0)/(max(meth,0) + max(unmeth,0) + p)
      B[j,i] <- beta


    } # for i

  } # for j


  B <- as.data.frame(B)
  colnames(B)<- colnames(clean.met2)
  rownames(B)<- rownames(clean.met2)
  cat("Beta matrix computed ...", "\n")

  # M values' matrix
  M <- matrix(0, nrow(clean.met2), ncol(clean.met2))

  for (j in 1:nrow(clean.met2)){

    for (i in 1:ncol(clean.met2)) {


      meth <- as.numeric(clean.met2[[j,i]])
      unmeth <-  unmeth <- as.numeric(clean.unmet2[[j,i]])

      m <- log2((max (meth,0) + p)/ (max(unmeth,0) +p))
      M[j,i] <- m
      M[j,i] <- m


    } # for i

  } # for j

  M <- as.data.frame(M)

  colnames(M)<- colnames(clean.met2)
  rownames(M)<- rownames(clean.met2)
  cat("M matrix computed ...", "\n")


  list.cleaned[[4]] <- B
  list.cleaned[[5]] <- M



  #dim(clean.m)

  clean.cov <- list.cleaned[[1]]
  clean.cov <- clean.cov[rownames(clean.cov) %in% rownames(clean.met2),]
  list.cleaned[[1]] <- clean.cov


  list.cleaned2 <- list(Coverage_matrix2 = list.cleaned[[1]],
                        Met_matrix2 = list.cleaned[[2]],
                        Unmet_matrix2 = list.cleaned[[3]],
                        Beta_matrix2 = as.matrix.data.frame(list.cleaned[[4]]),
                        M_matrix2 = as.matrix.data.frame(list.cleaned[[5]]))



  plots<- which_matrix2(list.cleaned2)


  results<- list(Cleaned2 = list.cleaned2, Plots= plots)

  return(results)

}

