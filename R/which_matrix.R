#' @title Calculating total amount of missing values in the methylated counts and unmethylated counts matrices
#'
#' @description
#' \code{\link{which_matrix}} computes the sum of missing values per dataset, highlighting which matrix should be use to compute partitions
#'
#' @param mat1 matrix of methylated counts after cleaning coverage
#' @param mat2 matrix of unmethylated counts after cleaning coverage
#'
#'
#' @name which_matrix
#'
#' @return
#' a list with two elements
#'  \item{NA_summary}{summary of the missing values in the two matrices}
#'  \item{NA_plot}{barplot of total amount of missing values}
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



which_matrix <- function(mat1=clean.coverage2[[2]], mat2=clean.coverage2[[3]]) {

  # Calculate missing value counts
  na_count_mat1 <- colSums(is.na(mat1))
  na_count_mat2 <- colSums(is.na(mat2))
  tot <- dim(mat1)[1]*dim(mat1)[2]

  # Compare missing value counts
  na_summary <- data.frame(
    Dataset = c("Methylated_counts", "Unmethylated_counts"),
    TotalMissingValues = as.numeric(c(sum(na_count_mat1), sum(na_count_mat2))),
    PercentageMissingValues = as.numeric(c(sum(na_count_mat1), sum(na_count_mat2)))/tot
  )


  # Or create a bar plot to visualize missing value counts
  library(ggplot2)

  na_plot <- ggplot(na_summary, aes(x = Dataset, y = PercentageMissingValues, fill = Dataset)) +
    geom_bar(stat = "identity") +
    labs(title = "Missing Value Comparison", y = "Pecrentage Missing Values") +
    theme_minimal() + scale_fill_manual(values = c("Methylated_counts" = "#800080", "Unmethylated_counts" = "#E8C32E"))

  res <- list(NA_summary = na_summary, NA_plot = na_plot)


  return(res)





}
