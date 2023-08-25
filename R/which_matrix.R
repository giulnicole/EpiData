#' @title Calculating total amount of missing values in the coverage, methylated and unmethylated counts matrices
#'
#' @description
#' \code{\link{which_matrix}} computes the sum of missing values per dataset, highlighting which matrix should be use to compute partitions
#'
#' @param mat1 matrix of coverage counts after cleaning coverage
#' @param mat2 matrix of methylated counts after cleaning coverage
#' @param mat3 matrix of unmethylated counts after cleaning coverage
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
which_matrix <- function(mat1, mat2, mat3) {

  # Calculate missing value counts
  na_count_mat1 <- colSums(is.na(mat1))
  na_count_mat2 <- colSums(is.na(mat2))
  na_count_mat3 <- colSums(is.na(mat3))

  tot <- dim(mat1)[1]*dim(mat1)[2]


  # Compare missing value counts
  na_summary <- data.frame(
    Dataset = c("Coverage_counts", "Methylated_counts", "Unmethylated_counts"),
    TotalMissingValues = c(sum(na_count_mat1), sum(na_count_mat2), sum(na_count_mat3)),
    PercentageNA = c(sum(na_count_mat1), sum(na_count_mat2), sum(na_count_mat3))/tot)



  # Or create a bar plot to visualize missing value counts
  library(ggplot2)

  na_plot <- ggplot(na_summary, aes(x = Dataset, y = PercentageNA, fill = Dataset)) +
    geom_bar(stat = "identity") +
    labs(title = "Missing Value Comparison", y = "Pecrentage Missing Values") +
    theme_minimal() + scale_fill_manual(values = c("Coverage_counts" = "#fdbf6f", "Methylated_counts" = "#c5b0d5", "Unmethylated_counts" = "#a1d99b"))

  res <- list(NA_summary = na_summary, NA_plot = na_plot)


  return(res)





}



