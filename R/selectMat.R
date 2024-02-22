#' @title selectMat
#' @description
#' \code{\link{selectMat}} internal function which helps in visualizing the percentage of missing values
#'
#'
#' @param list.cleaned list containing the cleaned datasets which have been cleaned.
#'
#'
#' @name selectMat
#'
#' @return
#' a list with 2 elements:
#'  \item{NA_summary}{summary statistics of missing data}
#'  \item{NA_plot}{barplot for percentages of missing data}
#'
#'
#' @examples
#' # data("matrices")
#' # C<- assay(dati, 1)
#' # M<- assay(dati, 2)
#' # U<- assay(dati, 3)
#' # C[C==0]<- NA
#' # M[M==0]<- NA
#' # U[U==0]<- NA
#' # list.cleaned <- list(coverage=C, meth=M, unmeth=U)
#' # statistics <- selectMat(list.cleaned)
#'
#' @noRd
#'
selectMat <- function(list.cleaned) {


  mat1  <- list.cleaned[[1]]
  mat2 <- list.cleaned[[2]]
  mat3 <-  list.cleaned[[3]]

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
