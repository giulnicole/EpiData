

#' @noRd
#'
whichMatrix2 <- function(list.cleaned) {

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


