#' @title cleanOutliers
#' @description
#' \code{\link{cleanOutliers}} helps in discarding those CpGs with high rate of missing values of methylated and unmethylated counts after
#' coverage cleaning and outlier values in the conversion.
#'
#'
#' @param filtered.obj SummarizedExperiment object deriving from cleanMatCov.
#' @param outlier_threshold Threshold for considering a value of the counts' matrices as an outlier.
#' @param remove_outliers TRUE/FALSE whether removing of outliers should be performed or not.
#'
#' @name cleanOutliers
#'
#' @return
#' SummarizedExperiment object with 2 elements:
#'  \item{Output_filtered}{SummarizedExperiment object with cleaned matrices (coverage, methylated and unmethylated counts) after filtering for outliers and Beta and M values matrices}
#'  \item{Plots}{Barplot with missing values rates after second cleaning}
#'
#'
#' @examples
#' # data("matrices")
#' # clean.coverage2 <- cleanCovMat(input.obj=dati, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)
#' # clean.out <- cleanOutliers(filtered.obj=clean.coverage2, outlier_threshold=5, remove_outliers = TRUE)
#'
#'
#' @export
#'
#'
cleanOutliers<- function(filtered.obj, outlier_threshold,  remove_outliers){

  mat1<- filtered.obj$Output_filtered@assays@data@listData$Met_matrix
  mat2<- filtered.obj$Output_filtered@assays@data@listData$Unmet_matrix

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

  # THIRD PART: calculating ratios

  # Calculating beta values matrix and M values on cleaned data
  cat("Calculating beta values and M values ...", "\n")
  p <- 1  # offset parameter


  # B values' matrix
  B <- matrix(0, nrow(clean.met2), ncol(clean.met2))

  # M values' matrix
  M <- matrix(0, nrow(clean.met2), ncol(clean.met2))

  for (j in 1:nrow(clean.met2)){

    for (i in 1:ncol(clean.met2)) {

      meth <- as.numeric(clean.met2[[j,i]])
      unmeth <- as.numeric(clean.unmet2[[j,i]])

      beta <- max(meth, 0)/(max(meth,0) + max(unmeth,0) + p)
      B[j,i] <- beta

      m <- log2((max (meth,0) + p)/ (max(unmeth,0) +p))
      M[j,i] <- m

    } # for i

  } # for j




  colnames(B)<- colnames(clean.met2)
  rownames(B)<- rownames(clean.met2)
  cat("Beta matrix computed ...", "\n")

  colnames(M)<- colnames(clean.met2)
  rownames(M)<- rownames(clean.met2)
  cat("M matrix computed ...", "\n")



  #list.cleaned[[4]] <- B
  #list.cleaned[[5]] <- M



  #dim(clean.m)

  clean.cov <- filtered.obj$Output_filtered@assays@data@listData$Coverage_matrix
  clean.cov <- as.data.frame(clean.cov[rownames(clean.cov) %in% rownames(clean.met2),])


  list.cleaned2 <- list(Coverage_matrix = as.data.frame(clean.cov),
                        Met_matrix = as.data.frame(clean.met2),
                        Unmet_matrix = as.data.frame(clean.unmet2),
                        Beta_matrix = as.data.frame(B),
                        M_matrix = as.data.frame(M))



  plots<- whichMatrix2(list.cleaned2)

  filtered.obj2 <-  filtered.obj
  filtered.obj2$Output_filtered@assays@data@listData <- list.cleaned2

  results<- list(Output_outliers = filtered.obj2$Output_filtered, Plots= plots)


  # objGR<- GRconversion(cleaned.list = results$Cleaned2)
  return(results)


}  # (main function)



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


} # function 1







