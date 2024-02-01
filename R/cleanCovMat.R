#' @title clenCovMat
#' @description
#' \code{\link{cleanCovMat}} helps in the conversion of missing values (0->NA) in the coverage matrix,
#'  removing rows and columns above pre-specified missingness threshold.
#'
#'
#' @param input.obj SummarizedExperiment object as input with list with coverage, methylated counts and unmethylated counts in assays/data/listData (rows are CpGs and individuals are columns).
#' @param max_na_cpg Threshold of missing values per each CpG in the coverage counts matrix.
#' @param max_na_ind Threshold of missing values per each individual.
#' @param cpg_removal_threshold Minimum threshold of coverage counts (across all individuals) for a CpG to be kept.
#'
#' @name cleanCovMat
#'
#' @return
#'  SummarizedExperiment object with 2 elements:
#'  \item{Output_filtered}{SummarizedExperiment object with cleaned matrices: coverage matrix, methylated counts matrix, unmethylated counts matrix, cleaned at the specified missingness thresholds removed}
#'  \item{Plots}{Barplot with missing values rates after first cleaning}
#'
#' @examples
#' data("matrices")
#' #clean.coverage2 <- cleanCovMat(input.obj=dati, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)
#'
#'
#' @export
#'
#'
cleanCovMat <- function(input.obj, max_na_cpg=0.5, max_na_ind=0.2, cpg_removal_threshold=10) {

  X<- input.obj@assays@data@listData$Coverage_matrix
  Y<- input.obj@assays@data@listData$Met_matrix
  Z<- input.obj@assays@data@listData$Unmet_matrix

  sites <- rownames(X)
  stopifnot("Input must be numeric dataframe" =is.data.frame(X), all(sapply(X, is.numeric)))
  cat("Converting 0s in NAs for coverage counts ...", "\n")
  X <- X %>% mutate_all(~na_if(., 0))

  # Convert the object to a matrix (if needed)
  is_matrix <- is.matrix(X)
  if (!is_matrix) {
    X <- as.matrix(X)
  }

  # Calculate the maximum allowed missing values for rows and columns
  max_missing_rows <- ncol(X) * max_na_cpg
  max_missing_cols <- nrow(X) * max_na_ind

  # Discard CpGs with > max_na_cpg % missing values
  row_indices_to_keep <- rowSums(is.na(X)) <= max_missing_rows

  X2 <- X[row_indices_to_keep, ]
  dim(X2)

  eliminated <- nrow(X) - nrow(X2)
  index1 <- which(row_indices_to_keep==TRUE)
  sites_filtered1 <- sites[index1]
  rownames(X2) <- sites_filtered1

  if (eliminated != 0) {
    message(paste(eliminated, " CpG(s) removed due to exceeding the pre-defined removal threshold (>",
                  max_na_cpg * 100, "%) for missingness.", sep = "")) }


  # Discard individuals with > max_na_ind % missing values
  col_indices_to_keep <- colSums(is.na(X2)) <= max_missing_cols
  X3 <- X2[, col_indices_to_keep]
  eliminated2 <- ncol(X2) - ncol(X3)
  index2 <- which(col_indices_to_keep==TRUE)


  if (eliminated2 != 0) {
    message(paste(eliminated2, " individual(s) removed due to exceeding the pre-defined removal threshold (>",
                  max_na_ind * 100, "%) for missingness.", sep = "")) }

  # Remove counts below 10 each
  index_above_threshold <- which(apply(X3, 1, median, na.rm=T)  > cpg_removal_threshold)
  X4 <- X3[index_above_threshold, ]
  eliminated3 <- nrow(X3) - nrow(X4)

  sites_filtered2 <- sites_filtered1[index_above_threshold]
  rownames(X4) <- sites_filtered2


  if (eliminated3 != 0) {
    message(paste(eliminated3, " Cpg(s) removed due to not exceeding the pre-defined removal threshold (<",
                  cpg_removal_threshold, ") for counts.", sep = "")) }



  # Filter methylated counts and unmethylated count matrices
  cat("Computing selection of CpGs in methylated counts matrix ...", "\n")

  # Return meth and unmeth counts
  Y <- Y[rownames(Y) %in% rownames(X4),]
  Y <- Y[, colnames(Y) %in% colnames(X4)]
  #Y <- Y %>% mutate_all(~na_if(., 0))

  cat("Computing selection of CpGs in unmethylated counts matrix ...", "\n")
  Z <- Z[rownames(Z) %in% rownames(X4),]
  Z <- Z[, colnames(Z) %in% colnames(X4)]
  #Z <- Z %>% mutate_all(~na_if(., 0))

  cat("Adjusting the NAs in methylated and unmethylated counts' matrices ...", "\n")
  # Adjusting the zeros in methylated and unmethylated counts
  res <- adjustCounts(coverage = X4, methylated = Y, unmethylated = Z)

  Y <- res[[1]]
  Z <- res[[2]]
  X4 <- as.data.frame(X4)
  rownames(X4)<- sites_filtered2

  # Results
  res<- list(Coverage_matrix = X4, Met_matrix = Y, Unmet_matrix = Z)

  cat("Coverting results to matrices ...", "\n")
  res2<- lapply(res, as.matrix.data.frame)

  # Plots
  plots<- whichMatrix(res2)

  obj<- GRconversion2(res2)

  objGR<- list(Output_filtered = obj, Plots= plots)

  #objGR<- GRconversion2(cleaned.list = results$Cleaned1)
  return(objGR)

}  # (main function)




#' @noRd
#'
whichMatrix <- function(list.cleaned) {


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


}  # function 1



#' @noRd
#'
adjustCounts <- function(coverage, methylated, unmethylated){


  # NA
  # Create a new matrix based on conditions

  # Initialize
  met2 <- methylated
  co <- coverage
  un <- unmethylated

  for (i in 1:nrow(met2)) {
    for (j in 1:ncol(met2)) {
      if (is.na(co[i, j]) && un[i, j] == 0) {
        met2[i, j] <- NA
      } else {
        met2[i, j] <-  methylated[i, j]
      }
    }
  }


  # Initialize
  unmet2 <- unmethylated

  for (i in 1:nrow(unmet2)) {
    for (j in 1:ncol(unmet2)) {
      if (is.na(co[i, j]) && is.na(met2[i, j])) {
        unmet2[i, j] <- NA
      } else {
        unmet2[i, j] <-  unmethylated[i, j]
      }
    }
  }


  # Initialize
  met3 <- met2

  for (i in 1:nrow(met3)) {

    for (j in 1:ncol(met3)) {

      if (!is.na(co[i, j]) && co[i, j] == unmet2[i, j]) {
        #if (met3[i, j] == 0) {

        cpg <- as.numeric(met3[i,])
        wil <- wilcox.test(cpg, mu = 0, alternative = "greater")
        pval <- wil$p.value


        if (pval< 0.01) {
          # cat("NA \n")
          met3[i, j] <- NA


        }   else {
          # cat("not NA \n")
          met3[i, j] <- 0

        }


      } else {
        met3[i, j] <-  met2[i, j]
      }
    }
  }


  # Initialize
  unmet3 <- unmet2

  for (i in 1:nrow(unmet3)) {

    for (j in 1:ncol(unmet3)) {

      if (!is.na(co[i, j]) && !is.na(met3[i, j]) && co[i, j] == met3[i, j]) {
        #if (met3[i, j] == 0) {

        cpg <- as.numeric(unmet3[i,])
        wil <- wilcox.test(cpg, mu = 0, alternative = "greater")
        pval <- wil$p.value


        if (pval< 0.01) {
          # cat("NA \n")
          unmet3[i, j] <- NA


        }   else {
          # cat("not NA \n")
          unmet3[i, j] <- 0

        }


      } else {
        unmet3[i, j] <-  unmet2[i, j]
      }
    }
  }

  methylated <- met3
  unmethylated <- unmet3

  res <- list(methylated, unmethylated)

  return(res)



}  # function 2


#' @noRd
#'
GRconversion2<- function(cleaned.list){


  nrows <- dim(cleaned.list[["Coverage_matrix"]])[1]
  ncols <- dim(cleaned.list[["Coverage_matrix"]])[2]
  counts <- cleaned.list[["Coverage_matrix"]]

  names_cpg<-rownames(counts)

  chr<- gsub("-.*", "", rownames(counts))
  chr<- unique(chr)
  freq <- length(chr)
  positions <-  gsub(".*-", "", rownames(counts))

  rowRanges <- GRanges(rep(chr,freq),
                       IRanges(positions),
                       feature_id=names_cpg)



  colData <- DataFrame(Matrix=rep("coverage", ncol(counts)))
  object <- SummarizedExperiment(assays=cleaned.list,
                                 rowRanges=rowRanges, colData=colData)
  return(object)



}  # function 3

