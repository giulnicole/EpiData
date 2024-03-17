#' @title consensusPart
#'
#' @description A function for computing the non-negative definite matrix of covariance matrix, by eigenvectors' decompisition and
#' and simulate new values drawn from the same distribution to impute remaining missing values in the original cleaned dataset.
#' \code{} computes the imputation based on correlation pattern between CpGs derived from eigen decomposition of the covariance matrix.
#'
#' @import cola
#'
#' @name consensusPart
#'
#' @param data SummarizedExperimet object cleaned and after statistics and NAs' pattern computation (statsCpG output).
#' @param dataset.simulated Simulated dataset with impute.corr.
#' @param partition_method Method to perform partition (see cola package for further details).
#' @param top_rows Method to use for finding top rows. (see cola package for further details).
#' @param k_to_test Number of partitions to perform to divide CpGs. Default is 2:8.
#' @param top_n Top rows (number) on which perform partition. Default value is 100.
#'
#'
#'
#' @return
#' dataset imputed:
#' \item{}{Matrix of imputed CpGs via eigenvalue decomposition of covariance matrix and generation of values from a distribution similar per each CpG.}
#'
#'
#' @examples
#' \dontrun{

#'  data("matrices")
#'  clean.coverage2 <- cleanCovMat(input.obj=dati, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)
#'  clean.out <- cleanOutliers(filtered.obj=clean.coverage2, outlier_threshold=5, remove_outliers = TRUE)
#'  stats <- statsCpG(cleaned.obj=clean.out, varselect = 5)
#'  simu3 <- imputeCorr(cleaned.obj=stats, matrix="M", varselect=5)
#'  cpart <- consensusPart(data = clean.out$M_matrix, dataset.simulated = simu3$Simulated, partition_method = "pam", top_rows = "ATC", k_to_test = 2:8, top_n = 100)
#'
#' }
#'
#'
#'
#'
#' @export
consensusPart <- function(data, dataset.simulated,
                          partition_method = "pam", top_rows = "ATC",
                          k_to_test = 2:8, top_n = 100){


  # list.raw for calculating the means per CpG and per partition (when the partition will be evaluated)
  # dataset.simulated for calculating the consensus partitions

  dataset.raw <- data

  # computing similarity matrix (default is euclidean distance)
  sim.euc <- negDistMat(dataset.simulated)

  cat("Partitioning\n")
  # pam (skmeans) for more stable partitions
  # users can provide SD or ATC method for top rows
  sim.clust <- consensus_partition(sim.euc, k=k_to_test, top_value_method = top_rows,
                                   partition_method= partition_method, top_n = top_n)

  cat("Getting results from partitioning\n")
  best.k <- as.vector(suggest_best_k(sim.clust))
  stats <- as.data.frame(get_stats(sim.clust))
  membership <- as.data.frame(cbind(rownames(dataset.raw), get_membership(sim.clust, k=best.k)))
  classes <- get_classes(sim.clust)



  cat("Computing KNN per cluster to the missing values\n")
  colnames(classes) <- gsub("k=", "",  colnames(classes))
  index<- which(colnames(classes) == best.k)
  df<- cbind(classes[,index], dataset.raw)
  colnames(df)[1] <- "class"
  df <- as.data.frame(df)
  df <-  rownames_to_column(df, var = "ID")

  # Knn on subgroups

  subgroups_df <- as_tibble(df) %>%
    group_split(class)

  subgroups_imputed <- list()

  for (i in seq_along(subgroups_df)) {

    assign(paste0("df_class_", i), as.data.frame(subgroups_df[[i]]))
    sub <- as.data.frame(subgroups_df[[i]])
    #rownames(sub) <- rownames(df)
    sub <- as.matrix(sub[,-c(1,2)])
    sub.imp <- impute.knn(sub, colmax = 1)
    sub.imp2 <- as.data.frame(sub.imp$data)


    # Convert all character columns to numeric
    sub.imp2<- as.matrix.data.frame(sub.imp2)
    subgroups_imputed[[i]] <- sub.imp2
    rownames(subgroups_imputed[[i]]) <- subgroups_df[[i]]$ID


  }


  processed_df <- list()


  for (i in 1:length(subgroups_imputed)) {
    # Creare un nuovo data frame e aggiungerlo alla lista
    dfs <- as.data.frame(subgroups_imputed[[i]])
    dfs <- rownames_to_column(dfs, var = "ID")
    processed_df[[i]] <- dfs
  }


  # Unire tutti i data frame in processed_dfs in un unico data frame
  row_joined_df <- do.call(bind_rows, processed_df)

  ordered_merged_df <- row_joined_df %>%
    arrange(ID)


  #ordered_merged_df$ID == df$ID
  ordered_merged_df <- ordered_merged_df[,-1]
  dataset.cluster <- as.data.frame(ordered_merged_df)

  res.part.knn <- list(best_k=best.k, dataset.cluster.knn = dataset.cluster,
                            membership_in_clusters = membership, partition_statistics= stats,
                            partition_classes = classes)


  return(res.part.knn)

}

