#' @title Data matrix cleaning for missing values after RRBS
#' and computation of beta and M coefficients
#'
#' @description
#' \code{\link{clean_matrix}} helps in the conversion of missing values (0->NA),
#' variable types and removes rows and columns above pre-specified missingness threshold
#'
#'
#' @param X Original Coverage dataframe with individuals in columns and CpGs in rows
#' @param Y Original Methylated counts dataframe with individuals in columns and CpGs in rows
#' @param Z Original Unmethylated counts dataframe with individuals in rows and CpG as columns
#' @param max_na_cpg threshold of missing values per each CpG in the coverage counts matrix
#' @param max_na_ind threshold of missing values per each individual
#' @param cpg_removal_threshold minimum threshold of coverage counts (across all individuals) for a CpG to be kept
#' @param p offset parameter to compute regularizatio of Beta matrix and M matrix
#'
#' @name clean_matrix
#'
#' @return
#' #' list with 7 datasets:
#'  \item{Coverage_matrix}{cleaned coverage dataset with NAs as missing values and rows/columns above the pre-specified missingness thresholds removed}
#'  \item{Met_matrix }{methylated counts dataset with the filtered CpGs}
#'  \item{Unmet_matrix}{unmethylated counts dataset with the filtered CpGs}
#'  \item{Beta_matrix}{beta values matrix}
#'  \item{M_matrix}{M values matrix}
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
clean_matrix <- function(X, Y, Z, max_na_cpg=0.5, max_na_ind=0.2,
                         cpg_removal_threshold=10, p=100) {


  sites <- rownames(X)

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
    Y <- Y %>% mutate_all(~na_if(., 0))

    cat("Computing selection of CpGs in unmethylated counts matrix ...", "\n")
    Z <- Z[rownames(Z) %in% rownames(X4),]
    Z <- Z[, colnames(Z) %in% colnames(X4)]
    Z <- Z %>% mutate_all(~na_if(., 0))

    # Calculating beta values matrix and M values on cleaned data
    cat("Calculating beta values and M values ...", "\n")
    #p <- 100  # offset parameter

    # B values' matrix
    B <- matrix(0, nrow(X4), ncol(X4))

    for (j in 1:nrow(X4)){

      for (i in 1:ncol(X4)) {

        meth <- as.numeric(Y[[j,i]])
        unmeth <- as.numeric(Z[[j,i]])

        beta <- max(meth, 0)/(max(meth,0) + max(unmeth,0) + p)
        B[j,i] <- beta


      } # for i

    } # for j


    B <- as.data.frame(B)
    colnames(B)<- colnames(X4)
    rownames(B)<- rownames(X4)
    cat("Beta matrix computed ...", "\n")

    # M values' matrix
    M <- matrix(0, nrow(X4), ncol(X4))

    for (j in 1:nrow(X4)){

      for (i in 1:ncol(X4)) {


        meth <- as.numeric(Y[[j,i]])
        unmeth <-  unmeth <- as.numeric(Z[[j,i]])

        m <- log2(max (meth,0) + p/ max(unmeth,0) +p)
        M[j,i] <- m


      } # for i

    } # for j

    M <- as.data.frame(M)

    colnames(M)<- colnames(X4)
    rownames(M)<- rownames(X4)
    cat("M matrix computed ...", "\n")

    # rownames were not assigned before this step
    X4 <- as.data.frame(X4)
    colnames(X4)<- colnames(M)
    rownames(X4)<- sites_filtered2


    # Results
    res<- list(Coverage_matrix = as.matrix(X4), Met_matrix = Y, Unmet_matrix = Z,
               Beta_matrix = B, M_matrix = M)

    cat("Coverting results to matrices ...", "\n")
    res2<- lapply(res, as.matrix.data.frame)



    return(res2)

}


