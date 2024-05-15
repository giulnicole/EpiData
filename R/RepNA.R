#' @title repNA
#'
#' @description A function for computing the comparisong between imputation methods, compute accuracy measure with replications and Kolmogorov-Smirnov test on CpG
#' \code{\link{repNA}} computes the imputation based on standard methods and correlation pattern between CpGs derived from eigen decomposition of the covariance matrix.
#'
#' @import missForest
#' @import cola
#' @import FNN
#' @import pcaMethods
#' @import softImpute
#' @import impute
#' @import apcluster
#' @import missMDA
#'
#' @name repNA
#'
#' @param cleaned.obj cleaned dataset from cleaning part 1 (coverage cleaning) and cleaning part 2 (outliers).
#' @param missing_prop missing proportion to be tested at each replication.
#' @param sel_method method for imputation to be selected.
#' @param matrix Matrix on which computing the correlation matrix; "M" as default.
#' @param varselect Index of the dataset to be used (numeric value 1-5). 1 = coverage counts, 2 = methylated counts, 3 = unmethylated counts, 4 = beta values, 5 = M values.
#' @param n.iter number of replication.
#' @param nb number of nearest neighbors for KNN method.
#' @param ncomp number of components for PCA methods.
#' @param trees number of tree for miss forest regression method.
#' @param plots plots for KS test on CpGs.
#'
#'
#' @return
#' list with datasets imputed from selected method and relative accuracy measures (root mean squared error - RMSE- and mean absolute error -MAE) and KS test results.
#'
#' @examples
#'  \dontrun{
#'  data("rangedObject")
#'  clean.coverage <- cleanCovMat(input.obj=data, max_na_cpg = 0.5,
#'                                  max_na_ind = 0.2,  cpg_removal_threshold = 10)
#'  clean.out <- cleanOutliers(filtered.obj=clean.coverage, outlier_threshold=5,
#'                              remove_outliers = TRUE)
#'  stats <- lapply(clean.out, na.omit)
#'  imp.M <- repNA(cleaned.obj=stats,
#'                 missing_prop= 0.2,
#'                 varselect=5,
#'                 n.iter= 10,
#'                 sel_method=c(1:10),
#'                 trees= 50, nb= 10, ncomp= 2,
#'                  matrix="M")
#' }
#'
#'
#'
#'
#' @export
repNA <- function(cleaned.obj, varselect=5,
                  missing_prop= 0.3,
                  n.iter= 2,
                  sel_method, matrix="M",
                  trees= 50, nb= 10, ncomp= 2, plots=FALSE) {



  # original cleaned data

  mat1<- cleaned.obj[[varselect]]

  # Define the proportion of missing values to introduce
  # missing_prop <- 0.2  # 20% missing values


# - - - - - - - - - - - - -  - - - - - - - - - - - - -
  #  # Method 1 = Mean imputation

  if (1 %in% sel_method) {

    origmat1 <- as.data.frame(mat1)

    cat("Mean imputation - in progress\n")

    mean.mat1 <- list()
    rmse.mean.mat1 <- NULL
    mae.mean.mat1 <- NULL
    time.mean <- NULL
    KStest.mean <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      start.time <- Sys.time()
      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- origmat1
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the mat1
      na_positions <- is.na(masked_data.mat1)


      # Impute the missing values
     t <- system.time(
        for (s in 1:ncol(masked_data.mat1)) {
          masked_data.mat1[is.na(masked_data.mat1[, s]), s] <- mean(masked_data.mat1[, s], na.rm = TRUE)

        })

     time.mean[i] <-t[[3]]

      #imputed_data.mat1 <- round(masked_data.mat1)
      #imputed_data.mat1[imputed_data.mat1<0] <- 0
      mean.mat1[[i]] <- masked_data.mat1


      err <- errors(origmat1, masked_data.mat1, na_positions)

      S3 <- err[1]
      M3 <- err[2]

      rmse.mean.mat1[i]  <- S3
      mae.mean.mat1[i]  <- M3

      KS_test <- KStesCpG(origmat1, masked_data.mat1)

      KStest.mean <- KS_test

      cat("Time elapsed at iteration:\n")
      print(t[[3]])
      cat("\n")



      results.mean <- list(imputed.mean.mat1 =  mean.mat1,
                           RMSE.mean.mat1 = rmse.mean.mat1,
                           MAE.mean.mat1 = mae.mean.mat1,
                           Comp_time = time.mean,
                           KS_statsitics = KStest.mean
                           )



      }# mean


    }  else results.mean <- NULL



# - - - - - - - - - - - - -  - - - - - - - - - - - - -
  #  # Method 2 = PPCA imputation

  if (2 %in% sel_method) {

    origmat1 <- as.data.frame(mat1)

    cat("PPCA imputation - in progress\n")

    ppca.mat1 <- list()
    rmse.ppca.mat1 <- NULL
    mae.ppca.mat1 <- NULL
    time.ppca <- NULL
    KStest.ppca <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      start.time <- Sys.time()
      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- origmat1
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the mat1
      na_positions <- is.na(masked_data.mat1)

      # Impute the missing values
      imputed_data.mat1<- pca(masked_data.mat1, nPcs=ncomp, method="ppca", center = F)
      imputed_data.mat1 <- as.data.frame(imputed_data.mat1@completeObs)
      masked_data.mat1<- imputed_data.mat1

      #time.ppca[i] <-t[[3]]

      #imputed_data.mat1 <- round(masked_data.mat1)
      #imputed_data.mat1[imputed_data.mat1<0] <- 0
      ppca.mat1[[i]] <- masked_data.mat1

      err <- errors(origmat1, masked_data.mat1, na_positions)

      S3 <- err[1]
      M3 <- err[2]

      rmse.ppca.mat1[i]  <- S3
      mae.ppca.mat1[i]  <- M3

      KS_test <- KStesCpG(origmat1, masked_data.mat1)

      KStest.ppca <- KS_test

      #cat("Time elapsed at iteration:\n")
      #print(t[[3]])
      #cat("\n")


      results.ppca <- list(imputed.ppca.mat1 =  ppca.mat1,
                           RMSE.ppca.mat1 = rmse.ppca.mat1,
                           MAE.ppca.mat1 = mae.ppca.mat1,
                           #Comp_time = time.ppca,
                           KS_statsitics = KStest.ppca
      )



    }# ppca




  }  else results.ppca <- NULL


# - - - - - - - - - - - - -  - - - - - - - - - - - - -
  #  # Method 3 = BPCA imputation

  if (3 %in% sel_method) {

    origmat1 <- as.data.frame(mat1)

    cat("BPCA imputation - in progress\n")

    bpca.mat1 <- list()
    rmse.bpca.mat1 <- NULL
    mae.bpca.mat1 <- NULL
    time.bpca <- NULL
    KStest.bpca <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      start.time <- Sys.time()
      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- origmat1
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the mat1
      na_positions <- is.na(masked_data.mat1)


      # Impute the missing values
      imputed_data.mat1<- pca(masked_data.mat1, nPcs=ncomp, method="bpca", center = F)
      imputed_data.mat1 <- as.data.frame(imputed_data.mat1@completeObs)
      masked_data.mat1 <- imputed_data.mat1


      #time.bpca[i] <-t[[3]]

      #imputed_data.mat1 <- round(masked_data.mat1)
      #imputed_data.mat1[imputed_data.mat1<0] <- 0
      bpca.mat1[[i]] <- masked_data.mat1

      err <- errors(origmat1, masked_data.mat1, na_positions)

      S3 <- err[1]
      M3 <- err[2]

      rmse.bpca.mat1[i]  <- S3
      mae.bpca.mat1[i]  <- M3

      KS_test <- KStesCpG(origmat1, masked_data.mat1)

      KStest.bpca <- KS_test

      #cat("Time elapsed at iteration:\n")
      #print(t[[3]])
      #cat("\n")


      results.bpca <- list(imputed.bpca.mat1 =  bpca.mat1,
                           RMSE.bpca.mat1 = rmse.bpca.mat1,
                           MAE.bpca.mat1 = mae.bpca.mat1,
                           #Comp_time = time.bpca,
                           KS_statsitics = KStest.bpca
      )



    }# bpca




  }  else results.bpca <- NULL





# - - - - - - - - - - - - -  - - - - - - - - - - - - -
  #  # Method 4 = NIPALS imputation

  if (4 %in% sel_method) {

    origmat1 <- as.data.frame(mat1)

    cat("NIPALS imputation - in progress\n")

    nipals.mat1 <- list()
    rmse.nipals.mat1 <- NULL
    mae.nipals.mat1 <- NULL
    time.nipals <- NULL
    KStest.nipals <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      start.time <- Sys.time()
      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- origmat1
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the mat1
      na_positions <- is.na(masked_data.mat1)


      # Impute the missing values
      imputed_data.mat1<- pca(masked_data.mat1, nPcs=ncomp, method="nipals", center = F)
      imputed_data.mat1 <- as.data.frame(imputed_data.mat1@completeObs)
      masked_data.mat1<- imputed_data.mat1

      #time.nipals[i] <-t[[3]]

      #imputed_data.mat1 <- round(masked_data.mat1)
      #imputed_data.mat1[imputed_data.mat1<0] <- 0
      nipals.mat1[[i]] <- masked_data.mat1

      err <- errors(origmat1, masked_data.mat1, na_positions)


      S3 <- err[1]
      M3 <- err[2]

      rmse.nipals.mat1[i]  <- S3
      mae.nipals.mat1[i]  <- M3

      KS_test <- KStesCpG(origmat1, masked_data.mat1)

      KStest.nipals <- KS_test

      #cat("Time elapsed at iteration:\n")
      #print(t[[3]])
      #cat("\n")


      results.nipals <- list(imputed.nipals.mat1 =  nipals.mat1,
                          RMSE.nipals.mat1 = rmse.nipals.mat1,
                          MAE.nipals.mat1 = mae.nipals.mat1,
                          #Comp_time = time.nipals,
                          KS_statsitics = KStest.nipals
      )



    }# nipals




  }  else results.nipals <- NULL

# - - - - - - - - - - - - -  - - - - - - - - - - - - -
  #  # Method 5 = MDA imputation

  if (5 %in% sel_method) {

    origmat1 <- as.data.frame(mat1)

    cat("MDA imputation - in progress\n")

    mda.mat1 <- list()
    rmse.mda.mat1 <- NULL
    mae.mda.mat1 <- NULL
    time.mda <- NULL
    KStest.mda <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      start.time <- Sys.time()
      cat("Iteration:", i, "\n")

      # Mask the missing values
      origmat2 <- t(origmat1)
      masked_data.mat1 <- origmat2
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the mat1
      na_positions <- is.na(t(masked_data.mat1))

      # Impute the missing values
      imputed_data.mat1 <- missMDA::imputePCA(masked_data.mat1, method = "Regularized")
      imputed_data.mat2 <- imputed_data.mat1$completeObs
      masked_data.mat1 <- as.data.frame(t(imputed_data.mat2))
      #time.mda[i] <-t[[3]]

      mda.mat1[[i]] <- masked_data.mat1

      err <- errors(origmat1, masked_data.mat1, na_positions)

      S3 <- err[1]
      M3 <- err[2]

      rmse.mda.mat1[i]  <- S3
      mae.mda.mat1[i]  <- M3

      KS_test <- KStesCpG(origmat1, masked_data.mat1)

      KStest.mda <- KS_test

      #cat("Time elapsed at iteration:\n")
      #print(t[[3]])
      #cat("\n")


      results.mda <- list(imputed.mda.mat1 =  mda.mat1,
                              RMSE.mda.mat1 = rmse.mda.mat1,
                              MAE.mda.mat1 = mae.mda.mat1,
                              #Comp_time = time.mda,
                              KS_statsitics = KStest.mda
      )



    }# mda



  }  else results.mda <- NULL




# - - - - - - - - - - - - -  - - - - - - - - - - - - -
  #  # Method 6 = EM imputation

  if (6 %in% sel_method) {

    origmat1 <- as.data.frame(mat1)

    cat("EM imputation - in progress\n")

    em.mat1 <- list()
    rmse.em.mat1 <- NULL
    mae.em.mat1 <- NULL
    time.em <- NULL
    KStest.em <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      start.time <- Sys.time()
      cat("Iteration:", i, "\n")

      # Mask the missing values
      origmat2 <- t(origmat1)
      masked_data.mat1 <- origmat2
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the mat1
      na_positions <- is.na(t(masked_data.mat1))

      # Impute the missing values
      imputed_data.mat1 <- missMDA::imputePCA(masked_data.mat1, method = "EM")
      imputed_data.mat2 <- imputed_data.mat1$completeObs
      masked_data.mat1 <- as.data.frame(t(imputed_data.mat2))
      #time.em[i] <-t[[3]]

      em.mat1[[i]] <- masked_data.mat1

      err <- errors(origmat1, masked_data.mat1, na_positions)


      S3 <- err[1]
      M3 <- err[2]

      rmse.em.mat1[i]  <- S3
      mae.em.mat1[i]  <- M3

      KS_test <- KStesCpG(origmat1, masked_data.mat1)

      KStest.em <- KS_test

      #cat("Time elapsed at iteration:\n")
      #print(t[[3]])
      #cat("\n")


      results.em <- list(imputed.em.mat1 =  em.mat1,
                             RMSE.em.mat1 = rmse.em.mat1,
                             MAE.em.mat1 = mae.em.mat1,
                             #Comp_time = time.em,
                             KS_statsitics = KStest.em
      )



    }# em



  }  else results.em <- NULL



# - - - - - - - - - - - - -  - - - - - - - - - - - - -
  #  # Method 7 = MissForest imputation

  if (7 %in% sel_method) {

    origmat1 <- as.data.frame(mat1)

    cat("MissForest imputation - in progress\n")

    mf.mat1 <- list()
    rmse.mf.mat1 <- NULL
    mae.mf.mat1 <- NULL
    time.mf <- NULL
    KStest.mf <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      start.time <- Sys.time()
      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- t(origmat1)
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the mat1
      na_positions <- is.na(t(masked_data.mat1))

      # Impute the missing values
      imputed_data.mat1 <-  missForest(masked_data.mat1, ntree = 10, replace = TRUE)
      masked_data.mat1  <- as.data.frame(t(imputed_data.mat1[[1]]))

      #time.mf[i] <-t[[3]]

      mf.mat1[[i]] <- masked_data.mat1

      err <- errors(origmat1, masked_data.mat1, na_positions)


      S3 <- err[1]
      M3 <- err[2]

      rmse.mf.mat1[i]  <- S3
      mae.mf.mat1[i]  <- M3

      KS_test <- KStesCpG(origmat1, masked_data.mat1)

      KStest.mf <- KS_test

      #cat("Time elapsed at iteration:\n")
      #print(t[[3]])
      #cat("\n")


      results.mf <- list(imputed.mf.mat1 =  mf.mat1,
                                 RMSE.mf.mat1 = rmse.mf.mat1,
                                 MAE.mf.mat1 = mae.mf.mat1,
                                 #Comp_time = time.mf,
                                 KS_statsitics = KStest.mf
      )



    }# mf



  }  else results.mf <- NULL




# - - - - - - - - - - - - -  - - - - - - - - - - - - -
  #  # Method 8 = knn imputation

  if (8 %in% sel_method) {

    origmat1 <- as.data.frame(mat1)

    cat("Knn imputation - in progress\n")

    knn.mat1 <- list()
    rmse.knn.mat1 <- NULL
    mae.knn.mat1 <- NULL
    time.knn <- NULL
    KStest.knn <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      start.time <- Sys.time()
      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- t(origmat1)
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the mat1
      na_positions <- is.na(t(masked_data.mat1))

      # Impute the missing values
      imputed_data.mat1<- impute::impute.knn(data = as.matrix(masked_data.mat1), k = nb, rowmax = 100)
      imputed_data.mat1 <- imputed_data.mat1$data

      #time.knn[i] <-t[[3]]

      knn.mat1[[i]] <- imputed_data.mat1
      masked_data.mat1 <- as.data.frame(t(imputed_data.mat1))

      err <- errors(origmat1, masked_data.mat1, na_positions)

      S3 <- err[1]
      M3 <- err[2]

      rmse.knn.mat1[i]  <- S3
      mae.knn.mat1[i]  <- M3

      KS_test <- KStesCpG(origmat1, masked_data.mat1, plots = F)

      KStest.knn <- KS_test

      #cat("Time elapsed at iteration:\n")
      #print(t[[3]])
      #cat("\n")


      results.knn <- list(imputed.knn.mat1 =  knn.mat1,
                           RMSE.knn.mat1 = rmse.knn.mat1,
                           MAE.knn.mat1 = mae.knn.mat1,
                           #Comp_time = time.knn,
                           KS_statsitics = KStest.knn
      )


   }# knn



  }  else results.knn <- NULL



# - - - - - - - - - - - - -  - - - - - - - - - - - - -
  # Method 9 = imputation based on cor

  if (9 %in% sel_method) {

    origmat1 <- as.data.frame(mat1)

    cat("Correlation-based imputation - in progress\n")

    corr.mat1 <- list()
    rmse.corr.mat1 <- NULL
    mae.corr.mat1 <- NULL
    time.corr <- NULL
    KStest.corr <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      start.time <- Sys.time()
      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- origmat1
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the mat1
      na_positions <- is.na(masked_data.mat1)

      cleaned.obj[[varselect]] <-  t(masked_data.mat1)

      stats.m <-statsCpG(cleaned.obj = cleaned.obj, varselect)

      simu3 <- imputeCorr(cleaned.obj=stats.m,
                             matrix="M",
                             varselect)


      imputed_data.mat1 <- simu3$Imputed
      corr.mat1[[i]] <- imputed_data.mat1
      masked_data.mat1 <- imputed_data.mat1

      err <- errors(origmat1, masked_data.mat1, na_positions)


      S3 <- err[1]
      M3 <- err[2]

      rmse.corr.mat1[i]  <- S3
      mae.corr.mat1[i]  <- M3

      KS_test <- KStesCpG(origmat1, masked_data.mat1, plots = F)

      KStest.corr <- KS_test

      #cat("Time elapsed at iteration:\n")
      #print(t[[3]])
      #cat("\n")


      results.corr <- list(imputed.corr.mat1 =  corr.mat1,
                            RMSE.corr.mat1 = rmse.corr.mat1,
                            MAE.corr.mat1 = mae.corr.mat1,
                            #Comp_time = time.corr,
                            KS_statsitics = KStest.corr
      )



    }# corr



  }  else results.corr <- NULL


  # - - - - - - - - - - - - -  - - - - - - - - - - - - -
  # Method 10 = imputation based on consensus partition

  if (10 %in% sel_method) {

    origmat1 <- as.data.frame(mat1)

    cat("Consensus partition imputation - in progress\n")

    cpart.mat1 <- list()
    rmse.cpart.mat1 <- NULL
    mae.cpart.mat1 <- NULL
    time.cpart <- NULL
    KStest.cpart <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      start.time <- Sys.time()
      cat("Iteration:", i, "\n")

      # Mask the missing values
      origmat2 <- (origmat1)
      masked_data.mat1<- prodNA(origmat2, missing_prop) # Mask some values

      # Assign NA values to the mat1
      na_positions <- is.na(masked_data.mat1)


      # Impute the missing values
      #time.cpart.knn[i] <- system.time(
      imputed_data.mat1<- consensusPart(data = masked_data.mat1,
                                                dataset.simulated=origmat1,
                                                partition_method = "pam", top_rows = "ATC",
                                                k_to_test = 2:8, top_n = 100)


      cpart_data.mat1 <-  t(imputed_data.mat1$dataset.cluster)
      cpart.mat1[[i]] <- cpart_data.mat1
      masked_data.mat1[na_positions] <- cpart_data.mat1[na_positions]

      err <- errors(origmat1, masked_data.mat1, na_positions)


      S3 <- err[1]
      M3 <- err[2]

      rmse.cpart.mat1[i]  <- S3
      mae.cpart.mat1[i]  <- M3

      KS_test <- KStesCpG(origmat1, masked_data.mat1, plots = F)

      KStest.cpart <- KS_test

      #cat("Time elapsed at iteration:\n")
      #print(t[[3]])
      #cat("\n")


      results.cpart <- list(imputed.cpart.mat1 = cpart.mat1,
                            RMSE.cpart.mat1 = rmse.cpart.mat1,
                            MAE.cpart.mat1 = mae.cpart.mat1,
                            #Comp_time = time.cpart,
                            KS_statsitics = KStest.cpart
      )



    }#cpart



  }  else results.cpart <- NULL






  return(list(MEAN_imputation = results.mean,
              PPCA_imputation = results.ppca,
              BPCA_imputation = results.bpca,
              MF_imputation = results.mf,
              NIPALS_imputation = results.nipals,
              MDA_imputation = results.mda,
              EM_imputation = results.em,
              KNN_imputation = results.knn,
              CORR_imputation = results.corr,
              CPART_imputation = results.cpart))






  }  # function



