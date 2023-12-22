#' @title NAreplicating
#' @description
#' \code{\link{NAreplicating}}
#'
#'
#' @param
#' @param
#' @param
#' @param
#' @param
#' @param
#'
#' @name NAreplicating
#'
#' @return
#'
#'  \item{}{}
#'  \item{}{}
#'
#'
#'
#'
#' @export
#'
NAreplicating <- function(list.cleaned, varselect, missing_prop= 0.3,
                         scale= TRUE, n.iter= 10, sel_method= 1,
                         trees= 50, nb= 10, ncomp= 2) {



  # original cleaned data
  mat1<- list.cleaned[[varselect]]
  origmat1 <- as.data.frame(t(mat1))

  # Define the proportion of missing values to introduce
  # missing_prop <- 0.2  # 20% missing values


  # 1 = Mean (run if selected)

  if (1 %in% sel_method) {

    cat("Mean imputation - in progress\n")

    mean.mat1 <- list()
    rmse.mean.mat1 <- NULL
    mae.mean.mat1 <- NULL
    mape.mean.mat1 <- NULL

    time.mean <- NULL

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
      time.mean[i] <- system.time(
        for (s in 1:ncol(masked_data.mat1)) {
          masked_data.mat1[is.na(masked_data.mat1[, s]), s] <- mean(masked_data.mat1[, s], na.rm = TRUE)

        })



      #imputed_data.mat1 <- round(masked_data.mat1)
      #imputed_data.mat1[imputed_data.mat1<0] <- 0
      mean.mat1[[i]] <- masked_data.mat1


      # Calculate RMSE, MAE, KL distance, KL divergence for mat1
      error.rmse.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mae.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mape.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))


      for (j in 1:nrow(origmat1)){

        for (l in 1:ncol(origmat1)) {


          if (na_positions[j,l] ==TRUE){
          # computing RMSE
          diff3 <- as.numeric(origmat1[j,l] - masked_data.mat1[j,l])
          sqdiff3 <- (diff3)^2
          error.rmse.mat1[j,l] <- sqdiff3

          #computing MAE
          abs3 <- abs(as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))
          error.mae.mat1[j,l] <- abs3

          #computing MAPE
          abs.rap <- abs((as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))/as.numeric(origmat1[j,l]))
          error.mape.mat1[j,l] <- abs.rap

          } # if

        } # for i


      } # for j

      x3 = mean(error.rmse.mat1)
      S3 = sqrt(x3)/(max(origmat1)-min(origmat1))
      rmse.mean.mat1[i]  <- S3

      M3 = mean(error.mae.mat1)/(max(origmat1)-min(origmat1))
      mae.mean.mat1[i] <- M3

      MP3 = mean(error.mape.mat1)
      mape.mean.mat1[i] <- MP3
      # Accuracy
      # Reg.eval.knn2[i] <- R2(origmat1, imputed_data.mat1)

      # ks.beta  <- stats::ks.test(origbeta, imputed_data.beta, exact=TRUE)$statistic

      # Print summary statistics
      # cat("Summary Statistics After Imputation:\n")
      # print(summary(data))
      # cat("\n")

      cat("RMSE at iteration:\n")
      print(S3)
      cat("\n")

      cat("MAE at iteration:\n")
      print(M3)
      cat("\n")


      # cat("Accuracy:\n")
      # print(c(Reg.eval.knn2[i]))
      # cat("\n")

      cat("time difference at iteration:\n")
      print(time.mean[i])
      cat("\n")


      results.mean <- list(imputed.mean.mat1 =  mean.mat1,
                           RMSE.mean.mat1 = rmse.mean.mat1,
                           MAE.mean.mat1 = mae.mean.mat1,
                           Comp_time = time.mean)


    }# iter


  }# mean

  else results.mean <- NULL



  origmat1 <- as.data.frame(t(mat1))

  if (2 %in% sel_method) {
    cat("PPCA in PCA Methods imputation - in progress\n")


    # original data (simulated)
    ppca.mat1 <- list()
    rmse.ppca.mat1 <- NULL
    mae.ppca.mat1 <- NULL
    mape.ppca.mat1 <- NULL

    time.ppca <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- origmat1
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the methylated counts
      na_positions <- is.na(masked_data.mat1)

      time.ppca[i]<- system.time(
        # Impute the missing values
        imputed_data.mat1<- pca(masked_data.mat1, nPcs=ncomp, method="ppca", center = F))

      imputed_data.mat1 <-imputed_data.mat1@completeObs
      #imputed_data.mat1 <- round(imputed_data.mat1@completeObs)
      #imputed_data.mat1[imputed_data.mat1<0] <- 0
      ppca.mat1[[i]] <- imputed_data.mat1

      # Update the original dataset with the imputed values
      masked_data.mat1[is.na(masked_data.mat1)] <- imputed_data.mat1[is.na(masked_data.mat1)]

      # Calculate RMSE, MAE, KL distance, KL divergence for met
      error.rmse.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mae.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mape.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))


      time.elapsed.ppca <- NULL

      for (j in 1:nrow(origmat1)){

        for (l in 1:ncol(origmat1)) {

          if (na_positions[j,l]==TRUE){

          # computing RMSE
          diff3 <- as.numeric(origmat1[j,l] - masked_data.mat1[j,l])
          sqdiff3 <- (diff3)^2
          error.rmse.mat1[j,l] <- sqdiff3

          #computing MAE
          abs3 <- abs(as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))
          error.mae.mat1[j,l] <- abs3

          #computing MAPE
          abs.rap <- abs((as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))/as.numeric(origmat1[j,l]))
          error.mape.mat1[j,l] <- abs.rap

          } # if

        } # for i


      } # for j


      x3 = mean(error.rmse.mat1)
      S3 = sqrt(x3)/(max(origmat1)-min(origmat1))
      rmse.ppca.mat1[i]  <- S3

      M3 = mean(error.mae.mat1)/(max(origmat1)-min(origmat1))
      mae.ppca.mat1[i] <- M3

      # Accuracy
      # Reg.eval.ppca[i] <- R2(origmat1, imputed_data.mat1)

      # ks.beta  <- stats::ks.test(origbeta, imputed_data.beta, exact=TRUE)$statistic

      # Print summary statistics
      # cat("Summary Statistics After Imputation:\n")
      # print(summary(data))
      # cat("\n")

      cat("RMSE at iteration:\n")
      print(S3)
      cat("\n")

      cat("MAE at iteration:\n")
      print(M3)
      cat("\n")

      # cat("Accuracy:\n")
      # print(c(Reg.eval.ppca[i]))
      # cat("\n")

      cat("time difference at iteration:\n")
      print(time.ppca[i])
      cat("\n")


      results.ppca <- list(imputed.ppca.mat1 =  ppca.mat1,
                           RMSE.ppca.mat1 = rmse.ppca.mat1,
                           MAE.ppca.mat1 = mae.ppca.mat1,
                           Comp_time = time.ppca)


    }# iter


  }# ppca

  else results.ppca <- NULL


  origmat1 <- as.data.frame(t(mat1))

  if (3 %in% sel_method) {
    cat("BPCA in PCA Methods imputation - in progress\n")


    # original data (simulated)
    bpca.mat1 <- list()
    rmse.bpca.mat1 <- NULL
    mae.bpca.mat1 <- NULL
    mape.bpca.mat1 <- NULL

    time.bpca <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- origmat1
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the methylated counts
      na_positions <- is.na(masked_data.mat1)


      # Impute the missing values
      time.bpca[i]<- system.time(
        imputed_data<- pca(masked_data.mat1, nPcs=ncomp, method="bpca", center = F))

      impued_data.mat1 <- as.data.frame(imputed_data@completeObs)
      #imputed_data.mat1 <- round(imputed_data.mat1@completeObs)
      #imputed_data.mat1[imputed_data.mat1<0] <- 0
      bpca.mat1[[i]] <- imputed_data.mat1

      # Update the original dataset with the imputed values
      masked_data.mat1[na_positions] <- imputed_data.mat1[na_positions]

      # Calculate RMSE, MAE, KL distance, KL divergence for met
      error.rmse.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mae.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mape.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))

      for (j in 1:nrow(origmat1)){

        for (l in 1:ncol(origmat1)) {

          if (na_positions[j,l]==TRUE) {

          # computing RMSE
          diff3 <- as.numeric(origmat1[j,l] - masked_data.mat1[j,l])
          sqdiff3 <- (diff3)^2
          error.rmse.mat1[j,l] <- sqdiff3

          #computing MAE
          abs3 <- abs(as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))
          error.mae.mat1[j,l] <- abs3

          #computing MAPE
          abs.rap <- abs((as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))/as.numeric(origmat1[j,l]))
          error.mape.mat1[j,l] <- abs.rap

          } # if

        } # for i


      } # for j


      x3 = mean(error.rmse.mat1)
      S3 = sqrt(x3)/(max(origmat1)-min(origmat1))
      rmse.bpca.mat1[i]  <- S3

      M3 = mean(error.mae.mat1)/(max(origmat1)-min(origmat1))
      mae.bpca.mat1[i] <- M3

      # Accuracy
      # Reg.eval.bpca[i] <- R2(origmat1, imputed_data.mat1)

      # ks.beta  <- stats::ks.test(origbeta, imputed_data.beta, exact=TRUE)$statistic

      # Print summary statistics
      # cat("Summary Statistics After Imputation:\n")
      # print(summary(data))
      # cat("\n")

      cat("RMSE at iteration:\n")
      print(S3)
      cat("\n")

      cat("MAE at iteration:\n")
      print(M3)
      cat("\n")

      # cat("Accuracy:\n")
      # print(c(Reg.eval.bpca[i]))
      # cat("\n")

      cat("time difference at iteration:\n")
      print(time.bpca[i])
      cat("\n")


      results.bpca <- list(imputed.bpca.mat1 =  bpca.mat1,
                           RMSE.bpca.mat1 = rmse.bpca.mat1,
                           MAE.bpca.mat1 = mae.bpca.mat1,
                           Comp_time = time.bpca)


    }# iter


  }# bpca

  else results.bpca <- NULL


  origmat1 <- as.data.frame(t(mat1))

  if  (4 %in% sel_method) {
    cat("missForest imputation - in progress\n")


    # original data (simulated)
    missf.mat1 <- list()
    rmse.missf.mat1 <- NULL
    mae.missf.mat1 <- NULL
    mape.missf.mat1 <- NULL

    time.missf <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {


      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- origmat1
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the methylated counts
      na_positions <- is.na(masked_data.mat1)


      # Impute the missing values
      time.missf[i] <- system.time(
        imputed_data.mat1 <-  missForest(masked_data.mat1, ntree = trees, replace = TRUE))
      imputed_data.mat1 <- round(imputed_data.mat1$ximp)
      imputed_data.mat1[imputed_data.mat1<0] <- 0
      missf.mat1[[i]] <- imputed_data.mat1

      # Update the original dataset with the imputed values
      masked_data.mat1[is.na(masked_data.mat1)] <- imputed_data.mat1[is.na(masked_data.mat1)]

      # Calculate RMSE, MAE, KL distance, KL divergence for met
      error.rmse.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mae.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mape.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))

      time.elapsed.missf <- NULL

      for (j in 1:nrow(origmat1)){

        for (l in 1:ncol(origmat1)) {

          if (na_positions[j,l]==TRUE){
          # computing RMSE
          diff3 <- as.numeric(origmat1[j,l] - masked_data.mat1[j,l])
          sqdiff3 <- (diff3)^2
          error.rmse.mat1[j,l] <- sqdiff3

          #computing MAE
          abs3 <- abs(as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))
          error.mae.mat1[j,l] <- abs3

          #computing MAPE
          abs.rap <- abs((as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))/as.numeric(origmat1[j,l]))
          error.mape.mat1[j,l] <- abs.rap

          } # if

        } # for i


      } # for j


      x3 = mean(error.rmse.mat1)
      S3 = sqrt(x3)/(max(origmat1)-min(origmat1))
      rmse.missf.mat1[i]  <- S3

      M3 = mean(error.mae.mat1)/(max(origmat1)-min(origmat1))
      mae.missf.mat1[i] <- M3

      # Accuracy
      # Reg.eval.missf[i] <- R2(origmat1, imputed_data.mat1)

      # ks.beta  <- stats::ks.test(origbeta, imputed_data.beta, exact=TRUE)$statistic

      # Print summary statistics
      # cat("Summary Statistics After Imputation:\n")
      # print(summary(data))
      # cat("\n")

      cat("RMSE at iteration:\n")
      print(S3)
      cat("\n")

      cat("MAE at iteration:\n")
      print(M3)
      cat("\n")

      # cat("Accuracy:\n")
      # print(c(Reg.eval.missf[i]))
      # cat("\n")

      cat("time difference at iteration:\n")
      print(time.missf[i])
      cat("\n")


      results.missf <- list(imputed.missf.mat1 =  missf.mat1,
                            RMSE.missf.mat1 = rmse.missf.mat1,
                            MAE.missf.mat1 = mae.missf.mat1,
                            Comp_time = time.missf)


    }# iter


  }# missf

  else results.missf <- NULL


  origmat1 <- as.data.frame(t(mat1))

  if (5 %in% sel_method) {
    cat("SVD in PCA Methods imputation - in progress\n")


    # original data (simulated)
    svd.mat1 <- list()
    rmse.svd.mat1 <- NULL
    mae.svd.mat1 <- NULL
    mape.svd.mat1 <- NULL

    time.svd <- NULL

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {


      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- origmat1
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the methylated counts
      na_positions <- is.na(masked_data.mat1)

      # Impute the missing values
      time.svd[i] <- system.time(
        imputed_data.mat1<- pca(masked_data.mat1, nPcs=ncomp, method="svd", center = F))
      imputed_data.mat1 <- imputed_data.mat1@completeObs
      #imputed_data.mat1 <- round(imputed_data.mat1@completeObs)
      #imputed_data.mat1[imputed_data.mat1<0] <- 0
      svd.mat1[[i]] <- imputed_data.mat1


      # Update the original dataset with the imputed values
      masked_data.mat1[is.na(masked_data.mat1)] <- imputed_data.mat1[is.na(masked_data.mat1)]

      # Calculate RMSE, MAE, KL distance, KL divergence for met
      error.rmse.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mae.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mape.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      time.elapsed.svd <- NULL

      for (j in 1:nrow(origmat1)){

        for (l in 1:ncol(origmat1)) {

          if (na_positions[j,l]==TRUE) {
          # computing RMSE
          diff3 <- as.numeric(origmat1[j,l] - masked_data.mat1[j,l])
          sqdiff3 <- (diff3)^2
          error.rmse.mat1[j,l] <- sqdiff3


          #computing MAE
          abs3 <- abs(as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))
          error.mae.mat1[j,l] <- abs3

          #computing MAPE
          abs.rap <- abs((as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))/as.numeric(origmat1[j,l]))
          error.mape.mat1[j,l] <- abs.rap


          } # if

        } # for i


      } # for j


      x3 = mean(error.rmse.mat1)
      S3 = sqrt(x3)/(max(origmat1)-min(origmat1))
      rmse.svd.mat1[i]  <- S3

      M3 = mean(error.mae.mat1)/(max(origmat1)-min(origmat1))
      mae.svd.mat1[i] <- M3

      # Accuracy
      # Reg.eval.svd[i] <- R2(origmat1, imputed_data.mat1)

      # ks.beta  <- stats::ks.test(origbeta, imputed_data.beta, exact=TRUE)$statistic

      # Print summary statistics
      # cat("Summary Statistics After Imputation:\n")
      # print(summary(data))
      # cat("\n")

      cat("RMSE at iteration:\n")
      print(S3)
      cat("\n")

      cat("MAE at iteration:\n")
      print(M3)
      cat("\n")

      # cat("Accuracy:\n")
      # print(c(Reg.eval.svd[i]))
      # cat("\n")

      cat("time difference at iteration:\n")
      print(time.svd[i])
      cat("\n")


      results.svd <- list(imputed.svd.mat1 =  svd.mat1,
                          RMSE.svd.mat1 = rmse.svd.mat1,
                          MAE.svd.mat1 = mae.svd.mat1,
                          Comp_time = time.svd)



    }# iter


  }# svd

  else results.svd <- NULL



  origmat1 <- as.data.frame(t(mat1))

  if (6 %in% sel_method) {
    cat("missMDA imputation - in progress\n")


    # original data (simulated)
    missmda.mat1 <- list()
    rmse.missmda.mat1 <- NULL
    mae.missmda.mat1 <- NULL
    mape.missmda.mat1 <- NULL

    time.missmda <- NULL

    #origmat1 <- as.data.frame(t(mat1))
    origmat2 <- t(origmat1)

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      start.time <- Sys.time()
      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- origmat2
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the methylated counts
      na_positions <- is.na(masked_data.mat1)

      # Impute the missing values
      time.missmda[i]<- system.time(
        imputed_data.mat1 <- missMDA::imputePCA(masked_data.mat1, ncp = ncomp, method = "Regularized"))
      imputed_data.mat1 <- t(imputed_data.mat1$completeObs)
      #imputed_data.mat1 <- t(round(imputed_data.mat1$completeObs))
      #imputed_data.mat1[imputed_data.mat1<0] <- 0
      missmda.mat1[[i]] <- imputed_data.mat1


      # Update the original dataset with the imputed values
      masked_data.mat1[is.na(masked_data.mat1)] <- imputed_data.mat1[is.na(masked_data.mat1)]

      # Calculate RMSE, MAE, KL distance, KL divergence for met
      error.rmse.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mae.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mape.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      masked_data.mat1 <- t(masked_data.mat1)


      time.elapsed.missmda <- NULL

      for (j in 1:nrow(origmat1)){

        for (l in 1:ncol(origmat1)) {

          if (na_positions[l,j]==TRUE){
          # computing RMSE
          diff3 <- as.numeric(origmat1[j,l] - masked_data.mat1[j,l])
          sqdiff3 <- (diff3)^2
          error.rmse.mat1[j,l] <- sqdiff3


          #computing MAE
          abs3 <- abs(as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))
          error.mae.mat1[j,l] <- abs3


          #computing MAPE
          abs.rap <- abs((as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))/as.numeric(origmat1[j,l]))
          error.mape.mat1[j,l] <- abs.rap

          } # if

        } # for i


      } # for j


      x3 = mean(error.rmse.mat1)
      S3 = sqrt(x3)/(max(origmat1)-min(origmat1))
      rmse.missmda.mat1[i]  <- S3

      M3 = mean(error.mae.mat1)/(max(origmat1)-min(origmat1))
      mae.missmda.mat1[i] <- M3

      # Accuracy
      # Reg.eval.missMDA[i] <- R2(origmat1, imputed_data.mat1)

      # ks.beta  <- stats::ks.test(origbeta, imputed_data.beta, exact=TRUE)$statistic

      # Print summary statistics
      # cat("Summary Statistics After Imputation:\n")
      # print(summary(data))
      # cat("\n")

      cat("RMSE at iteration:\n")
      print(S3)
      cat("\n")

      cat("MAE at iteration:\n")
      print(M3)
      cat("\n")

      # cat("Accuracy:\n")
      # print(c(Reg.eval.missMDA[i]))
      # cat("\n")

      cat("time difference at iteration:\n")
      print(time.missmda[i])
      cat("\n")


      results.missmda <- list(imputed.missmda.mat1 =  missmda.mat1,
                              RMSE.missmda.mat1 = rmse.missmda.mat1,
                              MAE.missmda.mat1 = mae.missmda.mat1,
                              Comp_time = time.missmda)


    }# iter


  }# missMDA

  else results.missmda <- NULL


  origmat1 <- as.data.frame(t(mat1))

  if (7 %in% sel_method) {
    cat("missEM imputation - in progress\n")


    # original data (simulated)
    missem.mat1 <- list()
    rmse.missem.mat1 <- NULL
    mae.missem.mat1 <- NULL
    mape.missem.mat1 <- NULL

    time.missem <- NULL

    # origmat1 <- as.data.frame(t(mat1))
    origmat3 <- t(origmat1)
    # Repeat the process n_iterations times
    for (i in 1:n.iter) {

      start.time <- Sys.time()
      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- origmat3
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the methylated counts
      na_positions <- is.na(masked_data.mat1)


      # Impute the missing values
      time.missem[i] <- system.time(
        imputed_data.mat1 <- missMDA::imputePCA(masked_data.mat1, ncomp , method = "EM"))
      imputed_data.mat1 <- t(imputed_data.mat1$completeObs)
      #imputed_data.mat1 <- t(round(imputed_data.mat1$completeObs))
      #imputed_data.mat1[imputed_data.mat1<0] <- 0
      missem.mat1[[i]] <- t(imputed_data.mat1)


      # Update the original dataset with the imputed values

      masked_data.mat1[is.na(masked_data.mat1)] <- imputed_data.mat1[is.na(masked_data.mat1)]
      masked_data.mat1 <- t(masked_data.mat1)

      # Calculate RMSE, MAE, KL distance, KL divergence for met
      error.rmse.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mae.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mape.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))

      time.elapsed.missem <- NULL

      for (j in 1:nrow(origmat1)) {

        for (l in 1:ncol(origmat1)) {

          if (na_positions[l,j]==TRUE) {
          # computing RMSE
          diff3 <- as.numeric(origmat1[j,l] - masked_data.mat1[j,l])
          sqdiff3 <- (diff3)^2
          error.rmse.mat1[j,l] <- sqdiff3


          #computing MAE
          abs3 <- abs(as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))
          error.mae.mat1[j,l] <- abs3

          #computing MAPE
          abs.rap <- abs((as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))/as.numeric(origmat1[j,l]))
          error.mape.mat1[j,l] <- abs.rap

          } # if

        } # for i


      } # for j


      x3 = mean(error.rmse.mat1)
      S3 = sqrt(x3)/(max(origmat1)-min(origmat1))
      rmse.missem.mat1[i]  <- S3

      M3 = mean(error.mae.mat1)/(max(origmat1)-min(origmat1))
      mae.missem.mat1[i] <- M3

      # Accuracy
      # Reg.eval.missEM[i] <- R2(origmat1, imputed_data.mat1)

      # ks.beta  <- stats::ks.test(origbeta, imputed_data.beta, exact=TRUE)$statistic

      # Print summary statistics
      # cat("Summary Statistics After Imputation:\n")
      # print(summary(data))
      # cat("\n")

      cat("RMSE at iteration:\n")
      print(S3)
      cat("\n")

      cat("MAE at iteration:\n")
      print(M3)
      cat("\n")


      # cat("Accuracy:\n")
      # print(c(Reg.eval.missEM[i]))
      # cat("\n")

      cat("time difference at iteration:\n")
      print(time.missem[i])
      cat("\n")


      results.missem <- list(imputed.missem.mat1 =  missem.mat1,
                             RMSE.missem.mat1 = rmse.missem.mat1, MAE.missem.mat1 = mae.missem.mat1,
                             Comp_time = time.missem)


    }# iter


  }# missEM

  else results.missem <- NULL



  origmat1 <- as.data.frame(t(mat1))

  # 10 = kNN impute (run if selected)
  if (10 %in% sel_method) {
    cat("KNN from impute imputation - in progress\n")

    knn2.mat1 <- list()
    rmse.knn2.mat1 <- NULL
    mae.knn2.mat1 <- NULL
    mape.knn2.mat1 <- NULL

    time.knn2 <- NULL


    #origmat1 <- as.data.frame(t(origmat1))
    origmat4 <- t(origmat1)   # NOTE: this method uses the transpose matrix

    # Repeat the process n_iterations times
    for (i in 1:n.iter) {


      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- origmat4
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the mat1
      na_positions <- is.na(masked_data.mat1)


      # Impute the missing values
      time.knn2[i] <- system.time(
        imputed_data.mat1<- impute::impute.knn(data = as.matrix(masked_data.mat1), k = nb))
      imputed_data.mat1 <- imputed_data.mat1$data

      #imputed_data.mat1 <- round(imputed_data.mat1$data)
      #imputed_data.mat1[imputed_data.mat1<0] <- 0
      knn2.mat1[[i]] <- imputed_data.mat1


      # Update the original dataset with the imputed values
      masked_data.mat1[is.na(masked_data.mat1)] <- imputed_data.mat1[is.na(masked_data.mat1)]
      masked_data.mat1 <- t(masked_data.mat1)

      # Calculate RMSE, MAE, KL distance, KL divergence for mat1
      error.rmse.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mae.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mape.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))


      for (j in 1:nrow(origmat1)){

        for (l in 1:ncol(origmat1)) {

          if (na_positions[l,j]==TRUE) {

          # computing RMSE
          diff3 <- as.numeric(origmat1[j,l] - masked_data.mat1[j,l])
          sqdiff3 <- (diff3)^2
          error.rmse.mat1[j,l] <- sqdiff3

          #computing MAE
          abs3 <- abs(as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))
          error.mae.mat1[j,l] <- abs3

          #computing MAPE
          abs.rap <- abs((as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))/as.numeric(origmat1[j,l]))
          error.mape.mat1[j,l] <- abs.rap

          } # if

        } # for i


      } # for j

      x3 = mean(error.rmse.mat1)
      S3 = sqrt(x3)/(max(origmat1)-min(origmat1))
      rmse.knn2.mat1[i]  <- S3

      M3 = mean(error.mae.mat1)/(max(origmat1)-min(origmat1))
      mae.knn2.mat1[i] <- M3

      # Accuracy
      # Reg.eval.knn2[i] <- R2(origmat1, imputed_data.mat1)

      # ks.beta  <- stats::ks.test(origbeta, imputed_data.beta, exact=TRUE)$statistic

      # Print summary statistics
      # cat("Summary Statistics After Imputation:\n")
      # print(summary(data))
      # cat("\n")

      cat("RMSE at iteration:\n")
      print(S3)
      cat("\n")

      cat("MAE at iteration:\n")
      print(M3)
      cat("\n")

      # cat("Accuracy:\n")
      # print(c(Reg.eval.knn2[i]))
      # cat("\n")

      cat("time difference at iteration:\n")
      print(time.knn2[i])
      cat("\n")


      results.knn2 <- list(imputed.knn2.mat1 =  knn2.mat1,
                           RMSE.knn2.mat1 = rmse.knn2.mat1,
                           MAE.knn2.mat1 = mae.knn2.mat1,
                           Comp_time = time.knn2)


    }# iter


  }# knn2

  else results.knn2 <- NULL




  origmat1 <- as.data.frame(t(mat1))



  #  origmat1 <- as.data.frame(t(mat1))

  # 12 = cluster-partition based on knn impute (run if selected)
  if (12 %in% sel_method) {

    cat("Simulation-correlation imputation - in progress\n")

    corr.mat1 <- list()
    rmse.corr.mat1 <- NULL
    mae.corr.mat1 <- NULL
    mape.corr.mat1 <- NULL

    time.corr <- NULL



    # Repeat the process n_iterations times
    for (i in 1:n.iter) {


      cat("Iteration:", i, "\n")

      # Mask the missing values
      masked_data.mat1 <- origmat1
      masked_data.mat1<- prodNA(masked_data.mat1, missing_prop) # Mask some values

      # Assign NA values to the mat1
      na_positions <- is.na(masked_data.mat1)


      stats.m <- stats_cpg2(X=list.cleaned, varselect)
      mat1<- list.cleaned[[varselect]]
      means <- rowMeans(mat1, na.rm=TRUE)
      sds <- apply(mat1, 1, sd, na.rm=T)


      simu3 <- patternCpG(rownum= stats.m$Rows, colnum= stats.m$Columns,
                           cormat= stats.m$Linear_correlation, meanval = means, sdval = sds,
                           dataset = mat1, matrix = "M")

      imputed_data.mat1 <- simu3$Simulated_matrix

      corr.mat1[[i]] <- imputed_data.mat1

      # Impute the missing values
      #time.simu[i] <- system.time(

      # Update the original dataset with the imputed values
      #masked_data.mat1[is.na(masked_data.mat1)] <- imputed_data.mat1[is.na(masked_data.mat1)
      masked_data.mat1[is.na(masked_data.mat1)] <- imputed_data.mat1[is.na(masked_data.mat1)]
      # masked_data.mat1 <- t(masked_data.mat1)


      # Calculate RMSE, MAE, KL distance, KL divergence for mat1
      error.rmse.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mae.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))
      error.mape.mat1 <- matrix(0,nrow(origmat1),ncol(origmat1))


      for (j in 1:nrow(origmat1)){

        for (l in 1:ncol(origmat1)) {

          if (na_positions[j,i]==TRUE){
          # computing RMSE
          diff3 <- as.numeric(origmat1[j,l] - masked_data.mat1[j,l])
          sqdiff3 <- (diff3)^2
          error.rmse.mat1[j,l] <- sqdiff3

          #computing MAE
          abs3 <- abs(as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))
          error.mae.mat1[j,l] <- abs3

          #computing MAPE
          abs.rap <- abs((as.numeric(origmat1[j,l] - masked_data.mat1[j,l]))/as.numeric(origmat1[j,l]))
          error.mape.mat1[j,l] <- abs.rap

          } # if

        } # for i


      } # for j

      x3 = mean(error.rmse.mat1)
      S3 = sqrt(x3)/(max(origmat1)-min(origmat1))
      rmse.corr.mat1[i]  <- S3

      M3 = mean(error.mae.mat1)/(max(origmat1)-min(origmat1))
      mae.corr.mat1[i] <- M3



      # Accuracy
      # Reg.eval.simu[i] <- R2(origmat1, imputed_data.mat1)

      # ks.beta  <- stats::ks.test(origbeta, imputed_data.beta, exact=TRUE)$statistic

      # Print summary statistics
      # cat("Summary Statistics After Imputation:\n")
      # print(summary(data))
      # cat("\n")

      cat("RMSE at iteration:\n")
      print(S3)
      cat("\n")

      cat("MAE at iteration:\n")
      print(M3)
      cat("\n")

      # cat("Accuracy:\n")
      # print(c(Reg.eval.simu[i]))
      # cat("\n")

      cat("time difference at iteration:\n")
      print(time.corr[i])
      cat("\n")


      results.corr <- list(imputed.corr.mat1 =  corr.mat1,
                                RMSE.corr.mat1 = rmse.corr.mat1,
                                MAE.corr.mat1 = mae.corr.mat1,
                                Comp_time = time.corr)


    }# iter


  }# simu

  else results.corr <- NULL

  # Final results to be returned

  return(list(MEAN_imputation = results.mean,
              PPCA_imputation = results.ppca,
              BPCA_imputation = results.bpca,
              MISSF_imputation = results.missf,
              SVD_imputation = results.svd,
              MISSMDA_imputation = results.missmda,
              MISSEM_imputation = results.missem,
              KNN2_imputation = results.knn2,
              CORR_imputation = results.corr))





}      # function




