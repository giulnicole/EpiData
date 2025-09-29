#' @title measureAccuracy
#'
#' @description A function for computing accuracy measure (RMSE and MAE computation) and perform the Kolmogorov-Smirnov test per each CpG and produce boxplots for comaprison.
#' \code{} computes the accuracy measure of imputation methods.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot theme element_text ggtitle xlab
#' @importFrom dplyr bind_rows
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom viridis scale_fill_viridis
#' @import stats
#'
#' @param list.imputed output of imputed datasets performed by repNA.
#'
#' @name measureAccuracy
#'
#' @return
#'  List with
#' \item{Accuracy_measure}{list with the available methods for imputation with respective errors calculated.}
#' \item{Dataset}{Dataset to create pots with all the errors computed.}
#' \item{Boxplot.rmse}{Boxplots comparing RNSEs.}
#' \item{Boxplot.mae}{Boxplots comparing MAEs.}
#'
#' @examples
#' \dontrun{
#' data("meth_data")
#' input.list<- assays(meth_data)
#' clean.coverage <- cleanCovMat(input.obj=input.list, max_na_cpg = 0.5,
#'                                  max_na_ind = 0.2,  cpg_removal_threshold = 10)
#'  clean.out <- cleanOutliers(filtered.obj=clean.coverage, outlier_threshold=5,
#'                              remove_outliers = TRUE)
#'  stats <- lapply(clean.out, na.omit)
#'
#'  imp.M <- repNA(list.cleaned=stats,
#'               missing_prop= 0.2, varselect=5,
#'               n.iter= 2, sel_method=c(1:10),
#'               trees= 50, nb= 10, ncomp= 2, matrix="M")
#'
#' measure_imp_m <- measureAccuracy(imp.M)
#' }
#'
#' @export
measureAccuracy <- function(list.imputed){


  name <- NULL
  Method <- NULL
  replication <- NULL
  RMSE <- NULL
  MAE <- NULL
  Method <- NULL
  Replication <- NULL

  RMSE.val <- list()
  MAE.val <- list()

  # Initialize lists
  RMSE.mean <- list()
  RMSE.sd <- list()
  MAE.mean <- list()
  MAE.sd <- list()
  Accuracy <- list()

  # Initialize an empty data frame
  df <- data.frame()


  for (k in 1:length(list.imputed)) {
    call <- names(list.imputed)[k]

    name <- gsub("_.*", "", call)
    Method <- name  # Assuming Method is a vector

    call2 <- tolower(name)

    if (length(list.imputed[[k]]) != 0) {
      RMSE.val <- list.imputed[[call]][[paste0("RMSE.", call2, ".mat1")]]
      MAE.val <- list.imputed[[call]][[paste0("MAE.", call2, ".mat1")]]

      # Calculate mean and sd
      RMSE.mean[[k]] <- mean(RMSE.val)
      RMSE.sd[[k]] <- sd(RMSE.val)
      MAE.mean[[k]] <- mean(MAE.val)
      MAE.sd[[k]] <- sd(MAE.val)

      n_rep <- length(RMSE.val)

      # Create a data frame for Accuracy
      data <- data.frame(
        RMSE = RMSE.val,
        MAE = MAE.val,
        Method = rep(call2, n_rep),
        Replication = seq(1:n_rep)
      )
      Accuracy[[k]] <- data

      # Combine data with df
      df <- bind_rows(df, data)
    } else {
      Accuracy[[k]] <- NULL
    }
  }


  boxplot1<- df %>%
    ggplot( aes(x=Method, y=RMSE, fill=Method)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("RMSE boxplot") +
    xlab("")


  boxplot2<- df %>%
    ggplot( aes(x=Method, y=MAE, fill=Method)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("MAE boxplot") +
    xlab("")




  results <- list(Accuracy_measures=Accuracy, Dataset=df,
             Boxplot.rmse=boxplot1, Boxplot.mae=boxplot2)


    return(results)

}




