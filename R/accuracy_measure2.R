

# Accuracy script



accuracy_measure2 <- function(list.imputed){
  
  library(ggplot2)
  library(tidyverse)
  library(tidyverse)
  library(hrbrthemes)
  library(viridis)
  
  
  name <- NULL
  Method <- NULL
  replication <- NULL
  Rmse <- NULL
  Mae <- NULL
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





