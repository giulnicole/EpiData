
adjust_counts <- function(coverage, methylated, unmethylated){
  
  
  # NA sicuri
  # Create a new matrix based on conditions
  met2 <- methylated # Initialize 
  
  for (i in 1:nrow(met2)) {
    for (j in 1:ncol(met2)) {
      if (is.na(co[i, j]) && un[i, j] == 0) {
        met2[i, j] <- NA
      } else {
        met2[i, j] <-  methylated[i, j]
      }
    }
  }
  
  
  unmet2 <- unmethylated # Initialize 
  
  for (i in 1:nrow(unmet2)) {
    for (j in 1:ncol(unmet2)) {
      if (is.na(co[i, j]) && is.na(met2[i, j])) {
        unmet2[i, j] <- NA
      } else {
        unmet2[i, j] <-  unmethylated[i, j]
      }
    }
  }
  
  
  met3 <- met2 # Initialize 
  

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
  
  
  
  unmet3 <- unmet2 # Initialize 
  

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
  
  
  
}  # function