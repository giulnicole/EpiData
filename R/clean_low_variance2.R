

clean_low_variance2 <- function(list.cleaned, sd_quantile=0.05, max_na = 0.25, varselect=2) {
  
  library(matrixStats)
  library(impute)
  
  n <- nrow(list.cleaned[[varselect]])
  
  # FIRST PART: cleaning low variance for meth counts 
  
  X <- list.cleaned[[varselect]]
  
  # Convert the object to a matrix (if needed)
  is_matrix <- is.matrix(X)
  if (!is_matrix) {
    X <- as.matrix(X)
  }
  
  row_sd = rowSds(X, na.rm = T)
  l = which(abs(row_sd) <= 1e-10)
  X2<- X[-l,]
  filtered1<- rownames(X2)
  
  #X2 = X[!l, , drop = T]
  #row_sd = row_sd[!l]
  
  l <- as.matrix(l)
  
  #dim(l)[1]
  
  if(dim(l)[1] !=0) {
    cat(paste(dim(l)[1], " rows have been removed with zero variance.\n"))
  }
  
  qa = quantile(unique(row_sd), sd_quantile, na.rm = TRUE)
  l2 = which(row_sd < qa)
  l2 <- as.matrix(l2)
  
  if(dim(l2)[1] !=0) {
    cat(paste(dim(l2)[1], " rows have been removed with too low variance (sd <= ", sd_quantile, " quantile).\n"))
  }
  
  
  X3 <- X2[-l2,]
  filtered2<- rownames(X3)
  
  throw<- which(is.na(rowMeans(X3, na.rm = T)))
  X4 <- X3[-throw,]
  
  throw <- as.matrix(throw)
  if(dim(throw)[1] !=0) {
    cat(paste(dim(throw)[1], " rows have been removed because all entries were NAs.\n"))
  }
  

  #rownames(X3) <- sites_filtered2
  #return(X4)
  
  clean.met <- X4
  
  
  # SECOND PART: Subsetting
  

  # Subset the same unmet and beta and m 
  clean.unmet <- list.cleaned[[3]]
  clean.unmet <- clean.unmet[rownames(clean.unmet) %in% rownames(clean.met),]
  #dim(clean.unmet)
  
  clean.beta <- list.cleaned[[4]]
  clean.beta <- clean.beta[rownames(clean.beta) %in% rownames(clean.met),]
  #dim(clean.beta)
  
  clean.m <- list.cleaned[[5]]
  clean.m <- clean.m[rownames(clean.m) %in% rownames(clean.met),]
  #dim(clean.m)
  
  
  list.cleaned2 <- list(met.mat2 = clean.met, unmet.mat2 = clean.unmet, 
                        beta.mat2 = clean.beta, m.mat2 = clean.m)
  
  

  
  return(list.cleaned2)
  
  
}
