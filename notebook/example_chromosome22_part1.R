
# setworking directory in the EpiData package
library(devtools)
library(magrittr)
library(dplyr)


#-------  Data
setwd("C:/Users/gnbal/Desktop/GitHub/EpiData")
load("data/chr22.RData")



# Sources until functions are not installed


#------  Example clean_matrix
source("R/clean_matrix.R")

# Cleaning the coverage matrix
clean.coverage2 <- clean_matrix(X, Y, Z, max_na_cpg = 0.5, max_na_ind = 0.2)



#------  Example clean_low_coverage
source("R/clean_low_variance2.R")

# Cleaning for low variance of the methylated counts (excessive na)
clean.new <- clean_low_variance2(list.cleaned=clean.coverage2, sd_quantile=0.05, max_na = 0.25, 2)


#which(is.na(rowMeans(clean.unmet, na.rm = T))) #ok

# match added to the function



