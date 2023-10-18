
# setworking directory in the EpiData package

# NOTE 1 for Ali: I need to connect DESCRIPTION file to each function I created
# some libraries are not read properly although functions are checked (Fix this)
# library(devtools)
# library(magrittr)
# library(dplyr)


#-------  Data
setwd("C:/Users/gnbal/Desktop/GitHub/EpiData")
load("data/chr22.RData")


# Sources until functions are not installed


#------  Example clean_matrix
source("R/clean_matrix_cov.R")

# Cleaning the coverage matrix
X2 <- X %>% replace(is.na(.), 0)
Y2 <- Y %>% replace(is.na(.), 0)
Z2 <- Z %>% replace(is.na(.), 0)
clean.coverage2 <- clean_matrix_cov(X2, Y2, Z2, max_na_cpg = 0.5, max_na_ind = 0.2)


#------  Example clean_low_coverage fopr methylated counts
source("R/clean_counts_outliers.R")

# Cleaning for low variance of the methylated counts (excessive NAs that occurr after clean low coverage)
clean.new <- clean_counts_outliers(list.cleaned=clean.coverage2$Cleaned1, outlier_threshold=5, remove_outliers = TRUE)


# match added to the function
# clean low variance cleans the excessive NAs on methylated, and then match the same CpGs in coverage, unmethylated counts, beta matrix and M matrix
#

# NOTE 2 the examples are then to be implemented in the function body, after fixing how data are read

# NOTE 3 After fixed NOTE 1 and 2, I will add the imputation + clustering according to similarity script
# NOTE 4 after NOTE3 is completed I will parallelize the steps for each chromosome
