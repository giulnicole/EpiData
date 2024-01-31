
# Setworking directory in the EpiData package

# some libraries are not read properly although functions are checked (Fix this)
# library(devtools)
# library(magrittr)
# library(dplyr)
# library(kableExtra)
# library(GenomicRanges)
# library(SummarizedExperiment)
# library(tidyverse)
# library(missForest)
# library(FNN)
# library(pcaMethods)
# library(softImpute)
# library(impute)

#-------           Data        ----------------------------

data("matrices")

#save(dati, file = "data/matrices.rda")

# Sources until functions are not installed
devtools::load_all()


#------  Example cleanMatCov   -----------------------------

# Cleaning low counts for coverage
clean.coverage2 <- cleanCovMat(input.obj = dati, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)

#saved to rda file
#save(clean.coverage2, file = "data/clean.coverage2.rda")

#------  Example cleanOutliers    --------------------------

# Cleaning for low variance of the methylated counts (excessive NAs that occurr after clean low coverage)
clean.out<- cleanOutliers(filtered.obj=clean.coverage2, outlier_threshold=5, remove_outliers = F)

#saved to rda file
#save(clean.new, file = "data/clean.new.rda")

# match added to the function
# clean low variance cleans the excessive NAs on methylated, and then match the same CpGs in coverage, unmethylated counts, beta matrix and M matrix
#
#------  Example statsCpG      -------------------------------
stats<- statsCpG(cleaned.obj= clean.out, varselect = 5)

#save to rda file
#save(stats, file = "data/stats.rda")

 #------  Example patternCpG     -----------------------------
# library(latentcor)
# library(propagate)

simu3 <- patternCpG(cleaned.obj=stats,
                    matrix="M",
                    varselect = 5)

#save to rda file
#save(simu3, file = "data/simulated.rda")

#---------------------------------------------------------------------------
