
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
input.list<- list(X, Y, Z)
names(input.list)<- c("Coverage_matrix", "Met_matrix", "Unmet_matrix")
dati2<- GRconversion2(input.list)


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

# correlation pattern via blocks adn submatrices (bigcor)
simu3 <- patternCpG(cleaned.obj=stats,
                    matrix="M",
                    varselect = 5)



#save to rda file
#save(simu3, file = "data/simulated.rda")

#---------------------------------------------------------------------------

############################ Parallelization

#save(clean.out, file="data/clean.out.rda")

data('clean.out')

# Data
list.cleaned1<- clean.out
list.cleaned2<- clean.out
list.cleaned<- list(list.cleaned1$Output_outliers@assays@data@listData, list.cleaned2$Output_outliers@assays@data@listData)

# Filtering coverage matrix
# Cleaning outliers (0->NA)

char <- c(1,2)

names(list.cleaned) <- char

input.stats<- list()

for (i in 1:length(list.cleaned)){

  input.stats[[i]] <- GRconversion2(list.cleaned[[i]])

  names(input.stats)[[i]] <- names(list.cleaned)[[i]]
}


# Statistics

stats2<- ParallStats(input.stats)

# Imputation

simu.mix2 <- ParallSimu(list.stats = stats2)










