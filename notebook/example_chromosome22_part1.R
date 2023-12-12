
# setworking directory in the EpiData package

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

#-------  Data

data("matrices")

# Sources until functions are not installed
devtools::load_all()

#------  Example cleanMatCov

# Cleaning low counts for coverage
clean.coverage2 <- cleanCovMat(X2, Y2, Z2, max_na_cpg = 0.5, max_na_ind = 0.2)


#------  Example cleanOutliers

# Cleaning for low variance of the methylated counts (excessive NAs that occurr after clean low coverage)
clean.new <- cleanOutliers(list.cleaned=clean.coverage2$Cleaned1, outlier_threshold=5, remove_outliers = TRUE)

# match added to the function
# clean low variance cleans the excessive NAs on methylated, and then match the same CpGs in coverage, unmethylated counts, beta matrix and M matrix
#

#------  Example statsCpG
stats<- stats_cpg2(clean.new$Cleaned2, varselect = 5)


#------  Example patternCpG
# library(latentcor)
# library(propagate)

mat.m <- clean.new$Cleaned2$M_matrix
means <- rowMeans(mat.m, na.rm=TRUE)
sds <- apply(mat.m, 1, sd, na.rm=T)

simu3 <- patternCpG(rownum= stats$Rows, colnum= stats$Columns,
                     cormat= stats$Linear_correlation, meanval = means, sdval = sds,
                     dataset = mat.m, matrix = "M")



# NA replication on different methods
mat2<- na.omit(mat.m)
list.new<- lapply(clean.new$Cleaned2, na.omit)

imp.m <- NAreplicating(list.cleaned=list.new, varselect=5, missing_prop= 0.3,
                       scale= TRUE, n.iter= 10, sel_method= c(1:7,10,12),
                       trees= 50, nb= 10, ncomp= 2)



# Accuracy meausre
measure_imp_met <- accuracy_measure2(imp.m)
measure_imp_met$Boxplot.rmse
measure_imp_met$Boxplot.mae



# Pheatmap of the inizitial distance between M values of CpGs
library(apcluster)
m.cleaned2 <- list.new$M_matrix
m.d0 <- negDistMat(m.cleaned2)
pheatmap::pheatmap(as.matrix(m.d0), cluster_rows = FALSE, cluster_cols = FALSE, na_col = "black",
                   display_numbers=F, show_colnames = F, show_rownames = F)


stats<- stats_cpg2(list.new, varselect = 5)

#------  Example patternCpG
mat.m <- list.new$M_matrix
means <- rowMeans(mat.m, na.rm=TRUE)
sds <- apply(mat.m, 1, sd, na.rm=T)

simu3 <- patternCpG(rownum= stats$Rows, colnum= stats$Columns,
                    cormat= stats$Linear_correlation, meanval = means, sdval = sds,
                    dataset = mat.m, matrix = "M")

m.cleaned3 <- list.new$M_matrix
m.d3 <- negDistMat(m.cleaned3)
pheatmap::pheatmap(as.matrix(m.d3), cluster_rows = FALSE, cluster_cols = FALSE, na_col = "black",
                   display_numbers=F, show_colnames = F, show_rownames = F)
