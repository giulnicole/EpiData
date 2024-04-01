## Description of the package
Welcome to the `EpiData` project, which is a package for enhancing reliability in DNA methylation analysis (contributing to normalize M value disstribution per CpG) and propose a novel approach for imputing missing data from bisulfite sequencing experiments. 
We hope you enjoy and we look forward to your contributions!


## Contributing
We welcome any and all contributions. Here are some ways you can get started:

**Report bugs**: please, feel free to report, if you encounter any bugs. Open up an issue and let us know the problem.

**Contribute code**: if you are a developer and want to contribute, follow the instructions below to get started!

**Suggestions**: if you don't want to code but have some awesome ideas, open up an issue explaining some updates or imporvements you would like to see!

**Documentation**: If you see the need for some additional documentation, feel free to add some!


## Installing the package
Please install devtools and Bioconductor if you haven't them yet.

Required Bioconductor/devtools packages:
```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

install.packages("devtools")
```

Installing EpiData package:
```{r}
BiocManager::install(c("Epidata"))
```



## Fork this repository
Clone the forked repository
Add your contributions (code or documentation)
Commit and push
Wait for pull request to be merged

## Get started with the analysis

**Epidata**
```{r}
library(EpiData)
```

**Required packages**
```{r}
library(devtools)
library(magrittr)
library(purrr)
library(dplyr)
library(psych)
library(GenomicRanges)
library(SummarizedExperiment)
library(radiant.data)
```

**Loading data** 
```{r}
data("matrices")
```
The data contains a list of assays from experiment obtained by bisulfite sequencing 43 individuals. 
As an example data derives from chromosome 22 (the same pipeline is applied to the whole chromosome matrix if it s provided): 

- coverage counts' matrix (Coverage_matrix)

- methylated counts' matrix (Met_matrix)

- unmenthylated counts' matrix (Unmet_matrix)

In this objects exact locations and CpGs' names are present (labeled as chr:bp). 


## Pipeline

### Part 1: cleaning the whole matrix of the CpGs
**Cleaning step 1**

```{r}
input.list<- assays(dati)

clean.coverage <- cleanCovMat(input.obj = input.list, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)
```

**Cleaning step 2**

```{r}
clean.out<- cleanOutliers(filtered.obj=clean.coverage, outlier_threshold=0, remove_outliers = T)
```
After cleaning the data, the pipeline follows with studying each chromosome seprately.the following function helps in separating the information of BS experiment per each chromosome. 




## Part2: parallelization per chromosome for imputation

Here is shown the example pipeline through each one of the chromosomes. 


```{r}
library(BiocParallel)
library(matrixcalc)
```

### Calculating statistics per CpG and imputing data
```{r}
stats<- ParallStats(clean.out)
```


```{r}
simu<- ParallSimu(stats)
```

Here is shown the example of `imputeCorr` function (used by ParallSimu) through each one of the chromosomes. Note that the input could be *pairwise* as well as *meanCpG*.



The default missing values' investigation is performed for M values. Note that percentage of missing values in M values are higher than counts' missing value percentages.


# Comparison of imputation methods

```{r}
data.missing <- lapply(clean.out, na.omit)

imp.M <- repNA(list.cleaned=data.missing ,
                       missing_prop= 0.2, varselect=5,
                       n.iter= 2, sel_method=c(1:10),
                       trees= 50, nb= 10, ncomp= 2, matrix="M")

```

## Measuring and visualizing the accuracy of the imputation methods

```{r}
measure_imp_m <- measureAccuracy(imp.M)
measure_imp_m$Boxplot.rmse
measure_imp_m$Boxplot.mae
```
