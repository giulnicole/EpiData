## Description of the package

Welcome to the `EpiData` project, which is a package for enhancing reliability in DNA methylation analysis (contributing to normalize M value disstribution per CpG) and propose a novel approach for imputing missing data from bisulfite sequencing experiments. We hope you enjoy and we look forward to your contributions!

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
BiocManager::install("Epidata")
```


## Fork this repository
Clone the forked repository
Add your contributions (code or documentation)
Commit and push
Wait for pull request to be merged


## Loading the package and useful libraries
```{r, warning=FALSE}
library(EpiData)
```

```{r, warning=FALSE, message=FALSE}
library(dplyr)
library(tidyverse)
library(purrr)
library(SummarizedExperiment)
library(GenomicRanges)
library(GenomicFeatures)
library(MASS)
library(stats)
```

# Get started with the analysis
## Loading data 
```{r}
data("meth_data")
```

The data contains a list of assays from experiment obtained by bisulfite sequencing 43 individuals. 

As an example data derives from chromosome 22 (the same pipeline is applied to the whole chromosome matrix if it is provided): 



The dataset contains the LargeSummarizedExperiment object with the following list of assays from experiment obtained by bisulfite sequencing 43 individuals. As an example data derives from chromosome 22: 


- coverage counts matrix (Coverage_matrix)

- methylated counts matrix (Met_matrix)

- unmenthylated counts matrix (Unmet_matrix)

In this objects exact locations and CpGs' names are present (labeled as chr:bp). 

Assume that from the experiment we endend up with a dataset not only contanining chromosome 22, the optimal situation would be split dataset per each chromosome and perform a parallelized analysis. 



## Pipeline

### Part 1: cleaning the whole matrix of the CpGs
**Cleaning step 1: coverage**

```{r}

input.list<- assays(meth_data)

clean.coverage <- cleanCovMat(input.obj = input.list, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)
```

**Cleaning step 2: outliers**

```{r}
clean.out<- cleanOutliers(filtered.obj=clean.coverage, outlier_threshold=0, remove_outliers = T)
```

## Part 2: methylation values

```{r}
mat.values <- methValues(clean.coverage)
```


## Part 3: imputation


### Calculating statistics per CpG and imputing data

```{r}
stat.result <- missingStats(list(mat.values))
```

### Imputing data
```{r}
imputed <- batchImputeCorr(stat.result)
```




