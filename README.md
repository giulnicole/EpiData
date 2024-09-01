## Description of the package
Welcome to the `EpiData` project, which is a package for enhancing reliability in DNA methylation analysis (contributing to normalize M value disstribution per CpG) and propose a novel approach for imputing missing data from bisulfite sequencing experiments. 
We hope you enjoy and we look forward to your contributions!


## Contributing
We welcome any and all contributions. Here are some ways you can get started:

**Report bugs**: please, feel free to report, if you encounter any bugs. Open up an issue and let us know the problem.

**Contribute code**: if you are a developer and want to contribute, follow the instructions below to get started!

**Suggestions**: if you don't want to code but have some awesome ideas, open up an issue explaining some updates or improvements you would like to see!

**Documentation**: if you see the need for some additional documentation, feel free to add some!


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


## Loading the package and useful libraries
```{r, warning=FALSE}
library(EpiData)
```
```{r, warning=FALSE, message=FALSE}
library(dplyr)
library(tidyverse)
library(purrr)
library(dplyr)
library(SummarizedExperiment)
library(GenomicRanges)
library(GenomicFeatures)
library(MASS)
library(stats)
```

# Get started with the analysis
## Loading data 
```{r}
data("rangedObject")
```
<<<<<<< HEAD
The data contains a list of assays from experiment obtained by bisulfite sequencing 43 individuals. 

As an example data derives from chromosome 22 (the same pipeline is applied to the whole chromosome matrix if it s provided): 
=======


The dataset contains the LargeSummarizedExperiment object with the following list of assays from experiment obtained by bisulfite sequencing 43 individuals. As an example data derives from chromosome 22: 
>>>>>>> dac400e184f2dba6db4c1cb845c97535c0a9ee6d

- coverage counts matrix (Coverage_matrix)

- methylated counts matrix (Met_matrix)

- unmenthylated counts matrix (Unmet_matrix)

In this objects exact locations and CpGs' names are present (labeled as chr:bp). 

Assume that from the experiment we endend up with a dataset not only contanining chromosome 22, the optimal situation would be split dataset per each chromosome and perform a parallelized analysis. 


Note that if data is not in this format and sequencing gave the three matrices, this operation can be done with: 
```{r}
n_sites <- 100
n_indiv <- 10

# Generate random data for methylated counts
set.seed(123)  # For reproducibility
methylated <- as.data.frame(matrix(sample(0:100, n_sites * n_indiv, replace = TRUE), nrow = n_sites, ncol = n_indiv))

# Generate random data for unmethylated counts
unmethylated <- as.data.frame(matrix(sample(0:100, n_sites * n_indiv, replace = TRUE), nrow = n_sites, ncol = n_indiv))

# Total coverage is the sum of methylated and unmethylated counts
coverage <- as.data.frame(methylated + unmethylated)

counts.list <- list(Coverage_matrix=coverage, Met_matrix=methylated, Unmet_matrix=unmethylated)
```


Without the ranged object the final input list can be obtained with `GRconvresion` function

```{r}
gr.object <- GRconversion(counts.list)

input.list <-  assays(gr.object)
```


Otherwise the ranged object can be used and the subset of 100 sites and 10 individuals can be performed:
```{r message=FALSE, warning = FALSE}

# with the ranged object
input.list <- assays(rangedObject)

subset_list <- lapply(input.list, function(df) df[1:100, 1:10])
```


## Pipeline

<<<<<<< HEAD
### Part 1: cleaning the whole matrix of the CpGs
**Cleaning step 1: coverage**

```{r}

input.list<- assays(rabgedObject)

clean.coverage <- cleanCovMat(input.obj = input.list, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)
```

**Cleaning step 2: outliers**

```{r}
clean.out<- cleanOutliers(filtered.obj=clean.coverage, outlier_threshold=0, remove_outliers = T)
```

## Part 2: methylation values

```{r}
values <- MethilationValues(clean.coverage)
```


## Part 3: imputation


### Calculating statistics per CpG and imputing data

```{r}
stats<- missingStats(values)
```

### Imputing data
```{r}
simu <- imputeCorr(stats)
```

=======
### Filtering 1
```{r message=FALSE, warning = FALSE}
clean.coverage <- cleanCovMat(input.obj = subset_list, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)
```


Function `cleanCovMat` filters coverage matrix, consideriign the following assumptions as default: 

1) if the coverage counts are zero, then NAs are put, because coverage is the sum given by unmethylated counts and methylated counts, hence it should not be equal to zero;

2)	CpGs (rows) with missing values greater than 50% are removed;

3)	individuals (columns) with a missing percentage greater than 20% are removed.

### Filtering 2
```{r message=FALSE, warning = FALSE}
clean.out <- cleanOutliers(filtered.obj=clean.coverage, outlier_threshold=0, remove_outliers = TRUE)
```

The function returns filtered counts' matrices for excessive NAs in methylated and unmenhylated counts' matrices. 

- `cleanOutliers` filters per rows (CpGs) and columns (individuals). First the function removes entries which may still have undesirable high percentages of NAs in methylated or unmethylated matrices, not yet cleaned from cleaning step 1. Secondly, the function clean outliers, since they can significantly impact the results of statistical analyses. 

### Calculating Beta and M values

```{r}
values <- MethilationValues(clean.coverage)
```

`MethilationValues`  function computes the matrices of beta values and M values. 

### Computing statistics and imputation

The input should be the assay method of the SummarizedExperiment object, if present. Otherwise a list with the three matrices described before.

The result is a list where we have the total coverage matrix, methylated maytrix and unmethylated matrix with all chromoeomse, but also a list called *final* which contains all these three matrices already divided per each chromosome. 


## Second part: imputation
After cleaning steps, CpGs' matrices are divided per each chromosome with the function

```{r warning=FALSE, message=FALSE}
library(BiocParallel)
library(matrixcalc)
```

### Calculating statistics per CpG 
```{r warning=FALSE}
stats <- statsCpG(values)
```

```{r}
values2 <- list(values, values)
meta.stat <- lapply(values2, statsCpG)
```

### Imputing data based on correlation
```{r  warning=FALSE}
simu <- imputeCorr(stats)
```

Here is shown the example of `imputeCorr` function (used by ParallSimu) through each one of the chromosomes. Note that the input could be *pairwise* as well as *meanCpG*.

When computing correlation matrix in presence of missing data, the correlation pattern behind the data may result affected by NAs presence. Since out aim is to impute the M values matrix, the  of CpGs' correlation can be reconstructed via computing the correlation matrix once imputed the mean vM value per each CpG, not to change the similarity between CpGs or via pairwise correlation. 


Here is shown the example of `imputeCorr` function (used by `eigenImpute`) through each one of the chromosomes. Note that the input could be *pairwise* as well as *meanCpG*.

The default missing values' investigation is performed for M values. Note that percentage of missing values in M values are higher than counts' missing value percentages.

>>>>>>> dac400e184f2dba6db4c1cb845c97535c0a9ee6d


