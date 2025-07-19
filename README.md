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
data("meth_data")
print(meth_data)
```

The data contains a list of assays from experiment obtained by bisulfite sequencing 43 individuals. 

As an example data derives from chromosome 22 (the same pipeline is applied to the whole chromosome matrix if it s provided): 
=======


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

input.list<- assays(rangedObject)

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




