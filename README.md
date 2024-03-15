## Description
Welcome to our EpiData project! This is a package for enhancing reliability in DNA methylation analysis (contributing to normalize M value disstribution per CpG) and propose a novel approach for imputing missing data from bisulfite sequencing experiments. 
We hope you enjoy and we look forward to your contributions!

## Contributing
We welcome any and all contributions. Here are some ways you can get started:

**Report bugs**: If you encounter any bugs, please let us know. Open up an issue and let us know the problem.

**Contribute code**: If you are a developer and want to contribute, follow the instructions below to get started!

**Suggestions**: If you don't want to code but have some awesome ideas, open up an issue explaining some updates or imporvements you would like to see!

**Documentation**: If you see the need for some additional documentation, feel free to add some!
Instructions

## Fork this repository
Clone the forked repository
Add your contributions (code or documentation)
Commit and push
Wait for pull request to be merged

## Get started

**Epidata**
```{r setup}
library(EpiData)
```

**Required packages**
```{r, warning=FALSE, message=FALSE}
library(devtools)
library(magrittr)
library(dplyr)
library(psych)
library(kableExtra)
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

**Cleaning step 1**

```{r message=FALSE, warning = FALSE}
input.list<- assays(dati)

clean.coverage2 <- cleanCovMat(input.obj = input.list, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)
```


