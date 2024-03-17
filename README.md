## Description
Welcome to the `EpiData` project, which is a package for enhancing reliability in DNA methylation analysis (contributing to normalize M value disstribution per CpG) and propose a novel approach for imputing missing data from bisulfite sequencing experiments. 
Welcome to our EpiData project! This is a package for enhancing reliability in DNA methylation analysis (contributing to normalize M value disstribution per CpG) and propose a novel approach for imputing missing data from bisulfite sequencing experiments. 
We hope you enjoy and we look forward to your contributions!


##Installing
Please install devtools if you haven't yet.
```{r setup}
install.packages("devtools")
```

Required Bioconductor packages:
```{r setup}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```
BiocManager::install(c("qvalue", "rain", "limma"))

## Contributing
We welcome any and all contributions. Here are some ways you can get started:

**Report bugs**: please, feel free to report, if you encounter any bugs. Open up an issue and let us know the problem.

**Contribute code**: if you are a developer and want to contribute, follow the instructions below to get started!

**Suggestions**: if you don't want to code but have some awesome ideas, open up an issue explaining some updates or imporvements you would like to see!

**Documentation**: If you see the need for some additional documentation, feel free to add some!

## Fork this repository
Clone the forked repository
Add your contributions (code or documentation)
Commit and push
Wait for pull request to be merged

## Get started with the analysis

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

clean.coverage <- cleanCovMat(input.obj = input.list, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)
```

**Cleaning step 2**

```{r message=FALSE, warning = FALSE}
clean.out<- cleanOutliers(filtered.obj=clean.coverage, outlier_threshold=0, remove_outliers = T)
```
After cleaning the data, the pipeline follows with studying each chromosome seprately.the following function helps in separating the information of BS experiment per each chromosome. 

**Splitting dataset per chromosomes**

The input should be the assay method of the SummarizedExperiment object, if present. Otherwise a list with the three matrices described before.
```{r}
splitted<- split5MatXChrom(clean.out)
```
The result is a list where we have the total coverage matrix, methylated maytrix and unmethylated matrix with all chromoeomse, but also a list called *final* which contains all these three matrices already divided per each chromosome. 

Here is shown the example pipeline through each one of the chromosomes. 

## Parallelization
After cleaning steps, CpGs' matrices are divided per each chromosome with the function

```{r warning=FALSE, message=FALSE}
library(BiocParallel)
library(matrixcalc)
```


### Calculating statistics per CpG 
```{r warning=FALSE}
input.stats<- splitted$final

# Statistics
stats2<- ParallStats(input.stats)
```


Here is shown the example of `statsCpG` function (used by ParallStats) through each one of the chromosomes. 

```{r message=FALSE, warning = FALSE}
final2<- splitted$final

stats<- statsCpG(cleaned.obj= final2[[1]], varselect = 5, plot = TRUE)
stats$Barplot
```



### Imputing data based on correlation
```{r  warning=FALSE}
simu2<- ParallSimu(stats2)
```

Here is shown the example of `imputeCorr` function (used by ParallSimu) through each one of the chromosomes. Note that the input could be *pairwise* as well as *meanCpG*.


Imputation based on correlation structure 
```{r message=FALSE, warning = FALSE}
obj.imp.pairw <- imputeCorr(cleaned.obj=stats,
                    matrix="M",
                    varselect = 5,
                    correlation.type="pairwise")


obj.imp.mean <- imputeCorr(cleaned.obj=stats,
                    matrix="M",
                    varselect = 5,
                    correlation.type="meanCpG")


```

```{r}
imputed.mat <- extractFinalMat(list.imputed=simu2)
```


The default missing values' investigation is performed for M values. Note that percentage of missing values in M values are higher than counts' missing value percentages.


# Comparison of imputation methods

```{r, message=FALSE, warning=FALSE}
stats3<- lapply(clean.out, na.omit)

imp.M <- repNA(list.cleaned=stats3,
                       missing_prop= 0.2, varselect=5,
                       n.iter= 2, sel_method=c(1:10),
                       trees= 50, nb= 10, ncomp= 2, matrix="M")

```


```{r, warning=FALSE, message=FALSE}
measure_imp_m <- measureAccuracy(imp.M)
measure_imp_m$Boxplot.rmse
measure_imp_m$Boxplot.mae
```
