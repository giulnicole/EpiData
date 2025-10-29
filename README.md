
# Description of the package

Welcome to the `EpiData` project, which is a package for enhancing reliability in DNA methylation analysis (contributing to normalize M value disstribution per CpG) and propose a novel approach for imputing missing data from bisulfite sequencing experiments. 
Accurate imputation of missing methylation values is a crucial step in bisulfite sequencing (BS-seq) analysis, as incomplete data can bias downstream association or differential methylation studies.
EpiData provides a modular and reproducible framework for preprocessing and imputation of BS-seq data.

The main objectives of EpiData are:

- to ensure data integrity and quality before statistical modeling

- to provide parallelized and scalable computation for large methylome datasets

- to allow flexible integration with downstream pipelines (e.g., differential methylation, EWAS, or multi-omics integration)


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


## Loading the package 

```{r, warning=FALSE}
library(EpiData)
```

# Get started with the analysis

## Loading data 
```{r}
data("bs_list")
```

The data contains a list of assays from experiment obtained by bisulfite sequencing 43 individuals. 

As an example data derives from chromosome 22 (the same pipeline is applied to the whole chromosome matrix if it is provided): 


The dataset contains the LargeSummarizedExperiment object with the following list of assays from experiment obtained by bisulfite sequencing 43 individuals. As an example data derives from chromosome 22: 

- coverage counts matrix (Coverage_matrix)

- methylated counts matrix (Met_matrix)

- unmenthylated counts matrix (Unmet_matrix)

In this objects exact locations and CpGs' names are present (labeled as chr:bp). 

Assume that from the experiment we endend up with a dataset not only contanining chromosome 22, the optimal situation would be split dataset per each chromosome and perform a parallelized analysis. 



## Analysis part 1: data set and cleaning the whole matrix of the CpGs


**NOTE** that if you start from three raw matrices (coverage matrix, methylated and unmethylated counts matrices) -not arranged as GRanged object- you can put them in a list of three dataframe as shown: 

```{r, warning=FALSE}
head(input.list[[1]]) # coverage matrix
head(input.list[[2]]) # methylated counts matrix
head(input.list[[3]]) # unmethylated counts matrix
```
And then run directly the passage of the cleaning step 1 (Pipeline Part 1).


**Cleaning step 1: coverage**

In this step, we clean the input data by applying coverage filters. Specifically, we remove CpGs and individuals with excessive missing values and exclude CpGs with very low coverage (fewer than 10 reads by default).

Coverage filtering ensures that downstream analyses are based on reliable and comparable data. CpGs with low coverage or individuals with too many missing values can introduce noise and bias, potentially leading to spurious results. By applying these thresholds, we retain high-quality CpG sites and individuals, improving the robustness of subsequent statistical analyses.

```{r}
clean.coverage <- cleanCovMat(input.obj = bs_list, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)
```

**Cleaning step 2: outliers**

```{r}
clean.out<- cleanOutliers(filtered.obj=clean.coverage, outlier_threshold=0, remove_outliers = T)
```

## Analysis part 2: methylation values

Here are calculated the values that are used for testing methylation data, i.e., the beta values and the M values from the matrices of counts data.

```{r}
meth <- methValues(clean.out)
```

## Analysis part 3: imputation

### Calculating statistics per CpG and imputing data

```{r}
stats <- statsCpG_parallel(meth)
```

### Imputing data

Our method imputes missing values in methylation datasets by computing the correlation structure between CpG sites. It first computes a non-negative definite covariance matrix via eigen-decomposition, then simulates new values from the same distribution to fill in missing entries. The approach is parallelized for efficiency and preserves the natural correlation patterns in the data, improving accuracy over standard imputation methods.

The method that we introduce, firstly fills missing values with row means. Secondly, it calculates a correlation matrix, it constructs a covariance matrix, ensuring it is positive semi-definite via eigen decomposition and finally simulates new values from a multivariate normal distribution to impute missing entries.

```{r}
imputed <- imputeCorr_parallel(stats)

combine_imputed <- function(chr_list) {
  imputed_list <- lapply(chr_list, function(x) x[["imputed"]]) 
  do.call(rbind, imputed_list)                                  
}

# Example use:
combined_imputed_df <- combine_imputed(imputed)
```



## Contributing  

Contributions are welcome!  
You can help by reporting bugs, suggesting improvements, or contributing code/documentation through pull requests.  


