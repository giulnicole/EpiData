#' @title statsCpG
#'
#' @description A function for calculating statistics of missing values of one of the dataset collected in the SummarizedExperiment object.
#'
#' \code{\link{statsCpG}} computes the statistics on missing values per dataset, highlighting the
#' pattern of missing values, mean and standard deviation per each CpG.
#'
#' @importFrom dplyr summarise_all
#' @importFrom tibble as_tibble
#' @import stats
#'
#' @param cleaned.obj SummarizedExperiment list object from filtering coverage (cleanCovMat) + cleaning outliers (cleanOutliers).
#' @param varselect Index of the dataset to be used (numeric value 1-5) to extract the information related to the CpGs' pattern. 1 = coverage counts, 2 = methylated counts, 3 = unmethylated counts, 4 = beta values, 5 = M values.
#' @param plot Barplot of missing values. Default = FALSE.
#'
#'
#' @name statsCpG
#'
#' @return
#'  Object and with statistics on missing data:
#'  \item{Matrix}{Matrix on which missing data pattern is explored.}
#'  \item{Individuals}{Subject passing filters.}
#'  \item{CpGs}{Filtered CpGs.}
#'  \item{Table_missing}{Table with NAs statistics per CpG with 5 columns: CpG's id, nObs, percObs, nNA, percNA.}
#'  \item{Fraction_missingness}{Fraction of NAs in total.}
#'  \item{K_table}{Table K table extra for table_missing.}
#'  \item{Linear_correlation}{Matrix with liean correlation.}
#'  \item{long_correlation_matrix}{Matrix  with linear correlation in longitudinal format.}
#'  \item{means}{Mean per CpG.}
#'  \item{sds}{Standard deviation per CpG.}
#'  \item{Barplot}{Barplot with proportion of missing data per matrix.}
#'
#' @examples
#' \dontrun{
#' # Data
#'  data("meth_data")
#'
#'  input.list<- assays(meth_data)
#'
#'  # Cleaning low counts for coverage
#'  clean.coverage <- cleanCovMat(input.obj = input.list,
#'                               max_na_cpg = 0.5, max_na_ind = 0.2,
#'                              cpg_removal_threshold = 10)
#'
#' # Cleaning outliers in methylated and unmethylated counts
#'  clean.out <- cleanOutliers(filtered.obj=clean.coverage,
#'                             outlier_threshold=5,
#'                             remove_outliers = TRUE)
#'
#'  clean.values <- methValues(clean.out)
#'
#'  stats <- statsCpG(cleaned.obj=clean.values, varselect = 5)
#' }
#'
#'
#' @export
#'
#'
#'
#'
statsCpG <- function(cleaned.obj, varselect = 5, plot=FALSE) {

  options(error = expression(NULL))

  cpg <- rownames(cleaned.obj[[varselect]])


  # Transposing selected variable

  X<-t(cleaned.obj[[varselect]])
  Y <- X

  # Transposing selected variable


  # Extraction of dimensions of selected variable matrix
  rows <- nrow(X)  # individuals
  cols <- ncol(X)  # CpG


  colnames(X) <- cpg
  names <- rownames(X)


  # Applying as.numeric function to all X elements
  vars_non_num <- (X)[!sapply(X, is.numeric)]
  if (length(vars_non_num) != 0){
    stop(paste("Warning! CpG(s) ", (paste(vars_non_num, collapse = ", ")),
               " is/are not numeric. Convert this/these variables to numeric using clean_cpg() and repeat the get_data() function until no warnings are shown.",
               sep = ""))

  }
  ### MODULE 1: table of statistics on observed values and NAs and statistics on missing values

  # Statistics on selected variable values
  # Extraction of statistics (complete observations and NAs)

  # Number and percentage of complete individuals
  no.complete <- data.frame(apply(X, 2, function(x) sum(!is.na(x))))

  perc.complete <- no.complete / rows * 100L

  # Number and percentage of incomplete individuals
  cat("Calculating missing values per individual ...\n")
  no.incomplete <- data.frame(apply(X, 2, function(x) sum(!is.na(x))))

  perc.incomplete <- no.incomplete / rows * 100L


  ###

  # Frequency table
  table_missing <- data.frame(CpG = cpg,
                              matrix(c(no.complete, perc.complete, no.incomplete,
                                       perc.incomplete), ncol = 4L,
                                      dimnames = list(NULL, c("nObs", "percObs", "nNA", "percNA"))),
                                      stringsAsFactors = FALSE)

# Kable extra object
# table_k <- kable(table_missing, row.names = F)  %>%
#   kable_styling(bootstrap_options = c("striped", "hover", "bordered"),
#                  full_width = T)

  ###

  # Missingness per variable
  cat("Calculating missing values per variable ...\n")
  missfrac_per_X <- sum(is.na(X))/(nrow(X) * ncol(X))
  missfrac_per_var <- colMeans(is.na(X))
  na_per_X <- sum(is.na(X))
  na_per_var <- sapply(X, function(x) sum(length(which(is.na(x)))))

  # Missing data pattern
  #cat("Computing missing data pattern per variable ...\n")
  # mdpat <- mice::md.pattern(X, plot = FALSE)     # this funciton returns
  # a matrix with \code{ncol(x)+1} columns, in which each row corresponds
  # to a missing data pattern (1=observed, 0=missing).

  #mdpat <- mdpat[, colnames(mdpat) %in% cpg]

  # Removing counts for leaving only the missing pattern coded as 1,0
  # mdpat_count <- mdpat[-c(1, ncol(mdpat)), ] # removing
  #rownames(mdpat_count) <- names


  ### MODULE 2: correlation between NAs and observed selected variable

  # Complete cases to calculate correlation
  #comp <- sum(stats::complete.cases(Y))
  #mdpat_count <- as.matrix(mdpat_count)

  # Linear Pearson correlation calculation (between the obeserved values)
  cat("Computing the linear correlation matrix ...\n")
  cormat.lin<- matrix(nrow = cols, ncol = cols)   # I put rows because from function stat_cpg cpg are the rows
  colnames(cormat.lin) <- colnames(Y)
  rownames(cormat.lin) <- colnames(Y)

  cormat.lin <- stats::cor(Y, use = "pairwise.complete.obs", method = "pearson")
  colnames(cormat.lin) <- cpg
  rownames(cormat.lin) <- cpg


  cat("Computing means and standard deviations per CpG ...\n")
  Y<- t(Y)
  means <- rowMeans(Y, na.rm=TRUE)
  sds <- apply(Y, 1, sd, na.rm=T)

  long_cormat.lin <- cbind(expand.grid(dimnames(cormat.lin)), value = as.vector(cormat.lin))


   if (plot == TRUE){

  ### MODULE 3: plots
  plots<- plotStats(cleaned.obj)

  cat("Converting results ...\n")
  res<- list(Matrix=as.data.frame(Y), Individuals = rows, CpGs = cols, Table_missing = table_missing,
             Fraction_missingnessp =  missfrac_per_X,  #K_table = table_k,
             Linear_correlation = cormat.lin,
             long_correlation_matrix = long_cormat.lin,
             means=means, sds= sds, Barplot = plots)

     }  else{


  # Results
  # Results
  cat("Converting results ...\n")
  res<- list(Matrix=as.data.frame(Y), Individuals = rows, CpGs = cols, Table_missing = table_missing,
             Fraction_missingnessp =  missfrac_per_X,  #K_table = table_k,
             Linear_correlation = cormat.lin,
             long_correlation_matrix = long_cormat.lin,
             means=means, sds= sds)


     }

  #cleaned.obj@metadata$statistics <- res

  return(res)

}
