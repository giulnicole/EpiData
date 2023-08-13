#' @title Calculating statistics and missing values patterns
#'
#' @description
#' \code{\link{stats_cpg2}} computes the statistics on missing values per dataset, highlighting the
#' pattern of missing values per each CpG
#'
#'
#' @param X Cleaned matrix list from coverage cleaning (clean_matrix)
#' @param varselect is the index of the dataset to be used in the output from clean_matrix (numeric value 1-5)

#' @name stats_cpg2
#'
#' @return
#' #' list with 7 datasets:
#'  \item{Rows}{}
#'  \item{Columns}{}
#'  \item{Table_missing}{}
#'  \item{Total_NA}{}
#'  \item{NA_per_variable}{}
#'  \item{Fraction_missingness}{}
#'  \item{K_table}{}
#'  \item{MD_Pattern}{}
#'  \item{MD_Pattern_count}{}
#'  \item{Linear_correlation}{}
#'  \item{long_correlation_matrix}{}
#'
#'@examples
#'
#'
#'
#'
#'
#' @export
#'
#'
stats_cpg2 <- function(X, varselect = 2) {


  options(error = expression(NULL))

  cpg <- rownames(X[[varselect]])


  # Transposing selected variable

  X[[varselect]] <-t(X[[varselect]])
  Y <- X[[varselect]]

  # Extraction of dimensions of selected variable matrix
  rows <- nrow(X[[varselect]])  # individuals
  cols <- ncol(X[[varselect]])  # CpG


  colnames(X[[varselect]]) <- cpg
  names <- rownames(X[[varselect]])


  # Applying as.numeric function to all X elements
  vars_non_num <- (X[[varselect]])[!sapply(X[[varselect]], is.numeric)]
  if (length(vars_non_num) != 0){
    stop(paste("Warning! CpG(s) ", (paste(vars_non_num, collapse = ", ")),
               " is/are not numeric. Convert this/these variables to numeric using clean_cpg() and repeat the get_data() function until no warnings are shown.",
               sep = ""))

  }
  ### MODULE 1: table of statistics on observed values and NAs and statistics on msising values

  # Statistics on selected variable values
  # Extraction of statistics (complete observations and NAs)

  # Number and percentage of complete individuals
  no.complete <- X[[varselect]] %>% as_tibble(X[[varselect]])%>%
    summarise_all(funs(sum(!is.na(.))))

  no.complete <- as.numeric(no.complete)
  perc.complete <- no.complete / rows * 100L

  # Number and percentage of incomplete individuals
  cat("Calculating missing values per individual ...\n")
  no.incomplete <- X[[varselect]] %>% as_tibble(X[[varselect]])%>%
    summarise_all(funs(sum(is.na(.))))

  no.incomplete <- as.numeric(no.incomplete)
  perc.incomplete <- no.incomplete / rows * 100L


  ###

  # Frequency table
  table_missing <- data.frame(CpG = cpg,
                              matrix(c(no.complete, perc.complete, no.incomplete,
                                       perc.incomplete), ncol = 4L,
                                      dimnames = list(NULL, c("nObs", "percObs", "nNA", "percNA"))),
                                      stringsAsFactors = FALSE)


  # Kable extra object
  table_k <- kable(table_missing, row.names = F)  %>%
    kable_styling(bootstrap_options = c("striped", "hover", "bordered"),
                  full_width = T)


  ###

  # Missingness per variable
  cat("Calculating missing values per variable ...\n")
  missfrac_per_X <- sum(is.na(X[[varselect]]))/(nrow(X[[varselect]]) * ncol(X[[varselect]]))
  missfrac_per_var <- colMeans(is.na(X[[varselect]]))
  na_per_X <- sum(is.na(X[[varselect]]))
  na_per_var <- sapply(X[[varselect]], function(x) sum(length(which(is.na(x)))))

  # Missing data pattern
  cat("Computing missing data pattern per variable ...\n")
  mdpat <- mice::md.pattern(X[[varselect]], plot = FALSE)     # this funciton returns
  # a matrix with \code{ncol(x)+1} columns, in which each row corresponds
  # to a missing data pattern (1=observed, 0=missing).

  mdpat <- mdpat[, colnames(mdpat) %in% cpg]

  # Removing counts for leaving only the missing pattern coded as 1,0
  mdpat_count <- mdpat[-c(1, ncol(mdpat)), ] # removing
  #rownames(mdpat_count) <- names


  ### MODULE 2: correlation between NAs and observed selected variable

  # Complete cases to calculate correlation
  comp <- sum(stats::complete.cases(Y))
  mdpat_count <- as.matrix(mdpat_count)

  # Linear Pearson correlation calculation (between the obeserved values)
  cat("Computing the linear correlation matrix ...\n")
  cormat.lin<- matrix(nrow = cols, ncol = cols)   # I put rows because from function stat_cpg cpg are the rows
  colnames(cormat.lin) <- colnames(Y)
  rownames(cormat.lin) <- colnames(Y)

  cormat.lin <- stats::cor(Y, use = "pairwise.complete.obs", method = "pearson")

  cat("Converting results ...\n")
  long_cormat.lin <- cbind(expand.grid(dimnames(cormat.lin)), value = as.vector(cormat.lin))


  # Results
  res<- list(Rows = rows, Columns = cols, Table_missing = table_missing,
             Total_NA = na_per_X, NA_per_variable = na_per_var,
             Fraction_missingness =  missfrac_per_X,  K_table = table_k,
             MD_Pattern = mdpat, MD_Pattern_count = mdpat_count,
             Linear_correlation = cormat.lin,
             long_correlation_matrix = long_cormat.lin)


  return(res)

}
