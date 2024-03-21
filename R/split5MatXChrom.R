#' @title split5MatChrom
#'
#' @description A function for splitting the 5 matrices per each chromosome.
#'
#' @import purrr
#' @name split5MatChrom
#'
#' @param list.cleaned List cleaned by filtering steps, with all the 3 matrices with counts (coverage, methylated, unmethylated counts) and matrices of beta and M values.
#'
#' @return
#' List with all the matrices divided per each chromosome.
#'
#'
#' @examples
#' \dontrun{
#'  data("matrices")
#'  clean.coverage2 <- cleanCovMat(input.obj=dati, max_na_cpg = 0.5, max_na_ind = 0.2,  cpg_removal_threshold = 10)
#'  clean.out <- cleanOutliers(filtered.obj=clean.coverage2, outlier_threshold=5, remove_outliers = TRUE)
#'  splitted.list <- split5MatXChrom(clean.out)
#'  }
#'
#' @export
split5MatXChrom <- function(list.cleaned)   {

  # split the chromosome

  cov.cleaned <- list.cleaned[[1]]
  met.cleaned <- list.cleaned[[2]]
  unmet.cleaned <- list.cleaned[[3]]
  beta.cleaned <- list.cleaned[[4]]
  m.cleaned <- list.cleaned[[5]]

  chr <- rownames(met.cleaned)
  chr <- gsub("-.*", "", chr)
  n <- ncol(met.cleaned)

  met.cleaned <- as.data.frame(cbind(met.cleaned, chr))
  colnames(met.cleaned)[n+1] <- "chr"

  unmet.cleaned <- as.data.frame(cbind(unmet.cleaned, chr))
  colnames(unmet.cleaned)[n+1] <- "chr"

  cov.cleaned <- as.data.frame(cbind(cov.cleaned, chr))
  colnames(cov.cleaned)[n+1] <- "chr"

  beta.cleaned <- as.data.frame(cbind(beta.cleaned, chr))
  colnames(beta.cleaned)[n+1] <- "chr"

  m.cleaned <- as.data.frame(cbind(m.cleaned, chr))
  colnames(m.cleaned)[n+1] <- "chr"

  label.chr <- unique(chr)


  split_met <- list()
  split_unmet <- list()
  split_cov <- list()
  split_beta <- list()
  split_m <- list()



  # Met
  split_met <- label.chr %>%
    map(~ met.cleaned %>%
          filter(chr == .x) %>%
          select(-n) %>%  # Remove the last column
          as.matrix() %>%
          as.numeric() %>%
          matrix(ncol = n))

  #
  split_met <- split_met %>%
    map(~ .[,-(n+1)])

  # Unmet
  split_unmet <- label.chr %>%
    map(~ unmet.cleaned %>%
          filter(chr == .x) %>%
          select(-n) %>%  # Remove the last column
          as.matrix() %>%
          as.numeric() %>%
          matrix(ncol = n))  # Assuming n is the desired number of columns


  split_unmet <- split_unmet %>%
    map(~ .[,-(n+1)])


  # Coverage
  split_cov <- label.chr %>%
    map(~ cov.cleaned %>%
          filter(chr == .x) %>%
          select(-n) %>%  # Remove the last column
          as.matrix() %>%
          as.numeric() %>%
          matrix(ncol = n))  # Assuming n is the desired number of columns

  #
  split_cov <- split_cov %>%
    map(~ .[,-(n+1)])



  # Beta
  split_beta <- label.chr %>%
    map(~ beta.cleaned %>%
          filter(chr == .x) %>%
          select(-n) %>%  # Remove the last column
          as.matrix() %>%
          as.numeric() %>%
          matrix(ncol = n))  # Assuming n is the desired number of columns

  #
  split_beta <- split_beta %>%
    map(~ .[,-(n+1)])



  # M
  split_m <- label.chr %>%
    map(~ m.cleaned %>%
          filter(chr == .x) %>%
          select(-n) %>%  # Remove the last column
          as.matrix() %>%
          as.numeric() %>%
          matrix(ncol = n))  # Assuming n is the desired number of columns

  #
  split_m <- split_m %>%
    map(~ .[,-(n+1)])


 # unifying the 5 lists per each chromosome

  final <- list()

  matrices <- list(Coverage = split_cov, Methylated = split_met, Unmethylated = split_unmet,
                   Beta = split_beta, M = split_m)


  for (i in 1:length(label.chr)){

    a <- as.numeric(i)
    aa <- list(Coverage_matrix = as.data.frame(matrices[["Coverage"]][[a]]),
               Met_matrix =  as.data.frame(matrices[["Methylated"]][[a]]),
               Unmet_matrix =  as.data.frame(matrices[["Unmethylated"]][[a]]),
               Beta_matrix = as.data.frame(matrices[["Beta"]][[a]]),
               M_matrix = as.data.frame(matrices[["M"]][[a]]))

    final[[i]]<- aa

  }

  res <- list(final = final, coverage = split_cov, methylated = split_met,
              unmethylated = split_unmet, beta = split_beta, m = split_m)

  return(res)

}
