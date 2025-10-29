#' @title split3MatChrom
#'
#' @description A function for splitting the 5 matrices per each chromosome.
#'
#' @import purrr
#' @name split3MatChrom
#'
#' @param list.cleaned List cleaned by filtering steps, with all the 3 matrices with counts (coverage, methylated, unmethylated counts).
#'
#' @return
#' List with all the matrices divided per each chromosome.
#'
#'
#' @examples
#' \dontrun{
#'  data("matrices")
#'  clean.coverage2 <- cleanCovMat(input.obj=dati, max_na_cpg = 0.5,
#'                                  max_na_ind = 0.2,  cpg_removal_threshold = 10)
#'  clean.out <- cleanOutliers(filtered.obj=clean.coverage2, outlier_threshold=5,
#'                              remove_outliers = TRUE)
#'  splitted.list <- split3MatXChrom(clean.out)
#'  }
#'
#' @export
split3MatXChrom <- function(list.cleaned)   {

  # split the chromosome
  cov.cleaned <- list.cleaned[[1]]
  met.cleaned <- list.cleaned[[2]]
  unmet.cleaned <- list.cleaned[[3]]

  cpg <- rownames(met.cleaned)
  chr <- gsub("-.*", "", cpg)
  n <- ncol(met.cleaned)

  met.cleaned <- as.data.frame(cbind(met.cleaned, chr))
  colnames(met.cleaned)[n+1] <- "chr"

  unmet.cleaned <- as.data.frame(cbind(unmet.cleaned, chr))
  colnames(unmet.cleaned)[n+1] <- "chr"

  cov.cleaned <- as.data.frame(cbind(cov.cleaned, chr))
  colnames(cov.cleaned)[n+1] <- "chr"
  label.chr <- unique(chr)


  split_met <- list()
  split_unmet <- list()
  split_cov <- list()


  # Met
  split_met <- met.cleaned %>%
    group_split(chr)%>%
    purrr::map(~ as.data.frame(.))

  #
  split_met <- split_met %>%
    purrr::map(~ .[,-(n+1)])

 split_met<- as.data.frame(split_met)
  rownames(split_met)<- cpg

  # Unmet
  split_unmet <- unmet.cleaned %>%
    group_split(chr)%>%
    purrr::map(~ as.data.frame(.))

  #
  split_unmet <- split_unmet %>%
    purrr::map(~ .[,-(n+1)])

  split_unmet<- as.data.frame(split_unmet)
  rownames(split_unmet)<- cpg

  # Coverage
  split_cov <- cov.cleaned %>%
    group_split(chr)%>%
    purrr::map(~ as.data.frame(.))

  #
  split_cov <- split_cov %>%
    map(~ .[,-(n+1)])

  split_cov<- as.data.frame(split_cov)
  rownames(split_cov)<- cpg

  final <- list()

  matrices <- list(Coverage = split_cov, Methylated = split_met, Unmethylated = split_unmet)


  if (length(label.chr)>1){
  for (i in 1:length(label.chr)){

    a <- as.numeric(i)

    aa <- list(Coverage_matrix = as.data.frame(matrices[[a]][["Coverage"]]),
               Met_matrix =  as.data.frame(matrices[[a]][["Methylated"]]),
               Unmet_matrix =  as.data.frame(matrices[[a]][["Unmethylated"]]) )


    final[[a]]<- aa

    }

    res <- list(final = final,
                coverage = split_cov,
                methylated = split_met,
                unmethylated = split_unmet)


  } else {

    res <- list(Coverage_matrix = split_cov, Met_matrix = split_met,
                Unmet_matrix = split_unmet)


    }


  return(res)

}
