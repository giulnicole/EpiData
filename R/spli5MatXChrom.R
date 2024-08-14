
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
#'
#'
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

  cpg <- rownames(met.cleaned)
  chr <- gsub("-.*", "", cpg)
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
  split_met <- met.cleaned %>%
    group_split(chr)%>%
    purrr::map(~ as.data.frame(.))

  #
  split_met <- split_met %>%
    map(~ .[,-(n+1)])

  split_met<- as.data.frame(split_met)
  rownames(split_met)<- cpg

  # Unmet
  split_unmet <- unmet.cleaned %>%
    group_split(chr)%>%
    purrr::map(~ as.data.frame(.))

  #
  split_unmet <- split_unmet %>%
    map(~ .[,-(n+1)])

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


  # Beta
   split_beta <- beta.cleaned %>%
     group_split(chr)%>%
     purrr::map(~ as.data.frame(.))

  #
   split_beta <- split_beta %>%
     map(~ .[,-(n+1)])

   split_beta<- as.data.frame(split_beta)
     rownames(split_beta)<- cpg


  # M
   split_m <- m.cleaned %>%
    group_split(chr)%>%
    purrr::map(~ as.data.frame(.))

  #
  split_m <- split_m %>%
    map(~ .[,-(n+1)])

  split_m<- as.data.frame(split_m)
   rownames(split_m)<- cpg

  # unifying the 5 lists per each chromosome

  final <- list()

  matrices <- list(Coverage = split_cov, Methylated = split_met, Unmethylated = split_unmet,
                    Beta = split_beta, M = split_m)


  if (length(label.chr)>1){
    for (i in 1:length(label.chr)){

      a <- as.numeric(i)

      aa <- list(Coverage_matrix = as.data.frame(matrices[[a]][["Coverage"]]),
                 Met_matrix =  as.data.frame(matrices[[a]][["Methylated"]]),
                 Unmet_matrix =  as.data.frame(matrices[[a]][["Unmethylated"]]),
                 Beta_matrix = as.data.frame(matrices[[a]][["Beta"]]),
                 M_matrix = as.data.frame(matrices[[a]][["M"]]))

      final[[a]]<- aa

    }

    res <- list(final = final,
                coverage = split_cov,
                methylated = split_met,
                unmethylated = split_unmet,
                beta = split_beta, m = split_m)

  } else {

    res <- list(Coverage_matrix = split_cov, Met_matrix = split_met,
                Unmet_matrix = split_unmet,
                Beta = split_beta, M = split_m)

  }


  return(res)

}
