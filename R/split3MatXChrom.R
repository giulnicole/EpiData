
#' @title split3MatXChrom
#' @description
#' \code{}A function for splitting the cleaned list per each chromosome.
#'
#'
#' @param list.input list with the three counts' matrices filtered from the `cleanCovMat` and `cleanOutliers` functions.
#'
#'
#' @name split3MatXChrom
#'
#' @return
#' List with 4 elements:
#'  \item{final}{coverage, methylated and unmethylated counts splitted per each chromosome.}
#'  \item{coverage}{coverage with all chromosomes.}
#'  \item{methylated}{methylated counts matrix with all chromosomes.}
#'  \item{unmethylated}{unmethylated counts matrix with all chromosomes.}
#'
#'
#' @examples
#' \dontrun{
#'  data("matrices")
#'  splitted<- split3MatXChrom(assays(dati))
#' }
#'
#' @export
#'
split3MatXChrom <- function(list.input)   {

  # split the chromosome

  cov.cleaned <- list.input[[1]]
  met.cleaned <- list.input[[2]]
  unmet.cleaned <- list.input[[3]]

  chr <- rownames(met.cleaned)
  chr <- gsub("-.*", "", chr)
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
  for (i in label.chr) {

    split_data <- met.cleaned[met.cleaned$chr== i, ]
    #split_data <- data.frame(sapply(split_data, as.numeric))
    is_matrix <- is.matrix(split_data)

    if (!is_matrix) {
      split_data <- as.matrix(split_data)

    }

    split_data2 <- apply(split_data, c(1, 2), as.numeric)
    split_met[[i]] <- split_data2
    #split_met2 <- lapply(split_met, as.matrix.data.frame)

    names(split_met)
  }


  # remove chr col
  for (i in label.chr) {

    split_met[[i]] <- split_met[[i]][,-(n+1)]


  }


  # Unmet
  for (i in label.chr) {

    split_data <- unmet.cleaned[unmet.cleaned$chr== i, ]
    #split_data <- data.frame(sapply(split_data, as.numeric))
    is_matrix <- is.matrix(split_data)

    if (!is_matrix) {
      split_data <- as.matrix(split_data)

    }

    split_data2 <- apply(split_data, c(1, 2), as.numeric)
    split_unmet[[i]] <- split_data2
    #split_met2 <- lapply(split_met, as.matrix.data.frame)

    names(split_unmet)
  }

  # remove chr col
  for (i in label.chr) {

    split_unmet[[i]] <- split_unmet[[i]][,-(n+1)]


  }



  # Coverage
  for (i in label.chr) {

    split_data <- cov.cleaned[cov.cleaned$chr== i, ]
    #split_data <- data.frame(sapply(split_data, as.numeric))
    is_matrix <- is.matrix(split_data)

    if (!is_matrix) {
      split_data <- as.matrix(split_data)

    }

    split_data2 <- apply(split_data, c(1, 2), as.numeric)
    split_cov[[i]] <- split_data2
    #split_met2 <- lapply(split_met, as.matrix.data.frame)

    names(split_cov)
  }

  # remove chr col
  for (i in label.chr) {

    split_cov[[i]] <- split_cov[[i]][,-(n+1)]


  }


  # unifying the 3 lists per each chromosome

  final <- list()

  matrices <- list(Coverage = split_cov, Methylated = split_met, Unmethylated = split_unmet)


  for (i in label.chr){

    a <- as.character(i)
    aa <- list(Coverage_matrix = as.data.frame(matrices[["Coverage"]][[a]]), Met_matrix =  as.data.frame(matrices[["Methylated"]][[a]]),
               Unmet_matrix =  as.data.frame(matrices[["Unmethylated"]][[a]]))

    final[[i]]<- aa

  }


  res <- list(final = final, coverage = split_cov,
              methylated = split_met,
              unmethylated = split_unmet)

  return(res)

}
