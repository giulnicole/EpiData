#' @title split3MatChrom_compact
#'
#' @description A function for splitting the 5 matrices per each chromosome.
#'
#' @import purrr
#' @name split3MatChrom_compact
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
#'  splitted.list <- split3MatXChrom_compact(clean.out)
#'  }
#'
#' @export

split3MatXChrom_compact <- function(list.cleaned) {

  coverage <- list.cleaned$Coverage_matrix
  methylated <- list.cleaned$Met_matrix
  unmethylated <- list.cleaned$Unmet_matrix

  chroms <- sub("-.*", "", rownames(coverage))

  # Creo la lista compatta
  setNames(
    lapply(unique(chroms), function(chr) {
      rows <- chroms == chr
      list(
        coverage = coverage[rows, , drop = FALSE],
        methylated = methylated[rows, , drop = FALSE],
        unmethylated = unmethylated[rows, , drop = FALSE]
      )
    }),
    paste0("chr", unique(chroms))
  )
}



