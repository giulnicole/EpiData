

split5MatXChrom_compact <- function(list.cleaned) {

  # Extract matrices
  cov <- list.cleaned$Coverage_matrix
  met <- list.cleaned$Met_matrix
  unmet <- list.cleaned$Unmet_matrix
  beta <- list.cleaned$Beta_matrix
  M <- list.cleaned$M_matrix

  # Extract chromosome from rownames
  chroms <- sub("-.*", "", rownames(cov))
  chroms_unique <- unique(chroms)

  # Build list per chromosome
  result <- setNames(
    lapply(chroms_unique, function(chr) {
      rows <- chroms == chr
      list(
        Coverage_matrix = cov[rows, , drop = FALSE],
        Met_matrix = met[rows, , drop = FALSE],
        Unmet_matrix = unmet[rows, , drop = FALSE],
        Beta_matrix = beta[rows, , drop = FALSE],
        M_matrix = M[rows, , drop = FALSE]
      )
    }),
    paste0("chr", chroms_unique)
  )

  return(result)
}
