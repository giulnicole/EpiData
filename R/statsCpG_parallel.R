#' @title statsCpG_parallel
#' @description
#' \code{} Parallelized version of statsCpG that processes multiple chromosomes. It takes a list where each element contains the cleaned object data.
#'
#' @importFrom dplyr summarise_all
#' @importFrom tibble as_tibble
#' @importFrom parallel makeCluster stopCluster parLapply detectCores clusterExport clusterEvalQ
#' @import stats
#'
#'
#' @param cleaned.obj List of lists. Each element contains the cleaned object with multiple assays.
#' @param varselect Index of the dataset to be used (1-5).
#' @param plot Barplot of missing values. Default = FALSE.
#' @param n.cores Number of cores to use. Default is NULL (uses detectCores() - 1).
#'
#' @name statsCpG_parallel
#'
#' @return List of statistics objects, one per chromosome
#'
#' @examples
#' \dontrun{
#'  data("bs_list")
#'  clean.coverage <- cleanCovMat(input.obj = bs_list, max_na_cpg = 0.5,
#'  max_na_ind = 0.2, cpg_removal_threshold = 10)
#'  clean.out <- cleanOutliers(filtered.obj=clean.coverage,
#'  outlier_threshold=0, remove_outliers = TRUE)
#'  mat.values <- methValues(clean.coverage)
#'  stat.result <- statsCpG_parallel(mat.values)
#' }
#'
#' @export
statsCpG_parallel <- function(cleaned.obj, varselect = 5, plot = FALSE, n.cores = NULL) {

  # Detect vignette or non-interactive mode
  if (!interactive() || !is.null(knitr::opts_knit$get("rmarkdown.pandoc.to"))) {
    message("Non-interactive or vignette build detected, running in serial mode.")
    n.cores <- 1
  } else {
    if (is.null(n.cores)) {
      n.cores <- parallel::detectCores() - 1
    }
    n.cores <- max(1, min(n.cores, length(cleaned.obj)))
  }

  cat(paste("Processing", length(cleaned.obj), "chromosomes using", n.cores, "cores (parLapply)\n"))

  var_select <- varselect
  plot_flag <- plot

  # --- SERIAL EXECUTION (for vignette / non-interactive) ---
  if (n.cores == 1) {
    results <- lapply(cleaned.obj, function(chr_data) {
      options(error = expression(NULL))

      cpg <- rownames(chr_data[[var_select]])
      X <- t(chr_data[[var_select]])
      Y <- X

      rows <- nrow(X)
      cols <- ncol(X)

      colnames(X) <- cpg
      names <- rownames(X)

      vars_non_num <- (X)[!sapply(X, is.numeric)]
      if (length(vars_non_num) != 0) {
        stop(paste("Warning! CpG(s) ", (paste(vars_non_num, collapse = ", ")),
                   " is/are not numeric.", sep = ""))
      }

      no.complete <- data.frame(apply(X, 2, function(x) sum(!is.na(x))))
      perc.complete <- no.complete / rows * 100L

      no.incomplete <- data.frame(apply(X, 2, function(x) sum(is.na(x))))
      perc.incomplete <- no.incomplete / rows * 100L

      table_missing <- data.frame(
        CpG = cpg,
        matrix(c(no.complete, perc.complete, no.incomplete, perc.incomplete),
               ncol = 4L,
               dimnames = list(NULL, c("nObs", "percObs", "nNA", "percNA"))),
        stringsAsFactors = FALSE
      )

      missfrac_per_X <- sum(is.na(X)) / (nrow(X) * ncol(X))
      missfrac_per_var <- colMeans(is.na(X))
      na_per_X <- sum(is.na(X))
      na_per_var <- sapply(X, function(x) sum(length(which(is.na(x)))))

      cormat.lin <- matrix(nrow = cols, ncol = cols)
      colnames(cormat.lin) <- colnames(Y)
      rownames(cormat.lin) <- colnames(Y)

      cormat.lin <- stats::cor(Y, use = "pairwise.complete.obs", method = "pearson")
      colnames(cormat.lin) <- cpg
      rownames(cormat.lin) <- cpg

      Y <- t(Y)
      means <- rowMeans(Y, na.rm = TRUE)
      sds <- apply(Y, 1, sd, na.rm = TRUE)

      long_cormat.lin <- cbind(expand.grid(dimnames(cormat.lin)), value = as.vector(cormat.lin))

      if (plot_flag == TRUE) {
        plots <- plotStats(chr_data)

        res <- list(
          Matrix = as.data.frame(Y),
          Individuals = rows,
          CpGs = cols,
          Table_missing = table_missing,
          Fraction_missingnessp = missfrac_per_X,
          Linear_correlation = cormat.lin,
          long_correlation_matrix = long_cormat.lin,
          means = means,
          sds = sds,
          Barplot = plots
        )
      } else {
        res <- list(
          Matrix = as.data.frame(Y),
          Individuals = rows,
          CpGs = cols,
          Table_missing = table_missing,
          Fraction_missingnessp = missfrac_per_X,
          Linear_correlation = cormat.lin,
          long_correlation_matrix = long_cormat.lin,
          means = means,
          sds = sds
        )
      }

      return(res)
    })

    return(results)
  }

  # --- PARALLEL EXECUTION (interactive mode only) ---
  cl <- parallel::makeCluster(n.cores)

  tryCatch({
    parallel::clusterEvalQ(cl, {
      library(stats)
      library(dplyr)
      library(tibble)
    })

    parallel::clusterExport(cl, varlist = c("var_select", "plot_flag"), envir = environment())

    cat("Starting parallel processing...\n")

    results <- parallel::parLapply(cl, cleaned.obj, function(chr_data) {

      options(error = expression(NULL))

      cpg <- rownames(chr_data[[var_select]])

      X <- t(chr_data[[var_select]])
      Y <- X

      rows <- nrow(X)
      cols <- ncol(X)

      colnames(X) <- cpg
      names <- rownames(X)

      vars_non_num <- (X)[!sapply(X, is.numeric)]
      if (length(vars_non_num) != 0) {
        stop(paste("Warning! CpG(s) ", (paste(vars_non_num, collapse = ", ")),
                   " is/are not numeric.", sep = ""))
      }

      no.complete <- data.frame(apply(X, 2, function(x) sum(!is.na(x))))
      perc.complete <- no.complete / rows * 100L

      no.incomplete <- data.frame(apply(X, 2, function(x) sum(is.na(x))))
      perc.incomplete <- no.incomplete / rows * 100L

      table_missing <- data.frame(
        CpG = cpg,
        matrix(c(no.complete, perc.complete, no.incomplete, perc.incomplete),
               ncol = 4L,
               dimnames = list(NULL, c("nObs", "percObs", "nNA", "percNA"))),
        stringsAsFactors = FALSE
      )

      missfrac_per_X <- sum(is.na(X)) / (nrow(X) * ncol(X))
      missfrac_per_var <- colMeans(is.na(X))
      na_per_X <- sum(is.na(X))
      na_per_var <- sapply(X, function(x) sum(length(which(is.na(x)))))

      cormat.lin <- matrix(nrow = cols, ncol = cols)
      colnames(cormat.lin) <- colnames(Y)
      rownames(cormat.lin) <- colnames(Y)

      cormat.lin <- stats::cor(Y, use = "pairwise.complete.obs", method = "pearson")
      colnames(cormat.lin) <- cpg
      rownames(cormat.lin) <- cpg

      Y <- t(Y)
      means <- rowMeans(Y, na.rm = TRUE)
      sds <- apply(Y, 1, sd, na.rm = TRUE)

      long_cormat.lin <- cbind(expand.grid(dimnames(cormat.lin)), value = as.vector(cormat.lin))

      if (plot_flag == TRUE) {
        plots <- plotStats(chr_data)

        res <- list(
          Matrix = as.data.frame(Y),
          Individuals = rows,
          CpGs = cols,
          Table_missing = table_missing,
          Fraction_missingnessp = missfrac_per_X,
          Linear_correlation = cormat.lin,
          long_correlation_matrix = long_cormat.lin,
          means = means,
          sds = sds,
          Barplot = plots
        )
      } else {
        res <- list(
          Matrix = as.data.frame(Y),
          Individuals = rows,
          CpGs = cols,
          Table_missing = table_missing,
          Fraction_missingnessp = missfrac_per_X,
          Linear_correlation = cormat.lin,
          long_correlation_matrix = long_cormat.lin,
          means = means,
          sds = sds
        )
      }

      return(res)
    })

    cat("Parallel processing completed!\n")

  }, finally = {
    parallel::stopCluster(cl)
  })

  return(results)
}
