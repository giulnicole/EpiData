#' @title BSImputeScore
#'
#' @description
#' Compute I-Scores (Näf et al.) to evaluate imputation quality
#' for BS-seq M-values under real missingness.
#'
#' @param bs BSseq object after imputation.
#' @param assay_original Name of assay with original M-values.
#' @param assay_imputed Name of assay with imputed M-values.
#' @param conditioning One of "auto", "cluster", or "sample".
#' @param clusters Optional integer vector of CpG cluster IDs.
#' @param min_obs Minimum observed values per conditioning unit.
#'
#' @return A list with overall I-score and per-unit I-scores.
#'
#'
#'
#' @export
BSImputeScore <- function(bs,
                           assay_original = "M",
                           assay_imputed = "M_values_imputed",
                           conditioning = c("auto", "cluster", "sample"),
                           clusters = NULL,
                           min_obs = 5) {

  conditioning <- match.arg(conditioning)

  if (!inherits(bs, "BSseq"))
    stop("bs must be a BSseq object")

  ## Extract assays
  M_orig <- SummarizedExperiment::assay(bs, assay_original)
  M_imp  <- SummarizedExperiment::assay(bs, assay_imputed)

  if (!identical(dim(M_orig), dim(M_imp)))
    stop("Original and imputed assays must have identical dimensions")

  # Save original missingness mask
  missing_mask <- is.na(M_orig)

  ## Determine conditioning automatically
  if (conditioning == "auto") {
    if ("cluster" %in% colnames(SummarizedExperiment::rowData(bs))) {
      conditioning <- "cluster"
    } else {
      conditioning <- "sample"
    }
  }

  ## Build conditioning groups
  if (conditioning == "cluster") {

    if (is.null(clusters)) {
      if (!"cluster" %in% colnames(SummarizedExperiment::rowData(bs))) {
        stop("Cluster-based conditioning requested but no clusters found.")
      }
      clusters <- SummarizedExperiment::rowData(bs)$cluster
    }

    groups <- as.integer(clusters)
    group_ids <- sort(unique(groups))

  } else {  # sample-wise conditioning

    # Each CpG is its own group for sample-wise evaluation
    groups <- seq_len(nrow(M_orig))
    group_ids <- groups
  }

  ## Helper: predictive distribution estimator
  estimate_predictive <- function(x) {
    mu <- mean(x, na.rm = TRUE)
    sdv <- stats::sd(x, na.rm = TRUE)
    if (is.na(sdv) || sdv == 0) sdv <- 1e-6
    list(mu = mu, sd = sdv)
  }

  ## Storage
  unit_scores <- numeric()
  total_scores <- numeric()

  ## I-Score computation
  for (g in group_ids) {

    idx <- which(groups == g)
    g_scores <- numeric()

    for (j in seq_len(ncol(M_orig))) {

      # Observed values used to estimate predictive distribution
      obs <- M_orig[idx, j]
      imp <- M_imp[idx, j]

      # Use original missingness mask
      miss_idx <- which(missing_mask[idx, j])
      obs_idx  <- which(!missing_mask[idx, j])

      # Skip if not enough observed values or no missing to evaluate
      if (length(obs_idx) < min_obs || length(miss_idx) == 0)
        next

      # Estimate predictive distribution
      pred <- estimate_predictive(obs[obs_idx])

      # Log-likelihood for originally missing positions
      ll <- stats::dnorm(
        imp[miss_idx],
        mean = pred$mu,
        sd = pred$sd,
        log = TRUE
      )

      g_scores <- c(g_scores, ll)
    }

    if (length(g_scores) > 0) {
      unit_scores[as.character(g)] <- mean(g_scores)
      total_scores <- c(total_scores, g_scores)
    }
  }

  ## Handle non-identifiable I-Score
  if (length(total_scores) == 0) {
    warning(
      "I-Score could not be computed: no conditioning unit had both ",
      "sufficient observed and missing values (min_obs = ", min_obs, ")."
    )

    return(list(
      I_score = NA_real_,
      unit_I_scores = numeric(),
      conditioning = conditioning,
      n_evaluated = 0
    ))
  }

  ## Final result
  list(
    I_score = mean(total_scores),
    unit_I_scores = unit_scores,
    conditioning = conditioning,
    n_evaluated = length(total_scores)
  )
}
