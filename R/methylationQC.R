#' @title methylationQC
#'
#' @description This function performs comprehensive quality control on methylation M-values
#' after imputation, including checks for outliers, data range, missing values,
#' and sample/probe quality metrics.
#'
#' @param m_values Data frame or matrix of M-values (rows = CpGs, columns = samples).
#' @param beta_values Optional. Data frame or matrix of beta-values for additional checks.
#' @param positions Optional. Data frame with columns: CpG, chr, pos for genomic context.
#' @param sample_metadata Optional. Data frame with sample information.
#' @param max_na_per_probe Maximum proportion of NA allowed per probe (default 0.2).
#' @param max_na_per_sample Maximum proportion of NA allowed per sample (default 0.1).
#' @param m_range Expected M-value range as c(min, max). Default c(-10, 10).
#' @param detect_outliers Logical. Perform outlier detection (default TRUE).
#' @param outlier_sd Number of SD from median for outlier detection (default 3).
#' @param verbose Logical. Print detailed report (default TRUE).
#'
#' @return A list containing:
#' \itemize{
#'   \item qc_summary: data frame with QC metrics
#'   \item failed_probes: Vector of probe IDs failing QC
#'   \item failed_samples: vector of sample IDs failing QC
#'   \item outlier_info: data frame with outlier statistics
#'   \item cleaned_data: M-values with failed probes/samples removed
#'   \item qc_passed: logical indicating if data passed QC
#' }
#'
#' @export
methylationQC <- function(m_values,
                          beta_values = NULL,
                          positions = NULL,
                          sample_metadata = NULL,
                          max_na_per_probe = 0.2,
                          max_na_per_sample = 0.1,
                          m_range = c(-10, 10),
                          detect_outliers = TRUE,
                          outlier_sd = 3,
                          verbose = TRUE) {

  # Convert to matrix
  m_values <- as.matrix(m_values)
  n_probes <- nrow(m_values)
  n_samples <- ncol(m_values)

  if (verbose) {
    cat("METHYLATION DATA QUALITY CONTROL \n")
    cat(sprintf("Dataset dimensions: %d probes x %d samples\n\n", n_probes, n_samples))
  }

  # Initialize results
  qc_results <- list()
  failed_probes <- character(0)
  failed_samples <- character(0)
  warnings_list <- character(0)

 # Missing values assesstment
  if (verbose) cat("1. Checking missing values...\n")

  na_per_probe <- rowSums(is.na(m_values)) / n_samples
  na_per_sample <- colSums(is.na(m_values)) / n_probes

  probes_high_na <- names(na_per_probe)[na_per_probe > max_na_per_probe]
  samples_high_na <- names(na_per_sample)[na_per_sample > max_na_per_sample]

  if (length(probes_high_na) > 0) {
    failed_probes <- c(failed_probes, probes_high_na)
    warnings_list <- c(warnings_list,
                       sprintf("%d probes exceed %.1f%% missing values",
                               length(probes_high_na), max_na_per_probe * 100))
  }

  if (length(samples_high_na) > 0) {
    failed_samples <- c(failed_samples, samples_high_na)
    warnings_list <- c(warnings_list,
                       sprintf("%d samples exceed %.1f%% missing values",
                               length(samples_high_na), max_na_per_sample * 100))
  }

  if (verbose) {
    cat(sprintf("  - Probes with >%.1f%% NA: %d\n", max_na_per_probe * 100, length(probes_high_na)))
    cat(sprintf("  - Samples with >%.1f%% NA: %d\n", max_na_per_sample * 100, length(samples_high_na)))
    cat(sprintf("  - Total remaining NAs: %d (%.2f%%)\n",
                sum(is.na(m_values)),
                sum(is.na(m_values)) / length(m_values) * 100))
  }


  # M value check
  if (verbose) cat("\n2. Checking M-value ranges...\n")

  m_min <- min(m_values, na.rm = TRUE)
  m_max <- max(m_values, na.rm = TRUE)
  out_of_range <- sum(m_values < m_range[1] | m_values > m_range[2], na.rm = TRUE)

  if (out_of_range > 0) {
    warnings_list <- c(warnings_list,
                       sprintf("%d values outside expected range [%.1f, %.1f]",
                               out_of_range, m_range[1], m_range[2]))
  }

  if (verbose) {
    cat(sprintf("  - M-value range: [%.2f, %.2f]\n", m_min, m_max))
    cat(sprintf("  - Expected range: [%.1f, %.1f]\n", m_range[1], m_range[2]))
    cat(sprintf("  - Values out of range: %d\n", out_of_range))
  }



  if (!is.null(beta_values)) {
    if (verbose) cat("\n3. Validating beta-values...\n")

    beta_values <- as.matrix(beta_values)
    beta_invalid <- sum(beta_values < 0 | beta_values > 1, na.rm = TRUE)

    if (beta_invalid > 0) {
      warnings_list <- c(warnings_list,
                         sprintf("%d beta-values outside [0, 1] range", beta_invalid))
    }

    if (verbose) {
      cat(sprintf("  - Beta range: [%.4f, %.4f]\n",
                  min(beta_values, na.rm = TRUE),
                  max(beta_values, na.rm = TRUE)))
      cat(sprintf("  - Invalid beta-values: %d\n", beta_invalid))
    }
  }


  outlier_info <- NULL
  if (detect_outliers) {
    if (verbose) cat("\n4. Detecting outliers...\n")

    probe_medians <- apply(m_values, 1, median, na.rm = TRUE)
    probe_mads <- apply(m_values, 1, mad, na.rm = TRUE)

    sample_medians <- apply(m_values, 2, median, na.rm = TRUE)
    sample_mads <- apply(m_values, 2, mad, na.rm = TRUE)

    # Sample-level outliers
    median_of_medians <- median(sample_medians, na.rm = TRUE)
    mad_of_medians <- mad(sample_medians, na.rm = TRUE)

    outlier_samples <- abs(sample_medians - median_of_medians) > outlier_sd * mad_of_medians
    outlier_sample_names <- names(sample_medians)[outlier_samples]

    if (length(outlier_sample_names) > 0) {
      failed_samples <- c(failed_samples, outlier_sample_names)
      warnings_list <- c(warnings_list,
                         sprintf("%d samples are statistical outliers (>%d SD)",
                                 length(outlier_sample_names), outlier_sd))
    }

    outlier_info <- data.frame(
      sample = colnames(m_values),
      median_m = sample_medians,
      mad_m = sample_mads,
      is_outlier = outlier_samples,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      cat(sprintf("  - Outlier samples detected: %d\n", sum(outlier_samples)))
      cat(sprintf("  - Median M-value range across samples: [%.2f, %.2f]\n",
                  min(sample_medians), max(sample_medians)))
    }
  }


  if (verbose) cat("\n5. Assessing probe variability...\n")

  probe_vars <- apply(m_values, 1, var, na.rm = TRUE)
  zero_var_probes <- sum(probe_vars == 0, na.rm = TRUE)
  low_var_probes <- sum(probe_vars < 0.01, na.rm = TRUE)

  if (verbose) {
    cat(sprintf("  - Zero variance probes: %d\n", zero_var_probes))
    cat(sprintf("  - Low variance probes (var < 0.01): %d\n", low_var_probes))
    cat(sprintf("  - Median probe variance: %.4f\n", median(probe_vars, na.rm = TRUE)))
  }


  if (!is.null(positions)) {
    if (verbose) cat("\n6. Checking genomic distribution...\n")

    if (all(c("CpG", "chr") %in% colnames(positions))) {
      positions_subset <- positions[positions$CpG %in% rownames(m_values), ]
      chr_counts <- table(positions_subset$chr)

      if (verbose) {
        cat(sprintf("  - Chromosomes covered: %d\n", length(chr_counts)))
        cat(sprintf("  - Probes per chromosome (median): %.0f\n", median(chr_counts)))
      }
    }
  }


  if (verbose) cat("\n7. Creating cleaned dataset...\n")

  failed_probes <- unique(failed_probes)
  failed_samples <- unique(failed_samples)

  keep_probes <- setdiff(rownames(m_values), failed_probes)
  keep_samples <- setdiff(colnames(m_values), failed_samples)

  cleaned_data <- m_values[keep_probes, keep_samples, drop = FALSE]

  if (verbose) {
    cat(sprintf("  - Removed probes: %d\n", length(failed_probes)))
    cat(sprintf("  - Removed samples: %d\n", length(failed_samples)))
    cat(sprintf("  - Cleaned dataset: %d probes x %d samples\n",
                nrow(cleaned_data), ncol(cleaned_data)))
  }


  qc_summary <- data.frame(
    metric = c("Total probes", "Total samples", "Failed probes", "Failed samples",
               "Probes remaining", "Samples remaining", "Missing values (%)",
               "M-value range", "Outlier samples", "Zero variance probes"),
    value = c(n_probes, n_samples, length(failed_probes), length(failed_samples),
              nrow(cleaned_data), ncol(cleaned_data),
              sprintf("%.2f", sum(is.na(m_values)) / length(m_values) * 100),
              sprintf("[%.2f, %.2f]", m_min, m_max),
              ifelse(is.null(outlier_info), 0, sum(outlier_info$is_outlier)),
              zero_var_probes),
    stringsAsFactors = FALSE
  )

  # Determine if QC passed
  qc_passed <- length(failed_probes) < n_probes * 0.5 &&
    length(failed_samples) < n_samples * 0.5 &&
    sum(is.na(cleaned_data)) / length(cleaned_data) < 0.05



  # Final report
  if (verbose) {
    cat("\nQC SUMMARY\n")
    print(qc_summary, row.names = FALSE)

    if (length(warnings_list) > 0) {
      cat("\nWARNINGS\n")
      for (w in warnings_list) cat(paste0(" ", w, "\n"))
    }

    cat(sprintf("\n QC STATUS: %s \n",
                ifelse(qc_passed, "PASSED", "FAILED")))
  }

  # Return comprehensive results
  return(list(
    qc_summary = qc_summary,
    failed_probes = failed_probes,
    failed_samples = failed_samples,
    outlier_info = outlier_info,
    cleaned_data = as.data.frame(cleaned_data),
    qc_passed = qc_passed,
    warnings = warnings_list,
    metrics = list(
      na_per_probe = na_per_probe,
      na_per_sample = na_per_sample,
      probe_variance = probe_vars
    )
  ))
}

