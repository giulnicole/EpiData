#' @title readBS
#' @description
#' Function for importing raw data matrices from text files.
#' Each file should represent a different matrix (coverage, methylated counts, and unmethylated counts).
#'
#' @importFrom utils read.table
#'
#' @param file_paths Character vector with paths to three `.txt` files in the following order:
#' coverage matrix, methylated count matrix, and unmethylated count matrix.
#'
#' @details
#' The function reads each file as a numeric data frame, using the first column as row names
#' and the first row as column names. All columns are converted to numeric values.
#'
#' @return
#' A list with three elements:
#' \describe{
#'   \item{Coverage_matrix}{Numeric matrix of coverage values.}
#'   \item{Met_matrix}{Numeric matrix of methylated counts.}
#'   \item{Unmet_matrix}{Numeric matrix of unmethylated counts.}
#' }
#'
#' @examples
#' coverage_file <- system.file("extdata", "coverage.txt", package = "EpiData")
#' met_file <- system.file("extdata", "methylated_counts.txt", package = "EpiData")
#' unmet_file <- system.file("extdata", "unmethylated_counts.txt", package = "EpiData")
#' file_paths <- c(coverage_file, met_file, unmet_file)
#' bs_list <- readBS(file_paths)
#'
#' @export
readBS <- function(file_paths) {
  read_single_file <- function(file) {
    # Read table with first column as row names
    df <- read.table(file, header = TRUE, row.names = 1, sep = "", check.names = FALSE)
    # Convert all columns to numeric
    df[] <- lapply(df, function(x) as.numeric(as.character(x)))

    return(df)
  }

  # Apply to all files
  df_list <- lapply(file_paths, read_single_file)

  # Assign fixed names to the three elements
  if (length(df_list) != 3) stop("You must provide exactly three files.")
  names(df_list) <- c("Coverage_matrix", "Met_matrix", "Unmet_matrix")

  return(df_list)
}

