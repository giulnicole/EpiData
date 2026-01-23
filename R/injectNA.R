#' @title injectNA
#'
#' @importFrom mice ampute
#'
#' @param data A numeric dataframe or matrix
#' @param prop Proportion of missing values to introduce (0 to 1)
#' @param mech Missing data mechanism: "MCAR", "MAR", or "MNAR"
#' @param patterns Optional: specific patterns for ampute (see mice::ampute documentation)
#' @param ... Additional arguments passed to mice::ampute
#'
#' @return Dataframe with injected NA values
#'
#' @examples
#' # Create sample data
#' df <- data.frame(
#'   x1 = rnorm(100),
#'   x2 = rnorm(100),
#'   x3 = rnorm(100)
#' )
#'
#' # MCAR: Missing Completely At Random
#' df_mcar <- injectNA(df, prop = 0.2, mech = "MCAR")
#'
#' # MAR: Missing At Random
#' df_mar <- injectNA(df, prop = 0.3, mech = "MAR")
#'
#' # MNAR: Missing Not At Random
#' df_mnar <- injectNA(df, prop = 0.25, mech = "MNAR")
#'
#' @export
injectNA <- function(data, prop = 0.2, mech = "MCAR", patterns = NULL, ...) {

  # Check if mice package is available
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required. Install it with: install.packages('mice')")
  }

  # Validate inputs
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a dataframe or matrix")
  }

  if (!all(sapply(data, is.numeric))) {
    stop("All columns must be numeric")
  }

  if (prop < 0 || prop > 1) {
    stop("prop must be between 0 and 1")
  }

  mech <- toupper(mech)
  if (!mech %in% c("MCAR", "MAR", "MNAR")) {
    stop("mech must be one of: 'MCAR', 'MAR', or 'MNAR'")
  }

  # Convert to dataframe if matrix
  is_matrix <- is.matrix(data)
  if (is_matrix) {
    col_names <- colnames(data)
    data <- as.data.frame(data)
    if (!is.null(col_names)) {
      colnames(data) <- col_names
    }
  }

  # Set mechanism for ampute
  mech_ampute <- switch(mech,
                        "MCAR" = "MCAR",
                        "MAR" = "MAR",
                        "MNAR" = "MNAR")

  # Use mice::ampute to introduce missing values
  if (is.null(patterns)) {
    result <- mice::ampute(data, prop = prop, mech = mech_ampute, ...)
  } else {
    result <- mice::ampute(data, prop = prop, patterns = patterns, mech = mech_ampute, ...)
  }

  # Extract the amputed data
  data_with_na <- result$amp

  # Convert back to matrix if input was matrix
  if (is_matrix) {
    data_with_na <- as.matrix(data_with_na)
  }

  return(data_with_na)
}


