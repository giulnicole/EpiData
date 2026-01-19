#' @title evaluateImputation
#'
#' @param original A numeric dataframe or matrix (original complete data)
#' @param imputed A numeric dataframe or matrix (imputed data, same dimensions as original)
#'
#' @return A list containing:
#'   - metrics: Dataframe with RMSE and MAE for each variable
#'   - overall: Overall RMSE and MAE across all variables
#' @export
#'
#' @examples
#' # Original complete data
#' original <- data.frame(
#'   x1 = rnorm(100),
#'   x2 = rnorm(100),
#'   x3 = rnorm(100)
#' )
#'
#' # Data after imputation
#' imputed <- original
#' imputed[sample(1:100, 20), 1] <- imputed[sample(1:100, 20), 1] + rnorm(20, 0, 0.5)
#'
#' # Evaluate
#' results <- evaluateImputation(original, imputed)
#' print(results$metrics)
#' print(results$overall)
#' @export
evaluateImputation <- function(original, imputed) {

  # Validate inputs
  if (!is.data.frame(original) && !is.matrix(original)) {
    stop("original must be a dataframe or matrix")
  }

  if (!is.data.frame(imputed) && !is.matrix(imputed)) {
    stop("imputed must be a dataframe or matrix")
  }

  if (!all(dim(original) == dim(imputed))) {
    stop("original and imputed must have the same dimensions")
  }

  # Convert to dataframe if needed
  if (is.matrix(original)) {
    original <- as.data.frame(original)
  }
  if (is.matrix(imputed)) {
    imputed <- as.data.frame(imputed)
  }

  # Initialize results dataframe
  var_metrics <- data.frame(
    variable = character(),
    RMSE = numeric(),
    MAE = numeric(),
    stringsAsFactors = FALSE
  )

  # Calculate metrics for each variable
  for (col in 1:ncol(original)) {

    # Get values
    true_vals <- as.numeric(original[, col])
    imputed_vals <- as.numeric(imputed[, col])

    # Calculate errors
    errors <- true_vals - imputed_vals

    # RMSE: Root Mean Squared Error
    rmse <- sqrt(mean(errors^2))

    # MAE: Mean Absolute Error
    mae <- mean(abs(errors))

    # Add to results
    var_metrics <- rbind(var_metrics, data.frame(
      variable = colnames(original)[col],
      RMSE = rmse,
      MAE = mae,
      stringsAsFactors = FALSE
    ))
  }

  # Calculate overall metrics (mean across all variables)
  overall_metrics <- data.frame(
    RMSE = mean(var_metrics$RMSE),
    MAE = mean(var_metrics$MAE)
  )

  # Return results
  return(list(
    metrics = var_metrics,
    overall = overall_metrics
  ))
}


# Example usage
if (interactive()) {
  # Create original complete data
  set.seed(123)
  original <- data.frame(
    age = rnorm(100, mean = 50, sd = 10),
    income = rnorm(100, mean = 50000, sd = 15000),
    score = rnorm(100, mean = 75, sd = 10)
  )

  # Simulate imputed data (slightly different from original)
  imputed <- original
  imputed$age <- imputed$age + rnorm(100, 0, 2)
  imputed$income <- imputed$income + rnorm(100, 0, 1000)
  imputed$score <- imputed$score + rnorm(100, 0, 3)

  # Evaluate imputation quality
  results <- evaluateImputation(original, imputed)

  cat("Per-variable metrics:\n")
  print(results$metrics)

  cat("\nOverall metrics:\n")
  print(results$overall)
}
