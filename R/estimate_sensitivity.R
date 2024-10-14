# R cmd check throws a warning if S3generic and S3methods do not accept identical arguments.
# This may cause problems when submitting the package to CRAN later.

estimate_sensitivity <- function(x){
  UseMethod("estimate_sensitivity")
}

#' @export
estimate_sensitivity.data.frame <- function(x){
  counts_table <- get_counts_table(x)
  sensitivity  <- estimate_sensitivity(counts_table)
  sensitivity
}

#' @export
estimate_sensitivity.matrix <- function(x){
  # Normalize to avoid infinite values
  counts_table <- x

  # Count number of predictions and correct predictions
  n_responses <- ncol(counts_table)
  idx_predicting_first_row  <- 1:(n_responses/2)
  idx_predicting_second_row <- (n_responses/2 + 1):n_responses

  n_correctly_predict_first_row <- sum(counts_table[1, idx_predicting_first_row])
  n_predict_first_row           <- sum(counts_table[ , idx_predicting_first_row])

  n_correctly_predict_second_row <- sum(counts_table[2, idx_predicting_second_row])
  n_predict_second_row           <- sum(counts_table[ , idx_predicting_second_row])

  hit_rate         <- n_correctly_predict_first_row / n_predict_first_row
  false_alarm_rate <- 1 - n_correctly_predict_second_row / n_predict_second_row

  sensitivity <- qnorm(hit_rate) - qnorm(false_alarm_rate)

  sensitivity
}

#' @export
estimate_sensitivity.table <- function(x){
  sensitivity <- estimate_sensitivity.matrix(x)
  sensitivity
}
