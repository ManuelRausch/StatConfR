#' Estimate Meta Sensitivity
#' 
#' Estimate regular sensitivity from SDT. Here, we assume that the first half
#' of the confidence-binned responses equals to predictions for the first
#' label (Y = +1) and the second half of the responses to predictions for the
#' second label (Y = -1).
#' 
#' Data can be input in three ways:

#' - A data frame with variables "y" for true labels and "r" for
#'   confidence-binned responses. "y" needs to contain values -1 and +1 while
#'   r needs to be a factor with ordered levels such that the first half of
#'   the levels correspond to predictions for y=-1 and the second half to
#'   predictions for y=+1.
#' - A counts table with joint absolute frequencies. Rows correspond to true
#'   labels (stimulus categories) and columns correspond to responses.
#' - A contingency matrix with joint relative frequencies (as before but
#'   normalized to sum up to 1).
#' 
#' @param x Data
#' 
#' @return Sensitivity value
estimate_sensitivity <- function(x)
{
  UseMethod("estimate_sensitivity")
}

estimate_sensitivity.data.frame <- function(msd)
{
  counts_table <- get_counts_table(msd)
  sensitivity  <- estimate_sensitivity(counts_table)
  sensitivity
}

estimate_sensitivity.matrix <- function(counts_table)
{
  # Normalize to avoid infinite values
  counts_table <- counts_table

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

estimate_sensitivity.table <- function(counts_table)
{  
  sensitivity <- estimate_sensitivity.matrix(counts_table)
  sensitivity
}
