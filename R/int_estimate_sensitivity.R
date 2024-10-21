estimate_sensitivity <- function(estimated_classifier){

  # Count number of predictions
  n_responses <- ncol(estimated_classifier)
  idx_predicting_first_row  <- 1:(n_responses/2)
  idx_predicting_second_row <- (n_responses/2 + 1):n_responses

  # Count number of correct predictions for first stimulus
  n_correctly_predict_first_row <- sum(estimated_classifier[1, idx_predicting_first_row])
  n_predict_first_row           <- sum(estimated_classifier[ , idx_predicting_first_row])

  # Count number of correct predictions for second stimulus
  n_correctly_predict_second_row <- sum(estimated_classifier[2, idx_predicting_second_row])
  n_predict_second_row           <- sum(estimated_classifier[ , idx_predicting_second_row])

  # Calculcate sensitivity
  hit_rate         <- n_correctly_predict_first_row / n_predict_first_row
  false_alarm_rate <- 1 - ( n_correctly_predict_second_row / n_predict_second_row )
  sensitivity      <- qnorm(hit_rate) - qnorm(false_alarm_rate)
  sensitivity
}
