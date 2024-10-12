#' Estimate Classifier
#' 
#' Construct a contingency table as an estimate for the classifier's joint
#' probability distribution between true labels and responses.
#' 
#' @param dat Data frame with true labels representing stimulus categories
#'  (variable "y") and confidence-binned responses (variable "r")
#' 
#' @return Estimated classifier represented as a contingency table with true
#'  labels as rows and responses as columns
estimate_classifier <- function(dat)
{
  counts_table         <- get_counts_table(dat)
  estimated_classifier <- counts_table / nrow(dat)
  estimated_classifier
}
