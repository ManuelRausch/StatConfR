#' Compute accuracy of a classifier
#'
#' The accuracy is the sum of highest entries for each response in the
#' classifier's contingency matrix. 
#' 
#' @param classifier A contingency matrix with joint probabilities. Rows
#'  correspond to true labels (stimulus categories) and columns correspond to
#'  responses.
#'
#' @return Accuracy of the classifier
get_accuracy <- function(classifier)
{
  accuracy <- sum( apply(classifier, 2, max) )
  accuracy
}
