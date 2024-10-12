#' Compute transmitted information of a classifier
#'
#' Compute the transmitted information for a given classifier.
#' 
#' @param classifier A contingency matrix with joint probabilities. Rows
#'  correspond to true labels (stimulus categories) and columns correspond to
#'  responses.
#'
#' @return Transmitted information of the classifier's output about the label
get_information = function(classifier)
{
  stopifnot( round( sum(classifier) - 1, 6) == 0 )
  info <- H(colSums(classifier)) + H(rowSums(classifier)) - H(classifier)
  info
}
