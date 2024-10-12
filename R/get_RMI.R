#' Get Relative Meta-Information
#' 
#' Compute the meta-information as a relative measure in the possible range
#' given the accuracy: Highest information for a given accuracy produces 
#' RMI = 1, lowest information for a given accuracy produces RMI = 0.
#' 
#' @param classifier A contingency matrix with joint probabilities. Rows
#'  correspond to true labels (stimulus categories) and columns correspond to
#'  responses.
#' 
#' @return Relative meta-information value
get_RMI <- function(classifier)
{
  prior       <- rowSums(classifier)
  accuracy    <- get_accuracy(classifier)
  information <- get_information(classifier)
  
  RMI <- normalize_information(prior, accuracy, information)

  RMI
}
