# This function is not getting documented because I think that the use case to ignore
# the mapping intended by the classifier so rare that the chance is much higher that
# someone will mess up using this function! The function is available internally, though.

# Estimate Accuracy
#
# Estimate the probability of correct predictions. Note that this ignores the
# mapping intended by the classifier and instead uses the optimal
# prediction-mapping for each response.
#
# Data can be input in three ways:
# - A data frame with variables "y" for true labels and "r" for
#   confidence-binned responses. "y" needs to contain values -1 and +1 while
#   r needs to be factor with levels ordered such that the first half of the
#   factor levels are predictions for -1 and the second half for
#   predictions +1. For example with 4 factor levels: 1: high confidence
#   prediction for y=-1, 2: low confidence prediction for y=-1, 3: low
#   confidence prediction for y=+1, 4: high confidence prediction for y=+1.
# - A counts table with joint absolute frequencies. Rows correspond to true
#   labels (stimulus categories) and columns correspond to responses.
# - A contingency matrix with joint relative frequencies (as before but
#   normalized to sum up to 1).
#
# @param x Data
#
# @return Estimated accuracy

estimate_accuracy <- function(x)
{
  UseMethod("estimate_accuracy")
}

#' @export
estimate_accuracy.matrix <- function(estimated_classifier)
{
  estimated_classifier <- estimated_classifier/sum(estimated_classifier)

  acc <- get_accuracy(estimated_classifier)
  acc
}
#'@export
estimate_accuracy.data.frame <- function(msd)
{
  estimated_classifier <- estimate_classifier(msd)

  acc <- estimate_accuracy.matrix(estimated_classifier)
  acc
}
#'@export
estimate_accuracy.table <- function(counts_table)
{
  estimated_classifier <- counts_table/sum(counts_table)

  acc <- estimate_accuracy.matrix(estimated_classifier)
  acc
}
