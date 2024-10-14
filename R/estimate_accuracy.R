# This function is not getting documented because I think that the use case to ignore
# the mapping intended by the classifier so rare that the chance is much higher that
# someone will mess up using this function! The function is available internally, though.
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
