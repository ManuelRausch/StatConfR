# R cmd check throws a warning if S3generic and S3methods do not accept identical arguments.
# This may cause problems when submitting the package to CRAN later.

estimate_accuracy <- function(x){
  UseMethod("estimate_accuracy")
}

#' @export
estimate_accuracy.matrix <- function(x)
{
  estimated_classifier <- x/sum(x)

  acc <- get_accuracy(estimated_classifier)
  acc
}
#'@export
estimate_accuracy.data.frame <- function(x)
{
  estimated_classifier <- estimate_classifier(x)

  acc <- estimate_accuracy.matrix(estimated_classifier)
  acc
}
#'@export
estimate_accuracy.table <- function(x)
{
  estimated_classifier <- x/sum(x)

  acc <- estimate_accuracy.matrix(estimated_classifier)
  acc
}
