get_accuracy <- function(classifier)
{
  accuracy <- sum( apply(classifier, 2, max) )
  accuracy
}
