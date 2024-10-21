get_information <- function(classifier)
{
  stopifnot( round( sum(classifier) - 1, 6) == 0 )
  info <- H(colSums(classifier)) + H(rowSums(classifier)) - H(classifier)
  info
}
