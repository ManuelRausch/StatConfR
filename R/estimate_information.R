estimate_information <- function(x)
{
  UseMethod("estimate_information")
}


#'@export
estimate_information.data.frame <- function(msd)
{
  estimated_classifier <- estimate_classifier(msd)

  info <- get_information(estimated_classifier)
  info
}

#'@export
estimate_information.matrix <- function(tab)
{
  estimated_classifier <- tab/sum(tab)

  info <- get_information(estimated_classifier)
  info
}

#'@export
estimate_information.table <- function(tab)
{
  estimated_classifier <- tab/sum(tab)

  info <- get_information(estimated_classifier)
  info
}
