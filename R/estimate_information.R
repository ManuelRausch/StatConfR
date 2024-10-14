# R cmd check throws a warning if S3generic and S3methods do not accept identical arguments.
# This may cause problems when submitting the package to CRAN later.

estimate_information <- function(x){
  UseMethod("estimate_information")
}


#'@export
estimate_information.data.frame <- function(x){
  estimated_classifier <- estimate_classifier(x)

  info <- get_information(estimated_classifier)
  info
}

#'@export
estimate_information.matrix <- function(x){
  estimated_classifier <- x/sum(x)

  info <- get_information(estimated_classifier)
  info
}

#'@export
estimate_information.table <- function(x){
  estimated_classifier <- x/sum(x)

  info <- get_information(estimated_classifier)
  info
}
