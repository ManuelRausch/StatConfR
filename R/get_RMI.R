get_RMI <- function(classifier){
  prior       <- rowSums(classifier)
  accuracy    <- get_accuracy(classifier)
  information <- get_information(classifier)

  RMI <- normalize_information(prior, accuracy, information)

  RMI
}
