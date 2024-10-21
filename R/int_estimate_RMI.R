estimate_RMI <- function(estimated_classifier) {
  prior       <- rowSums(estimated_classifier)         # Prior
  accuracy    <- get_accuracy(estimated_classifier)    # Accuracy
  information <- get_information(estimated_classifier) # Information

  information_bounds <- get_analytic_information_bounds(prior, accuracy)
  lower_bound <- information_bounds$lowest
  upper_bound <- information_bounds$highest

  RMI <- ( information - lower_bound ) /
         ( upper_bound - lower_bound )

  # Where bounds collapse because of accuracy edge cases, return NaN
  RMI[lower_bound == upper_bound] <- NaN

  RMI
}

  