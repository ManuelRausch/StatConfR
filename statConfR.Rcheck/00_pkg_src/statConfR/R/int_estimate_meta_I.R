estimate_meta_I <- function(estimated_table) {
  
  # Normalize
  estimated_classifier <- estimated_table/sum(estimated_table)

  p <- rowSums(estimated_classifier)         # Prior
  a <- get_accuracy(estimated_classifier)    # Accuracy
  
  info               <- get_information(estimated_classifier) # Transmitted information
  information_bounds <- get_analytic_information_bounds(p, a)
  lower_bound        <- information_bounds$lowest # Lower bound on trans. information

  meta_I <- info - lower_bound
  meta_I
}
