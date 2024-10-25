estimate_meta_Ir1_acc <- function(estimated_table) {
  
  # Normalize
  estimated_classifier <- estimated_table/sum(estimated_table)

  p <- rowSums(estimated_classifier)      # Prior
  a <- get_accuracy(estimated_classifier) # Sensitivity
  d <- 2*qnorm(a)                         # Accuracy

  meta_I <- estimate_meta_I(estimated_classifier) # Meta-I

  information_bounds <- get_analytic_information_bounds(p, a)
  lower_bound <- information_bounds$lowest
  info_normal   <- get_normal_noise_information(sensitivities = d)$info
  meta_I_normal <- info_normal - lower_bound # Gaussian Meta-I

  # Normalization
  meta_Ir1_acc <- meta_I / meta_I_normal
  meta_Ir1_acc
}
