estimate_meta_Ir1 <- function(estimated_table) {

  # Normalize
  estimated_classifier <- estimated_table/sum(estimated_table)

  p <- rowSums(estimated_classifier)         # Prior
  d <- estimate_sensitivity(estimated_table) # Sensitivity
  a <- pnorm(abs(d)/2)                       # Accuracy

  meta_I <- estimate_meta_I(estimated_classifier) # Meta-I
  
  information_bounds <- get_analytic_information_bounds(p, a)
  lower_bound        <- information_bounds$lowest # Lower bound on trans. information
  info_normal        <- get_normal_noise_information(sensitivities = d)$info
  meta_I_normal      <- info_normal - lower_bound # Gaussian Meta-I

  # Normalization
  meta_Ir1 <- meta_I / meta_I_normal
  meta_Ir1
}