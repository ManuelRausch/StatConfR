estimate_meta_Ir2 <- function(estimated_classifier){

  meta_I <- estimate_meta_I(estimated_classifier) # Meta-I

  a          <- get_accuracy(estimated_classifier) # Accuracy
  H_accuracy <- H2(a)                              # Remaining entropy

  meta_I_r2 <- meta_I / H_accuracy
  meta_I_r2
}
