get_pointwise_transmitted_information <- pti <- function(c, prior){
  H(prior) - H(c)
}
